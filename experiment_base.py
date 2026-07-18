# experiment_base.py
"""
DiodeExperiment — works with the new ExperimentConfig dataclass.

Config access: uses config.scans, config.mw, config.seq, etc. directly.
No backward-compat dict path — new config format only.

Loop order (outer → Nruns → inner):
    for each outer_combo (e.g. mw_power):
        for i_run in range(Nruns):
            for each inner_pt (e.g. freq):
                acquire

Data shape: (outer_combos, Nruns, inner_pts, daq_Nsamples * n_channels)

Execution modes (controlled by RuntimeFlags):
    Normal   (reload_pb=True,  load_all_params=False): Per-point PB reprogram + DAQ read.
    Fixed-PB (reload_pb=False, load_all_params=False): PB once per outer combo; setter only per point.
    Batch    (load_all_params=True):                   All inner pts compiled into one PB program;
                                                        single DAQ read per inner sweep.

Validation:
    exp.validate() checks the full Cartesian product of (outer × inner) PB timing
    parameters before execution.  Invalid cells are NaN-filled; the grid stays rectangular.

Example — sweep t_AOM as a config (shape-changing) parameter:

    from experiment_base import DiodeExperiment
    from spinapi import ms

    exp = DiodeExperiment(instruments, config)

    def apply_t_aom(val):
        config.times.t_AOM = val
        config.seq.args_values[0] = val
        new_N = int((val - 40*ms) / ms * 1e-3 * 10e3)
        config.seq.Nsamples = max(1, new_N)
        config.daq_ai.ai_samps_per_chan = config.seq.Nsamples

    results = exp.config_sweep(
        config_name='t_AOM',
        config_values=[10*ms, 30*ms, 70*ms, 100*ms],
        config_applier=apply_t_aom,
    )
"""

import numpy as np, time, logging, yaml, json, h5py
from pathlib import Path
from itertools import product
from typing import Optional, Callable, TYPE_CHECKING
from dataclasses import dataclass, field

from PBcontrol import PulseBlaster
from DAQcontrol import AnalogInputTask, AnalogOutputTask, CounterInputTask
from SGcontrol import SignalGenerator, SignalGenerator_sim
import sequencecontrol as seqctrl
from spinapi import ns, us, ms, Inst
from sweep_utils import Sweep, SweepIndex, make_pb_setter
from experiment_config import ExperimentConfig, DAQ_INFO

def printt(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")

# ═══════════════════════════════════════════════════════════════════════
#  Validation result container
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class ValidationResult:
    """Result of validate(): per-cell validity over the full sweep grid.

    Attributes
    ----------
    valid_mask : ndarray, shape (outer_count, inner_count), dtype bool
        True  = PB program compiled successfully for that cell.
        False = timing violation; cell will be NaN-filled at run-time.
    pb_errors : list of (outer_idx, inner_idx, detail_dict) tuples
    hw_errors : list of (param_name, value_str, reason) tuples
    n_total   : int — total cells checked
    n_invalid : int — cells that failed PB compilation
    """
    valid_mask: np.ndarray
    pb_errors: list = field(default_factory=list)
    hw_errors: list = field(default_factory=list)
    n_total: int = 0
    n_invalid: int = 0

    @property
    def all_valid(self) -> bool:
        return self.n_invalid == 0 and len(self.hw_errors) == 0

    def summary(self) -> str:
        lines = [f"Validation: {self.n_total} cells checked, "
                 f"{self.n_invalid} PB failures"]
        if self.hw_errors:
            lines.append(f"  {len(self.hw_errors)} hardware-bound violation(s):")
            for name, val, reason in self.hw_errors:
                lines.append(f"    {name}={val}: {reason}")
        if self.n_invalid > 0:
            for oc in range(self.valid_mask.shape[0]):
                bad = np.where(~self.valid_mask[oc])[0]
                if len(bad) > 0:
                    lines.append(
                        f"  outer[{oc}]: {len(bad)} bad inner pts "
                        f"(indices {bad[:8].tolist()}"
                        f"{'...' if len(bad) > 8 else ''})")
        return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
#  Shared helpers
# ═══════════════════════════════════════════════════════════════════════

def calc_contrast(signal, reference, op='s/r'):
    ops = {
        '+-':  lambda s, r: s - r,  '-+':  lambda s, r: r - s,
        's/r': lambda s, r: s / r,  'r/s': lambda s, r: r / s,
    }
    return ops.get(op, lambda s, r: (s - r) / (s + r))(signal, reference)


def read_details(sequence: str, n_channels: int, Nsamples_cfg: int):
    if sequence in ['modesr']:
        rpc = [4, 1] if n_channels == 2 else [4]
    elif sequence in ['drift_seq', 'T1ms0_train']:
        rpc = [1, 1] if n_channels == 2 else [1]
    elif 'lia' in sequence:
        rpc = [1, 1] if n_channels == 2 else [1]
    else:
        rpc = [2, 2] if n_channels == 2 else [2]
    return rpc, sum(rpc) * Nsamples_cfg


def process_data(i_max, data_row, Nsamples, sequence):
    # sig = data_row[:i_max, 0:int(Nsamples/2)]
    sig = data_row[:i_max, 0::Nsamples]
    mean_sig = sig.mean(axis=-1)
    has_ref = sequence.lower() not in ['t1ms0_train']
    if has_ref:
        # ref = data_row[:i_max, int(Nsamples/2):]
        ref = data_row[:i_max, 1::Nsamples]
        mean_ref = ref.mean(axis=-1)
        mean_contrast = calc_contrast(mean_sig, mean_ref)
    
    return (mean_sig,) if not has_ref else (mean_sig, mean_ref, mean_contrast)


def savefile(filename: Path, data) -> bool:
    def _convert(var):
        if isinstance(var, np.ndarray):  return var.tolist()
        if isinstance(var, np.generic):  return var.item()
        if isinstance(var, dict):        return {k: _convert(v) for k, v in var.items()}
        if isinstance(var, list):        return [_convert(v) for v in var]
        if isinstance(var, Path):        return str(var)
        return var
    try:
        s = filename.suffix.lower()
        if s == '.npy':
            np.save(filename, data, allow_pickle=False)
        elif s == '.yaml':
            with open(filename, 'w') as f:
                yaml.dump(_convert(data), f, indent=4, default_flow_style=False)
        elif s in ('.h5', '.hdf5'):
            with h5py.File(filename, 'a') as f:
                g = f.require_group("data")
                for idx in range(data.shape[0]):
                    grp = g.require_group(f"outer_{idx}")
                    for r in range(data.shape[1]):
                        dk = f"run_{r}"
                        if dk in grp: del grp[dk]
                        grp.create_dataset(dk, data=data[idx, r])
        else:
            np.save(filename.with_suffix('.npy'), data, allow_pickle=False)
        return True
    except Exception as e:
        logging.exception(f"Save failed: {filename}: {e}")
        return False


# ═══════════════════════════════════════════════════════════════════════
#  Hardware bounds for instrument-setting parameters
# ═══════════════════════════════════════════════════════════════════════

# (min_val, max_val, description) — extend as instruments are added
HW_BOUNDS = {
    'freq':      (0.0,        8.1e9,  'SG384 RF frequency (Hz)'),
    'mw_freq':   (0.0,        8.1e9,  'SG384 RF frequency (Hz)'),
    'mw_power':  (-110.0,     8.99,   'SG384 RF power (dBm) — amp limit'),
    'mod_rate':   (0.0,       2.0e6,  'SG384 modulation rate (Hz)'),
    'mod_dev':    (0.0,       50e6,   'SG384 modulation deviation (Hz)'),
}


# ═══════════════════════════════════════════════════════════════════════
#  DiodeExperiment
# ═══════════════════════════════════════════════════════════════════════

class DiodeExperiment:
    """
    Stateful experiment for diode-detector NV measurements.

    Accepts an ExperimentConfig dataclass (from experiment_config.py).
    All intermediate state exposed as attributes.
    """

    # ── setter registry: name → lambda(self) → callable ────────────────
    SETTER_REGISTRY = {
        'freq':      lambda self: self.sg.set_freq,
        'mw_freq':   lambda self: self.sg.set_freq,
        'mw_power':  lambda self: self.sg.set_amp_rf,
        'mod_rate':  lambda self: self.sg.set_mod_rate,
        'mod_dev':   lambda self: self.sg.set_mod_dev,
    }

    FREQ_ONLY_SEQUENCES = frozenset(['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr',
                                     'drift_seq'])

    # Parameters that change instrument settings (not PB timing)
    INSTRUMENT_PARAMS = frozenset(['freq', 'mw_freq', 'mw_power', 'mod_rate', 'mod_dev'])

    def __init__(self, instruments: dict, config: 'ExperimentConfig'):
        self.sg: SignalGenerator = instruments['sg']
        self.pb: PulseBlaster = instruments['pb']
        self.ao_task: AnalogOutputTask|None = instruments.get('ao_task')

        self.cfg = config
        self.instr = 'diode'

        # State (populated by setup / execute)
        self.sweep: Sweep|None = None
        self.input_task: AnalogInputTask|CounterInputTask|None = None
        self.data_array: np.ndarray|None = None
        self.daq_Nsamples: int = 0
        if self.cfg.runtime.detector == 'counter':
            # assert config.daq_ci is not None
            self.n_channels = len(config.daq_ci.channel)
        else:
            # assert config.daq_ai is not None
            self.n_channels = len(config.daq_ai.channels)
        self.scan_times: list[float] = []
        self.elapsed: float = 0.0
        self.i_outer: int = 0
        self.i_run: int = 0
        self.i_scanpt: int = 0

        # Derived from config
        # assert config.seq.name and config.seq.Nsamples and config.pb.channels and config.seq.Ncycles
        self.sequence: str = config.seq.name
        self.Nsamples_cfg: int = config.seq.Nsamples
        self.Ncycles: int = config.seq.Ncycles
        self.Nruns: int = config.runtime.Nruns
        self.pb_channels: list = config.pb.channels

        # Execution mode flags
        self._reload_pb: bool = config.runtime.reload_pb
        self._load_all_params: bool = config.runtime.load_all_params

        # Primary scan
        scan_names = config.scan_names
        self.scan_name: str = scan_names[0] if scan_names else ''
        self.primary_scan = config.scans.get(self.scan_name)
        self.param_values: np.ndarray = (
            self.primary_scan.values if self.primary_scan else np.array([]))
        self.Nscanpts: int = len(self.param_values)

        self._is_freq_sweep = self.sequence in self.FREQ_ONLY_SEQUENCES

        # Validation state
        self._validation: Optional[ValidationResult] = None
        self._valid_mask: Optional[np.ndarray] = None

    # ══════════════════════════════════════════════════════════════════
    #  VALIDATE — full Cartesian product check
    # ══════════════════════════════════════════════════════════════════

    def validate(self, plot_errors: bool = False) -> ValidationResult:
        """Check every (outer_combo, inner_pt) cell before execution.

        1. Hardware bounds for instrument-setting params (freq, power, …)
        2. PB timing validity via PB_program(…, err_check=True) over
           the full outer×inner Cartesian product.

        The grid stays rectangular.  Invalid cells are marked False in
        valid_mask and NaN-filled during execute().

        Parameters
        ----------
        plot_errors : bool
            If True, matplotlib-plot every failing PB sequence
            (useful in notebooks, disable for GUI runs).

        Returns
        -------
        ValidationResult
            Also stored as self._validation / self._valid_mask.
        """
        cfg = self.cfg
        hw_errors = []

        # ── 1. Hardware bounds for every scan axis ────────────────────
        for name, axis in cfg.scans.items():
            if name in HW_BOUNDS:
                lo, hi, desc = HW_BOUNDS[name]
                vals = axis.values
                below = vals[vals < lo]
                above = vals[vals > hi]
                if len(below) > 0:
                    hw_errors.append((
                        name, f"min={below.min():.6g}",
                        f"Below lower bound {lo:.6g} ({desc})"))
                if len(above) > 0:
                    hw_errors.append((
                        name, f"max={above.max():.6g}",
                        f"Above upper bound {hi:.6g} ({desc})"))

        # Static MW power check
        if self.sequence not in ['aom_timing', 'T1ms0_train', 'drift_seq']:
            # assert cfg.mw
            if cfg.mw.power >= 9.0:
                hw_errors.append((
                    'mw.power', f"{cfg.mw.power} dBm",
                    "MW power >= 9 dBm — risk of damaging RF amplifier"))

        # ── 2. Determine grid dimensions ─────────────────────────────
        inner_count = len(self.param_values)
        outer_count = self._calc_outer_count()

        # assert cfg.seq.args_values
        # ── 3. Freq-only sequences: PB is fixed ─────────────────────
        if self._is_freq_sweep:
            valid_mask = np.ones((outer_count, inner_count), dtype=bool)
            # Compile once to verify
            sal = list(cfg.seq.args_values)
            _, the_list = PulseBlaster.PB_program(
                self.instr, self.sequence,
                sal + [self.pb_channels], err_check=True)
            for i in range(len(the_list)):
                if the_list[i][1] > 0:
                    valid_mask[:, :] = False
                    hw_errors.append((
                        'pb_static', 'N/A',
                        f"PB program has {the_list[i][1]} timing errors "
                        f"at instructions {the_list[i][2]}"))
            return self._store_result(valid_mask, [], hw_errors)

        # ── 4. PB timing sweep: full Cartesian product ───────────────
        args_names = cfg.seq.args_names; # assert args_names
        seq_args_base = list(cfg.seq.args_values)
        scan_names = cfg.scan_names
        inner_name = scan_names[0] if scan_names else ''

        # Collect outer axes with their PB/instrument classification
        outer_axes = []  # (name, values, is_pb_param)
        for oname in scan_names[1:]:
            oaxis = cfg.scans.get(oname)
            if oaxis is None or len(oaxis.values) <= 1:
                continue
            outer_axes.append((
                oname, oaxis.values,
                oname not in self.INSTRUMENT_PARAMS))

        # Build outer combo index grid
        if outer_axes:
            outer_combos = list(product(
                *(range(len(v)) for _, v, _ in outer_axes)))
        else:
            outer_combos = [()]

        outer_count = len(outer_combos)
        valid_mask = np.ones((outer_count, inner_count), dtype=bool)
        pb_errors = []

        inner_is_pb = inner_name not in self.INSTRUMENT_PARAMS

        # Inner param's index in seq_args (-1 = legacy prepend)
        inner_pidx = (args_names.index(inner_name)
                      if inner_is_pb and inner_name in args_names
                      else -1 if inner_is_pb else None)
        # assert inner_pidx
        printt(f"🔃 Validating {outer_count} × {inner_count} = "
              f"{outer_count * inner_count} cells...")

        for oc_flat, oc_indices in enumerate(outer_combos):
            # Build seq_args template with this outer combo's PB values
            template = list(seq_args_base)
            for k, idx_val in enumerate(oc_indices):
                oname, ovals, is_pb = outer_axes[k]
                if is_pb and oname in args_names:
                    template[args_names.index(oname)] = ovals[idx_val]

            if not inner_is_pb:
                # Inner is instrument param — PB same for all inner pts.
                # Check once for this outer combo.
                _, the_list = PulseBlaster.PB_program(
                    self.instr, self.sequence,
                    list(template) + [self.pb_channels], err_check=True)
                for i in range(len(the_list)):
                    if the_list[i][1] > 0:
                        valid_mask[oc_flat, :] = False
                        pb_errors.append((oc_flat, -1, {
                            'error_count': the_list[i][1],
                            'inst_errors': the_list[i][2],
                            'error_times': the_list[i][3]}))
                continue

            # Inner IS a PB timing param — check every cell
            for j, inner_val in enumerate(self.param_values):
                test_args = list(template)
                if inner_pidx == -1:
                    test_args = [inner_val] + test_args
                else:
                    test_args[inner_pidx] = inner_val

                _, the_list = PulseBlaster.PB_program(
                    self.instr, self.sequence,
                    test_args + [self.pb_channels], err_check=True)

                for i in range(len(the_list)):
                    if the_list[i][1] > 0:
                        valid_mask[oc_flat, j] = False
                        pb_errors.append((oc_flat, j, {
                            'inner_val': inner_val,
                            'error_count': the_list[i][1],
                            'inst_errors': the_list[i][2],
                            'error_times': the_list[i][3]}))

                        if plot_errors:
                            import matplotlib.pyplot as plt
                            inst_list = the_list[i][0]
                            plt.figure()
                            t_us, chP, yT = seqctrl.plot_sequence(
                                inst_list, self.pb_channels)
                            for ch in chP:
                                plt.plot(t_us, list(ch))
                            plt.xlabel('Time (us)')
                            plt.ylabel('Channel')
                            plt.title(
                                f'Err: outer[{oc_flat}] '
                                f'inner={inner_val:.4g}\n'
                                f'Bad inst: {the_list[i][2]}',
                                color='r', fontsize=10)

        return self._store_result(valid_mask, pb_errors, hw_errors)

    def _store_result(self, valid_mask, pb_errors, hw_errors):
        """Package and store a ValidationResult."""
        n_total = valid_mask.size
        n_invalid = int(np.sum(~valid_mask))
        result = ValidationResult(
            valid_mask=valid_mask, pb_errors=pb_errors,
            hw_errors=hw_errors, n_total=n_total, n_invalid=n_invalid)
        self._validation = result
        self._valid_mask = valid_mask
        printt(result.summary())
        return result

    def _calc_outer_count(self) -> int:
        """Count outer combos from config (before Sweep is built)."""
        n = 1
        for name in self.cfg.scan_names[1:]:
            axis = self.cfg.scans.get(name)
            if axis is not None and len(axis.values) > 1:
                n *= len(axis.values)
        return max(n, 1)

    # ══════════════════════════════════════════════════════════════════
    #  SETUP
    # ══════════════════════════════════════════════════════════════════
    def setup(self):
        cfg = self.cfg

        # assert cfg.times
        # Per-sequence total time (used for DAQ timeout)
        self._t_total_s = (cfg.times.t_tot / 1e9 if cfg.times.t_tot > 0
                           else 2 * cfg.times.t_AOM / 1e9)

        # assert cfg.seq.args_values
        self.seq_args = list(cfg.seq.args_values)

        # ── Validate if not already done ─────────────────────────────
        if self._validation is None:
            self.validate(plot_errors=True)

        # Fatal HW errors
        if self._validation and self._validation.hw_errors:
            for name, val, reason in self._validation.hw_errors:
                if 'MW power >= 9' in reason:
                    printt(f"❌ FATAL: {reason}. Aborting.")
                    raise ValueError(reason)
                printt(f"⚠ HW bound: {name}={val}: {reason}")

        # DAQ sizing
        self.daq_Nsamples = self.Nsamples_cfg * self.Ncycles

        # ── Initial PB program ───────────────────────────────────────
        # assert cfg.seq.args_names
        if self._is_freq_sweep:
            sal = list(self.seq_args)
        elif self.scan_name in cfg.seq.args_names:
            sal = list(self.seq_args)
        else:
            sal = [self.param_values[0]] + self.seq_args

        _, the_list = PulseBlaster.PB_program(
            self.instr, self.sequence, sal + [self.pb_channels])
        self.pb.run_sequence_for_diode(
            [the_list[i][0] for i in range(len(the_list))])

        # ── Build Sweep ──────────────────────────────────────────────
        self.sweep = Sweep()
        self._build_sweeps(cfg)

        printt(self.sweep.info())

        # SG
        if self.sequence not in ['aom_timing', 'T1ms0_train', 'drift_seq']:
            # assert cfg.mw
            self.sg.set_freq(cfg.mw.freq)
            self.sg.set_amp_rf(cfg.mw.power)
            self.sg.enable_ntype(1)
            self.sg.setup_ext_pulse_mod()

        # AO test field
        tf = cfg.extra.get('test_field', [0, 0, 0])
        if self.ao_task is not None:
            self.ao_task.set_outputs_to_constant(output_field_in_gauss=tf)

        # AI / CI task (per-point read size; batch mode re-creates in execute)
        if cfg.runtime.detector == 'counter':
            ci = cfg.daq_ci
            # assert ci
            # assert ci.start_trigger_source and ci.pause_trigger_source
            self.input_task = CounterInputTask(
                dev = 'P6363',
                channel = ci.channel,
                counter = ci.counter,
                sample_source = ci.sample_source,
                sample_rate = ci.sample_rate,
                sample_mode = ci.sample_mode,
                samps_per_chan = int(self.daq_Nsamples),
                start_trigger_source = ci.start_trigger_source,
                pause_trigger_source = ci.pause_trigger_source,
            )
            self.n_channels = 1  # counter is always single-channel
            printt(f"▶ Detector: COUNTER ({ci.counter}, input={ci.channel})")
        else:
            ai = cfg.daq_ai
            # assert ai
            # assert ai.voltage_ranges and ai.sample_source and ai.start_trigger_source and ai.pause_trigger_source
            self.input_task = AnalogInputTask(
                channels=ai.channels,
                voltage_ranges=ai.voltage_ranges,
                sample_source=ai.sample_source,
                sample_rate=ai.sample_rate,
                sample_mode=ai.sample_mode,
                samps_per_chan=int(self.daq_Nsamples),
                start_trigger_source=ai.start_trigger_source,
                pause_trigger_source=ai.pause_trigger_source,
            )
            print(f"▶ Detector: ANALOG (ai{ai.channels})")

        # Data array
        self.data_array = np.zeros((self.sweep.outer_count, self.Nruns,
                                    self.sweep.inner_count, self.daq_Nsamples * self.n_channels,
                                    ))

        # Shuffle
        shuffle = self.primary_scan.shuffle if self.primary_scan else False
        self._shuffle = shuffle
        if shuffle:
            self._inner_order = np.random.permutation(self.sweep.inner_count)
            print("⚡ Inner sweep shuffled")
        else:
            self._inner_order = np.arange(self.sweep.inner_count)

        # Verify valid_mask shape matches sweep
        if self._valid_mask is not None:
            expected = (self.sweep.outer_count, self.sweep.inner_count)
            if self._valid_mask.shape != expected:
                print("⚠ valid_mask shape mismatch — re-validating")
                self.validate()

        # Print mode
        if self._load_all_params:
            mode_str = "BATCH (single PB + single DAQ read per sweep)"
        elif not self._reload_pb:
            mode_str = "FIXED-PB (PB once per outer; setter per point)"
        else:
            mode_str = "NORMAL (per-point PB reprogram + DAQ read)"
        print(f"▶ Mode: {mode_str}")
        print(f"▶ data {self.data_array.shape}  "
              f"({self.sweep.outer_count} outer × {self.Nruns} runs × "
              f"{self.sweep.inner_count} inner)")

    def _build_sweeps(self, cfg):
        """Add inner + outer axes to self.sweep.

        KEY DESIGN: All PB-timing setters share ONE template list.

        When the inner scan param uses legacy-prepend mode (scan name is
        NOT in args_names, e.g. Rabi 'tau' prepended at index 0), the
        template is [inner_val, *args_values].  All outer PB setters
        must use this same extended template so their param_index is
        shifted by +1 to account for the prepend slot.

        When the inner IS in args_names, the template is just
        list(args_values) — no shift needed.

        Instrument-setting axes (freq, power) don't touch PB at all;
        they get their own setter from SETTER_REGISTRY.
        """
        args_names = cfg.seq.args_names

        # ── Determine template layout ────────────────────────────────
        # Two cases:
        #   (a) Inner scan param IS in args_names → template = args_values
        #       Inner modifies template[inner_pidx].
        #   (b) Inner scan param is NOT in args_names (legacy prepend)
        #       → template = [inner_val] + args_values
        #       Inner modifies template[0], all args_names indices shift +1.
        inner_in_args = self.scan_name in args_names
        if inner_in_args:
            self._pb_template = list(self.seq_args)
            self._pb_prepend = False
            inner_pidx = args_names.index(self.scan_name)
            self._pb_template[inner_pidx] = self.param_values[0]
            args_offset = 0   # no shift
        else:
            self._pb_template = [self.param_values[0]] + self.seq_args
            self._pb_prepend = True
            inner_pidx = 0    # prepend slot
            args_offset = 1   # everything in args_names shifts by 1

        # ── Inner axis ───────────────────────────────────────────────
        # assert self.sweep
        if self._is_freq_sweep:
            self.sweep.add(self.scan_name, self.param_values,
                           setter=self.sg.set_freq)
        elif self._load_all_params:
            self.sweep.add(self.scan_name, self.param_values, setter=None)
        elif not self._reload_pb:
            if self.scan_name in self.INSTRUMENT_PARAMS:
                self.sweep.add(self.scan_name, self.param_values,
                               setter=self._resolve_setter(self.scan_name))
            else:
                printt("⚠ Fixed-PB with PB-timing inner param: PB programmed once with first value only.")
                self.sweep.add(self.scan_name, self.param_values,
                               setter=None)
        else:
            # Normal mode: per-point PB reprogram
            pb_setter = make_pb_setter(
                self.pb, param_index=inner_pidx, instr=self.instr,
                sequence=self.sequence,
                seq_args_template=self._pb_template,
                pb_channels=self.pb_channels)
            self.sweep.add(self.scan_name, self.param_values,
                           setter=pb_setter)

        # ── Outer axes ───────────────────────────────────────────────
        for name in cfg.scan_names[1:]:
            axis = cfg.scans.get(name)
            if axis is None or len(axis.values) <= 1:
                continue

            if name not in self.INSTRUMENT_PARAMS and name in args_names:
                # Outer PB timing param — use the SAME shared template
                pidx = args_names.index(name) + args_offset
                pb_setter = make_pb_setter(
                    self.pb, param_index=pidx, instr=self.instr,
                    sequence=self.sequence,
                    seq_args_template=self._pb_template,
                    pb_channels=self.pb_channels)
                self.sweep.add(name, axis.values, setter=pb_setter)
            else:
                self.sweep.add(name, axis.values,
                               setter=self._resolve_setter(name))

    def _resolve_setter(self, name: str) -> Optional[Callable]:
        factory = self.SETTER_REGISTRY.get(name)
        if factory is not None:
            return factory(self)
        printt(f"⚠ No setter for '{name}'")
        return None

    # ══════════════════════════════════════════════════════════════════
    #  VIEW SEQUENCES — plot PB timing at selected scan points
    # ══════════════════════════════════════════════════════════════════

    def view_sequences(self, indices: list[int] = [], dpi: int = 100) -> list:
        """Plot PB pulse sequences at selected inner-param values.

        Builds the correct seqArgList from config (handling legacy
        prepend vs in-args), then delegates to sequencecontrol.view_sequence.

        Parameters
        ----------
        indices : list[int], optional
            Which indices into the inner param_values array to plot.
            Supports negative indexing ([-1] = last point).
            Default: from config.runtime.seq_plot_indices, or [0, -1].
        dpi : int
            Matplotlib figure DPI.

        Can be called before or after setup().  Does NOT mutate any
        experiment state — uses copies of all args.
        """
        cfg = self.cfg

        if indices is []:
            indices = getattr(cfg.runtime, 'seq_plot_indices', [0, -1])

        # Build the seqArgList in the same layout PB_program expects
        args_names = cfg.seq.args_names
        # assert args_names and cfg.seq.args_values
        seq_args_base = list(cfg.seq.args_values)

        if self._is_freq_sweep:
            # PB program is fixed for freq sweeps — just plot once
            full_args = seq_args_base + [self.pb_channels]
            return seqctrl.view_sequence(
                instr=self.instr, sequence=self.sequence, seq_args=full_args,
                param_values=None, indices=[0], dpi=dpi)

        # Determine inner param position
        if self.scan_name in args_names:
            # In-args: param lives at its args_names index
            pidx = args_names.index(self.scan_name)
            full_args = seq_args_base + [self.pb_channels]
        else:
            # Legacy prepend: param goes at index 0
            pidx = 0
            full_args = [self.param_values[0]] + seq_args_base + [self.pb_channels]

        return seqctrl.view_sequence(
            instr=self.instr, sequence=self.sequence, seq_args=full_args,
            param_values=self.param_values, indices=indices, param_index=pidx, dpi=dpi)

    # ══════════════════════════════════════════════════════════════════
    #  PLOT + VIEW SEQUENCES — play sequence at selected scan point
    # ══════════════════════════════════════════════════════════════════

    def play_sequence(self, indices: list[int] = []):
        """Play PB pulse sequence at given inner-param value.

        Builds the correct seqArgList from config (handling legacy
        prepend vs in-args), then delegates to sequencecontrol.view_sequence.

        Parameters
        ----------
        indices : list[int]
            Which indices into the inner param_values array to plot.
            Supports negative indexing ([-1] = last point).
            Default: from config.runtime.seq_plot_indices, or [0, -1].

        Can be called before or after setup().  Does NOT mutate any
        experiment state — uses copies of all args.

        """
        printt("Playing sequence...")
        cfg = self.cfg
        if indices is []:
            indices = getattr(cfg.runtime, 'seq_plot_indices', [0, -1])

        # Build the seqArgList in the same layout PB_program expects
        args_names = cfg.seq.args_names
        # assert args_names and cfg.seq.args_values
        seq_args_base = list(cfg.seq.args_values)
        if self._is_freq_sweep:
            # PB program is fixed for freq sweeps — just plot once
            full_args = seq_args_base + [self.pb_channels]
            the_lists = seqctrl.view_sequence(
                instr=self.instr, sequence=self.sequence, seq_args=full_args,
                param_values=self.param_values, indices=indices)
        else:        
            # Determine inner param position
            if self.scan_name in args_names:
                # In-args: param lives at its args_names index
                pidx = args_names.index(self.scan_name)
                full_args = seq_args_base + [self.pb_channels]
            else:
                # Legacy prepend: param goes at index 0
                pidx = 0
                full_args = [self.param_values[0]] + seq_args_base + [self.pb_channels]

            the_lists = seqctrl.view_sequence(
                instr=self.instr, sequence=self.sequence, seq_args=full_args,
                param_values=self.param_values, indices=indices, param_index=pidx)
        
        self.pb.run_sequence_for_diode(
                [the_lists[-1][i][0] for i in range(len(the_lists[-1]))])
        
        return (self.sequence, self.scan_name, self.param_values[indices[-1]])

    # ══════════════════════════════════════════════════════════════════
    #  EXECUTE — dispatch to appropriate mode
    # ══════════════════════════════════════════════════════════════════

    def execute(
        self,
        callback: Optional[Callable] = None,
        stop_check: Optional[Callable] = None,
        save_path: Path|None = None,
        folder_number: str|int|None = None,
    ) -> np.ndarray:
        """
        Acquisition loop: outer → Nruns → inner (optionally shuffled).

        Dispatches to _execute_normal, _execute_fixed_pb, or _execute_batch
        based on RuntimeFlags.

        `ExperimentThread._callback(inner_idx, param_value, i_run, oc_idx, oc_vals, processed, raw)` 
        """
        if self._load_all_params:
            return self._execute_batch(callback, stop_check, save_path, folder_number)
        elif not self._reload_pb:
            return self._execute_fixed_pb(callback, stop_check, save_path, folder_number)
        else:
            return self._execute_normal(callback, stop_check, save_path, folder_number)

    # ── Normal mode ────────────────────────────────────────────────────

    def _execute_normal(self, callback: Callable|None, stop_check: Callable|None,
                        save_path: Path|None, folder_number: str|int|None):
        """Per-point PB reprogram (via inner setter) + DAQ read.
        Real-time plot per point."""
        return self._execute_per_point(
            callback, stop_check, save_path, folder_number)

    # ── Fixed-PB mode ──────────────────────────────────────────────────

    def _execute_fixed_pb(self, callback: Callable|None, stop_check: Callable|None,
                          save_path: Path|None, folder_number: str|int|None):
        """PB programmed once per outer combo (outer setter handles it).
        Only instrument setter (e.g. set_freq) per inner point.
        DAQ read per point.  Real-time plot per point."""
        return self._execute_per_point(
            callback, stop_check, save_path, folder_number)

    def _execute_per_point(self, callback: Callable|None, stop_check: Callable|None,
                           save_path: Path|None, folder_number: str|int|None):
        """Shared implementation for Normal and Fixed-PB modes.

        Both iterate inner points one-by-one with per-point DAQ reads.
        The only difference is what the inner setter does (PB reprogram
        vs instrument-only), which is already baked into the Sweep object.

        IMPORTANT: We do NOT use sweep.inner_iter() here because it calls
        the setter before yielding — but we need to check _valid_mask
        BEFORE the setter fires (invalid cells would cause PB_program to
        crash with bad timing values).  Instead, we iterate inner_values
        directly and call the setter ourselves only for valid cells.
        """
        # assert self.data_array and self.sweep
        self.scan_times = []
        t_start = time.perf_counter()
        t_total_s = self._t_total_s
        # timeout = t_total_s * self.Nsamples_cfg + 5
        timeout = 60
        stopped = False

        inner_vals = self.sweep.inner_values
        inner_setter = self.sweep._inner[2] if self.sweep._inner else None

        try:
            for oc_idx, oc_vals in self.sweep.outer_combos():
                if stopped: break
                self.i_outer = oc_idx
                if oc_vals:
                    lbl = "  ".join(f"{k}={v:.6g}"
                                    for k, v in oc_vals.items())
                    printt(f"▶ Outer [{oc_idx+1}/"
                          f"{self.sweep.outer_count}]: {lbl}")

                for self.i_run in range(self.Nruns):
                    if stop_check and stop_check():
                        stopped = True; break
                    printt(f"  Run {self.i_run+1}/{self.Nruns}")

                    for order_pos in range(self.sweep.inner_count):
                        if stop_check and stop_check():
                            stopped = True; break

                        actual_i = int(self._inner_order[order_pos])
                        actual_val = inner_vals[actual_i]

                        # Skip invalid cells → NaN (BEFORE setter fires)
                        if (self._valid_mask is not None and not self._valid_mask[oc_idx, actual_i]):
                            self.data_array[oc_idx, self.i_run, actual_i, :] = np.nan
                            continue

                        # Shuffled: manually call inner setter
                        if self._shuffle and self.sweep._inner[2] is not None:
                            self.sweep._inner[2](actual_val)

                        # Call inner setter only for valid cells
                        if inner_setter is not None:
                            inner_setter(actual_val)

                        self.i_scanpt = actual_i

                        # if order_pos % max(1, self.sweep.inner_count // 10) == 0:
                        print(f"    \x1b[38;2;250;250;0m{order_pos+1}/ {self.sweep.inner_count}: "
                              f"{actual_val:.6g}\x1b[0m")

                        t0 = time.perf_counter()
                        # assert self.input_task is not None
                        self.input_task._task.stop()
                        self.input_task._task.start()
                        self.input_task._task_state = "running"

                        cts = self.input_task.read_daq(self.daq_Nsamples, timeout=timeout)
                        self.data_array[oc_idx, self.i_run, actual_i, :] = np.ravel(np.array(cts))
                        self.scan_times.append(time.perf_counter() - t0)

                        # Callback: to process and plot per point
                        if callback is not None:
                            self._fire_callback(callback, actual_i, actual_val,
                                                self.i_run, oc_idx, oc_vals)

                    self._autosave_run(save_path, folder_number)

        except KeyboardInterrupt:
            printt("🛑 Interrupted")

        self._finish(t_start, save_path, folder_number)
        return self.data_array

    # ── Batch mode ─────────────────────────────────────────────────────

    def _execute_batch(self, callback: Callable|None, stop_check: Callable|None,
                       save_path: Path|None, folder_number:str|int|None):
        """All inner scan points compiled into a single PB program.
        Single DAQ read captures the entire inner sweep.
        Real-time plot updates once per run (after the read).
        Fastest mode — eliminates per-point PB + DAQ overhead.

        The sequence function must accept a list/array of parameter
        values in the appropriate seqArgList slot.  This mirrors the
        existing acquire_data_all_params() pattern.
        """
        self.scan_times = []
        t_start = time.perf_counter()
        stopped = False

        cfg = self.cfg
        args_names = cfg.seq.args_names

        # Inner param's position in seq_args
        # assert args_names and self.sweep
        if self.scan_name in args_names:
            inner_pidx = args_names.index(self.scan_name)
        else:
            inner_pidx = -1  # legacy prepend

        # Batch DAQ: all inner points × Nsamples at once
        n_inner = self.sweep.inner_count
        per_pt = self.daq_Nsamples * self.n_channels
        batch_Nsamples = self.daq_Nsamples * n_inner
        batch_timeout = 60 * 10

        # Re-create AI/CI task for the larger batch read
        # assert self.input_task is not None
        self.input_task._task.stop()
        self.input_task._task.close()
        if cfg.runtime.detector == 'counter':
            ci = cfg.daq_ci
            # assert ci
            # assert ci.start_trigger_source and ci.pause_trigger_source
            self.input_task = CounterInputTask(dev='P6363',
                                               channel=ci.channel,
                                               counter=ci.counter,
                                               sample_source=ci.sample_source,
                                               sample_rate=ci.sample_rate,
                                               sample_mode=ci.sample_mode,
                                               samps_per_chan=int(batch_Nsamples),
                                               start_trigger_source=ci.start_trigger_source,
                                               pause_trigger_source=ci.pause_trigger_source,
                                               )
            self.n_channels = 1  # counter is always single-channel
            printt(f"▶ Detector: COUNTER ({ci.counter}, input={ci.channel})")
        else:
            ai = cfg.daq_ai
            # assert ai
            # assert ai.voltage_ranges and ai.sample_source and ai.pause_trigger_source and ai.start_trigger_source
            self.input_task = AnalogInputTask(channels=ai.channels,
                                              voltage_ranges=ai.voltage_ranges,
                                              sample_source=ai.sample_source,
                                              sample_rate=ai.sample_rate,
                                              sample_mode=ai.sample_mode,
                                              samps_per_chan=int(batch_Nsamples),
                                              start_trigger_source=ai.start_trigger_source,
                                              pause_trigger_source=ai.pause_trigger_source,
                                              )
            printt(f"▶ Detector: ANALOG (ai{ai.channels})")

        # assert self.data_array
        try:
            for oc_idx, oc_vals in self.sweep.outer_combos():
                if stopped: break
                self.i_outer = oc_idx
                if oc_vals:
                    lbl = "  ".join(f"{k}={v:.6g}"
                                    for k, v in oc_vals.items())
                    printt(f"▶ Outer [{oc_idx+1}/{self.sweep.outer_count}]: {lbl}")

                # Build batch PB program with all inner values
                batch_args = list(self.seq_args)
                if inner_pidx == -1:
                    batch_args = [list(self.param_values)] + batch_args
                else:
                    batch_args[inner_pidx] = list(self.param_values)

                _, the_list = PulseBlaster.PB_program(
                    self.instr, self.sequence,
                    batch_args + [self.pb_channels])
                
                # instructionList = [element[0] for element in the_list]

                self.pb.run_sequence_for_diode(
                    [the_list[i][0] for i in range(len(the_list))])

                for self.i_run in range(self.Nruns):
                    if stop_check and stop_check():
                        stopped = True; break
                    printt(f"  Run {self.i_run+1}/{self.Nruns} (batch)")

                    t0 = time.perf_counter()
                    self.input_task._task.stop()
                    self.input_task._task.start()
                    self.input_task._task_state = "running"

                    cts = self.input_task.read_daq(batch_Nsamples, timeout=batch_timeout)

                    dt = time.perf_counter() - t0
                    self.scan_times.append(dt)

                    # Reshape into (inner_pts, daq_per_pt)
                    raw_flat = np.ravel(np.array(cts))
                    n_got = min(len(raw_flat) // per_pt, n_inner)

                    if n_got == n_inner:
                        self.data_array[oc_idx, self.i_run] = raw_flat[:n_inner * per_pt].reshape(n_inner, per_pt)
                    elif n_got > 0:
                        printt(f"  ⚠ Batch read: got {n_got}/{n_inner} pts")
                        self.data_array[
                            oc_idx, self.i_run, :n_got] = raw_flat[:n_got * per_pt].reshape(n_got, per_pt)

                    # NaN-fill invalid cells
                    if self._valid_mask is not None:
                        bad = ~self._valid_mask[oc_idx]
                        self.data_array[oc_idx, self.i_run, bad, :] = np.nan

                    # Callback: once per run - to process and plot
                    if callback is not None:
                        self._fire_batch_callback(callback, oc_idx, oc_vals, raw_flat, per_pt)

                    printt(f"    Batch: {dt:.2f}s ({dt/max(n_got,1)*1e3:.1f} ms/pt effective)")

                    self._autosave_run(save_path, folder_number)

        except KeyboardInterrupt:
            printt("🛑 Interrupted")

        self._finish(t_start, save_path, folder_number)
        return self.data_array

    # ── shared helpers ─────────────────────────────────────────────────

    def _fire_callback(self, callback: Callable, actual_i, actual_val, i_run, oc_idx, oc_vals):
        """Process current data slice and fire per-point callback to plot data and update status bar.
        
        Parameters
        ----------
        callback : Callable
            callback after each parameter: `ExperimentThread._callback()`
        actual_i : int
            Inner loop / parameter index
        actual_val : float
            Inner loop / parameter value
        i_run : int
            Run index
        oc_idx : int
            Outer sweep indices
        oc_vals : float
            Outer sweep values
        """
        # assert self.data_array
        filled = self.data_array[oc_idx, i_run, :, :]
        valid_rows = ~np.isnan(filled).any(axis=1) & filled.any(axis=1)
        n_filled = int(np.count_nonzero(valid_rows))
        proc = process_data(n_filled,
                            filled[valid_rows],
                            self.Nsamples_cfg, self.sequence,
                            ) if n_filled > 0 else None
        raw_row = self.data_array[oc_idx, i_run, actual_i, :]
        callback(actual_i, actual_val, i_run, oc_idx, oc_vals, proc, raw_row)

    def _fire_batch_callback(self, callback: Callable, oc_idx, oc_vals, raw_flat, per_pt):
        """Process full-run data and fire callback once for batch mode to plot data and update status bar.

        Parameters
        ----------
        callback : Callable
            callback after each run: `ExperimentThread._callback()`
        oc_idx : int
            Outer sweep indices
        oc_vals : float
            Outer sweep values
        raw_flat : Array-like
            Flattened raw data
        per_pt : int

        """
        # assert self.data_array and self.sweep
        filled = self.data_array[oc_idx, self.i_run]
        valid_rows = ~np.isnan(filled).any(axis=1) & filled.any(axis=1)
        n_filled = int(np.count_nonzero(valid_rows))
        proc = process_data(n_filled,
                            filled[valid_rows],
                            self.Nsamples_cfg, self.sequence,
                            ) if n_filled > 0 else None
        
        n_inner = self.sweep.inner_count
        callback(n_inner - 1, self.sweep.inner_values[-1], self.i_run, oc_idx, oc_vals, proc,
                 raw_flat[-per_pt:] if len(raw_flat) >= per_pt else None)

    def _autosave_run(self, save_path: Path|None, folder_number: str|int|None):
        if save_path and folder_number:
            df = save_path / f"data_{folder_number}.npy"
            savefile(df, self.data_array)
            printt(f"    💾 Run {self.i_run+1} → {df.name}")

    def _finish(self, t_start: float, save_path: Path|None, folder_number:str|int|None):
        self.elapsed = time.perf_counter() - t_start
        if self.scan_times:
            printt(f"✔ {self.elapsed:.1f}s  "
                  f"({np.mean(self.scan_times)*1e3:.1f} ± "
                  f"{np.std(self.scan_times)*1e3:.1f} ms/pt)")
        if self._valid_mask is not None:
            n_nan = int(np.sum(~self._valid_mask))
            if n_nan > 0:
                printt(f"  ℹ {n_nan} cells NaN-filled (invalid PB timing)")
        if save_path and folder_number:
            self._save_params(save_path, folder_number)

    def _save_params(self, save_path: Path, folder_number: int|str):
        # TODO: Confirm: update params from DAQ class and experiment threads
        """Save config + sweep metadata + validation + mode info."""
        # assert self.data_array and self.sweep
        cfg_dict = self.cfg.to_dict()
        cfg_dict['experiment_type'] = 'sweep'
        cfg_dict['_sweep_info'] = {
            'inner': self.sweep.inner_name,
            'inner_count': self.sweep.inner_count,
            'outer_names': self.sweep.outer_names,
            'outer_shape': list(self.sweep.outer_shape),
            'data_shape': list(self.data_array.shape),
            'shuffled': self._shuffle,
        }
        cfg_dict['_timing'] = {
            'elapsed_s': self.elapsed,
            'mean_per_pt_ms': (float(np.mean(self.scan_times) * 1e3)
                               if self.scan_times else 0),
            'std_per_pt_ms': (float(np.std(self.scan_times) * 1e3)
                              if self.scan_times else 0),
        }
        cfg_dict['_execution_mode'] = {
            'reload_pb': self._reload_pb,
            'load_all_params': self._load_all_params,
            'mode': ('batch' if self._load_all_params
                     else 'fixed_pb' if not self._reload_pb
                     else 'normal'),
        }
        if self._validation is not None:
            cfg_dict['_validation'] = {
                'n_total': self._validation.n_total,
                'n_invalid': self._validation.n_invalid,
                'hw_errors': [(n, v, r) for n, v, r in self._validation.hw_errors],
            }
        pf = save_path / f"params_{folder_number}.yaml"
        savefile(pf, cfg_dict)

        vf = save_path / f"scan_values_{folder_number}.npz"
        arrays = {name: axis.values
                  for name, axis in self.cfg.scans.items()}
        if self._valid_mask is not None:
            arrays['_valid_mask'] = self._valid_mask.astype(np.uint8)
        np.savez_compressed(vf, **arrays)
        printt(f"  💾 Params → {pf.name}, values → {vf.name}")

    # ── teardown ────────────────────────────────────────────────────────

    def teardown(self):
        if self.input_task is not None:
            try:
                self.input_task.stop()
                self.input_task._task.close()
            except Exception:
                pass
            self.input_task = None

    # ── config_sweep ────────────────────────────────────────────────────

    def config_sweep(self, config_name: str, config_values,
                     config_applier: Callable,
                     **execute_kwargs) -> list[np.ndarray]:
        """
        Outer loop over a CONFIGURATION parameter that requires
        AI task rebuild (e.g. t_AOM, daq_Nsamples, sample_rate).
        """
        results = []
        for i, val in enumerate(config_values):
            printt(f"\n{'='*60}")
            printt(f"Config sweep [{i+1}/{len(config_values)}]: "
                  f"{config_name} = {val}")
            printt(f"{'='*60}")

            config_applier(val)

            self.teardown()
            # assert self.cfg.seq.Nsamples and self.cfg.seq.args_values
            self.Nsamples_cfg = self.cfg.seq.Nsamples
            self.seq_args = list(self.cfg.seq.args_values)
            self.param_values = self.cfg.primary_scan.values
            self.Nscanpts = len(self.param_values)
            self._validation = None  # force re-validation
            self.setup()
            data = self.execute(**execute_kwargs)
            results.append(data)

        return results
