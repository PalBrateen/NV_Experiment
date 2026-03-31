"""
experiment_config.py — Base config framework for NV-center experiments.

Design goals:
  - Dataclasses for type safety, IDE autocomplete, and trivial serialization
  - ScanAxis owns its own value generation (linear / log / nonlinear / explicit)
  - Lazy: scan values computed on first access, not at import time
  - Clean separation: sweep defs | sequence timing | hardware | runtime flags
  - Fully serializable to dict/JSON/YAML via .to_dict() / .save_json() / .save_yaml()
  - diff_config() for comparing experiment configs (single pair or batch)

Usage:
    from esr_config_new import config

    config.scans['freq'].values        # lazily generates if needed
    config.mw.power               # attribute access (IDE autocomplete)
    config.pb.channels                 # dict
    config.to_dict()                   # full serializable snapshot
    config.save_json('esr_run1.json')  # save to disk
    config.save_yaml('esr_run1.yaml')  # save to disk (needs pyyaml)

    # Compare configs:
    from experiment_config import diff_config
    diffs = diff_config(config_old, config_new)
    for path, (v1, v2) in diffs.items():
        print(f"  {path}: {v1} -> {v2}")
"""

from __future__ import annotations
from dataclasses import dataclass, field, fields
from typing import Optional, Union, Any
from enum import Enum
import numpy as np
from pathlib import Path
import json

# ── Unit imports (unchanged from your codebase) ────────────────────────
try:
    from spinapi import ns, us, ms
    from SGcontrol import Hz, kHz, MHz, GHz
    from connectionConfig import (PBclk, laser, samp_clk, start_trig,
                                  MW, camera, lia1)
except ImportError:
    ns, us, ms = 1, 1e3, 1e6
    Hz, kHz, MHz, GHz = 1, 1e3, 1e6, 1e9
    PBclk = 500
    laser = samp_clk = start_trig = MW = camera = lia1 = 0


# ═══════════════════════════════════════════════════════════════════════
# 1. SCAN AXIS — one per swept parameter
# ═══════════════════════════════════════════════════════════════════════

class ScanMode(Enum):
    """How sweep values are generated."""
    LINEAR = 'linear'
    LOG = 'log'
    SEGMENTS = 'segments'     # piecewise-linear
    EXPLICIT = 'explicit'     # user-provided list


@dataclass
class NonlinearSegment:
    """One segment of a piecewise-linear sweep.

    Args:
        start:  segment start (in the segment's natural units, e.g. GHz)
        stop:   segment stop
        step:   step size in *milli*-units of the segment spacing
                (matching your existing f-dict convention: step=10 in a
                GHz-defined segment → 10 MHz spacing)
    """
    start: float
    stop: float
    step: float


@dataclass
class ScanAxis:
    """A single swept parameter (frequency, tau, mw_power, etc.).

    Supports four generation modes controlled by what you provide:

    1. LINEAR (default):  start + stop + (step OR Nscanpts)
    2. LOG:               start + stop + Nscanpts + mode=ScanMode.LOG
    3. SEGMENTS:          segments=[NonlinearSegment(...), ...] + unit_scale
    4. EXPLICIT:          explicit_values=[v1, v2, v3, ...]

    All modes support shuffle=True for randomized order.

    Examples
    --------
    # Linear: 2.82 GHz to 2.92 GHz in 1 MHz steps
    ScanAxis('freq', start=2.82*GHz, stop=2.92*GHz, step=1*MHz)

    # Logarithmic: 10 ns to 1 ms, 50 points
    ScanAxis('tau', start=10*ns, stop=1*ms, Nscanpts=50, mode=ScanMode.LOG)

    # Explicit list of arbitrary values
    ScanAxis('freq', explicit_values=[2.85*GHz, 2.87*GHz, 2.89*GHz])

    # Piecewise (nonlinear) with unit scaling
    ScanAxis('freq', segments=[...], unit_scale=GHz)
    """
    name: str
    start: float = 0.0
    stop: float = 0.0
    step: Optional[float] = None
    Nscanpts: Optional[int] = None
    mode: ScanMode = ScanMode.LINEAR
    segments: Optional[list[NonlinearSegment]] = None
    unit_scale: float = 1.0
    explicit_values: Optional[list[float]] = None
    shuffle: bool = False

    # ── internal cache ──
    _values: Optional[np.ndarray] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        # Auto-detect mode from what was provided
        if self.explicit_values is not None:
            self.mode = ScanMode.EXPLICIT
            self.Nscanpts = len(self.explicit_values)
            return
        if self.segments is not None:
            self.mode = ScanMode.SEGMENTS
            return

        # For linear/log: resolve step ↔ Nscanpts
        if self.mode in (ScanMode.LINEAR, ScanMode.LOG):
            if self.step is not None and self.Nscanpts is None:
                if self.step != 0:
                    self.Nscanpts = round((self.stop - self.start) / self.step) + 1
                else:
                    self.Nscanpts = 1
            elif self.Nscanpts is not None and self.step is None:
                if self.Nscanpts > 1:
                    self.step = (self.stop - self.start) / (self.Nscanpts - 1)
                else:
                    self.step = 0

    @property
    def values(self) -> np.ndarray:
        """Lazily generate and cache the sweep array."""
        if self._values is None:
            self._values = self._generate()
        return self._values

    @values.setter
    def values(self, arr):
        """Allow direct assignment: scan.values = my_array."""
        self._values = np.asarray(arr)
        self.Nscanpts = len(self._values)
        self.mode = ScanMode.EXPLICIT

    def invalidate(self):
        """Clear cached values and re-resolve Nscanpts/step.
        Call after modifying start/stop/step/Nscanpts."""
        self._values = None
        # Re-run the step ↔ Nscanpts resolution
        if self.mode in (ScanMode.LINEAR, ScanMode.LOG):
            if self.step is not None and self.Nscanpts is None:
                if self.step != 0:
                    self.Nscanpts = round((self.stop - self.start) / self.step) + 1
                else:
                    self.Nscanpts = 1
            elif self.Nscanpts is not None and self.step is None:
                if self.Nscanpts > 1:
                    self.step = (self.stop - self.start) / (self.Nscanpts - 1)
                else:
                    self.step = 0

    def _generate(self) -> np.ndarray:
        if self.mode == ScanMode.EXPLICIT:
            vals = np.asarray(self.explicit_values, dtype=float)
        elif self.mode == ScanMode.SEGMENTS:
            vals = self._generate_segments()
        elif self.mode == ScanMode.LOG:
            vals = np.geomspace(self.start, self.stop, self.Nscanpts, endpoint=True)
        else:  # LINEAR
            vals = np.linspace(self.start, self.stop, self.Nscanpts, endpoint=True)

        if self.shuffle:
            vals = vals.copy()
            np.random.shuffle(vals)
        return vals

    def _generate_segments(self) -> np.ndarray:
        """Piecewise-linear sweep from segment list."""
        pieces = []
        for i, seg in enumerate(self.segments):
            n = round((seg.stop - seg.start) / (seg.step / 1e3))
            endpoint = (i == len(self.segments) - 1)
            arr = np.linspace(seg.start, seg.stop,
                              n + (1 if endpoint else 0),
                              endpoint=endpoint)
            pieces.append(arr)
        raw = np.concatenate(pieces)
        return raw * self.unit_scale

    def to_dict(self, include_values: bool = False) -> dict:
        """Serializable recipe (no numpy, no cache).

        Args:
            include_values: If True, append the full values list.
                Default False — values go in the companion .npz instead.
        """
        d = {
            'name': self.name,
            'mode': self.mode.value,
            'start': self.start,
            'stop': self.stop,
            'step': self.step,
            'Nscanpts': self.Nscanpts if self.Nscanpts else len(self.values),
            'shuffle': self.shuffle,
        }
        if self.mode == ScanMode.SEGMENTS and self.segments is not None:
            d['segments'] = [
                {'start': s.start, 'stop': s.stop, 'step': s.step}
                for s in self.segments
            ]
            d['unit_scale'] = self.unit_scale
        if self.mode == ScanMode.EXPLICIT and self.explicit_values is not None:
            # For explicit mode the recipe IS the values — no way to regenerate
            d['explicit_values'] = list(self.explicit_values)
        if include_values:
            d['values'] = self.values.tolist()
        return d


# ═══════════════════════════════════════════════════════════════════════
# 2. INSTRUMENT / HARDWARE CONFIG BLOCKS
#    Each is a plain dataclass. To add a new instrument or parameter,
#    just add a field — existing configs that don't set it get the default.
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class TimeConfig:
    """Time variables to Pulse Blaster"""
    t_AOM: float = 10 *ms
    AOM_lag: float = 700 *ns
    MW_lag: float = 160 *ns

    t_pi: float = 150 *ns
    # t_piby2: float = t_pi/2
    ro_delay: float = 400 *ns
    t_tot: float = 0.

@dataclass
class MicrowaveConfig:
    """SG384 signal generator settings."""
    power: float = -24.0            # dBm
    freq: float = 2.87e9            # Hz
    mod_freq: float = 10e3          # Hz
    # mod_type: str = 

# @dataclass
# class LaserConfig:
#     t_AOM: float = 0.0     # ns
#     power: float = 0.0     # placeholder (Verdi watts, or AOM drive V)

@dataclass
class PulseBlasterConfig:
    """PulseBlaster channel map + clock."""
    clock_MHz: float = PBclk
    channels: dict[str, int] = field(default_factory=dict)

    @property
    def clk_cyc_ns(self) -> float:
        return 1e3 / self.clock_MHz

@dataclass
class SequenceConfig:
    """Pulse sequence identity + timing args."""
    name: str = ''
    args_names: list[str] = field(default_factory=list)
    args_values: list[float] = field(default_factory=list)  # ns
    Nsamples: int = 0

@dataclass
class DAQAIConfig:
    """NI-DAQ acquisition — extend fields when ready."""
    ai_channels: list[str] = field(default_factory=list)
    ai_voltage_ranges: list[tuple[float,float]] = field(default_factory=list)

    ai_sample_source: str = ''          # '' means internal
    ai_sample_rate: float = 10e3        # Sa/s
    ai_sample_mode: str = 'finite'      # finite sampling mode
    ai_samps_per_chan: int = 1
    
    ai_start_trigger_type: str = 'digital'
    ai_start_trigger_source: str = ''
    ai_start_trigger_edge: str = 'rising'

    ai_pause_trigger_type: str = 'digital'
    ai_pause_trigger_source: str = ''
    ai_pause_trigger_edge: str = ''


@dataclass
class DAQAOConfig:
    """NI-DAQ acquisition — extend fields when ready."""

    ao_channels: list[str] = field(default_factory=list)
    ao_voltage_ranges: list[tuple[float,float]] = field(default_factory=list)

    ao_sample_source: str = ''          # '' means internal
    ao_sample_rate: float = 10e3        # Sa/s
    ao_sample_mode: str = 'finite'      # finite sampling mode
    ao_samps_per_chan: int = 1
    
    ao_start_trigger_type: str = 'digital'
    ao_start_trigger_source: str = ''
    ao_start_trigger_edge: str = 'rising'


@dataclass
class DAQCIConfig:
    """NI-DAQ counter input configuration.

    Counters read photon counts from a single-photon detector (SPD).
    The counter accumulates edges; read_daq() returns np.diff (counts per bin).

    Fields
    ------
    ci_counter : str
        Physical counter channel, e.g. 'P6363/ctr0'.
    ci_input_terminal : str
        PFI terminal where SPD pulses arrive, e.g. '/P6363/PFI2'.
    ci_sample_source : str
        External sample clock source (PB), e.g. 'PFI14'.  '' = internal.
    ci_sample_rate : float
        Nominal rate passed to cfg_samp_clk_timing (Hz).
    ci_sample_mode : str
        'finite' for sweep, 'continuous' for timeseries.
        (Typically set by the experiment class, not the user.)
    ci_samps_per_chan : int
        Buffer / read size per acquisition.
    ci_start_trigger_source : str
        Arm-start trigger from PB, e.g. 'PFI15'.  '' = no trigger.
    ci_start_trigger_edge : str
        Edge type for the arm-start trigger.
    """
    ci_counter: str = 'P6363/ctr0'
    ci_input_terminal: str = '/P6363/PFI2'

    ci_sample_source: str = ''          # '' = internal clock
    ci_sample_rate: float = 10e3        # Sa/s
    ci_sample_mode: str = 'finite'
    ci_samps_per_chan: int = 1

    ci_start_trigger_source: str = ''
    ci_start_trigger_edge: str = 'rising'


@dataclass
class UHFLIConfig:
    """Zurich UHFLI lock-in settings — for experiments that use it.

    This is a *parameter snapshot* stored alongside the experiment config,
    complementary to the profile-based save/load in uhfli_control.py.
    Populate it manually here, or attach a raw profile dict from
    UHFLIConfigManager.get_config().
    """
    demod_index: int = 0
    oscillator_freq: float = 0.0
    timeconstant: float = 1e-3
    filter_order: int = 4
    input_range: float = 1.0
    output_amplitude: float = 0.0

    # Attach a raw profile dict from UHFLIConfigManager for full fidelity
    _raw_profile: Optional[dict] = field(default=None, repr=False)

@dataclass
class PlotConfig:
    x_units: float = GHz
    x_label: str = 'Frequency (GHz)'

@dataclass
class SaveConfig:
    prefix: str = 'Experiment'


# ═══════════════════════════════════════════════════════════════════════
# 3. RUNTIME FLAGS (separated from scan definitions)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class RuntimeFlags:
    Nruns: int = 1
    reload_pb: bool = True
    load_all_params: bool = False
    seq_plot_indices: list[int] = field(default_factory=lambda: [0, -1])

    # Detector type: 'analog' (photodiode via AI) or 'counter' (SPD via CI)
    detector: str = 'analog'

    # Timeseries
    ts_display_seconds: float = 30.0
    ts_max_duration: float = 0.0      # 0 = manual stop only
    ts_buffer_size: int = 0           # 0 = auto
    ts_contrast_op: str = 's/r'

# ═══════════════════════════════════════════════════════════════════════
# 4. EXPERIMENT CONFIG — the top-level container
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class ExperimentConfig:
    """Top-level experiment configuration.

    Attributes
    ----------
    scans : dict[str, ScanAxis]
        Keyed by parameter name. First key = primary sweep axis.
    mw, times, pb, seq, daq_ai, daq_ao, uhfli, plot, save, runtime :
        Typed sub-configs for each instrument / concern.
    extra : dict
        Escape hatch for one-off parameters that don't fit elsewhere.
        Saved as-is to JSON/YAML (must be JSON-serializable).

    Access patterns
    ---------------
        config.scans['freq'].values        # scan array
        config.mw.power               # attribute access
        config.pb.channels['laser']        # dict inside dataclass
        config.to_dict()                   # full serializable snapshot
        config.save_json('run1.json')      # save to disk
        config.save_yaml('run1.yaml')      # save to disk (needs pyyaml)

    Adding a new instrument
    -----------------------
    1. Define a new dataclass (e.g. MagnetConfig) with default values.
    2. Add it as a field here with default_factory.
    3. Add it to to_dict() serialization.
    That's it — existing config files that don't set it just get defaults.
    """
    scans: dict[str, ScanAxis] = field(default_factory=dict)
    mw: MicrowaveConfig = field(default_factory=MicrowaveConfig)
    # laser: LaserConfig = field(default_factory=LaserConfig)
    times: TimeConfig = field(default_factory=TimeConfig)
    pb: PulseBlasterConfig = field(default_factory=PulseBlasterConfig)
    seq: SequenceConfig = field(default_factory=SequenceConfig)
    daq_ai: DAQAIConfig = field(default_factory=DAQAIConfig)
    daq_ao: DAQAOConfig = field(default_factory=DAQAOConfig)
    daq_ci: DAQCIConfig = field(default_factory=DAQCIConfig)
    uhfli: UHFLIConfig = field(default_factory=UHFLIConfig)
    plot: PlotConfig = field(default_factory=PlotConfig)
    save_opts: SaveConfig = field(default_factory=SaveConfig)
    runtime: RuntimeFlags = field(default_factory=RuntimeFlags)
    extra: dict = field(default_factory=dict)

    # Which optional instrument blocks were explicitly provided.
    # Populated automatically by __init_subclass__ / __post_init__.
    _active_instruments: set = field(default_factory=set, repr=False)

    # Optional instrument fields — only serialized/printed when active.
    _OPTIONAL_FIELDS = {'mw', 'times', 'daq_ai', 'daq_ao', 'daq_ci', 'uhfli'}

    def __post_init__(self):
        """Detect which optional instruments were explicitly passed.

        Works by checking if __init__ received a non-default-factory value.
        We do this by comparing object identity: if the user passed their
        own instance, it won't be the same object as a fresh default.
        
        Since dataclass __init__ always calls default_factory for unset 
        fields, we instead rely on a class-method constructor pattern:
        the _active_instruments set is filled by __init__ automatically 
        because we intercept via __init_subclass__ ... 
        
        Actually, the simplest reliable approach: compare each optional 
        field against a freshly constructed default. If ANY public field 
        differs, the instrument is active. If the user intentionally set 
        values that happen to match defaults (e.g., power=-24 which 
        IS the default), we still mark it active by having the user 
        explicitly include it in _active_instruments if desired.
        
        In practice: we auto-detect by value comparison, AND provide
        .use() for explicit opt-in when defaults happen to match.
        """
        for fname in self._OPTIONAL_FIELDS:
            obj = getattr(self, fname)
            default_cls = type(obj)
            default_obj = default_cls()
            for f in fields(obj):
                if f.name.startswith('_'):
                    continue
                if getattr(obj, f.name) != getattr(default_obj, f.name):
                    self._active_instruments.add(fname)
                    break

    def use(self, *instrument_names: str) -> 'ExperimentConfig':
        """Explicitly mark instruments as active (even if they match defaults).
        
        Use this when you intentionally set an instrument to its default 
        values and still want it serialized/printed.
        
        Example:
            config = ExperimentConfig(
                mw=MicrowaveConfig(power=-24),  # happens to match default
                ...
            ).use('mw', 'daq_ai')
        
        Returns self for chaining.
        """
        for name in instrument_names:
            if name not in self._OPTIONAL_FIELDS:
                raise ValueError(f"'{name}' is not an optional instrument. "
                                 f"Valid: {self._OPTIONAL_FIELDS}")
            self._active_instruments.add(name)
        return self

    def _is_active(self, field_name: str) -> bool:
        """Check if an instrument block should be serialized/printed."""
        if field_name not in self._OPTIONAL_FIELDS:
            return True  # core fields always active
        return field_name in self._active_instruments

    # ── Convenience properties ──

    @property
    def primary_scan(self) -> ScanAxis:
        """First scan axis = the primary swept parameter."""
        return next(iter(self.scans.values()))

    @property
    def scan_names(self) -> list[str]:
        return list(self.scans.keys())

    # ── Serialization ──

    # Instrument fields that are only serialized/printed when active.
    # Detected automatically in __post_init__, or forced via .use().

    def to_dict(self, include_values: bool = False, include_defaults: bool = False) -> dict:
        """Full serializable snapshot for JSON / HDF5 / YAML.

        Args:
            include_values: If True, embed scan value arrays in the dict.
                Default False — values go in the companion .npz via save().
            include_defaults: If True, include all instrument blocks even
                if they weren't configured. Default False — only active
                instruments appear in the output.
        """
        def _safe_asdict(obj) -> dict:
            """asdict but skip private/underscore fields."""
            return {
                f.name: getattr(obj, f.name)
                for f in fields(obj)
                if not f.name.startswith('_')
            }

        def _include(name: str) -> bool:
            if include_defaults:
                return True
            return self._is_active(name)

        d = {
            'experiment': self.seq.name,
            'scans': {k: v.to_dict(include_values=include_values)
                      for k, v in self.scans.items()},
        }
        d['times'] = _safe_asdict(self.times)
        if _include('mw'):
            d['mw'] = _safe_asdict(self.mw)
        # if _include('laser'):
        #     d['laser'] = _safe_asdict(self.laser)
        # if _include('pb'):
        d['pb'] = {
            'clock_MHz': self.pb.clock_MHz,
            'clk_cyc_ns': self.pb.clk_cyc_ns,
            'channels': self.pb.channels,
        }

        # seq, runtime, plot, save are always included (core experiment identity)
        d['seq'] = {
            'name': self.seq.name,
            'args_names': self.seq.args_names,
            'args_values': [float(v) for v in self.seq.args_values],
            'Nsamples': self.seq.Nsamples,
        }

        if _include('daq_ai'):
            d['daq_ai'] = _safe_asdict(self.daq_ai)

        if _include('daq_ao'):
            d['daq_ao'] = _safe_asdict(self.daq_ao)

        if _include('daq_ci'):
            d['daq_ci'] = _safe_asdict(self.daq_ci)

        if _include('uhfli'):
            d['uhfli'] = _safe_asdict(self.uhfli)

        d['plot'] = {'x_units': self.plot.x_units, 'x_label': self.plot.x_label}
        d['save'] = _safe_asdict(self.save_opts)
        d['runtime'] = _safe_asdict(self.runtime)

        if self.extra:
            d['extra'] = self.extra
        return d

    # ── Save: recipe (.json) + values (.npz) ──

    def save(self, filepath: Union[str, Path], fmt: str = 'yaml') -> tuple[Path, Path]:
        """Save config recipe + scan values as a pair of files.

        Produces:
            <filepath>.yaml  (or .json)  — human-readable recipe (small, inspectable)
            <filepath>.npz               — scan value arrays (one per axis)

        Args:
            filepath: Base path WITHOUT extension (extensions added automatically).
                      E.g. 'data/ESR_2026-03-24_run1'
            fmt: 'json' or 'yaml' for the recipe file.

        Returns:
            (recipe_path, values_path) tuple.

        Example:
            recipe, values = config.save('data/ESR_run1')
            # creates: data/ESR_run1.json + data/ESR_run1.npz

            # To reload later:
            recipe_dict, scan_values = ExperimentConfig.load('data/ESR_run1')
        """
        base = Path(filepath)

        # Recipe file (clean, small, human-readable)
        if fmt == 'json':
            recipe_path = base.with_suffix('.json')
            self.save_json(recipe_path)
        else:
            recipe_path = base.with_suffix('.yaml')
            self.save_yaml(recipe_path)

        # Values file (compact binary, one array per scan axis)
        values_path = base.with_suffix('.npz')
        arrays = {name: axis.values for name, axis in self.scans.items()}
        np.savez_compressed(values_path, **arrays)

        return recipe_path, values_path

    @staticmethod
    def load(filepath: Union[str, Path]) -> tuple[dict, dict[str, np.ndarray]]:
        """Load a saved config pair (recipe + values).

        Args:
            filepath: Base path without extension, or path to the .json/.yaml.
                      The companion .npz is found automatically.

        Returns:
            (recipe_dict, scan_values) where scan_values maps
            axis name → numpy array.

        Example:
            recipe, vals = ExperimentConfig.load('data/ESR_run1')
            print(recipe['mw']['power'])
            print(vals['freq'])  # numpy array
        """
        base = Path(filepath)
        # Find recipe file
        if base.suffix in ('.json', '.yaml', '.yml'):
            recipe_path = base
            npz_path = base.with_suffix('.npz')
        else:
            # Try json first, then yaml
            if base.with_suffix('.json').exists():
                recipe_path = base.with_suffix('.json')
            elif base.with_suffix('.yaml').exists():
                recipe_path = base.with_suffix('.yaml')
            else:
                raise FileNotFoundError(f"No .json or .yaml found at {base}")
            npz_path = base.with_suffix('.npz')

        # Load recipe
        if recipe_path.suffix == '.json':
            recipe = ExperimentConfig.load_json(recipe_path)
        else:
            recipe = ExperimentConfig.load_yaml(recipe_path)

        # Load values
        scan_values = {}
        if npz_path.exists():
            with np.load(npz_path) as data:
                for key in data.files:
                    scan_values[key] = data[key]
        else:
            import warnings
            warnings.warn(f"No .npz found at {npz_path}; scan values not loaded.")

        return recipe, scan_values

    # ── Individual format saves (still available) ──

    def save_json(self, filepath: Union[str, Path], indent: int = 2) -> Path:
        """Save recipe-only to JSON file (no scan values)."""
        filepath = Path(filepath)
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=indent, default=str)
        return filepath

    def save_yaml(self, filepath: Union[str, Path]) -> Path:
        """Save recipe-only to YAML file (no scan values). Requires pyyaml."""
        import yaml
        filepath = Path(filepath)
        with open(filepath, 'w') as f:
            yaml.dump(self.to_dict(), f,
                      indent=4, default_flow_style=False, sort_keys=False)
        return filepath

    @staticmethod
    def load_json(filepath: Union[str, Path]) -> dict:
        """Load a saved config dict from JSON. Returns raw dict."""
        with open(filepath) as f:
            return json.load(f)

    @staticmethod
    def load_yaml(filepath: Union[str, Path]) -> dict:
        """Load a saved config dict from YAML. Returns raw dict."""
        import yaml
        with open(filepath) as f:
            return yaml.safe_load(f)

    # ── Pretty-print ──

    def print_config(self, show_values: bool = False, show_all: bool = False) -> None:
        """Pretty-print the config to terminal.

        Args:
            show_values: If True, print first/last few scan values.
            show_all: If True, print all instrument blocks including
                unconfigured ones.
        """
        def _show(name: str) -> bool:
            if show_all:
                return True
            return self._is_active(name)

        w = 60
        print("=" * w)
        print(f"  {self.seq.name.upper()}  |  {self.save_opts.prefix}")
        print("=" * w)

        # Scans
        print("\n  SCANS:")
        for name, ax in self.scans.items():
            v = ax.values
            primary = " (primary)" if name == self.scan_names[0] else ""
            print(f"    {name}{primary}:")
            print(f"      mode={ax.mode.value}  |  {len(v)} pts  "
                  f"|  [{v[0]:.6g} → {v[-1]:.6g}]")
            if ax.mode in (ScanMode.LINEAR, ScanMode.LOG):
                print(f"      start={ax.start:.6g}  stop={ax.stop:.6g}  step={ax.step}")
            if ax.shuffle:
                print(f"      shuffle=True")
            if show_values and len(v) > 0:
                if len(v) <= 8:
                    print(f"      values: {np.array2string(v, precision=4, separator=', ')}")
                else:
                    print(f"      values: [{v[0]:.4g}, {v[1]:.4g}, {v[2]:.4g}, "
                          f"... {v[-2]:.4g}, {v[-1]:.4g}]")

        # Instruments — only if configured (or show_all)
        if _show('mw'):
            print(f"\n  MICROWAVE:  {self.mw.power} dBm  |  {self.mw.freq:.6g} Hz")

        # if _show('laser'):
        #     print(f"\n  LASER:  t_AOM={self.laser.t_AOM:.6g} ns  |  power={self.laser.power}")
        
        print(f"\n  {self.times}")

        # if _show('pb'):
        print(f"\n  PULSEBLASTER:  {self.pb.clock_MHz} MHz  "
                f"(clk_cyc={self.pb.clk_cyc_ns:.1f} ns)")
        if self.pb.channels:
            ch_str = ', '.join(f"{k}={v}" for k, v in self.pb.channels.items())
            print(f"    channels: {ch_str}")

        print(f"\n  SEQUENCE:  {self.seq.name}")
        if self.seq.args_names:
            args = ', '.join(f"{n}={v:.6g}"
                             for n, v in zip(self.seq.args_names, self.seq.args_values))
            print(f"    args: {args}")
        print(f"    Nsamples={self.seq.Nsamples}")

        if _show('daq_ai'):
            print(f"\n  DAQ:  Nsamples={self.daq_ai.ai_samps_per_chan}  "
                  f"|  rate={self.daq_ai.ai_sample_rate:.0f} Sa/s")
        
        if _show('daq_ao'):
            print(f"\n  DAQ AO:  Nsamples={self.daq_ao.ao_samps_per_chan}  "
                  f"|  rate={self.daq_ao.ao_sample_rate:.0f} Sa/s")

        if _show('daq_ci'):
            print(f"\n  DAQ CI:  counter={self.daq_ci.ci_counter}  "
                  f"|  input={self.daq_ci.ci_input_terminal}  "
                  f"|  rate={self.daq_ci.ci_sample_rate:.0f} Sa/s")

        if _show('uhfli'):
            print(f"\n  UHFLI:  demod={self.uhfli.demod_index}  "
                  f"|  osc={self.uhfli.oscillator_freq:.6g} Hz  "
                  f"|  TC={self.uhfli.timeconstant:.2g} s  "
                  f"|  order={self.uhfli.filter_order}")

        det = self.runtime.detector
        print(f"\n  RUNTIME:  Nruns={self.runtime.Nruns}  "
              f"|  reload_pb={self.runtime.reload_pb}  "
              f"|  detector={det}")

        print(f"\n  PLOT:  {self.plot.x_label}")
        if self.extra:
            print(f"\n  EXTRA:  {self.extra}")
        print("=" * w)

    def summary(self) -> str:
        """One-line summary for logging."""
        parts = []
        for name, scan in self.scans.items():
            v = scan.values
            parts.append(f"{name}: {len(v)} pts [{v[0]:.4g} → {v[-1]:.4g}]")
        return f"{self.seq.name} | " + " | ".join(parts) + f" | Nruns={self.runtime.Nruns}"


# ═══════════════════════════════════════════════════════════════════════
# 5. CONFIG DIFF — compare two (or more) experiment configs
# ═══════════════════════════════════════════════════════════════════════

def _flatten(d: dict, prefix: str = '') -> dict[str, Any]:
    """Flatten nested dict to {'path.to.key': value}.

    Skips keys starting with '_'.
    Converts lists to tuples for hashability in set comparisons.
    """
    items = {}
    for k, v in d.items():
        if isinstance(k, str) and k.startswith('_'):
            continue
        path = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict):
            items.update(_flatten(v, path))
        elif isinstance(v, list):
            items[path] = tuple(v)
        else:
            items[path] = v
    return items


def diff_config(
    config_a: Union[ExperimentConfig, dict],
    config_b: Union[ExperimentConfig, dict],
    ignore_keys: Optional[set[str]] = None,
) -> dict[str, tuple[Any, Any]]:
    """Compare two configs and return changed paths.

    Args:
        config_a, config_b: ExperimentConfig instances or raw dicts
            (e.g. from .to_dict(), .load_json(), or UHFLIConfigManager).
        ignore_keys: Set of dotted-path prefixes to skip, e.g.
            {'scans.freq.values', 'save'} to ignore value arrays
            and save prefix.

    Returns:
        Dict mapping 'dotted.path' → (value_in_a, value_in_b)
        for every key that differs. Missing keys show as None.

    Example:
        diffs = diff_config(esr_config, esr_config_modified)
        for path, (old, new) in diffs.items():
            print(f"  {path}: {old} -> {new}")
    """
    da = config_a.to_dict() if isinstance(config_a, ExperimentConfig) else config_a
    db = config_b.to_dict() if isinstance(config_b, ExperimentConfig) else config_b

    flat_a = _flatten(da)
    flat_b = _flatten(db)

    all_keys = set(flat_a.keys()) | set(flat_b.keys())
    if ignore_keys:
        all_keys = {k for k in all_keys
                    if not any(k.startswith(ik) for ik in ignore_keys)}

    diffs = {}
    for k in sorted(all_keys):
        va = flat_a.get(k)
        vb = flat_b.get(k)
        if va != vb:
            diffs[k] = (va, vb)

    return diffs


def diff_configs_batch(
    configs: dict[str, Union[ExperimentConfig, dict]],
    reference: Optional[str] = None,
    ignore_keys: Optional[set[str]] = None,
) -> dict[str, dict[str, tuple[Any, Any]]]:
    """Compare multiple configs against a reference.

    Args:
        configs: Dict mapping label → config. E.g.:
            {'monday': config1, 'tuesday': config2, 'wednesday': config3}
        reference: Label of the reference config. If None, uses first key.
        ignore_keys: Passed to diff_config.

    Returns:
        Dict mapping label → diff_dict.

    Example:
        all_diffs = diff_configs_batch({
            'run1': ExperimentConfig.load_json('run1.json'),
            'run2': ExperimentConfig.load_json('run2.json'),
            'run3': ExperimentConfig.load_json('run3.json'),
        }, reference='run1')

        for label, diffs in all_diffs.items():
            print(f"\\n--- {label} vs run1 ---")
            for path, (old, new) in diffs.items():
                print(f"  {path}: {old} -> {new}")
    """
    labels = list(configs.keys())
    ref = reference or labels[0]
    ref_config = configs[ref]
    return {
        label: diff_config(ref_config, configs[label], ignore_keys=ignore_keys)
        for label in labels if label != ref
    }


def print_diff(
    diffs: dict[str, tuple[Any, Any]],
    label_a: str = 'A',
    label_b: str = 'B',
) -> None:
    """Pretty-print a diff_config result."""
    if not diffs:
        print(f"  (no differences between {label_a} and {label_b})")
        return
    max_path = max(len(k) for k in diffs)
    for path, (va, vb) in diffs.items():
        print(f"  {path:<{max_path}}  {va}  →  {vb}")
