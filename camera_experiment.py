# camera_experiment.py
"""
CameraExperiment — widefield NV diamond imaging using Hamamatsu ORCA.

Separate from DiodeExperiment: different detector (CameraWorker), different
threading model (Event/Condition), different data shape (uint16 images),
different PB runners (run_sequence_for_camera_*).

Shared with DiodeExperiment: Sweep, make_pb_setter, calc_contrast, savefile.

Data shape: (outer_combos, Nruns, inner_pts, Nframes, vsize, hsize)  uint16
    Nframes typically 2 (signal + reference per scan point)

Loop order (outer → Nruns → inner):
    for each outer_combo (e.g. mw_power):
        for i_run in range(Nruns):
            focus_adjust(i_run)  — if focus_interval matches
            for each inner_pt (e.g. freq):
                PBThread programs + triggers → camera reads Nframes
                callback(inner_idx, val, i_run, oc_idx, oc_vals,
                         processed, last_frame)

Phase 1: cam_levelm mode, static B field, no AO patterns.
"""

import numpy as np, time, threading, logging
from pathlib import Path
from typing import Optional, Callable, TYPE_CHECKING

from PBcontrol import PulseBlaster
from sweep_utils import Sweep, make_pb_setter
from experiment_base import calc_contrast, savefile, HW_BOUNDS

from experiment_config import ExperimentConfig
from Camcontrol import CameraWorker


# ═══════════════════════════════════════════════════════════════════════
#  Camera-specific helpers
# ═══════════════════════════════════════════════════════════════════════

def process_camera_data(i_max, data_slice):
    """
    Process camera frames → intensity + contrast.

    Parameters
    ----------
    i_max : int
        Number of valid inner points to process.
    data_slice : ndarray, shape (inner_pts, Nframes, vsize, hsize), uint16
        Nframes typically 2: index 0::2 = signal, index 1::2 = reference.

    Returns
    -------
    (mean_sig, mean_ref, contrast) : each shape (i_max,)
        mean_sig/ref = ROI pixel sum of averaged signal/reference frames.
        contrast = signal / reference ratio.
    """
    sig = data_slice[:i_max, 0::2, :, :]     # signal frames
    ref = data_slice[:i_max, 1::2, :, :]     # reference frames

    # Average over repeat frames, then sum over pixels (ROI already applied)
    mean_sig = np.mean(sig, axis=1).sum(axis=(1, 2)).astype(float)
    mean_ref = np.mean(ref, axis=1).sum(axis=(1, 2)).astype(float)

    contrast = np.where(mean_ref > 0, mean_sig / mean_ref, np.nan)
    return mean_sig, mean_ref, contrast


# ═══════════════════════════════════════════════════════════════════════
#  PB Thread — programs PB per scan point, signals camera via Event
# ═══════════════════════════════════════════════════════════════════════

class _CameraPBThread(threading.Thread):
    """Scan-control thread for camera experiments.

    For each inner scan point:
        1. Program PB sequence (or set instrument)
        2. Call appropriate PB runner (e.g. run_sequence_for_camera_level_trigger_many)
        3. Set trigger_event → camera thread reads frames
        4. Wait on condition until camera is done
    """

    def __init__(self, experiment, trigger_event, condition):
        super().__init__(daemon=True)
        self.exp = experiment
        self.trigger_event = trigger_event
        self.condition = condition
        self.inst_set_times = []
        self.error = None

    def run(self):
        try:
            self._scan_loop()
        except Exception as e:
            self.error = e
            logging.exception("PBThread error")
            # Ensure camera thread doesn't hang
            self.trigger_event.set()

    def _scan_loop(self):
        exp = self.exp
        inner_vals = exp.sweep.inner_values
        n_inner = exp.sweep.inner_count

        for order_pos in range(n_inner):
            if exp._stop_requested:
                break

            actual_i = int(exp._inner_order[order_pos])
            actual_val = inner_vals[actual_i]

            t0 = time.perf_counter()

            # ESR-like: set SG frequency
            if exp._is_freq_sweep:
                exp.sg.set_freq(actual_val)
            else:
                # PB-timing sweep: update seqArgList and reprogram
                exp._seq_arg_list[exp._inner_pidx] = actual_val

            # Program PB and run the camera sequence
            _, the_list = PulseBlaster.PB_program(
                exp.instr_mode, exp.sequence,
                exp._seq_arg_list + [exp.pb_channels])
            instruction_list = [the_list[i][0] for i in range(len(the_list))]

            # Dispatch to the correct PB runner
            exp._run_pb_camera(instruction_list, actual_i)

            self.inst_set_times.append(time.perf_counter() - t0)

            # Signal camera to start reading frames
            self.trigger_event.set()

            # Wait for camera to finish reading frames at this scan point
            with self.condition:
                while self.trigger_event.is_set():
                    self.condition.wait()


# ═══════════════════════════════════════════════════════════════════════
#  CameraExperiment
# ═══════════════════════════════════════════════════════════════════════

class CameraExperiment:
    """
    Camera-based NV experiment.  Separate from DiodeExperiment.

    Shared: Sweep, make_pb_setter, calc_contrast, savefile
    Different: detector (CameraWorker), threading (Event/Condition),
               data shape (uint16 images), PB runner methods
    """

    # ESR / frequency-only sequences — PB fixed, sweep via SG
    FREQ_ONLY_SEQUENCES = frozenset([
        'esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq'])

    # Instrument-setting parameters (not PB timing)
    INSTRUMENT_PARAMS = frozenset([
        'freq', 'mw_freq', 'mw_power', 'mod_rate', 'mod_dev'])

    SETTER_REGISTRY = {
        'freq':      lambda self: self.sg.set_freq,
        'mw_freq':   lambda self: self.sg.set_freq,
        'mw_power':  lambda self: self.sg.set_amp_rf,
        'mod_rate':  lambda self: self.sg.set_mod_rate,
        'mod_dev':   lambda self: self.sg.set_mod_dev,
    }

    def __init__(self, instruments: dict, config: 'ExperimentConfig'):
        self.sg = instruments['sg']
        self.pb = instruments['pb']
        self.ao_task = instruments.get('ao_task')
        self.camera_worker: 'CameraWorker' = instruments['camera_worker']

        self.cfg = config
        self.cw = self.camera_worker  # short alias

        # Camera config
        cam_cfg = config.camera
        # assert cam_cfg
        self.instr_mode: str = cam_cfg.instr_mode
        self.exposure: float = cam_cfg.exposure_s            # seconds
        self.roi: list = cam_cfg.roi
        self.Nframes: int = cam_cfg.frames_per_cycle
        self.focus_interval: int = cam_cfg.focus_interval
        self.focus_time: float = cam_cfg.min_focus_time_s    # seconds

        # Field config
        field_cfg = config.field_cfg
        # assert field_cfg
        self.align_field = field_cfg.align_field
        self.test_field = field_cfg.test_field
        self.t_align_dc_ms = field_cfg.t_align_dc_ms

        # Sequence
        # assert config.seq.name and config.seq.Nsamples and config.pb.channels
        self.sequence: str = config.seq.name
        self.Nsamples_cfg: int = config.seq.Nsamples  # = frames_per_cycle
        self.Nruns: int = config.runtime.Nruns
        self.pb_channels: list = config.pb.channels

        # Primary scan
        scan_names = config.scan_names
        self.scan_name: str = scan_names[0] if scan_names else ''
        self.primary_scan = config.scans.get(self.scan_name)
        self.param_values = (
            self.primary_scan.values if self.primary_scan else np.array([]))
        self.Nscanpts = len(self.param_values)

        self._is_freq_sweep = self.sequence in self.FREQ_ONLY_SEQUENCES

        # State
        self.sweep: Optional[Sweep] = None
        self.data_array: Optional[np.ndarray] = None
        self.scan_times: list[float] = []
        self.elapsed: float = 0.0
        self.i_outer: int = 0
        self.i_run: int = 0
        self._stop_requested: bool = False
        self._inner_order: np.ndarray = np.array([])

        # PB sequence arg list — built during setup
        self._seq_arg_list: list = []
        self._inner_pidx: int = 0

    # ══════════════════════════════════════════════════════════════════
    #  SETUP
    # ══════════════════════════════════════════════════════════════════

    def setup(self):
        """Configure camera, build sweep, allocate data array."""
        cfg = self.cfg
        cw = self.camera_worker

        # 1. Configure camera for triggered acquisition
        cw.configure_camera(self.instr_mode)
        cw.exposure = self.exposure
        if not cw.simulate:
            import dcamcon
            cw.hdcamcon.set_propertyvalue(
                dcamcon.DCAM_IDPROP.EXPOSURETIME, self.exposure)

        # Set ROI if provided
        if self.roi and self.roi != [0, 0, 2048, 2048]:
            import dcamcon
            cw.roi = self.roi
            cw.subarray_mode = dcamcon.DCAMPROP.MODE.ON
            cw.set_roi()
        else:
            self.roi = cw.roi  # read back whatever is set

        # Frame size from ROI: [x0, y0, width, height]
        self.hsize = self.roi[2] if len(self.roi) >= 4 else 2048
        self.vsize = self.roi[3] if len(self.roi) >= 4 else 2048

        # 2. Build seq arg list for PB programming
        # assert cfg.seq.args_names and cfg.seq.args_values
        args_names = cfg.seq.args_names
        self._seq_arg_list = list(cfg.seq.args_values)

        # Determine inner param index in seq_args
        if self._is_freq_sweep:
            self._inner_pidx = -1  # not used
        elif self.scan_name in args_names:
            self._inner_pidx = args_names.index(self.scan_name)
        else:
            # Legacy prepend
            self._seq_arg_list = [self.param_values[0]] + self._seq_arg_list
            self._inner_pidx = 0

        # 3. Compute per-scan-point timing for PB runner
        # t_seq_total: array of total PB sequence times per scan point
        # (needed by run_sequence_for_camera_level_trigger_many)
        self._compute_seq_totals()

        # 4. Build Sweep
        self.sweep = Sweep()
        self._build_sweep(cfg)
        print(self.sweep.info())

        # 5. Set SG
        if self.sequence not in ['aom_timing', 'T1ms0_train', 'drift_seq']:
            # assert cfg.mw
            self.sg.set_freq(cfg.mw.freq)
            self.sg.set_amp_rf(cfg.mw.power)

        # 6. Set AO field
        if self.ao_task is not None:
            self.ao_task.set_outputs_to_constant(
                output_field_in_gauss=self.test_field)

        # 7. Allocate data array: (outer, Nruns, inner, Nframes, vsize, hsize)
        self.data_array = np.zeros((
            self.sweep.outer_count,
            self.Nruns,
            self.sweep.inner_count,
            self.Nframes,
            self.vsize,
            self.hsize,
        ), dtype=np.uint16)

        # 8. Shuffle order
        shuffle = self.primary_scan.shuffle if self.primary_scan else False
        self._shuffle = shuffle
        if shuffle:
            self._inner_order = np.random.permutation(self.sweep.inner_count)
            print("⚡ Inner sweep shuffled")
        else:
            self._inner_order = np.arange(self.sweep.inner_count)

        # 9. Start camera capture (prep_only — allocate buffer, start capture)
        buffer_size = self.Nframes + 2  # a few extra for safety
        cw.start_capture(buffer_size=buffer_size, prep_only=True)

        # 10. Prepare AO pattern for trigger_ao modes
        self._ao_pattern = None
        if 'trigger_ao' in self.instr_mode and self.ao_task is not None:
            self._ao_pattern = self._prepare_ao_pattern()
            if self._ao_pattern is not None:
                print(f"▶ AO pattern: {self._ao_pattern.shape}")

        print(f"▶ Camera: {self.instr_mode}  "
              f"exp={self.exposure*1e3:.1f}ms  "
              f"ROI={self.roi}  "
              f"frames/cyc={self.Nframes}")
        print(f"▶ data {self.data_array.shape}  "
              f"({self.sweep.outer_count} outer × {self.Nruns} runs × "
              f"{self.sweep.inner_count} inner × "
              f"{self.Nframes} frames × {self.vsize}×{self.hsize})")

    def _compute_seq_totals(self):
        """Compute per-scan-point total PB sequence time.

        Used by PB runners to calculate N_total (loop counts) and buffer.
        For constant-length sequences (ESR), this is a single value repeated.
        For variable-length (Rabi), each scan point has its own total.
        """
        cfg = self.cfg
        # Use PB_program to get instruction times for each scan point
        # For simplicity in Phase 1, compute from the sequence:
        # t_seq_total = sum of all instruction durations in one cycle
        # We get this from PB_program's the_list[i][4] (inst_times dict)

        t_totals = []
        for val in self.param_values:
            sal = list(self._seq_arg_list)
            if not self._is_freq_sweep:
                sal[self._inner_pidx] = val
            _, the_list = PulseBlaster.PB_program(
                self.instr_mode, self.sequence, sal + [self.pb_channels])
            # the_list is [(inst_list, n_err, err_inst, err_times, inst_times), ...]
            # inst_times is a dict {inst_idx: time_ns}
            # Total time = sum of all instruction times for both signal + ref
            total = 0
            for part in the_list:
                if len(part) > 4 and isinstance(part[4], dict):
                    total = sum(part[4].values())
                else:
                    # Fallback: estimate from instruction list
                    for inst in part[0]:
                        total += inst[3]  # inst[3] = duration
            t_totals.append(total)

        # PB runners expect t_seq_total as a list (one per sub-sequence, sig+ref)
        # For cam_levelm: two sub-sequences → total for each
        self._t_seq_total_per_point = t_totals

        # Also build the per-point [sig_total, ref_total] arrays
        # For simplicity: the PB_program returns one list entry per sub-seq
        self._t_seq_total_pairs = []
        for val in self.param_values:
            sal = list(self._seq_arg_list)
            if not self._is_freq_sweep:
                sal[self._inner_pidx] = val
            _, the_list = PulseBlaster.PB_program(
                self.instr_mode, self.sequence, sal + [self.pb_channels])
            pair = []
            for part in the_list:
                t = 0
                if len(part) > 4 and isinstance(part[4], dict):
                    t = sum(part[4].values())
                else:
                    for inst in part[0]:
                        t += inst[3]
                pair.append(t)
            self._t_seq_total_pairs.append(pair)

    def _prepare_ao_pattern(self):
        """Build the AO voltage pattern for trigger_ao rotating field modes.

        Uses FieldConfig align_field and t_align_dc_ms to create the
        retriggerable AO waveform that the DAQ outputs synchronously with PB.

        Returns ndarray (3, n_samples) or None if not applicable.
        """
        cfg = self.cfg
        field_cfg = cfg.field_cfg
        # assert field_cfg
        align_field = np.array(field_cfg.align_field, dtype=float)

        if np.allclose(align_field, 0):
            print("  ℹ AO pattern: align_field is zero — skipping")
            return None

        try:
            from DAQcontrol import AnalogOutputTask
            # Build a simple DC alignment pattern
            # The pattern length is determined by the alignment time
            samp_rate = 10e3  # default AO sample rate
            t_align_s = field_cfg.t_align_dc_ms / 1e3
            n_samples = max(int(samp_rate * t_align_s), 10)

            # Simple pattern: field ON for half, OFF for half
            pattern = np.zeros((3, n_samples))
            half = n_samples // 2
            # First half: Bz + MW
            pattern[:, :half] = (align_field / AnalogOutputTask.coil_calibration
                                 ).reshape(3, 1)
            # Second half: Bx + By
            # (simplified — the full rotating pattern is experiment-specific)
            pattern[:, half:] = (align_field / AnalogOutputTask.coil_calibration
                                 ).reshape(3, 1)

            return AnalogOutputTask.prepare_data_for_write(pattern)
        except Exception as e:
            print(f"⚠ AO pattern prep: {e}")
            return None

    def _build_sweep(self, cfg):
        """Add inner + outer axes to self.sweep."""
        # Inner axis — for camera experiments, sweep only calls setter for freq;
        # PB reprogramming is done inside the PBThread, not via Sweep setter
        # assert self.sweep
        if self._is_freq_sweep:
            self.sweep.add(self.scan_name, self.param_values,
                           setter=self.sg.set_freq)
        else:
            # PB-timing sweep: no setter in Sweep — PBThread handles it
            self.sweep.add(self.scan_name, self.param_values, setter=None)

        # Outer axes
        for name in cfg.scan_names[1:]:
            axis = cfg.scans.get(name)
            if axis is None or len(axis.values) <= 1:
                continue
            self.sweep.add(name, axis.values,
                           setter=self._resolve_setter(name))

    def _resolve_setter(self, name: str) -> Optional[Callable]:
        factory = self.SETTER_REGISTRY.get(name)
        if factory is not None:
            return factory(self)
        print(f"⚠ No setter for '{name}'")
        return None

    # ══════════════════════════════════════════════════════════════════
    #  EXECUTE
    # ══════════════════════════════════════════════════════════════════

    def execute(
        self,
        callback: Optional[Callable] = None,
        stop_check: Optional[Callable] = None,
        save_path: Optional[Path] = None,
        folder_number: Optional[str] = None,
    ) -> np.ndarray:
        """
        Two-thread acquisition: PBThread + camera read in main thread.

        callback(inner_idx, param_val, i_run, oc_idx, oc_vals,
                 processed, last_frame)
            processed = (mean_sig, mean_ref, contrast) arrays
            last_frame = uint16 ndarray (vsize, hsize) — the signal frame
        """
        self._stop_requested = False
        self.scan_times = []
        t_start = time.perf_counter()

        cw = self.camera_worker
        timeout_ms = int(self.exposure * 1e3 + 500)
        # assert self.sweep and self.data_array

        try:
            for oc_idx, oc_vals in self.sweep.outer_combos():
                if self._stop_requested:
                    break
                self.i_outer = oc_idx
                if oc_vals:
                    lbl = "  ".join(f"{k}={v:.6g}"
                                    for k, v in oc_vals.items())
                    print(f"▶ Outer [{oc_idx+1}/"
                          f"{self.sweep.outer_count}]: {lbl}")

                for self.i_run in range(self.Nruns):
                    if stop_check and stop_check():
                        self._stop_requested = True
                        break
                    print(f"  Run {self.i_run+1}/{self.Nruns}")

                    # Focus adjustment
                    if (self.focus_interval > 0 and
                            self.i_run % self.focus_interval == 0):
                        self._focus_adjust(
                            self.i_run,
                            callback=callback)

                    # Start retriggerable AO for trigger_ao modes
                    if (self._ao_pattern is not None
                            and 'trigger_ao' in self.instr_mode
                            and self.ao_task is not None):
                        try:
                            self.ao_task.create_retriggerable_ao_task(
                                self._ao_pattern.shape)
                            self.ao_task.start_retriggerable_ao_task(
                                self._ao_pattern)
                        except Exception as e:
                            print(f"⚠ AO retrig start: {e}")

                    # Create synchronization primitives
                    trigger_event = threading.Event()
                    condition = threading.Condition()

                    # Start PB thread
                    pb_thread = _CameraPBThread(self, trigger_event, condition)
                    pb_thread.start()

                    # Camera read loop (runs in this thread)
                    inner_vals = self.sweep.inner_values
                    n_inner = self.sweep.inner_count

                    for order_pos in range(n_inner):
                        if self._stop_requested:
                            break

                        actual_i = int(self._inner_order[order_pos])
                        actual_val = inner_vals[actual_i]

                        # Wait for PB thread to signal
                        trigger_event.wait()

                        t0 = time.perf_counter()

                        # Start the PB sequence (PB already programmed by PBThread)
                        self.pb.start_sequence()

                        # Read Nframes from camera
                        for f_idx in range(self.Nframes):
                            frame_ready = cw.hdcamcon.wait_capevent_frameready(
                                timeout_ms)
                            if frame_ready is True:
                                if f_idx == self.Nframes - 1:
                                    # Last frame: stop PB first, then read
                                    self.pb.stop_sequence()

                                frame = cw.hdcamcon.get_lastframedata()
                                if frame is not False:
                                    self.data_array[
                                        oc_idx, self.i_run, actual_i,
                                        f_idx, :, :] = frame
                                else:
                                    print(f"⚠ No frame at pt {actual_i} "
                                          f"frame {f_idx}")
                            else:
                                print(f"⚠ Frame timeout at pt {actual_i} "
                                      f"frame {f_idx}")
                                if f_idx == self.Nframes - 1:
                                    self.pb.stop_sequence()

                        self.scan_times.append(time.perf_counter() - t0)

                        # Notify PBThread that camera is done
                        with condition:
                            trigger_event.clear()
                            condition.notify_all()

                        # Progress
                        if order_pos % max(1, n_inner // 10) == 0:
                            print(f"    {order_pos+1}/{n_inner}: "
                                  f"{actual_val:.6g}")

                        # Callback: per-point update
                        if callback is not None:
                            self._fire_callback(
                                callback, actual_i, actual_val,
                                self.i_run, oc_idx, oc_vals)

                    # Wait for PB thread to finish
                    pb_thread.join(timeout=10)
                    if pb_thread.error is not None:
                        print(f"⚠ PBThread error: {pb_thread.error}")

                    # Auto-save after each run
                    self._autosave_run(save_path, folder_number)

        except KeyboardInterrupt:
            print("🛑 Interrupted")
            self._stop_requested = True

        self._finish(t_start, save_path, folder_number)
        return self.data_array

    def _run_pb_camera(self, instruction_list, scan_idx):
        """Dispatch to the correct PB runner method based on instr_mode.

        Called by _CameraPBThread for each scan point.
        The PB sequence is programmed but NOT started — that happens
        in the camera-read loop after trigger_event.wait().
        """
        t_exp_ns = self.exposure * 1e9
        N_total = []  # empty = auto-calculate in PB runner
        Nsamples = self.Nsamples_cfg

        # Get per-point seq totals
        t_seq_total_pair = self._t_seq_total_pairs[scan_idx]

        if self.instr_mode == 'cam_levelm':
            self.pb.run_sequence_for_camera_level_trigger_many(
                instruction_list, t_exp_ns, t_seq_total_pair,
                N_total, Nsamples)

        elif self.instr_mode == 'cam_syncm':
            t_align_ns = self.t_align_dc_ms * 1e6  # ms → ns
            self.pb.run_sequence_for_camera_sync_trigger_many(
                instruction_list, t_exp_ns, t_align_ns,
                t_seq_total_pair, N_total, Nsamples)

        elif 'levelm_trigger' in self.instr_mode:
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.run_sequence_for_camera_level_trigger_many_bac(
                instruction_list, t_exp_ns, t_align_ns,
                t_seq_total_pair, N_total, Nsamples)

        elif 'syncm_trigger_ao' in self.instr_mode:
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.run_sequence_for_camera_sync_trigger_many_bac(
                instruction_list, t_exp_ns, t_align_ns,
                t_seq_total_pair, N_total, Nsamples)

        elif self.instr_mode == 'cam_timeseries':
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.custom_trigger(t_exp_ns, t_align_ns)

        elif self.instr_mode == 'cam_timeseries_trigger_ao':
            t_align_ns = self.t_align_dc_ms * 1e6
            # custom_trigger_rot_field expects:
            #   t_exposure, t_align_rot, t_measurement, t_align_rot_extended
            # For Phase 2: use config.extra for t_meas and extended alignment
            t_meas_ns = self.cfg.extra.get('t_meas_ns', 50e6)  # 50 ms default
            t_align_ext_ns = self.cfg.extra.get('t_align_extended_ns', t_align_ns)
            self.pb.custom_trigger_rot_field(
                t_exp_ns, t_align_ns, t_meas_ns, t_align_ext_ns)

        else:
            raise ValueError(f"Unknown camera mode: {self.instr_mode}")

    def _fire_callback(self, callback, actual_i, actual_val,
                       i_run, oc_idx, oc_vals):
        """Process current run's data and fire per-point callback."""
        # assert self.data_array
        # Process accumulated data for current (oc, run)
        data_slice = self.data_array[oc_idx, i_run]  # (inner, Nframes, v, h)
        # Count valid (non-zero) inner points
        sums = data_slice.sum(axis=(1, 2, 3))
        n_filled = int(np.count_nonzero(sums))

        processed = None
        if n_filled > 0:
            processed = process_camera_data(n_filled, data_slice)

        # Last signal frame at current scan point
        last_frame = self.data_array[oc_idx, i_run, actual_i, 0, :, :]

        callback(actual_i, actual_val, i_run, oc_idx, oc_vals,
                 processed, last_frame)

    # ══════════════════════════════════════════════════════════════════
    #  FOCUS ADJUSTMENT
    # ══════════════════════════════════════════════════════════════════

    def _focus_adjust(self, i_run, callback=None):
        """Switch to internal trigger, show live frames, switch back.

        During focus, the callback receives (frame,) for pyqtgraph display.
        """
        cw = self.camera_worker
        print(f"  🔍 Focus adjustment (run {i_run+1})...")

        try:
            import dcamcon

            # Switch to internal trigger for live view
            cw.hdcamcon.set_propertyvalue(
                dcamcon.DCAM_IDPROP.TRIGGERSOURCE,
                dcamcon.DCAMPROP.TRIGGERSOURCE.INTERNAL)

            # PB runs focus sequence (laser ON, field channels open)
            self.pb.focus_adjustment_sequence(
                self.instr_mode,
                t_exposure=self.exposure * 1e9,
                Nsamples=self.Nframes,
                focus_time=self.focus_time * 1e9)
            self.pb.start_sequence()

            # Live loop — callback sends frames to session for pyqtgraph display
            t_start = time.perf_counter()
            timeout_focus_ms = int(self.exposure * 1e3 + 50)

            while (time.perf_counter() - t_start) < self.focus_time:
                res = cw.hdcamcon.wait_capevent_frameready(
                    timeout_millisec=timeout_focus_ms)
                if res is True:
                    frame = cw.hdcamcon.get_lastframedata()
                    if frame is not False and callback is not None:
                        # Send frame through callback for display
                        # Convention: when processed=None, frame is focus frame
                        callback(-1, 0.0, i_run, self.i_outer, {},
                                 None, frame)

            self.pb.stop_sequence()

            # Switch back to external trigger
            cw.configure_camera(self.instr_mode)

            # Clear stale frame from buffer (level trigger artifact)
            if 'level' in self.instr_mode:
                cw.hdcamcon.wait_capevent_frameready(
                    timeout_millisec=int(self.exposure * 1e3 + 50))
                cw.hdcamcon.get_lastframedata()

            print(f"  ✔ Focus done ({time.perf_counter() - t_start:.1f}s)")

        except ImportError:
            print("  ⚠ dcamcon not available — skipping focus")

    # ══════════════════════════════════════════════════════════════════
    #  SAVE
    # ══════════════════════════════════════════════════════════════════

    def _autosave_run(self, save_path, folder_number):
        """Save TIFF per (outer, run) + cumulative npy after each run."""
        if not save_path or not folder_number:
            return

        try:
            import tifffile
        except ImportError:
            print("  ⚠ tifffile not installed — skipping TIFF save")
            return
        # assert self.sweep and self.data_array
        for oc in range(self.sweep.outer_count):
            stack = self.data_array[oc, self.i_run]  # (inner, Nframes, v, h)
            stack_flat = stack.reshape(-1, *stack.shape[-2:])  # flatten inner×frames
            tiff_name = (save_path /
                f"frames_{folder_number}_outer{oc}_run{self.i_run:02d}.tiff")
            tifffile.imwrite(str(tiff_name), stack_flat, ome=True)

        # Cumulative npy (crash-safe)
        npy_path = save_path / f"data_{folder_number}.npy"
        np.save(str(npy_path), self.data_array)
        print(f"    💾 Run {self.i_run+1} → {npy_path.name}")

    def _finish(self, t_start, save_path, folder_number):
        self.elapsed = time.perf_counter() - t_start
        if self.scan_times:
            print(f"✔ {self.elapsed:.1f}s  "
                  f"({np.mean(self.scan_times)*1e3:.1f} ± "
                  f"{np.std(self.scan_times)*1e3:.1f} ms/pt)")
        if save_path and folder_number:
            self._save_params(save_path, folder_number)

    def _save_params(self, save_path: Path, folder_number: str):
        """Save config + sweep metadata."""
        d = self.cfg.to_dict()
        # assert self.sweep and self.data_array        
        d['_sweep_info'] = {
            'inner': self.sweep.inner_name,
            'inner_count': self.sweep.inner_count,
            'outer_names': self.sweep.outer_names,
            'outer_shape': list(self.sweep.outer_shape),
            'data_shape': list(self.data_array.shape),
            'shuffled': self._shuffle,
        }
        d['_timing'] = {
            'elapsed_s': self.elapsed,
            'mean_per_pt_ms': (float(np.mean(self.scan_times) * 1e3)
                               if self.scan_times else 0),
            'std_per_pt_ms': (float(np.std(self.scan_times) * 1e3)
                              if self.scan_times else 0),
        }
        d['_camera'] = {
            'instr_mode': self.instr_mode,
            'exposure_s': self.exposure,
            'roi': self.roi,
            'frames_per_cycle': self.Nframes,
            'frame_shape': [self.vsize, self.hsize],
        }
        pf = save_path / f"params_{folder_number}.yaml"
        savefile(pf, d)

        vf = save_path / f"scan_values_{folder_number}.npz"
        arrays = {name: axis.values
                  for name, axis in self.cfg.scans.items()}
        np.savez_compressed(str(vf), **arrays)
        print(f"  💾 Params → {pf.name}, values → {vf.name}")

    # ══════════════════════════════════════════════════════════════════
    #  TEARDOWN
    # ══════════════════════════════════════════════════════════════════

    def teardown(self):
        """Stop capture, release buffer. Do NOT uninit camera (session owns it)."""
        try:
            self.camera_worker.stop_capture(free_buffer=True)
        except Exception as e:
            print(f"⚠ Camera teardown: {e}")


# ═══════════════════════════════════════════════════════════════════════
#  CameraTimeseriesExperiment — continuous frame monitoring
# ═══════════════════════════════════════════════════════════════════════

class CameraTimeseriesExperiment:
    """Continuous PB-triggered camera acquisition at a fixed parameter point.

    Camera analogue of the DAQ-based TimeseriesExperiment:
    - PB programs the same pulse sequence as sweep (ESR, Rabi, etc.) at a
      fixed parameter value, then runs it in a continuous loop (BRANCH)
    - Camera is in EXTERNAL trigger mode — PB level/sync triggers each frame
    - Each PB cycle produces Nframes (signal + reference), exactly like sweep
    - ROI pixel-mean of each frame → mean_sig, mean_ref → contrast per cycle
    - Ring buffers for display: contrast, sig intensity, ref intensity
    - 30 Hz timer in session reads the rings and updates scrolling plots

    This is NOT free-running camera capture.  The PB sequence determines
    when each frame is acquired, what MW/laser state it captures, and
    the signal vs reference assignment.  The camera just responds to
    external triggers and delivers frames.

    Usage (from session):
        exp = CameraTimeseriesExperiment(instruments, config)
        exp.setup()
        exp.start()       # PB starts, polling thread reads frames
        ...               # session timer reads exp.contrast_ring etc.
        exp.stop()        # PB stops, thread exits
        exp.save('cam_ts')
        exp.teardown()
    """

    FREQ_ONLY = frozenset([
        'esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq'])

    def __init__(self, instruments: dict, config: 'ExperimentConfig'):
        self.sg = instruments['sg']
        self.pb = instruments['pb']
        self.ao_task = instruments.get('ao_task')
        self.camera_worker = instruments['camera_worker']

        self.cfg = config
        self.cw = self.camera_worker

        # Camera config
        cam_cfg = config.camera
        # assert cam_cfg
        self.instr_mode = cam_cfg.instr_mode
        self.exposure = cam_cfg.exposure_s
        self.roi = cam_cfg.roi
        self.Nframes = cam_cfg.frames_per_cycle

        # Field config
        field_cfg = config.field_cfg
        # assert field_cfg
        self.test_field = field_cfg.test_field
        self.t_align_dc_ms = field_cfg.t_align_dc_ms

        # Sequence
        self.sequence = config.seq.name
        self.Nsamples_cfg = config.seq.Nsamples
        self.pb_channels = config.pb.channels

        # Timeseries config
        rt = config.runtime
        self.display_seconds = getattr(rt, 'ts_display_seconds', 30.0)
        self.max_duration = getattr(rt, 'ts_max_duration', 0.0)

        # Fixed scan point (same logic as DAQ TimeseriesExperiment)
        scan_names = config.scan_names
        self.scan_name = scan_names[0] if scan_names else ''
        scan = config.scans.get(self.scan_name)
        indices = getattr(rt, 'seq_plot_indices', [0, -1])
        if scan is not None and len(scan.values) > 0:
            idx = indices[-1]
            if idx < 0:
                idx = len(scan.values) + idx
            self.fixed_value = scan.values[idx]
        else:
            self.fixed_value = 0.0

        # State
        self.contrast_ring = None
        self.sig_ring = None
        self.ref_ring = None
        self.full_contrasts = []
        self.full_sig = []
        self.full_ref = []
        self.last_frame = None
        self.cycle_count = 0
        self.elapsed = 0.0
        self._start_time = 0.0
        self._running = False
        self._stop_requested = False
        self._poll_thread = None

    # ── setup ───────────────────────────────────────────────────────

    def setup(self):
        """Program PB at fixed point, configure camera for external trigger."""
        cfg = self.cfg
        cw = self.camera_worker

        # ── SG ──
        if self.sequence not in ['aom_timing', 'T1ms0_train', 'drift_seq']:
            # assert cfg.mw
            self.sg.set_freq(cfg.mw.freq)
            self.sg.set_amp_rf(cfg.mw.power)

        is_freq = self.sequence in self.FREQ_ONLY
        if is_freq:
            self.sg.set_freq(self.fixed_value)

        # ── PB: program at fixed point ──
        args_names = cfg.seq.args_names
        seq_args = list(cfg.seq.args_values)
        if is_freq:
            sal = seq_args
        elif self.scan_name in args_names:
            sal = list(seq_args)
            sal[args_names.index(self.scan_name)] = self.fixed_value
        else:
            sal = [self.fixed_value] + seq_args

        _, the_list = PulseBlaster.PB_program(
            self.instr_mode, self.sequence, sal + [self.pb_channels])
        instruction_list = [the_list[i][0] for i in range(len(the_list))]

        # Compute per-sequence timing for PB runner
        t_seq_total_pair = []
        for part in the_list:
            t = 0
            if len(part) > 4 and isinstance(part[4], dict):
                t = sum(part[4].values())
            else:
                for inst in part[0]:
                    t += inst[3]
            t_seq_total_pair.append(t)

        # ── PB runner: programs the continuous trigger sequence ──
        # Same runner as sweep — PB generates external triggers that gate
        # camera exposures.  The Nsamples loop count makes PB run the
        # signal+reference cycle Nsamples times, then STOP.  Since we
        # want continuous, we pass a large Nsamples (PB loops internally).
        t_exp_ns = self.exposure * 1e9
        N_total = []  # auto-calculate

        if self.instr_mode == 'cam_levelm':
            self._pb_returns = self.pb.run_sequence_for_camera_level_trigger_many(
                instruction_list, t_exp_ns, t_seq_total_pair,
                N_total, self.Nsamples_cfg)
        elif self.instr_mode == 'cam_syncm':
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.run_sequence_for_camera_sync_trigger_many(
                instruction_list, t_exp_ns, t_align_ns,
                t_seq_total_pair, N_total, self.Nsamples_cfg)
        elif 'levelm_trigger' in self.instr_mode:
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.run_sequence_for_camera_level_trigger_many_bac(
                instruction_list, t_exp_ns, t_align_ns,
                t_seq_total_pair, N_total, self.Nsamples_cfg)
        elif 'syncm_trigger_ao' in self.instr_mode:
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.run_sequence_for_camera_sync_trigger_many_bac(
                instruction_list, t_exp_ns, t_align_ns,
                t_seq_total_pair, N_total, self.Nsamples_cfg)
        elif 'timeseries' in self.instr_mode:
            t_align_ns = self.t_align_dc_ms * 1e6
            self.pb.custom_trigger(t_exp_ns, t_align_ns)
        else:
            raise ValueError(f"Unknown camera mode: {self.instr_mode}")

        # ── Configure camera: EXTERNAL trigger (PB-governed) ──
        cw.configure_camera(self.instr_mode)
        cw.exposure = self.exposure
        if not cw.simulate:
            import dcamcon
            cw.hdcamcon.set_propertyvalue(
                dcamcon.DCAM_IDPROP.EXPOSURETIME, self.exposure)

        # ROI
        if self.roi and self.roi != [0, 0, 2048, 2048]:
            import dcamcon
            cw.roi = self.roi
            cw.subarray_mode = dcamcon.DCAMPROP.MODE.ON
            cw.set_roi()
        else:
            self.roi = cw.roi

        self.hsize = self.roi[2] if len(self.roi) >= 4 else 2048
        self.vsize = self.roi[3] if len(self.roi) >= 4 else 2048

        # AO field
        if self.ao_task is not None:
            self.ao_task.set_outputs_to_constant(
                output_field_in_gauss=self.test_field)

        # ── Buffers ──
        # Estimate contrast rate from PB cycle time
        # One contrast point per Nframes (signal+reference)
        # Cycle time ≈ Nframes × (exposure + cam_response + gap)
        t_cam_response = 87.7e-6 if 'level' in self.instr_mode else 38.96e-6
        cycle_time = self.Nframes * (self.exposure + t_cam_response + 8e-3)
        self.contrast_rate = 1.0 / cycle_time
        contrast_cap = int(self.display_seconds * self.contrast_rate)

        from ring_buffer import RingBuffer
        self.contrast_ring = RingBuffer(max(contrast_cap, 100))
        self.sig_ring = RingBuffer(max(contrast_cap, 100))
        self.ref_ring = RingBuffer(max(contrast_cap, 100))
        self.full_contrasts = []
        self.full_sig = []
        self.full_ref = []
        self.last_frame = None
        self.cycle_count = 0

        # Start camera capture (buffer allocated, capture started, waiting for triggers)
        buffer_size = self.Nframes + 4
        cw.start_capture(buffer_size=buffer_size, prep_only=True)

        print(f"▶ CamTS: {self.sequence} @ "
              f"{self.scan_name}={self.fixed_value:.4g}  "
              f"mode={self.instr_mode}  "
              f"contrast_rate≈{self.contrast_rate:.1f} pts/s  "
              f"display={self.display_seconds}s")

    # ── start / stop ────────────────────────────────────────────────

    def start(self):
        """Begin PB-triggered continuous camera acquisition."""
        if self._running:
            print("⚠ Already running")
            return

        self._stop_requested = False
        self._start_time = time.perf_counter()
        self._running = True

        # Start PB — begins external triggering of camera
        self.pb.start_sequence()

        # Start polling thread that reads externally-triggered frames
        self._poll_thread = _CameraPollingThread(self)
        self._poll_thread.start()

        # Auto-stop
        if self.max_duration > 0:
            self._auto_stop_time = self._start_time + self.max_duration
        else:
            self._auto_stop_time = None

        print("▶ CamTS: PB started → camera acquiring")

    def stop(self):
        """Stop PB and camera acquisition."""
        if not self._running:
            return

        self._stop_requested = True
        if self._poll_thread is not None:
            self._poll_thread.join(timeout=5)
            self._poll_thread = None

        self.pb.stop_sequence()
        self._running = False
        self.elapsed = time.perf_counter() - self._start_time

        n_c = self.contrast_ring.total_count if self.contrast_ring else 0
        print(f"⏹ CamTS: stopped — {self.elapsed:.1f}s, "
              f"{n_c} contrast pts, {self.cycle_count} cycles")

    @property
    def is_running(self) -> bool:
        return self._running

    # ── save ────────────────────────────────────────────────────────

    def save(self, filepath, fmt='npz'):
        """Save camera timeseries data (contrast + intensity traces)."""
        if not self.full_contrasts:
            print("⚠ No data to save")
            return

        base = Path(filepath)
        all_contrast = np.concatenate(self.full_contrasts)
        all_sig = np.concatenate(self.full_sig) if self.full_sig else np.array([])
        all_ref = np.concatenate(self.full_ref) if self.full_ref else np.array([])

        # assert self.cfg.mw
        meta = {
            'sequence': self.sequence,
            'scan_name': self.scan_name,
            'fixed_value': self.fixed_value,
            'instr_mode': self.instr_mode,
            'exposure_s': self.exposure,
            'roi': self.roi,
            'frames_per_cycle': self.Nframes,
            'elapsed_s': self.elapsed,
            'contrast_rate': self.contrast_rate,
            'n_cycles': self.cycle_count,
            'mw_freq': self.cfg.mw.freq,
            'mw_power': self.cfg.mw.power,
        }

        p = base.with_suffix('.npz')
        np.savez_compressed(str(p),
                            contrast=all_contrast,
                            mean_sig=all_sig,
                            mean_ref=all_ref,
                            **{f'meta_{k}': v for k, v in meta.items()})
        print(f"💾 CamTS → {p} ({len(all_contrast)} pts, "
              f"{self.cycle_count} cycles)")

        cfg_path = base.with_suffix('.yaml')
        savefile(cfg_path, self.cfg.to_dict())

    # ── teardown ───────────────────────────────────────────────────

    def teardown(self):
        if self._running:
            self.stop()
        try:
            self.camera_worker.stop_capture(free_buffer=True)
        except Exception as e:
            print(f"⚠ CamTS teardown: {e}")


class _CameraPollingThread(threading.Thread):
    """Background thread that reads PB-triggered camera frames continuously.

    PB generates external triggers → camera captures frames → this thread
    reads them via wait_capevent_frameready + get_lastframedata.

    Each cycle: reads Nframes (signal + reference), computes ROI means
    and contrast, updates the experiment's ring buffers.
    """

    def __init__(self, exp: CameraTimeseriesExperiment):
        super().__init__(daemon=True)
        self.exp = exp

    def run(self):
        exp = self.exp
        cw = exp.camera_worker
        timeout_ms = int(exp.exposure * 1e3 + 500)
        Nframes = exp.Nframes

        # assert exp.contrast_ring and exp.sig_ring and exp.ref_ring
        try:
            while not exp._stop_requested:
                # Check auto-stop
                if (exp._auto_stop_time is not None and
                        time.perf_counter() > exp._auto_stop_time):
                    break

                # Read one PB cycle of externally-triggered frames
                frames = []
                cycle_ok = True
                for f_idx in range(Nframes):
                    res = cw.hdcamcon.wait_capevent_frameready(timeout_ms)
                    if res is True:
                        frame = cw.hdcamcon.get_lastframedata()
                        if frame is not False:
                            frames.append(frame)
                        else:
                            cycle_ok = False
                            break
                    else:
                        # Timeout — PB might not be running yet, or cycle gap
                        cycle_ok = False
                        break

                if not cycle_ok or len(frames) < Nframes:
                    continue  # incomplete cycle, wait for next

                # Stop PB sequence after reading last frame of cycle
                # (PB will restart on next loop iteration via BRANCH)
                # Actually: PB runs continuously, so no stop needed.
                # Just process the frames.

                # Signal = even-indexed frames, Reference = odd-indexed
                sig_frames = frames[0::2]
                ref_frames = frames[1::2]

                # ROI pixel mean per frame, then average across repeats
                mean_sig = float(np.mean([np.mean(f.astype(np.float64))
                                          for f in sig_frames]))
                mean_ref = float(np.mean([np.mean(f.astype(np.float64))
                                          for f in ref_frames]))
                contrast = mean_sig / mean_ref if mean_ref > 0 else np.nan

                # Update ring buffers (single-writer from this thread)
                c_arr = np.array([contrast])
                s_arr = np.array([mean_sig])
                r_arr = np.array([mean_ref])

                exp.contrast_ring.append(c_arr)
                exp.sig_ring.append(s_arr)
                exp.ref_ring.append(r_arr)

                # Full data for save
                exp.full_contrasts.append(c_arr)
                exp.full_sig.append(s_arr)
                exp.full_ref.append(r_arr)

                # Latest signal frame for image display
                exp.last_frame = frames[0]
                exp.cycle_count += 1

        except Exception as e:
            logging.exception(f"CamTS polling error: {e}")

