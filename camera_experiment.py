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

import numpy as np
import time
import threading
import logging
from pathlib import Path
from typing import Optional, Callable, TYPE_CHECKING

from PBcontrol import PulseBlaster
from sequencecontrol import sequencecontrol
from sweep_utils import Sweep, make_pb_setter
from experiment_base import calc_contrast, savefile, HW_BOUNDS

if TYPE_CHECKING:
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
        self.instr_mode: str = cam_cfg.instr_mode
        self.exposure: float = cam_cfg.exposure_s            # seconds
        self.roi: list = cam_cfg.roi
        self.Nframes: int = cam_cfg.frames_per_cycle
        self.focus_interval: int = cam_cfg.focus_interval
        self.focus_time: float = cam_cfg.min_focus_time_s    # seconds

        # Field config
        field_cfg = config.field_cfg
        self.align_field = field_cfg.align_field
        self.test_field = field_cfg.test_field
        self.t_align_dc_ms = field_cfg.t_align_dc_ms

        # Sequence
        self.sequence: str = config.seq.name
        self.Nsamples_cfg: int = config.seq.Nsamples  # = frames_per_cycle
        self.Nruns: int = config.runtime.Nruns
        self.pb_channels: dict = config.pb.channels

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

    def _build_sweep(self, cfg):
        """Add inner + outer axes to self.sweep."""
        # Inner axis — for camera experiments, sweep only calls setter for freq;
        # PB reprogramming is done inside the PBThread, not via Sweep setter
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

        else:
            raise ValueError(f"Unknown camera mode: {self.instr_mode}")

    def _fire_callback(self, callback, actual_i, actual_val,
                       i_run, oc_idx, oc_vals):
        """Process current run's data and fire per-point callback."""
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

        oc_str = ("  ".join(f"{k}={v:.4g}" for k, v in oc_vals.items())
                  if oc_vals else "")
        callback(actual_i, actual_val, i_run, oc_idx, oc_str,
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
                        callback(-1, 0.0, i_run, self.i_outer, "",
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
