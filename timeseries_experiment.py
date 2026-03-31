# timeseries_experiment.py
"""
TimeseriesExperiment — continuous fluorescence / contrast monitoring
at a fixed parameter point.

Shares instrument setup with DiodeExperiment but has its own streaming
acquisition loop (CONTINUOUS DAQ mode, polling backend).

Data model: 1D streaming timeseries of contrast values (and raw voltages).
Display: scrolling strip chart (ring buffer → pyqtgraph 30 Hz refresh).
Save: growing list → HDF5 or npz on stop.

Usage
-----
From session GUI (Timeseries mode):
    exp = TimeseriesExperiment(instruments, config)
    exp.setup()       # PB programmed once at fixed param, SG set, DAQ ready
    exp.start()       # begins continuous acquisition + callbacks
    ...               # GUI receives per-chunk callbacks
    exp.stop()        # stops DAQ, saves data
    exp.teardown()

From a script (standalone):
    exp = TimeseriesExperiment(instruments, config)
    exp.setup()
    exp.start(max_duration=60.0)   # auto-stop after 60s
    exp.wait()                     # blocks until done
    exp.save('ts_data')
    exp.teardown()
"""

import numpy as np
import time
import logging
from pathlib import Path
from typing import Optional, Callable, TYPE_CHECKING
from dataclasses import dataclass

from PySide6.QtCore import QObject, Signal, QThread, QMutex, QMutexLocker, QTimer

import connectionConfig as concfg
from PBcontrol import PulseBlaster
from experiment_base import calc_contrast, read_details, savefile
from ring_buffer import RingBuffer

if TYPE_CHECKING:
    from experiment_config import ExperimentConfig


# ═══════════════════════════════════════════════════════════════════════
#  Contrast Processor (from ts_prototype, adapted)
# ═══════════════════════════════════════════════════════════════════════

class ContrastProcessor:
    """Splits raw DAQ chunks into signal/reference and computes contrast.

    Parameters
    ----------
    reads_per_cycle : int
        Gated reads per PB cycle (typically 2: sig + ref).
    samples_per_read : int
        DAQ samples per read window (= Nsamples from config).
    contrast_op : str
        Operation: 's/r', 'r/s', '+-', '-+'.
    """

    OPS = {
        's/r': lambda s, r: np.where(r != 0, s / r, np.nan),
        'r/s': lambda s, r: np.where(s != 0, r / s, np.nan),
        '+-':  lambda s, r: s - r,
        '-+':  lambda s, r: r - s,
    }

    def __init__(self, reads_per_cycle: int = 2,
                 samples_per_read: int = 50,
                 contrast_op: str = 's/r'):
        self.reads_per_cycle = reads_per_cycle
        self.samples_per_read = samples_per_read
        self.contrast_op = contrast_op
        self._op_fn = self.OPS[contrast_op]

    @property
    def chunk_samples(self) -> int:
        """Total DAQ samples per cycle."""
        return self.reads_per_cycle * self.samples_per_read

    def process_multi(self, raw: np.ndarray):
        """Process a chunk containing one or more complete cycles.

        Returns
        -------
        contrasts : ndarray, shape (n_cycles,)
        means : ndarray, shape (n_cycles, reads_per_cycle)
        """
        cs = self.chunk_samples
        n_cycles = len(raw) // cs
        if n_cycles == 0:
            return np.array([]), np.array([]).reshape(0, self.reads_per_cycle)

        usable = raw[:n_cycles * cs].reshape(
            n_cycles, self.reads_per_cycle, self.samples_per_read)
        means = usable.mean(axis=2)  # (n_cycles, reads_per_cycle)
        contrasts = self._op_fn(means[:, 0], means[:, 1])
        return contrasts, means


# ═══════════════════════════════════════════════════════════════════════
#  DAQ Worker — signal bridge
# ═══════════════════════════════════════════════════════════════════════

class TSWorker(QObject):
    """Qt signal bridge for timeseries DAQ → main thread."""
    chunk_ready = Signal(object)      # np.ndarray (raw DAQ samples)
    error = Signal(str)
    status = Signal(str)
    finished = Signal()               # acquisition ended (timeout or manual)


# ═══════════════════════════════════════════════════════════════════════
#  Polling DAQ Thread (from ts_prototype, adapted for PB sample clock)
# ═══════════════════════════════════════════════════════════════════════

class TSPollingThread(QThread):
    """Continuous DAQ acquisition using blocking task.read() in a loop.

    Configured for external sample clock (PB) + start trigger.
    The PB sequence runs in a loop (BRANCH), continuously generating
    sample clock pulses; each cycle produces reads_per_cycle × Nsamples
    DAQ samples.
    """

    def __init__(self, worker: TSWorker, channel: str,
                 sample_rate: float, chunk_size: int,
                 buffer_size: int,
                 voltage_range: tuple[float, float] = (-10, 10),
                 sampling_source: str = '',
                 start_trigger_source: str = '',
                 parent=None):
        super().__init__(parent)
        self.worker = worker
        self.channel = channel
        self.sample_rate = sample_rate
        self.chunk_size = chunk_size
        self.buffer_size = buffer_size
        self.voltage_range = voltage_range
        self.sampling_source = sampling_source
        self.start_trigger_source = start_trigger_source
        self._stop_requested = False
        self._mutex = QMutex()

    def request_stop(self):
        with QMutexLocker(self._mutex):
            self._stop_requested = True

    def run(self):
        import nidaqmx
        from nidaqmx.constants import AcquisitionType

        chunk_duration = self.chunk_size / max(self.sample_rate, 1)
        timeout = max(chunk_duration * 5, 5.0)

        self.worker.status.emit(
            f'TS Polling | {self.sample_rate/1e3:.1f} kSa/s '
            f'| chunk={self.chunk_size} | buf={self.buffer_size}')

        task = None
        try:
            task = nidaqmx.Task('ts_experiment_ai')
            task.ai_channels.add_ai_voltage_chan(
                self.channel,
                min_val=self.voltage_range[0],
                max_val=self.voltage_range[1],
            )
            task.timing.cfg_samp_clk_timing(
                rate=self.sample_rate,
                sample_mode=AcquisitionType.CONTINUOUS,
                samps_per_chan=self.buffer_size,
            )

            # External sample clock from PB (if configured)
            if self.sampling_source:
                task.timing.samp_clk_src = self.sampling_source

            # Start trigger from PB
            if self.start_trigger_source:
                task.triggers.start_trigger.cfg_dig_edge_start_trig(
                    self.start_trigger_source)

            task.start()

            while True:
                with QMutexLocker(self._mutex):
                    if self._stop_requested:
                        break
                try:
                    data = task.read(
                        number_of_samples_per_channel=self.chunk_size,
                        timeout=timeout,
                    )
                    self.worker.chunk_ready.emit(
                        np.array(data, dtype=np.float64))
                except Exception as e:
                    err_str = str(e).lower()
                    if 'timeout' in err_str:
                        self.worker.status.emit('Waiting for PB triggers...')
                        continue
                    if 'not yet been acquired' in err_str:
                        time.sleep(0.001)
                        continue
                    raise

        except Exception as e:
            self.worker.error.emit(f'TS polling: {e}')
        finally:
            if task:
                try: task.stop()
                except: pass
                try: task.close()
                except: pass
            self.worker.finished.emit()

    def start_acquisition(self):
        self._stop_requested = False
        self.start()  # QThread.start()

    def stop(self):
        self.request_stop()
        self.wait(5000)


# ═══════════════════════════════════════════════════════════════════════
#  Counter Polling Thread (CONTINUOUS mode, for SPD photon counting)
# ═══════════════════════════════════════════════════════════════════════

class CounterPollingThread(QThread):
    """Continuous counter acquisition using blocking task.read() in a loop.

    Mirrors TSPollingThread but uses a CI count-edges channel instead of
    an AI voltage channel.  Returns counts-per-bin (via np.diff) so the
    downstream data pipeline sees the same 1D numeric array shape as the
    analog path.

    The NI counter accumulates edges; each read returns monotonically
    increasing totals.  np.diff converts to per-bin counts.  Overflow
    (32-bit wrap) is handled by int64 cast + correction.
    """

    def __init__(self, worker: TSWorker,
                 counter: str,
                 input_terminal: str,
                 sample_rate: float,
                 chunk_size: int,
                 buffer_size: int,
                 sampling_source: str = '',
                 start_trigger_source: str = '',
                 start_trigger_edge: str = 'rising',
                 parent=None):
        super().__init__(parent)
        self.worker = worker
        self.counter = counter              # e.g. 'P6363/ctr0'
        self.input_terminal = input_terminal  # e.g. '/P6363/PFI2'
        self.sample_rate = sample_rate
        self.chunk_size = chunk_size
        self.buffer_size = buffer_size
        self.sampling_source = sampling_source
        self.start_trigger_source = start_trigger_source
        self.start_trigger_edge = start_trigger_edge
        self._stop_requested = False
        self._mutex = QMutex()

    def request_stop(self):
        with QMutexLocker(self._mutex):
            self._stop_requested = True

    def run(self):
        import nidaqmx
        from nidaqmx.constants import AcquisitionType, Edge

        chunk_duration = self.chunk_size / max(self.sample_rate, 1)
        timeout = max(chunk_duration * 5, 5.0)

        self.worker.status.emit(
            f'Counter Polling | {self.sample_rate/1e3:.1f} kSa/s '
            f'| chunk={self.chunk_size} | buf={self.buffer_size}')

        task = None
        try:
            task = nidaqmx.Task('ts_experiment_ci')

            # Counter channel
            ci_chan = task.ci_channels.add_ci_count_edges_chan(
                counter=self.counter,
                edge=Edge.RISING,
            )
            ci_chan.ci_count_edges_term = self.input_terminal

            # CONTINUOUS sample clock timing
            task.timing.cfg_samp_clk_timing(
                rate=self.sample_rate,
                source=self.sampling_source if self.sampling_source else '',
                sample_mode=AcquisitionType.CONTINUOUS,
                samps_per_chan=self.buffer_size,
            )

            # Arm-start trigger from PB (CI uses arm_start, not start_trigger)
            if self.start_trigger_source:
                edge_val = (Edge.RISING if self.start_trigger_edge == 'rising'
                            else Edge.FALLING)
                task.triggers.arm_start_trigger.dig_edge_src = (
                    self.start_trigger_source)
                task.triggers.arm_start_trigger.dig_edge_edge = edge_val

            task.start()

            while True:
                with QMutexLocker(self._mutex):
                    if self._stop_requested:
                        break
                try:
                    raw = task.read(
                        number_of_samples_per_channel=self.chunk_size,
                        timeout=timeout,
                    )
                    # Convert accumulated counts → counts per bin
                    raw_i64 = np.array(raw, dtype=np.int64)
                    counts = np.diff(raw_i64, prepend=raw_i64[0])
                    # Handle 32-bit counter overflow (wrap-around)
                    counts[counts < 0] += 2**32
                    self.worker.chunk_ready.emit(
                        counts.astype(np.float64))
                except Exception as e:
                    err_str = str(e).lower()
                    if 'timeout' in err_str:
                        self.worker.status.emit('Waiting for PB triggers...')
                        continue
                    if 'not yet been acquired' in err_str:
                        time.sleep(0.001)
                        continue
                    raise

        except Exception as e:
            self.worker.error.emit(f'Counter polling: {e}')
        finally:
            if task:
                try: task.stop()
                except: pass
                try: task.close()
                except: pass
            self.worker.finished.emit()

    def start_acquisition(self):
        self._stop_requested = False
        self.start()  # QThread.start()

    def stop(self):
        self.request_stop()
        self.wait(5000)


# ═══════════════════════════════════════════════════════════════════════
#  Buffer safety helpers (from ts_prototype)
# ═══════════════════════════════════════════════════════════════════════

MAX_CHUNK_FRACTION = 0.7

def clamp_chunk_to_buffer(chunk_size: int, buffer_size: int) -> int:
    max_chunk = int(buffer_size * MAX_CHUNK_FRACTION)
    return min(chunk_size, max(1, max_chunk))

def auto_buffer_size(chunk_size: int, sample_rate: float) -> int:
    return max(chunk_size * 20, int(sample_rate * 2), 10000)


# ═══════════════════════════════════════════════════════════════════════
#  TimeseriesExperiment
# ═══════════════════════════════════════════════════════════════════════

class TimeseriesExperiment:
    """Continuous acquisition at a fixed parameter point.

    Programs PB once (at the fixed scan value), configures DAQ in
    CONTINUOUS mode with external sample clock from PB, and streams
    data through a polling thread.

    Contrast is computed from signal/reference reads within each
    PB cycle (same structure as sweep experiments, but instead of
    sweeping a parameter, we repeat the same point and watch
    contrast evolve over time).

    Attributes
    ----------
    contrast_ring : RingBuffer — last N seconds of contrast values
    raw_ring : RingBuffer — last N seconds of raw DAQ samples
    full_chunks : list[ndarray] — all raw data for saving
    contrast_proc : ContrastProcessor — signal/ref → contrast
    elapsed : float — seconds since start
    """

    # Freq-only sequences where PB doesn't depend on scan param
    FREQ_ONLY = frozenset([
        'esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq'])

    def __init__(self, instruments: dict, config: 'ExperimentConfig'):
        self.sg = instruments['sg']
        self.pb = instruments['pb']
        self.ao_task = instruments.get('ao_task')
        self.cfg = config
        self.instr = 'diode'

        # Config
        self.sequence = config.seq.name
        self.pb_channels = config.pb.channels
        self.Nsamples_cfg = config.seq.Nsamples
        self.n_channels = len(concfg.input_terminals)

        # Timeseries config (from RuntimeFlags)
        rt = config.runtime
        self.display_seconds = getattr(rt, 'ts_display_seconds', 30.0)
        self.max_duration = getattr(rt, 'ts_max_duration', 0.0)
        self.contrast_op = getattr(rt, 'ts_contrast_op', 's/r')

        # Reads per cycle from sequence type
        self.reads_per_cyc, self.daq_Nsamples = read_details(
            self.sequence, self.n_channels, self.Nsamples_cfg)

        # Contrast processor
        rpc = self.reads_per_cyc[0] if self.reads_per_cyc else 2
        self.contrast_proc = ContrastProcessor(
            reads_per_cycle=rpc,
            samples_per_read=self.Nsamples_cfg,
            contrast_op=self.contrast_op,
        )
        self.chunk_samples = self.contrast_proc.chunk_samples

        # Which scan point to fix at
        scan_names = config.scan_names
        self.scan_name = scan_names[0] if scan_names else ''
        scan = config.scans.get(self.scan_name)
        # Use the last seq_plot_indices value, or midpoint of scan
        indices = getattr(rt, 'seq_plot_indices', [0, -1])
        if scan is not None and len(scan.values) > 0:
            idx = indices[-1]
            if idx < 0:
                idx = len(scan.values) + idx
            self.fixed_value = scan.values[idx]
        else:
            self.fixed_value = 0.0

        # State
        self._worker: Optional[TSWorker] = None
        self._thread: Optional[TSPollingThread | CounterPollingThread] = None
        self._auto_stop_timer: Optional[QTimer] = None
        self.contrast_ring: Optional[RingBuffer] = None
        self.raw_ring: Optional[RingBuffer] = None
        self.full_chunks: list = []
        self.elapsed: float = 0.0
        self._start_time: float = 0.0
        self._running = False

    # ── setup ───────────────────────────────────────────────────────

    def setup(self):
        """Program PB at fixed point, configure SG, prepare buffers."""
        cfg = self.cfg

        # ── SG ───────────────────────────────────────────────────────
        if self.sequence not in ['aom_timing', 'T1ms0_train', 'drift_seq']:
            self.sg.set_freq(cfg.mw.freq)
            self.sg.set_amp_rf(cfg.mw.power)

        # For freq-only sequences, also set the fixed frequency
        is_freq = self.sequence in self.FREQ_ONLY
        if is_freq:
            self.sg.set_freq(self.fixed_value)

        # ── PB program at fixed point ────────────────────────────────
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
            self.instr, self.sequence, sal + [self.pb_channels])
        self.pb.run_sequence_for_diode(
            [the_list[i][0] for i in range(len(the_list))])

        print(f"▶ TS: PB programmed — {self.sequence} @ "
              f"{self.scan_name}={self.fixed_value:.4g}")

        # ── AO ───────────────────────────────────────────────────────
        tf = cfg.extra.get('test_field', [0, 0, 0])
        if self.ao_task is not None:
            self.ao_task.set_outputs_to_constant(output_field_in_gauss=tf)

        # ── Buffers ──────────────────────────────────────────────────
        # Contrast rate: one contrast point per PB cycle
        sample_rate = cfg.daq_ai.ai_sample_rate
        contrast_rate = sample_rate / self.chunk_samples
        contrast_cap = int(self.display_seconds * contrast_rate)
        raw_cap = int(self.display_seconds * sample_rate)

        self.contrast_ring = RingBuffer(max(contrast_cap, 100))
        self.raw_ring = RingBuffer(max(raw_cap, 1000))
        self.full_chunks = []

        # ── Buffer sizing ────────────────────────────────────────────
        buf_size = getattr(cfg.runtime, 'ts_buffer_size', 0)
        if buf_size <= 0:
            buf_size = auto_buffer_size(self.chunk_samples, sample_rate)
        self._buffer_size = buf_size

        print(f"▶ TS: chunk={self.chunk_samples} samples/cycle, "
              f"buffer={buf_size}, contrast_rate={contrast_rate:.1f} pts/s, "
              f"display={self.display_seconds}s")

    # ── start / stop ────────────────────────────────────────────────

    def start(self, callback: Optional[Callable] = None,
              max_duration: Optional[float] = None):
        """Begin continuous acquisition.

        Parameters
        ----------
        callback : callable(contrasts, means, raw, elapsed)
            Called in main thread for each chunk.
        max_duration : float, optional
            Override config's ts_max_duration. 0 = manual stop only.
        """
        if self._running:
            print("⚠ Already running")
            return

        cfg = self.cfg
        sample_rate = cfg.daq_ai.ai_sample_rate

        # Worker + thread
        self._worker = TSWorker()
        self._callback = callback

        # Connect signals
        self._worker.chunk_ready.connect(
            self._on_chunk, type=2)   # Qt.QueuedConnection = 2
        self._worker.error.connect(
            lambda msg: print(f"❌ TS: {msg}"), type=2)
        self._worker.finished.connect(self._on_finished, type=2)

        channel = f"P6363/ai{concfg.input_terminals[0]}"

        # Dispatch: analog (AI) vs counter (CI) polling thread
        if cfg.runtime.detector == 'counter':
            ci = cfg.daq_ci
            self._thread = CounterPollingThread(
                worker=self._worker,
                counter=ci.ci_counter,
                input_terminal=ci.ci_input_terminal,
                sample_rate=sample_rate,
                chunk_size=self.chunk_samples,
                buffer_size=self._buffer_size,
                sampling_source=concfg.samp_clk_terminal,
                start_trigger_source=concfg.start_trig_terminal,
            )
            print(f"▶ TS: detector=COUNTER ({ci.ci_counter})")
        else:
            self._thread = TSPollingThread(
                worker=self._worker,
                channel=channel,
                sample_rate=sample_rate,
                chunk_size=self.chunk_samples,
                buffer_size=self._buffer_size,
                sampling_source=concfg.samp_clk_terminal,
                start_trigger_source=concfg.start_trig_terminal,
            )
            print(f"▶ TS: detector=ANALOG ({channel})")

        self._start_time = time.perf_counter()
        self._running = True
        self._thread.start_acquisition()

        # Auto-stop timer
        dur = max_duration if max_duration is not None else self.max_duration
        if dur > 0:
            self._auto_stop_timer = QTimer()
            self._auto_stop_timer.setSingleShot(True)
            self._auto_stop_timer.timeout.connect(self.stop)
            self._auto_stop_timer.start(int(dur * 1000))
            print(f"▶ TS: auto-stop in {dur:.1f}s")

        print("▶ TS: acquisition started")

    def stop(self):
        """Stop continuous acquisition."""
        if not self._running:
            return

        if self._auto_stop_timer:
            self._auto_stop_timer.stop()
            self._auto_stop_timer = None

        if self._thread:
            self._thread.stop()
            self._thread = None

        self._running = False
        self.elapsed = time.perf_counter() - self._start_time

        total = sum(len(c) for c in self.full_chunks)
        n_contrast = self.contrast_ring.total_count if self.contrast_ring else 0
        print(f"⏹ TS: stopped — {self.elapsed:.1f}s, "
              f"{total:,} raw samples, {n_contrast} contrast pts")

    def wait(self, timeout_ms: int = -1):
        """Block until acquisition finishes (for scripted use)."""
        if self._thread and self._thread.isRunning():
            self._thread.wait(timeout_ms)

    @property
    def is_running(self) -> bool:
        return self._running

    # ── data handlers ──────────────────────────────────────────────

    def _on_chunk(self, raw: np.ndarray):
        """Called in main thread for each DAQ chunk."""
        self.full_chunks.append(raw)
        self.raw_ring.append(raw)

        # Compute contrast
        contrasts, means = self.contrast_proc.process_multi(raw)
        if len(contrasts) > 0:
            self.contrast_ring.append(contrasts)

        # User callback
        if self._callback is not None:
            elapsed = time.perf_counter() - self._start_time
            self._callback(contrasts, means, raw, elapsed)

    def _on_finished(self):
        """Called when the polling thread exits."""
        if self._running:
            self.stop()

    # ── save ────────────────────────────────────────────────────────

    def save(self, filepath, fmt: str = 'npz'):
        """Save timeseries data.

        Parameters
        ----------
        filepath : str or Path
            Base path without extension. Extensions added automatically.
        fmt : str
            'npz' or 'h5' / 'hdf5'.
        """
        if not self.full_chunks:
            print("⚠ No data to save")
            return

        base = Path(filepath)
        all_raw = np.concatenate(self.full_chunks)

        # Also compute all contrast from raw
        all_contrast, all_means = self.contrast_proc.process_multi(all_raw)

        meta = {
            'sequence': self.sequence,
            'scan_name': self.scan_name,
            'fixed_value': self.fixed_value,
            'detector': self.cfg.runtime.detector,
            'sample_rate': self.cfg.daq_ai.ai_sample_rate,
            'Nsamples': self.Nsamples_cfg,
            'reads_per_cycle': self.contrast_proc.reads_per_cycle,
            'samples_per_read': self.contrast_proc.samples_per_read,
            'contrast_op': self.contrast_op,
            'elapsed_s': self.elapsed,
            'mw_freq': self.cfg.mw.freq,
            'mw_power': self.cfg.mw.power,
        }

        if fmt in ('h5', 'hdf5'):
            import h5py
            p = base.with_suffix('.h5')
            with h5py.File(p, 'w') as f:
                f.create_dataset('raw', data=all_raw)
                f.create_dataset('contrast', data=all_contrast)
                f.create_dataset('means', data=all_means)
                for k, v in meta.items():
                    f.attrs[k] = v
            print(f"💾 TS saved → {p} ({len(all_raw):,} samples)")
        else:
            p = base.with_suffix('.npz')
            np.savez_compressed(p,
                                raw=all_raw,
                                contrast=all_contrast,
                                means=all_means,
                                **{f'meta_{k}': v for k, v in meta.items()})
            print(f"💾 TS saved → {p} ({len(all_raw):,} samples)")

        # Also save config as YAML
        cfg_path = base.with_suffix('.yaml')
        savefile(cfg_path, self.cfg.to_dict())

    # ── teardown ───────────────────────────────────────────────────

    def teardown(self):
        """Stop acquisition if running, cleanup."""
        if self._running:
            self.stop()
        self._worker = None
