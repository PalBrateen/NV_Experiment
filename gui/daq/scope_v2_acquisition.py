"""
Scope Viewer v2 - Acquisition Module (Refactored)
=================================================

Architecture:
- Acquisition thread reads chunks of size = fft_size (continuous data)
- Each chunk is emitted directly for immediate Scope/FFT display
- Chunks also written to CircularBuffer for trends history
- Display throttled to max 30 Hz for fast sample rates

Author: BP Lab
Date: 2025
"""

import numpy as np, time, logging, nidaqmx, traceback, nidaqmx.errors as daq_error
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass
from nidaqmx.constants import AcquisitionType, TerminalConfiguration, Edge, Level, TriggerType

from PySide6.QtCore import QThread, Signal, QMutex, QMutexLocker, QTimer, QObject

from scope_v2_core import (
    ObservableConfig, AcquisitionConfig, ChannelConfig,
    CircularBuffer, HAS_NIDAQMX,
    TriggerConfig, ClockConfig, TriggerMode, TriggerEdge, PauseWhen, ClockSource
)


# =============================================================================
# ACQUISITION THREAD
# =============================================================================

class AcquisitionThread(QThread):
    """
    Background acquisition thread.
    
    Reads chunks of exactly fft_size samples (continuous, no gaps).
    Emits each chunk directly for immediate display.
    Also writes to CircularBuffer for trends history.
    """
    
    errorOccurred = Signal(str)
    statusChanged = Signal(str)
    
    # New: emit complete chunks directly for display
    # Dict[int, np.ndarray] = {channel_idx: data_array}
    chunkReady = Signal(object)
    
    def __init__(self, obs_config: ObservableConfig, buffers: Dict[int, CircularBuffer], 
                 chunk_size: int = 8192, parent=None):
        super().__init__(parent)
        self.obs_config = obs_config
        self.buffers = buffers  # For trends history
        self.chunk_size = chunk_size  # = fft_size
        
        self._stop_requested = False
        self._restart_requested = False
        self._mutex = QMutex()
        
        self._total_samples = 0
        self.input_type = 'analog'  # 'analog' or 'counter'
        
    def set_chunk_size(self, size: int):
        """Update chunk size (call before starting or will take effect on restart)"""
        self.chunk_size = size
    
    def set_input_type(self, input_type: str):
        """Set input type: 'analog' or 'counter'"""
        self.input_type = input_type
        
    def request_stop(self):
        with QMutexLocker(self._mutex):
            self._stop_requested = True
            
    def request_restart(self):
        with QMutexLocker(self._mutex):
            self._restart_requested = True
            
    def run(self):
        self.statusChanged.emit("Starting...")
        self._stop_requested = False
        self._total_samples = 0
        
        while not self._stop_requested:
            self._restart_requested = False
            
            if HAS_NIDAQMX:
                if self.input_type == 'counter':
                    self._run_hardware_counter()
                else:
                    self._run_hardware()
            else:
                self._run_simulation()
            
            if not self._restart_requested:
                break
                
        self.statusChanged.emit("Stopped")
        
    def _run_hardware(self):
        """Hardware acquisition - reads chunks of fft_size samples"""
        config = self.obs_config.get_config()
        sample_rate = config.sample_rate
        chunk_size = self.chunk_size
        
        # Calculate timing
        chunk_duration = chunk_size / sample_rate
        natural_rate = sample_rate / chunk_size
        timeout = max(chunk_duration * 2, 1.0)  # At least 1s or 2x chunk duration
        
        self.statusChanged.emit(f"Configuring... (chunk={chunk_size}, rate={natural_rate:.1f} Hz)")
        
        task = None
        try:
            task = nidaqmx.Task()
            
            # Add channels
            enabled_channels = []
            for i, ch_config in enumerate(config.channels):
                if not ch_config.enabled:
                    continue
                    
                phys_chan = f"{config.device}/{ch_config.physical_channel}"
                
                term_map = {
                    "RSE": TerminalConfiguration.RSE,
                    "NRSE": TerminalConfiguration.NRSE,
                    "DIFF": TerminalConfiguration.DIFF,
                    "PSEUDODIFF": TerminalConfiguration.PSEUDO_DIFF
                }
                term_cfg = term_map.get(ch_config.terminal_config, TerminalConfiguration.RSE)
                
                task.ai_channels.add_ai_voltage_chan(
                    physical_channel=phys_chan,
                    terminal_config=term_cfg,
                    min_val=ch_config.min_voltage,
                    max_val=ch_config.max_voltage
                )
                enabled_channels.append(i)
            
            if not enabled_channels:
                self.errorOccurred.emit("No enabled channels")
                return
            
            # Configure sample clock
            clock_cfg = config.clock
            trigger_cfg = config.trigger
            
            if clock_cfg.source == ClockSource.INTERNAL:
                clock_source = ""
                clock_rate = clock_cfg.rate
                clock_edge = Edge.RISING
            else:
                clock_source = f"/{config.device}/{clock_cfg.external_source}"
                clock_rate = clock_cfg.rate  # Expected rate for buffer sizing
                clock_edge = Edge.RISING if clock_cfg.external_edge == TriggerEdge.RISING else Edge.FALLING
            
            # Configure timing - DAQ buffer = 4x chunk size
            task.timing.cfg_samp_clk_timing(
                rate=clock_rate,
                source=clock_source,
                active_edge=clock_edge,
                sample_mode=AcquisitionType.CONTINUOUS,
                samps_per_chan=chunk_size * 4
            )
            
            # Configure triggers
            if trigger_cfg.mode != TriggerMode.NONE:
                self._configure_triggers(task, config.device, trigger_cfg)
            
            task.start()
            self.statusChanged.emit(f"Running | {sample_rate/1e3:.1f} kSa/s | {natural_rate:.1f} Hz update")
            
            # Main acquisition loop
            while not self._stop_requested and not self._restart_requested:
                try:
                    # Read exactly chunk_size samples (blocks until ready)
                    data = task.read(
                        number_of_samples_per_channel=chunk_size,
                        timeout=timeout
                    )
                    
                    data = np.array(data, dtype=np.float32)
                    
                    # Build chunk dict for emission
                    chunk_data = {}
                    
                    if data.ndim == 1:
                        # Single channel
                        ch_idx = enabled_channels[0]
                        chunk_data[ch_idx] = data
                        if ch_idx in self.buffers:
                            self.buffers[ch_idx].write(data)
                    else:
                        # Multiple channels
                        for buf_idx, ch_idx in enumerate(enabled_channels):
                            chunk_data[ch_idx] = data[buf_idx]
                            if ch_idx in self.buffers:
                                self.buffers[ch_idx].write(data[buf_idx])
                    
                    # Emit chunk for immediate display
                    self.chunkReady.emit(chunk_data)
                    
                    self._total_samples += chunk_size
                    
                except daq_error.DaqReadError as e:
                    err_str = str(e).lower()
                    if "timeout" in err_str:
                        self.statusChanged.emit("Waiting for data/trigger...")
                        continue
                    if "not yet been acquired" in err_str:
                        time.sleep(0.001)
                        continue
                    raise
                    
        except Exception as e:
            traceback.print_exc()
            self.errorOccurred.emit(str(e))
            
        finally:
            if task:
                try:
                    task.stop()
                    task.close()
                except:
                    pass
    
    def _run_hardware_counter(self):
        """Counter input acquisition - reads edge counts and converts to counts/second"""
        config = self.obs_config.get_config()
        sample_rate = config.sample_rate
        chunk_size = self.chunk_size
        
        # Calculate timing
        chunk_duration = chunk_size / sample_rate
        bin_time = 1.0 / sample_rate  # Time per sample bin
        natural_rate = sample_rate / chunk_size
        timeout = max(chunk_duration * 2, 1.0)
        
        self.statusChanged.emit(f"Configuring counter... (chunk={chunk_size}, rate={natural_rate:.1f} Hz)")
        
        task = None
        try:
            task = nidaqmx.Task()
            
            # Get counter channel from first enabled channel config
            # Counter channels use format "Dev/ctr0", "Dev/ctr1" etc.
            enabled_channels = []
            for i, ch_config in enumerate(config.channels):
                if not ch_config.enabled:
                    continue
                # For counter, physical_channel should be "ctr0", "ctr1", etc.
                counter_chan = f"{config.device}/{ch_config.physical_channel}"
                
                # Add counter edge counting channel
                task.ci_channels.add_ci_count_edges_chan(
                    counter=counter_chan,
                    edge=nidaqmx.constants.Edge.RISING
                )
                enabled_channels.append(i)
            
            if not enabled_channels:
                self.errorOccurred.emit("No enabled counter channels")
                return
            
            # Configure sample clock timing
            clock_cfg = config.clock
            trigger_cfg = config.trigger
            
            if clock_cfg.source == ClockSource.INTERNAL:
                clock_source = ""  # Use internal timebase (100MHz on X series)
                clock_rate = clock_cfg.rate
            else:
                clock_source = f"/{config.device}/{clock_cfg.external_source}"
                clock_rate = clock_cfg.rate
            
            # Configure timing for buffered counter acquisition
            task.timing.cfg_samp_clk_timing(
                rate=clock_rate,
                source=clock_source if clock_source else None,
                sample_mode=AcquisitionType.CONTINUOUS,
                samps_per_chan=chunk_size * 4
            )
            
            # Configure triggers if needed
            if trigger_cfg.mode != TriggerMode.NONE:
                self._configure_triggers(task, config.device, trigger_cfg)
            
            task.start()
            self.statusChanged.emit(f"Counter | {sample_rate/1e3:.1f} kSa/s | {natural_rate:.1f} Hz update")
            
            # Track last raw count for diff calculation
            last_counts = None
            
            # Main acquisition loop
            while not self._stop_requested and not self._restart_requested:
                try:
                    # Read accumulated counts
                    raw_counts = task.read(
                        number_of_samples_per_channel=chunk_size,
                        timeout=timeout
                    )
                    
                    raw_counts = np.array(raw_counts, dtype=np.int64)
                    
                    # Convert to counts per bin using diff
                    if raw_counts.ndim == 1:
                        # Single counter
                        if last_counts is None:
                            counts_per_bin = np.diff(raw_counts, prepend=raw_counts[0])
                        else:
                            counts_per_bin = np.diff(raw_counts, prepend=last_counts)
                        last_counts = raw_counts[-1]
                        
                        # Handle 32-bit overflow (unlikely but safe)
                        counts_per_bin[counts_per_bin < 0] += 2**32
                        
                        # Normalize to counts per second
                        counts_per_second = counts_per_bin.astype(np.float32) / bin_time
                        
                        ch_idx = enabled_channels[0]
                        chunk_data = {ch_idx: counts_per_second}
                        if ch_idx in self.buffers:
                            self.buffers[ch_idx].write(counts_per_second)
                    else:
                        # Multiple counters
                        chunk_data = {}
                        for buf_idx, ch_idx in enumerate(enabled_channels):
                            raw_ch = raw_counts[buf_idx]
                            if last_counts is None:
                                counts_per_bin = np.diff(raw_ch, prepend=raw_ch[0])
                            else:
                                counts_per_bin = np.diff(raw_ch, prepend=last_counts[buf_idx] if last_counts is not None else raw_ch[0])
                            
                            counts_per_bin[counts_per_bin < 0] += 2**32
                            counts_per_second = counts_per_bin.astype(np.float32) / bin_time
                            
                            chunk_data[ch_idx] = counts_per_second
                            if ch_idx in self.buffers:
                                self.buffers[ch_idx].write(counts_per_second)
                        
                        if raw_counts.ndim > 1:
                            last_counts = raw_counts[:, -1]
                        else:
                            last_counts = raw_counts[-1]
                    
                    # Emit chunk for display
                    self.chunkReady.emit(chunk_data)
                    self._total_samples += chunk_size
                    
                except daq_error.DaqReadError as e:
                    err_str = str(e).lower()
                    if "timeout" in err_str:
                        continue
                    elif self._stop_requested or self._restart_requested:
                        break
                    else:
                        self.errorOccurred.emit(f"Counter read error: {str(e)}")
                        break
                    
        except Exception as e:
            traceback.print_exc()
            self.errorOccurred.emit(str(e))
            
        finally:
            if task:
                try:
                    task.stop()
                    task.close()
                except:
                    pass
    
    def _configure_triggers(self, task, device: str, trigger_cfg: TriggerConfig):
        """Configure start and pause triggers"""
        
        # Start trigger
        if trigger_cfg.mode in [TriggerMode.START, TriggerMode.START_PAUSE]:
            start_source = f"/{device}/{trigger_cfg.start_source}"
            start_edge = Edge.RISING if trigger_cfg.start_edge == TriggerEdge.RISING else Edge.FALLING
            
            task.triggers.start_trigger.cfg_dig_edge_start_trig(
                trigger_source=start_source,
                trigger_edge=start_edge
            )
            
            if trigger_cfg.retriggerable:
                task.triggers.start_trigger.retriggerable = True
        
        # Pause trigger
        if trigger_cfg.mode == TriggerMode.START_PAUSE:
            pause_source = f"/{device}/{trigger_cfg.pause_source}"
            pause_level = Level.HIGH if trigger_cfg.pause_when == PauseWhen.HIGH else Level.LOW
            
            task.triggers.pause_trigger.trig_type = TriggerType.DIGITAL_LEVEL
            task.triggers.pause_trigger.dig_lvl_src = pause_source
            task.triggers.pause_trigger.dig_lvl_when = pause_level
                    
    def _run_simulation(self):
        """Simulated acquisition - same chunk-based approach"""
        config = self.obs_config.get_config()
        sample_rate = config.sample_rate
        chunk_size = self.chunk_size
        
        chunk_duration = chunk_size / sample_rate
        natural_rate = sample_rate / chunk_size
        
        self.statusChanged.emit(f"Simulation | {sample_rate/1e3:.1f} kSa/s | {natural_rate:.1f} Hz update")
        
        dt = 1.0 / sample_rate
        sim_time = 0.0
        
        # Find enabled channels
        enabled_channels = [i for i, ch in enumerate(config.channels) if ch.enabled]
        if not enabled_channels:
            enabled_channels = [0]
            
        while not self._stop_requested and not self._restart_requested:
            t = sim_time + np.arange(chunk_size, dtype=np.float32) * dt
            sim_time += chunk_duration
            
            chunk_data = {}
            
            for ch_idx in enabled_channels:
                freq = 1000 * (ch_idx + 1)  # 1kHz, 2kHz, etc.
                signal = (
                    0.5 * np.sin(2 * np.pi * freq * t) +
                    0.1 * np.sin(2 * np.pi * 50 * t) +  # 50 Hz noise
                    0.02 * np.random.randn(len(t))
                ).astype(np.float32)
                
                chunk_data[ch_idx] = signal
                
                if ch_idx in self.buffers:
                    self.buffers[ch_idx].write(signal)
            
            # Emit chunk
            self.chunkReady.emit(chunk_data)
            
            self._total_samples += chunk_size
            
            # Sleep to simulate real-time acquisition
            time.sleep(chunk_duration * 0.95)


# =============================================================================
# DISPLAY CONFIG
# =============================================================================

@dataclass
class DisplayConfig:
    """Configuration for display/processing"""
    scope_points: int = 10000       # Max points to show on scope (decimation target)
    fft_size: int = 8192            # FFT size = chunk size
    fft_window: str = "Hann"        # Window function
    max_display_rate: float = 30.0  # Maximum display update rate (Hz)
    trends_interval: float = 0.2    # Trends update interval (seconds)
    power_correction: bool = False  # True = Power correction, False = Amplitude correction


# Enum for spectrum modes (matches FFTDisplayMode in widgets)
class SpectrumMode:
    NORMAL = "NORMAL"              # V
    SPECTRAL_DENSITY = "SPECTRAL_DENSITY"  # V/√Hz (ASD)
    POWER = "POWER"                # V²
    PSD = "PSD"                    # V²/Hz


# =============================================================================
# DISPLAY MANAGER
# =============================================================================

class DisplayManager(QObject):
    """
    Manages display updates from acquisition chunks.
    
    - Receives chunks directly from AcquisitionThread
    - Throttles display to max 30 Hz for fast sample rates
    - Computes FFT on each chunk (continuous data)
    - Updates trends from CircularBuffer at slower rate
    """
    
    # Signals to UI
    scopeDataReady = Signal(int, object, object)  # ch_idx, time_array, data_array
    fftDataReady = Signal(int, object, object)    # ch_idx, freq_array, magnitude_array
    statsReady = Signal(object)                    # {mean, std, min, max} per channel
    displayRateUpdated = Signal(float)            # Actual display rate in Hz
    hangDetected = Signal()                        # Emitted when display hangs for >5s
    chunkProcessed = Signal(int)                   # Emitted with chunk count for single mode
    
    def __init__(self, obs_config: ObservableConfig, buffers: Dict[int, CircularBuffer], parent=None):
        super().__init__(parent)
        self.obs_config = obs_config
        self.buffers = buffers  # For trends
        
        self.config = DisplayConfig()
        self._fft_window_cache = {}
        
        # Display throttling
        self._last_display_time = 0.0
        self._min_display_interval = 1.0 / self.config.max_display_rate  # 33ms for 30 Hz
        self._display_count = 0
        self._display_rate_timer = time.time()
        self._actual_display_rate = 0.0
        
        # Chunk counting for single mode
        self._chunk_count = 0
        
        # Hang detection
        self._last_chunk_time = 0.0
        self._hang_timeout = 5.0  # seconds
        
        # Stats timer (for trends - runs at slower rate)
        self._stats_timer = QTimer()
        self._stats_timer.timeout.connect(self._update_stats)
        
        # Display rate reporting timer
        self._rate_timer = QTimer()
        self._rate_timer.timeout.connect(self._report_display_rate)
        
        # Hang detection timer
        self._hang_timer = QTimer()
        self._hang_timer.timeout.connect(self._check_for_hang)
        
        self._running = False
        self._fft_enabled = False
        
    def start(self):
        """Start display manager"""
        self._running = True
        self._last_display_time = 0.0  # Ensures first chunk is always displayed
        self._last_chunk_time = time.time()
        self._display_count = 0
        self._chunk_count = 0  # Reset for single mode
        self._display_rate_timer = time.time()
        
        # Start stats timer
        stats_interval_ms = int(self.config.trends_interval * 1000)
        self._stats_timer.start(stats_interval_ms)
        
        # Start rate reporting (every 1 second)
        self._rate_timer.start(1000)
        
        # Start hang detection (check every 1 second)
        self._hang_timer.start(1000)
        
    def stop(self):
        """Stop display manager"""
        self._running = False
        self._stats_timer.stop()
        self._rate_timer.stop()
        self._hang_timer.stop()
        
    def _check_for_hang(self):
        """Check if acquisition has hung (no data for >5s)"""
        if not self._running:
            return
        elapsed = time.time() - self._last_chunk_time
        if elapsed > self._hang_timeout:
            self.hangDetected.emit()
            # Reset timer to avoid repeated signals
            self._last_chunk_time = time.time()
        
    def on_chunk_ready(self, chunk_data: Dict[int, np.ndarray]):
        """
        Handle incoming chunk from acquisition thread.
        
        This is called for every chunk, but display may be throttled.
        """
        if not self._running:
            return
            
        now = time.time()
        
        # Reset hang detection timer
        self._last_chunk_time = now
        
        # Increment chunk count and emit for single mode
        self._chunk_count += 1
        self.chunkProcessed.emit(self._chunk_count)
        
        # Throttle display updates if coming too fast
        elapsed = now - self._last_display_time
        if elapsed < self._min_display_interval:
            # Skip this display update (data still goes to buffer via acq thread)
            return
        
        self._last_display_time = now
        self._display_count += 1
        
        # Process and display this chunk
        sample_rate = self.obs_config.sample_rate
        if not chunk_data:
            return
            
        chunk_size = len(next(iter(chunk_data.values())))  # Get size from first channel
        
        # Create time axis for scope (chunk_size points, ending at t=0)
        duration = chunk_size / sample_rate
        t = np.linspace(-duration, 0, chunk_size, dtype=np.float32)
        
        for ch_idx, data in chunk_data.items():
            # Scope display (with optional decimation)
            if len(data) > self.config.scope_points:
                decimate_factor = len(data) // self.config.scope_points
                data_display = data[::decimate_factor]
                t_display = t[::decimate_factor]
            else:
                data_display = data
                t_display = t
            
            self.scopeDataReady.emit(ch_idx, t_display, data_display)
            
            # FFT (on full chunk - continuous data!)
            if self._fft_enabled:
                freq, mag = self._compute_fft(data, sample_rate)
                self.fftDataReady.emit(ch_idx, freq, mag)
    
    def _compute_fft(self, data: np.ndarray, sample_rate: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute FFT on chunk data with proper normalization.
        
        Returns raw amplitude spectrum (V). Display mode conversion happens in FFTWidget.
        Window metrics are passed along for proper spectral density calculations.
        """
        N = len(data)
        fs = sample_rate
        
        # Apply window
        window = self._get_window(N)
        xw = data * window
        
        # Compute FFT
        X = np.fft.rfft(xw)
        frequencies = np.fft.rfftfreq(N, 1/fs).astype(np.float32)
        
        # Window metrics
        W_sum = np.sum(window)      # For amplitude normalization
        W_sq = np.sum(window**2)    # For power normalization
        
        # Frequency resolution
        df = fs / N
        
        # Store window metrics for spectral density calculations
        # These will be used by FFTWidget for mode conversions
        self._last_window_sum = W_sum
        self._last_window_sq = W_sq
        self._last_df = df
        
        if self.config.power_correction:
            # Power-normalized: compute PSD first, then derive amplitude
            # PSD = (2 / (fs * W_sq)) * |X|²   [V²/Hz]
            psd = (2.0 / (fs * W_sq)) * np.abs(X)**2
            # Amplitude from PSD: A = sqrt(PSD * df)
            magnitude = np.sqrt(psd * df).astype(np.float32)
        else:
            # Amplitude-normalized (standard)
            # amp = (2 / W_sum) * |X|   [V]
            magnitude = ((2.0 / W_sum) * np.abs(X)).astype(np.float32)
        
        return frequencies, magnitude
    
    def _get_window(self, size: int) -> np.ndarray:
        """Get cached window function"""
        key = (self.config.fft_window, size)
        if key not in self._fft_window_cache:
            if self.config.fft_window == "Hann":
                self._fft_window_cache[key] = np.hanning(size).astype(np.float32)
            elif self.config.fft_window == "Hamming":
                self._fft_window_cache[key] = np.hamming(size).astype(np.float32)
            elif self.config.fft_window == "Blackman":
                self._fft_window_cache[key] = np.blackman(size).astype(np.float32)
            elif self.config.fft_window == "Blackman-Harris":
                n = np.arange(size)
                a0, a1, a2, a3 = 0.35875, 0.48829, 0.14128, 0.01168
                w = a0 - a1*np.cos(2*np.pi*n/(size-1)) + a2*np.cos(4*np.pi*n/(size-1)) - a3*np.cos(6*np.pi*n/(size-1))
                self._fft_window_cache[key] = w.astype(np.float32)
            elif self.config.fft_window == "Flat Top":
                n = np.arange(size)
                a0, a1, a2, a3, a4 = 0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368
                w = a0 - a1*np.cos(2*np.pi*n/(size-1)) + a2*np.cos(4*np.pi*n/(size-1)) - a3*np.cos(6*np.pi*n/(size-1)) + a4*np.cos(8*np.pi*n/(size-1))
                self._fft_window_cache[key] = w.astype(np.float32)
            else:  # Rectangular
                self._fft_window_cache[key] = np.ones(size, dtype=np.float32)
        return self._fft_window_cache[key]
    
    def _update_stats(self):
        """Compute statistics for trends (from CircularBuffer history)"""
        if not self._running or not self.buffers:
            return
            
        # Use recent samples for stats
        chunk_size = int(self.obs_config.sample_rate * self.config.trends_interval)
        chunk_size = max(100, min(chunk_size, 10000))
        
        stats = {
            'mean': [],
            'std': [],
            'min': [],
            'max': [],
            'timestamp': time.time()
        }
        
        for ch_idx in sorted(self.buffers.keys()):
            buffer = self.buffers[ch_idx]
            available = buffer.available
            
            n = min(available, chunk_size) if available > 0 else 0
            
            if n > 0:
                data = buffer.read_latest(n)
                stats['mean'].append(float(np.mean(data)))
                stats['std'].append(float(np.std(data)))
                stats['min'].append(float(np.min(data)))
                stats['max'].append(float(np.max(data)))
            else:
                stats['mean'].append(0.0)
                stats['std'].append(0.0)
                stats['min'].append(0.0)
                stats['max'].append(0.0)
                
        self.statsReady.emit(stats)
    
    def _report_display_rate(self):
        """Report actual display rate"""
        now = time.time()
        elapsed = now - self._display_rate_timer
        
        if elapsed > 0:
            self._actual_display_rate = self._display_count / elapsed
            self.displayRateUpdated.emit(self._actual_display_rate)
        
        # Reset counters
        self._display_count = 0
        self._display_rate_timer = now
    
    # === Configuration Methods ===
    
    def set_fft_enabled(self, enabled: bool):
        """Enable/disable FFT computation"""
        self._fft_enabled = enabled
        
    def set_scope_points(self, points: int):
        """Set max points to display (decimation target)"""
        self.config.scope_points = points
        
    def set_fft_size(self, size: int):
        """Set FFT size (also updates chunk size expectation)"""
        self.config.fft_size = size
        self._fft_window_cache.clear()
        
    def set_fft_window(self, window: str):
        """Set FFT window function"""
        self.config.fft_window = window
        self._fft_window_cache.clear()
        
    def set_max_display_rate(self, rate: float):
        """Set maximum display update rate"""
        self.config.max_display_rate = rate
        self._min_display_interval = 1.0 / rate
    
    def set_power_correction(self, enabled: bool):
        """Set power correction mode (True=Power, False=Amplitude)"""
        self.config.power_correction = enabled
        
    def get_fft_resolution(self) -> float:
        """Get FFT frequency resolution in Hz"""
        sample_rate = self.obs_config.sample_rate
        return sample_rate / self.config.fft_size
    
    def get_actual_display_rate(self) -> float:
        """Get the actual display update rate"""
        return self._actual_display_rate
    
    def get_window_metrics(self) -> Tuple[float, float, float]:
        """Get last computed window metrics (W_sum, W_sq, df)"""
        return (
            getattr(self, '_last_window_sum', 1.0),
            getattr(self, '_last_window_sq', 1.0),
            getattr(self, '_last_df', 1.0)
        )
    
    def is_power_correction(self) -> bool:
        """Check if power correction mode is enabled"""
        return self.config.power_correction
