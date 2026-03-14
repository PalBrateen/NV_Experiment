"""
Real-Time DAQ Data Streaming with PyQtGraph

This module provides a high-performance real-time data visualization system
for DAQ devices streaming at up to 2 MSa/s.

Architecture (mirroring LabVIEW VI):
- DAQ Assistant → AnalogInputTask (nidaqmx)
- Waveform Chart → Scrolling time-series plot
- XY Graph → FL vs RF Frequency plot
- Statistics → Mean, Std, Variance calculations
- Load Resistance → Current calculation

Key performance strategies:
1. Ring buffer for efficient memory management
2. Downsampling for display (can't render 2M points/sec)
3. Separate acquisition thread to prevent GUI blocking
4. Chunk-based updates to reduce plot overhead

Author: Based on LabVIEW VI conversion
"""

import sys
import numpy as np
from collections import deque
import time
from typing import Optional, Callable
from dataclasses import dataclass, field

# Qt imports
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QSpinBox, QDoubleSpinBox,
    QGroupBox, QCheckBox, QComboBox, QLineEdit, QFrame
)
from PySide6.QtCore import QThread, Signal, QTimer, Qt, QMutex, QMutexLocker

# PyQtGraph for high-performance plotting
import pyqtgraph as pg

# DAQ imports (conditional for systems without hardware)
# try:
import nidaqmx
from nidaqmx.constants import AcquisitionType, Edge, TerminalConfiguration
HAS_NIDAQMX = True
# except ImportError:
#     HAS_NIDAQMX = False
#     print("⚠ nidaqmx not available - using simulation mode")


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class DAQConfig:
    """DAQ acquisition configuration (mirrors LabVIEW DAQ Assistant settings)"""
    device: str = "P6363"
    channel: str = "ai15"
    sample_rate: float = 2e6  # 2 MSa/s max
    samples_per_chunk: int = 10000  # Samples per read operation
    voltage_min: float = -0.2
    voltage_max: float = 0.2
    terminal_config: str = "RSE"  # RSE, NRSE, Diff, PseudoDiff
    
@dataclass
class DisplayConfig:
    """Display and processing configuration"""
    load_resistance: float = 1000.0  # Ohms (from LabVIEW VI)
    display_points: int = 10000  # Max points shown on waveform chart
    update_rate_hz: float = 30.0  # GUI update rate
    downsample_factor: int = 100  # Downsample for display
    history_seconds: float = 5.0  # Seconds of data to show


# =============================================================================
# RING BUFFER FOR EFFICIENT STREAMING
# =============================================================================

class RingBuffer:
    """
    Lock-free ring buffer for high-speed data streaming.
    
    Critical for 2 MSa/s - can't use Python lists or frequent allocations.
    Pre-allocates numpy array and uses modular indexing.
    """
    
    def __init__(self, capacity: int, dtype=np.float64):
        self.capacity = capacity
        self.buffer = np.zeros(capacity, dtype=dtype)
        self.write_idx = 0
        self.count = 0
        self._mutex = QMutex()
    
    def write(self, data: np.ndarray):
        """Write chunk of data to buffer (thread-safe)"""
        with QMutexLocker(self._mutex):
            n = len(data)
            if n >= self.capacity:
                # Data larger than buffer - keep only most recent
                self.buffer[:] = data[-self.capacity:]
                self.write_idx = 0
                self.count = self.capacity
            else:
                # Wrap-around write
                end_idx = self.write_idx + n
                if end_idx <= self.capacity:
                    self.buffer[self.write_idx:end_idx] = data
                else:
                    # Split write
                    first_part = self.capacity - self.write_idx
                    self.buffer[self.write_idx:] = data[:first_part]
                    self.buffer[:n - first_part] = data[first_part:]
                
                self.write_idx = end_idx % self.capacity
                self.count = min(self.count + n, self.capacity)
    
    def read_all(self) -> np.ndarray:
        """Read all valid data in order (thread-safe)"""
        with QMutexLocker(self._mutex):
            if self.count == 0:
                return np.array([])
            
            if self.count < self.capacity:
                # Buffer not full yet
                return self.buffer[:self.count].copy()
            else:
                # Buffer full - rearrange to chronological order
                return np.roll(self.buffer, -self.write_idx).copy()
    
    def clear(self):
        """Clear buffer"""
        with QMutexLocker(self._mutex):
            self.write_idx = 0
            self.count = 0


# =============================================================================
# DAQ ACQUISITION THREAD
# =============================================================================

class DAQAcquisitionThread(QThread):
    """
    Background thread for continuous DAQ acquisition.
    
    Emits chunks of data to avoid GUI thread bottleneck.
    This mirrors the LabVIEW continuous acquisition loop.
    """
    
    # Signals
    data_ready = Signal(np.ndarray)  # Raw voltage data chunk
    stats_ready = Signal(dict)  # Processed statistics
    error_occurred = Signal(str)
    acquisition_started = Signal()
    acquisition_stopped = Signal()
    
    def __init__(self, config: DAQConfig, display_config: DisplayConfig):
        super().__init__()
        self.config = config
        self.display_config = display_config
        self._running = False
        self._task = None
        self._simulation_mode = not HAS_NIDAQMX
        
        # For simulation mode
        self._sim_time = 0.0
        self._sim_rf_freq = 2.87  # GHz - ESR center
        
    def set_rf_frequency(self, freq_ghz: float):
        """Update RF frequency for XY plot (from signal generator)"""
        self._sim_rf_freq = freq_ghz
        
    def run(self):
        """Main acquisition loop (runs in background thread)"""
        self._running = True
        self.acquisition_started.emit()
        
        try:
            if self._simulation_mode:
                self._run_simulation()
            else:
                self._run_hardware()
        except Exception as e:
            self.error_occurred.emit(str(e))
        finally:
            self._running = False
            self.acquisition_stopped.emit()
    
    def _run_hardware(self):
        """Hardware acquisition using nidaqmx"""
        with nidaqmx.Task() as task:
            # Configure analog input channel
            channel_path = f"{self.config.device}/{self.config.channel}"
            
            # Terminal configuration mapping
            term_configs = {
                "RSE": TerminalConfiguration.RSE,
                "NRSE": TerminalConfiguration.NRSE,
                "Diff": TerminalConfiguration.DIFF,
                "PseudoDiff": TerminalConfiguration.PSEUDO_DIFF
            }
            term_config = term_configs.get(self.config.terminal_config, 
                                           TerminalConfiguration.RSE)
            
            task.ai_channels.add_ai_voltage_chan(
                channel_path,
                min_val=self.config.voltage_min,
                max_val=self.config.voltage_max,
                terminal_config=term_config
            )
            
            # Configure timing for continuous acquisition
            task.timing.cfg_samp_clk_timing(
                rate=self.config.sample_rate,
                sample_mode=AcquisitionType.CONTINUOUS,
                samps_per_chan=self.config.samples_per_chunk * 10  # Buffer size
            )
            
            task.start()
            self._task = task
            
            while self._running:
                try:
                    # Read chunk of samples
                    data = np.array(task.read(
                        number_of_samples_per_channel=self.config.samples_per_chunk,
                        timeout=1.0
                    ))
                    
                    # Emit raw data
                    self.data_ready.emit(data)
                    
                    # Calculate and emit statistics
                    self._emit_statistics(data)
                    
                except nidaqmx.errors.DaqReadError:
                    # Timeout - no data available yet
                    pass
    
    def _run_simulation(self):
        """Simulated acquisition for testing without hardware"""
        dt = 1.0 / self.config.sample_rate
        chunk_duration = self.config.samples_per_chunk / self.config.sample_rate
        
        while self._running:
            # Generate time array for this chunk
            t = self._sim_time + np.arange(self.config.samples_per_chunk) * dt
            self._sim_time += chunk_duration
            
            # Simulate fluorescence signal with ESR dip
            # Base signal: ~1V with small drift
            base_signal = 1.0 + 0.1 * np.sin(2 * np.pi * 0.1 * t)
            
            # Add ESR resonance dip (Lorentzian) at 2.87 GHz
            resonance_center = 2.87  # GHz
            linewidth = 0.01  # GHz
            dip_depth = 0.1  # 10% contrast
            
            # Calculate dip based on RF frequency
            detuning = self._sim_rf_freq - resonance_center
            lorentzian = 1.0 / (1.0 + (detuning / linewidth)**2)
            
            signal = base_signal * (1.0 - dip_depth * lorentzian)
            
            # Add realistic noise
            noise = np.random.normal(0, 0.01, len(t))
            signal += noise
            
            # Emit data
            self.data_ready.emit(signal)
            self._emit_statistics(signal)
            
            # Sleep to match real acquisition rate
            time.sleep(chunk_duration * 0.9)  # Slightly fast to prevent underrun
    
    def _emit_statistics(self, data: np.ndarray):
        """Calculate statistics (mirrors LabVIEW statistics blocks)"""
        R_load = self.display_config.load_resistance
        
        # Convert to mV (×1000 in LabVIEW)
        data_mv = data * 1000
        
        # Statistics
        mean_mv = np.mean(data_mv)
        std_mv = np.std(data_mv)
        variance = np.var(data_mv)
        
        # Current calculation: I = V / R_load (in mA)
        current_ma = mean_mv / R_load
        
        # Percent deviation
        if mean_mv != 0:
            percent_deviation = (std_mv / abs(mean_mv)) * 100
        else:
            percent_deviation = 0.0
        
        stats = {
            'mean_mv': mean_mv,
            'std_mv': std_mv,
            'variance': variance,
            'current_ma': current_ma,
            'percent_deviation': percent_deviation,
            'rf_freq_ghz': self._sim_rf_freq
        }
        
        self.stats_ready.emit(stats)
    
    def stop(self):
        """Stop acquisition"""
        self._running = False
        self.wait()  # Wait for thread to finish


# =============================================================================
# MAIN GUI WINDOW
# =============================================================================

class RealTimeDAQViewer(QMainWindow):
    """
    Main application window for real-time DAQ visualization.
    
    Layout mirrors LabVIEW VI:
    - Waveform Chart (time series)
    - XY Graph (FL vs RF Frequency)
    - Statistics displays
    - Control buttons
    """
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Real-Time DAQ Viewer (2 MSa/s)")
        self.setGeometry(100, 100, 1400, 800)
        
        # Configuration
        self.daq_config = DAQConfig()
        self.display_config = DisplayConfig()
        
        # Data storage
        buffer_size = int(self.display_config.history_seconds * self.daq_config.sample_rate)
        self.data_buffer = RingBuffer(buffer_size)
        self.time_buffer = RingBuffer(buffer_size)
        
        # XY plot data (FL mean vs RF freq)
        self.xy_data = {'freq': [], 'fl': []}
        
        # Acquisition thread
        self.acq_thread: Optional[DAQAcquisitionThread] = None
        
        # Current time reference
        self.start_time = 0.0
        self.current_time = 0.0
        
        # Setup UI
        self._setup_ui()
        self._setup_connections()
        
        # Update timer for plots
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self._update_plots)
        
        # Performance monitoring
        self.frame_count = 0
        self.fps_timer = QTimer()
        self.fps_timer.timeout.connect(self._update_fps)
        self.fps_timer.start(1000)
        
    def _setup_ui(self):
        """Create the user interface"""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QHBoxLayout(central_widget)
        
        # Left panel: Controls
        control_panel = self._create_control_panel()
        main_layout.addWidget(control_panel)
        
        # Right panel: Plots and displays
        plot_panel = self._create_plot_panel()
        main_layout.addWidget(plot_panel, stretch=3)
        
    def _create_control_panel(self) -> QWidget:
        """Create control panel (left side)"""
        panel = QGroupBox("Controls")
        layout = QVBoxLayout(panel)
        
        # DAQ Settings
        daq_group = QGroupBox("DAQ Settings")
        daq_layout = QGridLayout(daq_group)
        
        # Sample rate
        daq_layout.addWidget(QLabel("Sample Rate:"), 0, 0)
        self.rate_spin = QDoubleSpinBox()
        self.rate_spin.setRange(1, 2e6)
        self.rate_spin.setValue(self.daq_config.sample_rate)
        self.rate_spin.setSuffix(" Sa/s")
        self.rate_spin.setDecimals(0)
        daq_layout.addWidget(self.rate_spin, 0, 1)
        
        # Samples per chunk
        daq_layout.addWidget(QLabel("Chunk Size:"), 1, 0)
        self.chunk_spin = QSpinBox()
        self.chunk_spin.setRange(100, 100000)
        self.chunk_spin.setValue(self.daq_config.samples_per_chunk)
        daq_layout.addWidget(self.chunk_spin, 1, 1)
        
        # Device/Channel (mirrors LabVIEW "Analog Input Device")
        daq_layout.addWidget(QLabel("Device:"), 2, 0)
        self.device_edit = QLineEdit(self.daq_config.device)
        daq_layout.addWidget(self.device_edit, 2, 1)
        
        daq_layout.addWidget(QLabel("Channel:"), 3, 0)
        self.channel_edit = QLineEdit(self.daq_config.channel)
        daq_layout.addWidget(self.channel_edit, 3, 1)
        
        layout.addWidget(daq_group)
        
        # Processing Settings (mirrors LabVIEW constants)
        proc_group = QGroupBox("Processing")
        proc_layout = QGridLayout(proc_group)
        
        # Load resistance (from LabVIEW)
        proc_layout.addWidget(QLabel("Load R (Ω):"), 0, 0)
        self.load_r_spin = QDoubleSpinBox()
        self.load_r_spin.setRange(1, 1e6)
        self.load_r_spin.setValue(self.display_config.load_resistance)
        proc_layout.addWidget(self.load_r_spin, 0, 1)
        
        # RF Frequency input (for XY plot)
        proc_layout.addWidget(QLabel("RF Freq (GHz):"), 1, 0)
        self.rf_freq_spin = QDoubleSpinBox()
        self.rf_freq_spin.setRange(0, 10)
        self.rf_freq_spin.setValue(2.87)
        self.rf_freq_spin.setDecimals(4)
        self.rf_freq_spin.setSingleStep(0.001)
        proc_layout.addWidget(self.rf_freq_spin, 1, 1)
        
        layout.addWidget(proc_group)
        
        # Start/Stop buttons
        btn_layout = QHBoxLayout()
        
        self.start_btn = QPushButton("▶ Start")
        self.start_btn.setStyleSheet("background-color: #4CAF50; color: white;")
        btn_layout.addWidget(self.start_btn)
        
        self.stop_btn = QPushButton("■ Stop")
        self.stop_btn.setStyleSheet("background-color: #f44336; color: white;")
        self.stop_btn.setEnabled(False)
        btn_layout.addWidget(self.stop_btn)
        
        layout.addLayout(btn_layout)
        
        # Reset buttons (mirrors LabVIEW reset controls)
        reset_layout = QHBoxLayout()
        
        self.reset_fl_btn = QPushButton("Reset FL Plot")
        reset_layout.addWidget(self.reset_fl_btn)
        
        self.reset_xy_btn = QPushButton("Reset XY Plot")
        reset_layout.addWidget(self.reset_xy_btn)
        
        layout.addLayout(reset_layout)
        
        # Status
        self.status_label = QLabel("Status: Stopped")
        layout.addWidget(self.status_label)
        
        self.fps_label = QLabel("FPS: --")
        layout.addWidget(self.fps_label)
        
        layout.addStretch()
        
        return panel
        
    def _create_plot_panel(self) -> QWidget:
        """Create plot panel (right side)"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        # Configure PyQtGraph for performance
        pg.setConfigOptions(
            antialias=False,  # Disable antialiasing for speed
            useOpenGL=True,   # Use OpenGL acceleration
            enableExperimental=True
        )
        
        # Top row: Waveform charts (mirrors LabVIEW Waveform Chart)
        charts_layout = QHBoxLayout()
        
        # Main waveform chart
        self.waveform_widget = pg.PlotWidget(title="Fluorescence Signal (Waveform Chart)")
        self.waveform_widget.setLabel('left', 'Voltage', units='V')
        self.waveform_widget.setLabel('bottom', 'Time', units='s')
        self.waveform_widget.showGrid(x=True, y=True, alpha=0.3)
        self.waveform_curve = self.waveform_widget.plot(pen=pg.mkPen('y', width=1))
        charts_layout.addWidget(self.waveform_widget)
        
        # Secondary waveform chart (Waveform Chart 2 in LabVIEW)
        self.waveform2_widget = pg.PlotWidget(title="Processed Signal")
        self.waveform2_widget.setLabel('left', 'Current', units='mA')
        self.waveform2_widget.setLabel('bottom', 'Time', units='s')
        self.waveform2_widget.showGrid(x=True, y=True, alpha=0.3)
        self.waveform2_curve = self.waveform2_widget.plot(pen=pg.mkPen('c', width=1))
        charts_layout.addWidget(self.waveform2_widget)
        
        layout.addLayout(charts_layout)
        
        # Bottom row: XY Graph and Statistics
        bottom_layout = QHBoxLayout()
        
        # XY Graph (FL vs RF Frequency - mirrors Build XY Graph)
        self.xy_widget = pg.PlotWidget(title="FL vs RF Frequency (XY Graph)")
        self.xy_widget.setLabel('left', 'FL Mean', units='mV')
        self.xy_widget.setLabel('bottom', 'RF Frequency', units='GHz')
        self.xy_widget.showGrid(x=True, y=True, alpha=0.3)
        self.xy_curve = self.xy_widget.plot(
            pen=None,
            symbol='o',
            symbolSize=5,
            symbolBrush='m'
        )
        bottom_layout.addWidget(self.xy_widget)
        
        # Statistics panel (mirrors LabVIEW numeric indicators)
        stats_widget = self._create_stats_panel()
        bottom_layout.addWidget(stats_widget)
        
        layout.addLayout(bottom_layout)
        
        return panel
    
    def _create_stats_panel(self) -> QWidget:
        """Create statistics display panel"""
        panel = QGroupBox("Statistics")
        layout = QGridLayout(panel)
        
        # FL Mean (mV)
        layout.addWidget(QLabel("FL Mean (mV):"), 0, 0)
        self.fl_mean_display = QLineEdit("--")
        self.fl_mean_display.setReadOnly(True)
        self.fl_mean_display.setAlignment(Qt.AlignRight)
        layout.addWidget(self.fl_mean_display, 0, 1)
        
        # FL Std (mV)
        layout.addWidget(QLabel("FL Std (mV):"), 1, 0)
        self.fl_std_display = QLineEdit("--")
        self.fl_std_display.setReadOnly(True)
        self.fl_std_display.setAlignment(Qt.AlignRight)
        layout.addWidget(self.fl_std_display, 1, 1)
        
        # FL % Deviation
        layout.addWidget(QLabel("FL % Deviation:"), 2, 0)
        self.fl_dev_display = QLineEdit("--")
        self.fl_dev_display.setReadOnly(True)
        self.fl_dev_display.setAlignment(Qt.AlignRight)
        layout.addWidget(self.fl_dev_display, 2, 1)
        
        # Anode Current (mA)
        layout.addWidget(QLabel("Current (mA):"), 3, 0)
        self.current_display = QLineEdit("--")
        self.current_display.setReadOnly(True)
        self.current_display.setAlignment(Qt.AlignRight)
        layout.addWidget(self.current_display, 3, 1)
        
        # Excess current indicator (like LabVIEW boolean)
        layout.addWidget(QLabel("Excess Current:"), 4, 0)
        self.excess_current_indicator = QLabel("●")
        self.excess_current_indicator.setStyleSheet("color: green; font-size: 20px;")
        layout.addWidget(self.excess_current_indicator, 4, 1)
        
        # Load Resistance display
        layout.addWidget(QLabel("Load R (Ω):"), 5, 0)
        self.load_r_display = QLineEdit(f"{self.display_config.load_resistance:.1f}")
        self.load_r_display.setReadOnly(True)
        self.load_r_display.setAlignment(Qt.AlignRight)
        layout.addWidget(self.load_r_display, 5, 1)
        
        return panel
    
    def _setup_connections(self):
        """Connect signals and slots"""
        self.start_btn.clicked.connect(self._start_acquisition)
        self.stop_btn.clicked.connect(self._stop_acquisition)
        self.reset_fl_btn.clicked.connect(self._reset_fl_plot)
        self.reset_xy_btn.clicked.connect(self._reset_xy_plot)
        
        # Parameter changes
        self.load_r_spin.valueChanged.connect(self._update_load_resistance)
        self.rf_freq_spin.valueChanged.connect(self._update_rf_frequency)
        
    def _start_acquisition(self):
        """Start DAQ acquisition"""
        # Update config from UI
        self.daq_config.sample_rate = self.rate_spin.value()
        self.daq_config.samples_per_chunk = self.chunk_spin.value()
        self.daq_config.device = self.device_edit.text()
        self.daq_config.channel = self.channel_edit.text()
        
        # Clear buffers
        self.data_buffer.clear()
        self.time_buffer.clear()
        self.start_time = time.time()
        self.current_time = 0.0
        
        # Create and start acquisition thread
        self.acq_thread = DAQAcquisitionThread(self.daq_config, self.display_config)
        self.acq_thread.data_ready.connect(self._handle_data)
        self.acq_thread.stats_ready.connect(self._handle_stats)
        self.acq_thread.error_occurred.connect(self._handle_error)
        self.acq_thread.acquisition_started.connect(
            lambda: self.status_label.setText("Status: Running")
        )
        self.acq_thread.acquisition_stopped.connect(
            lambda: self.status_label.setText("Status: Stopped")
        )
        
        self.acq_thread.start()
        
        # Start update timer
        update_interval = int(1000 / self.display_config.update_rate_hz)
        self.update_timer.start(update_interval)
        
        # Update UI
        self.start_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        
    def _stop_acquisition(self):
        """Stop DAQ acquisition"""
        if self.acq_thread:
            self.acq_thread.stop()
            self.acq_thread = None
        
        self.update_timer.stop()
        
        # Update UI
        self.start_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        
    def _handle_data(self, data: np.ndarray):
        """Handle incoming data chunk (called from acquisition thread)"""
        # Update buffers
        self.data_buffer.write(data)
        
        # Generate time array for this chunk
        dt = 1.0 / self.daq_config.sample_rate
        chunk_duration = len(data) * dt
        t = self.current_time + np.arange(len(data)) * dt
        self.current_time += chunk_duration
        self.time_buffer.write(t)
        
    def _handle_stats(self, stats: dict):
        """Handle statistics update"""
        # Update displays (mirrors LabVIEW numeric indicators)
        self.fl_mean_display.setText(f"{stats['mean_mv']:.3f}")
        self.fl_std_display.setText(f"{stats['std_mv']:.3f}")
        self.fl_dev_display.setText(f"{stats['percent_deviation']:.2f}")
        self.current_display.setText(f"{stats['current_ma']:.4f}")
        
        # Update excess current indicator
        threshold = 0.1  # mA threshold - adjust as needed
        if abs(stats['current_ma']) > threshold:
            self.excess_current_indicator.setStyleSheet(
                "color: red; font-size: 20px;"
            )
        else:
            self.excess_current_indicator.setStyleSheet(
                "color: green; font-size: 20px;"
            )
        
        # Add point to XY plot
        self.xy_data['freq'].append(stats['rf_freq_ghz'])
        self.xy_data['fl'].append(stats['mean_mv'])
        
    def _handle_error(self, error_msg: str):
        """Handle acquisition error"""
        self.status_label.setText(f"Error: {error_msg}")
        self._stop_acquisition()
        
    def _update_plots(self):
        """Update plot displays (called by timer)"""
        self.frame_count += 1
        
        # Get data from buffers
        data = self.data_buffer.read_all()
        time_data = self.time_buffer.read_all()
        
        if len(data) == 0:
            return
        
        # Downsample for display (critical for performance at 2 MSa/s)
        ds = self.display_config.downsample_factor
        if len(data) > ds:
            # Use efficient decimation
            display_data = data[::ds]
            display_time = time_data[::ds]
        else:
            display_data = data
            display_time = time_data
        
        # Limit display points
        max_points = self.display_config.display_points
        if len(display_data) > max_points:
            display_data = display_data[-max_points:]
            display_time = display_time[-max_points:]
        
        # Update waveform chart
        self.waveform_curve.setData(display_time, display_data)
        
        # Update processed signal (current)
        R_load = self.display_config.load_resistance
        current_data = (display_data * 1000) / R_load  # mA
        self.waveform2_curve.setData(display_time, current_data)
        
        # Update XY plot
        if len(self.xy_data['freq']) > 0:
            self.xy_curve.setData(
                self.xy_data['freq'],
                self.xy_data['fl']
            )
    
    def _update_fps(self):
        """Update FPS display"""
        self.fps_label.setText(f"FPS: {self.frame_count}")
        self.frame_count = 0
        
    def _reset_fl_plot(self):
        """Reset FL waveform plot (mirrors LabVIEW Reset FL plot button)"""
        self.data_buffer.clear()
        self.time_buffer.clear()
        self.start_time = time.time()
        self.current_time = 0.0
        self.waveform_curve.setData([], [])
        self.waveform2_curve.setData([], [])
        
    def _reset_xy_plot(self):
        """Reset XY plot (mirrors LabVIEW Reset Plots button)"""
        self.xy_data = {'freq': [], 'fl': []}
        self.xy_curve.setData([], [])
        
    def _update_load_resistance(self, value: float):
        """Update load resistance"""
        self.display_config.load_resistance = value
        self.load_r_display.setText(f"{value:.1f}")
        
    def _update_rf_frequency(self, value: float):
        """Update RF frequency for acquisition thread"""
        if self.acq_thread:
            self.acq_thread.set_rf_frequency(value)
            
    def closeEvent(self, event):
        """Clean shutdown"""
        self._stop_acquisition()
        event.accept()


# =============================================================================
# INTEGRATION WITH YOUR EXISTING DAQ SYSTEM
# =============================================================================

class DAQStreamAdapter:
    """
    Adapter to integrate with your existing DAQcontrol.py AnalogInputTask.
    
    Usage:
        from DAQcontrol import AnalogInputTask
        
        # Create your existing task
        ai_task = AnalogInputTask(...)
        ai_task.configure_continuous(sample_rate=2e6)
        
        # Wrap with adapter
        adapter = DAQStreamAdapter(ai_task)
        adapter.start_streaming(callback=my_data_handler)
    """
    
    def __init__(self, ai_task):
        """
        Args:
            ai_task: Your existing AnalogInputTask instance
        """
        self.ai_task = ai_task
        self._streaming = False
        self._thread = None
        
    def start_streaming(self, callback: Callable[[np.ndarray], None], 
                       samples_per_chunk: int = 10000):
        """
        Start continuous data streaming with callback.
        
        Args:
            callback: Function called with each data chunk
            samples_per_chunk: Samples per callback invocation
        """
        self._streaming = True
        self._callback = callback
        self._samples_per_chunk = samples_per_chunk
        
        # Start background thread
        import threading
        self._thread = threading.Thread(target=self._stream_loop, daemon=True)
        self._thread.start()
        
    def _stream_loop(self):
        """Background streaming loop"""
        while self._streaming:
            try:
                data = self.ai_task.read_daq(self._samples_per_chunk)
                if data is not None and len(data) > 0:
                    self._callback(np.array(data))
            except Exception as e:
                print(f"Stream error: {e}")
                break
                
    def stop_streaming(self):
        """Stop streaming"""
        self._streaming = False
        if self._thread:
            self._thread.join(timeout=1.0)


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    """Run the real-time DAQ viewer application"""
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    # Create and show main window
    window = RealTimeDAQViewer()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()