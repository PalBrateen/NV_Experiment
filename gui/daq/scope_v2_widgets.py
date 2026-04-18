"""
Scope Viewer v2 - Widgets Module
================================

Custom Qt widgets and plot widgets.

Updates:
- Multi-channel FFT support
- Cursors for scope and FFT
- Averaging support
- Spectral density / Power modes
- Improved gridlines
"""

import numpy as np, pyqtgraph as pg, time
from typing import Optional, List, Dict, Tuple
from enum import Enum

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, 
    QPushButton, QCheckBox, QComboBox, QDoubleSpinBox, QGroupBox,
    QScrollArea
)
from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtGui import QColor, QPainter, QBrush, QPen, QFont

from scope_v2_core import (
    ChannelConfig, DAQDeviceInfo, CircularBuffer,
    ObservableConfig, get_device_info, enumerate_devices
)


# =============================================================================
# ENUMS
# =============================================================================

class FFTDisplayMode(Enum):
    """FFT display mode - what quantity to show"""
    NORMAL = "V"               # Voltage (V)
    SPECTRAL_DENSITY = "V/√Hz" # Voltage spectral density
    POWER = "V²"               # Power (voltage squared)
    PSD = "V²/Hz"              # Power spectral density


class FrequencyUnit(Enum):
    """Frequency axis units"""
    HZ = ("Hz", 1.0)
    KHZ = ("kHz", 1e3)
    MHZ = ("MHz", 1e6)
    
    def __init__(self, label: str, divisor: float):
        self.label = label
        self.divisor = divisor


# =============================================================================
# BASIC WIDGETS
# =============================================================================

class ToggleSwitch(QCheckBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedSize(38, 21)
        self._knob_position = 4  # starting position
        self.animation = QPropertyAnimation(self, b"knob_position")
        self.animation.setDuration(150)
        self.animation.setEasingCurve(QEasingCurve.Type.InOutCubic)
        self.stateChanged.connect(self._start_animation)
    
    def hitButton(self, pos):
        return self.rect().contains(pos)
    
    def _start_animation(self, state):
        self.animation.setStartValue(self._knob_position)
        self.animation.setEndValue(20 if state else 4)
        self.animation.start()

    def get_knob_position(self):
        return self._knob_position

    def set_knob_position(self, pos):
        self._knob_position = pos
        self.update()

    knob_position = Property(float, get_knob_position, set_knob_position)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Track
        track_color = QColor("#2196F3") if self.isChecked() else QColor("#CCCCCC")
        p.setBrush(QBrush(track_color))
        p.setPen(Qt.PenStyle.NoPen)
        p.drawRoundedRect(0, 0, self.width(), self.height(), 10, 10)
        
        # Knob
        p.setBrush(QBrush(QColor("white")))
        p.drawEllipse(int(self._knob_position), 3, 15, 15)


class SegmentedButton(QWidget):
    """Segmented button group"""
    valueChanged = Signal(int)
    
    def __init__(self, options: List[str], parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        
        self.buttons = []
        self._value = 0
        
        for i, opt in enumerate(options):
            btn = QPushButton(opt)
            btn.setCheckable(True)
            btn.setChecked(i == 0)
            btn.clicked.connect(lambda _, idx=i: self._on_clicked(idx))
            
            if i == 0:
                btn.setStyleSheet("border-top-right-radius:0;border-bottom-right-radius:0;")
            elif i == len(options) - 1:
                btn.setStyleSheet("border-top-left-radius:0;border-bottom-left-radius:0;")
            else:
                btn.setStyleSheet("border-radius:0;")
                
            layout.addWidget(btn)
            self.buttons.append(btn)
            
    def _on_clicked(self, idx: int):
        self._value = idx
        for i, btn in enumerate(self.buttons):
            btn.setChecked(i == idx)
        self.valueChanged.emit(idx)
        
    def value(self) -> int:
        return self._value
        
    def setValue(self, idx: int):
        if 0 <= idx < len(self.buttons):
            self._on_clicked(idx)


# =============================================================================
# CUSTOM AXIS - NO MINOR TICK LABELS
# =============================================================================

class CleanAxis(pg.AxisItem):
    """Axis that only shows labels for major ticks (hides minor tick labels).
    
    For pyqtgraph 0.14.0 - overrides drawPicture to filter tick labels.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def tickSpacing(self, minVal, maxVal, size):
        """Override to return only major tick spacing"""
        # Get default spacings
        spacings = super().tickSpacing(minVal, maxVal, size)
        
        # Return only the largest (major) spacing to reduce tick density
        if len(spacings) > 1:
            # Keep only major ticks (first one or two levels)
            return spacings[:1]
        return spacings


# =============================================================================
# PLOT WIDGETS
# =============================================================================

class ScopeWidget(pg.PlotWidget):
    """Main oscilloscope display with mouse crosshair cursor, averaging, and histogram"""
    
    cursorMoved = Signal(dict)  # Emits cursor positions
    
    def __init__(self, parent=None):
        # Create with custom axes
        super().__init__(parent, axisItems={
            'bottom': CleanAxis(orientation='bottom'),
            'left': CleanAxis(orientation='left')
        })
        
        self.setBackground('#1a1a1a')
        self.setLabel('left', 'Voltage', 'V')
        self.setLabel('bottom', 'Time', 's')
        
        # Configure grid
        self.showGrid(x=True, y=True, alpha=0.3)
        
        self.curves = {}
        self.colors = ['#00BFFF', '#FF6B6B', '#4ECDC4', '#FFE66D', '#95E1D3']
        
        # Mouse crosshair cursor
        self._cursors_enabled = False
        self._vline = pg.InfiniteLine(angle=90, movable=False, 
                                       pen=pg.mkPen('#FFAA00', width=1, style=Qt.PenStyle.DashLine))
        self._hline = pg.InfiniteLine(angle=0, movable=False, 
                                       pen=pg.mkPen('#FFAA00', width=1, style=Qt.PenStyle.DashLine))
        self._vline.setVisible(False)
        self._hline.setVisible(False)
        self.addItem(self._vline, ignoreBounds=True)
        self.addItem(self._hline, ignoreBounds=True)
        
        # Floating cursor label on plot
        self._cursor_label = pg.TextItem(color='#FFAA00', anchor=(0, 1))
        self._cursor_label.setFont(pg.QtGui.QFont('Consolas', 9))
        self._cursor_label.setVisible(False)
        self.addItem(self._cursor_label, ignoreBounds=True)
        
        # Enable mouse tracking
        self.setMouseTracking(True)
        self.scene().sigMouseMoved.connect(self._on_mouse_moved)
        
        # Averaging
        self._averaging_enabled = False
        self._avg_count = 1
        self._avg_buffers: Dict[int, List[np.ndarray]] = {}  # ch_idx -> list of traces
        
        # Histogram display (rotated, shown at left edge)
        self._histogram_enabled = False
        self._histogram_bins = 50
        self._histogram_curves: Dict[int, pg.PlotCurveItem] = {}  # Per-channel histograms
        self._histogram_fills: Dict[int, pg.FillBetweenItem] = {}
        self._last_data: Dict[int, np.ndarray] = {}  # Store last data for histogram
        self._histogram_width_fraction = 0.15  # Histogram takes 15% of plot width
        
    def update_data(self, ch_idx: int, time_data: np.ndarray, volt_data: np.ndarray):
        if ch_idx not in self.curves:
            color = self.colors[ch_idx % len(self.colors)]
            self.curves[ch_idx] = self.plot(pen=pg.mkPen(color, width=1))
            
        # Apply averaging if enabled
        if self._averaging_enabled and self._avg_count > 1:
            volt_data = self._apply_averaging(ch_idx, volt_data)
            
        self.curves[ch_idx].setData(time_data, volt_data)
        
        # Store data for histogram
        self._last_data[ch_idx] = volt_data.copy()
        
        # Update histogram if enabled
        if self._histogram_enabled:
            self._update_histogram(ch_idx, volt_data, time_data)
    
    def _update_histogram(self, ch_idx: int, data: np.ndarray, time_data: np.ndarray):
        """Update rotated histogram for channel, positioned at left edge"""
        if len(data) < 10:
            return
        
        color = self.colors[ch_idx % len(self.colors)]
        qcolor = pg.mkColor(color)
        
        # Create histogram items if needed
        if ch_idx not in self._histogram_curves:
            # Curve for histogram outline
            hist_curve = pg.PlotCurveItem(pen=pg.mkPen(color, width=1))
            self.addItem(hist_curve)
            
            # Zero line for fill reference
            zero_curve = pg.PlotCurveItem(pen=None)
            self.addItem(zero_curve)
            
            # Fill between histogram and zero
            fill_color = (qcolor.red(), qcolor.green(), qcolor.blue(), 60)
            fill = pg.FillBetweenItem(zero_curve, hist_curve, brush=fill_color)
            self.addItem(fill)
            
            self._histogram_curves[ch_idx] = (hist_curve, zero_curve)
            self._histogram_fills[ch_idx] = fill
        
        hist_curve, zero_curve = self._histogram_curves[ch_idx]
        
        # Compute histogram
        counts, bin_edges = np.histogram(data, bins=self._histogram_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Scale counts to fit in left portion of plot
        # Get current x-range from time data
        if len(time_data) > 0:
            x_min = time_data[0]
            x_max = time_data[-1]
            x_range = x_max - x_min
        else:
            x_min, x_max = -1, 0
            x_range = 1
        
        # Normalize and scale counts
        max_count = np.max(counts) if np.max(counts) > 0 else 1
        scaled_counts = (counts / max_count) * (x_range * self._histogram_width_fraction)
        
        # Position histogram at left edge (x_min), extending right
        hist_x = x_min + scaled_counts
        
        # Set data: X is position, Y is voltage (bin centers)
        hist_curve.setData(hist_x, bin_centers)
        
        # Zero line at x_min
        zero_x = np.full_like(bin_centers, x_min)
        zero_curve.setData(zero_x, bin_centers)
        
    def set_histogram_enabled(self, enabled: bool):
        """Enable/disable histogram display"""
        self._histogram_enabled = enabled
        
        # Show/hide histogram items
        for ch_idx in list(self._histogram_curves.keys()):
            hist_curve, zero_curve = self._histogram_curves[ch_idx]
            hist_curve.setVisible(enabled)
            zero_curve.setVisible(enabled)
            if ch_idx in self._histogram_fills:
                self._histogram_fills[ch_idx].setVisible(enabled)
        
        # If enabling and we have data, update histograms
        if enabled:
            for ch_idx, data in self._last_data.items():
                if ch_idx in self.curves:
                    xdata, _ = self.curves[ch_idx].getData()
                    if xdata is not None and len(xdata) > 0:
                        self._update_histogram(ch_idx, data, xdata)
    
    def set_histogram_bins(self, bins: int):
        """Set number of histogram bins"""
        self._histogram_bins = max(10, min(200, bins))
        
    def _apply_averaging(self, ch_idx: int, data: np.ndarray) -> np.ndarray:
        """Apply exponential moving average or block averaging"""
        if ch_idx not in self._avg_buffers:
            self._avg_buffers[ch_idx] = []
            
        buf = self._avg_buffers[ch_idx]
        
        # Store trace in buffer
        if len(buf) == 0 or len(buf[-1]) != len(data):
            # Size changed, reset buffer
            buf.clear()
            
        buf.append(data.copy())
        
        # Keep only avg_count traces
        while len(buf) > self._avg_count:
            buf.pop(0)
            
        # Average all stored traces
        if len(buf) == 1:
            return data
        return np.mean(buf, axis=0).astype(np.float32)
        
    def set_averaging(self, enabled: bool, count: int = 1):
        """Enable/disable averaging"""
        self._averaging_enabled = enabled
        self._avg_count = max(1, count)
        if not enabled:
            self._avg_buffers.clear()
            
    def clear_averaging(self):
        """Clear averaging buffers"""
        self._avg_buffers.clear()
        
    def clear_all(self):
        for curve in self.curves.values():
            curve.setData([], [])
        self._avg_buffers.clear()
        self._last_data.clear()
        # Clear histogram displays
        for ch_idx in list(self._histogram_curves.keys()):
            hist_curve, zero_curve = self._histogram_curves[ch_idx]
            hist_curve.setData([], [])
            zero_curve.setData([], [])
    
    def remove_channel(self, ch_idx: int):
        """Remove a channel's curve from the plot"""
        if ch_idx in self.curves:
            self.removeItem(self.curves[ch_idx])
            del self.curves[ch_idx]
        if ch_idx in self._avg_buffers:
            del self._avg_buffers[ch_idx]
        if ch_idx in self._last_data:
            del self._last_data[ch_idx]
        # Remove histogram items
        if ch_idx in self._histogram_curves:
            hist_curve, zero_curve = self._histogram_curves[ch_idx]
            self.removeItem(hist_curve)
            self.removeItem(zero_curve)
            del self._histogram_curves[ch_idx]
        if ch_idx in self._histogram_fills:
            self.removeItem(self._histogram_fills[ch_idx])
            del self._histogram_fills[ch_idx]
        
    # === Mouse Cursor Methods ===
    
    def set_cursors_enabled(self, enabled: bool):
        """Enable/disable crosshair cursor"""
        self._cursors_enabled = enabled
        self._vline.setVisible(enabled)
        self._hline.setVisible(enabled)
        self._cursor_label.setVisible(enabled)
        if not enabled:
            self.cursorMoved.emit({})
            
    def _on_mouse_moved(self, pos):
        """Handle mouse movement for crosshair"""
        if not self._cursors_enabled:
            return
        
        vb = self.getPlotItem().getViewBox()
        
        if self.sceneBoundingRect().contains(pos):
            mouse_point = vb.mapSceneToView(pos)
            x = mouse_point.x()
            y = mouse_point.y()
            
            self._vline.setPos(x)
            self._hline.setPos(y)
            
            # Find data values at cursor for all channels
            channel_values = {}
            for ch_idx, curve in self.curves.items():
                xdata, ydata = curve.getData()
                if xdata is not None and len(xdata) > 0:
                    idx = np.argmin(np.abs(xdata - x))
                    channel_values[ch_idx] = float(ydata[idx])
            
            # Format floating label
            if abs(x) >= 1:
                t_str = f"T: {x:.4f} s"
            elif abs(x) >= 1e-3:
                t_str = f"T: {x*1e3:.3f} ms"
            else:
                t_str = f"T: {x*1e6:.2f} µs"
            
            lines = [t_str]
            for ch_idx, val in sorted(channel_values.items()):
                lines.append(f"Ch{ch_idx+1}: {val*1e3:.2f} mV")
            
            self._cursor_label.setText("\n".join(lines))
            
            # Position label near cursor but offset to stay visible
            view_range = self.viewRange()
            x_range = view_range[0][1] - view_range[0][0]
            y_range = view_range[1][1] - view_range[1][0]
            
            # Offset label to upper-right of cursor
            label_x = x + x_range * 0.02
            label_y = y + y_range * 0.02
            self._cursor_label.setPos(label_x, label_y)
            
            cursor_data = {
                "x": x,
                "y": y,
                "channel_values": channel_values
            }
            self.cursorMoved.emit(cursor_data)
            
    def get_cursor_values(self) -> dict:
        return {}


class TrendsWidget(pg.PlotWidget):
    """Trends display with circular buffer - supports multiple channels"""
    
    COLORS = ['#00BFFF', '#FF6B6B', '#4ECDC4', '#FFE66D', '#95E1D3']
    
    def __init__(self, max_seconds: float = 60.0, update_interval: float = 0.1, parent=None):
        # Create with custom axes
        super().__init__(parent, axisItems={
            'bottom': CleanAxis(orientation='bottom'),
            'left': CleanAxis(orientation='left')
        })
        
        self.max_seconds = max_seconds
        self.update_interval = update_interval
        self.max_points = int(max_seconds / update_interval)
        
        self.setBackground('#1a1a1a')
        self.showGrid(x=True, y=True, alpha=0.3)
        self.setLabel('left', 'Value', 'µV')
        self.setLabel('bottom', 'Time', 's')
        self.addLegend()
        
        # Per-channel buffers and curves
        self._channel_data = {}
        self._channel_curves = {}
        
        self.start_time = None
        self._running = False
        
    def _ensure_channel(self, ch_idx: int):
        """Create buffers and curves for a channel if needed"""
        if ch_idx in self._channel_data:
            return
            
        self._channel_data[ch_idx] = {
            'time': CircularBuffer(self.max_points),
            'avg': CircularBuffer(self.max_points),
            'std': CircularBuffer(self.max_points)
        }
        
        color = self.COLORS[ch_idx % len(self.COLORS)]
        # Parse color for fill brush with alpha
        qcolor = pg.mkColor(color)
        fill_color = (qcolor.red(), qcolor.green(), qcolor.blue(), 50)  # 50 alpha for shading
        
        # Create curves for FillBetweenItem (hidden, no pen)
        upper_curve = pg.PlotCurveItem(pen=None)
        lower_curve = pg.PlotCurveItem(pen=None)
        self.addItem(upper_curve)
        self.addItem(lower_curve)
        
        # Create fill between for std deviation shading
        fill_between = pg.FillBetweenItem(upper_curve, lower_curve, brush=fill_color)
        self.addItem(fill_between)
        
        # Mean line (solid, visible in legend)
        mean_curve = self.plot(pen=pg.mkPen(color, width=2), name=f'Ch{ch_idx+1}')
        
        self._channel_curves[ch_idx] = {
            'avg': mean_curve,
            'std_upper': upper_curve,
            'std_lower': lower_curve,
            'fill': fill_between
        }
        
    def start(self):
        self._running = True
        if self.start_time is None:
            self.start_time = time.time()
        
    def stop(self):
        self._running = False
        
    def clear_and_restart(self):
        """Clear all channel data and restart"""
        for ch_data in self._channel_data.values():
            ch_data['time'].clear()
            ch_data['avg'].clear()
            ch_data['std'].clear()
            
        for ch_curves in self._channel_curves.values():
            ch_curves['avg'].setData([], [])
            ch_curves['std_upper'].setData([], [])
            ch_curves['std_lower'].setData([], [])
            
        self.start_time = time.time()
        self._running = True
        
    def set_max_history(self, seconds: float):
        """Update max history duration"""
        self.max_seconds = seconds
        self.max_points = int(seconds / self.update_interval)
        
    def add_point(self, ch_idx: int, mean_uv: float, std_uv: float):
        """Add a statistics point for a specific channel"""
        if not self._running or self.start_time is None:
            return
            
        self._ensure_channel(ch_idx)
        
        t = time.time() - self.start_time
        
        ch_data = self._channel_data[ch_idx]
        ch_data['time'].write(np.array([t], dtype=np.float32))
        ch_data['avg'].write(np.array([mean_uv], dtype=np.float32))
        ch_data['std'].write(np.array([std_uv], dtype=np.float32))
        
        self._update_channel_display(ch_idx)
        
    def add_stats(self, stats: dict):
        """Add statistics for all channels at once"""
        if not self._running or self.start_time is None:
            return
            
        if not stats.get('mean'):
            return
            
        for i, (mean, std) in enumerate(zip(stats['mean'], stats['std'])):
            self.add_point(i, mean * 1e6, std * 1e6)
        
    def _update_channel_display(self, ch_idx: int):
        """Update plot for a specific channel"""
        ch_data = self._channel_data[ch_idx]
        ch_curves = self._channel_curves[ch_idx]
        
        n = ch_data['time'].available
        if n < 2:
            return
            
        t = ch_data['time'].read_latest(n)
        avg = ch_data['avg'].read_latest(n)
        std = ch_data['std'].read_latest(n)
        
        # Update mean line
        ch_curves['avg'].setData(t, avg)
        
        # Update upper/lower curves for fill
        ch_curves['std_upper'].setData(t, avg + std)
        ch_curves['std_lower'].setData(t, avg - std)


class FFTWidget(pg.PlotWidget):
    """FFT spectrum display with multi-channel, cursors, and spectral modes"""
    
    COLORS = ['#00BFFF', '#FF6B6B', '#4ECDC4', '#FFE66D', '#95E1D3']
    
    cursorMoved = Signal(dict)  # Emits cursor positions
    
    def __init__(self, parent=None):
        # Create with custom axes
        super().__init__(parent, axisItems={
            'bottom': CleanAxis(orientation='bottom'),
            'left': CleanAxis(orientation='left')
        })
        
        self.setBackground('#1a1a1a')
        self.showGrid(x=True, y=True, alpha=0.3)
        self.setLabel('left', 'Magnitude', 'dBV')
        self.setLabel('bottom', 'Frequency', 'Hz')
        
        # Multi-channel FFT curves
        self._fft_curves: Dict[int, pg.PlotDataItem] = {}
        
        # Peak markers per channel
        self._peak_scatters: Dict[int, pg.ScatterPlotItem] = {}
        
        # State
        self._x_log = False  # Log scale for axis VIEW (not data)
        self._y_log = False  # Log scale for axis VIEW (not data)
        self._use_db = True  # Show in dB
        self._display_mode = FFTDisplayMode.NORMAL
        self._freq_unit = FrequencyUnit.HZ
        self._show_peaks = True
        self._power_correction = False  # False=Amplitude, True=Power correction
        
        # Store raw data per channel (always in Hz for frequency)
        self._channel_data: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}  # ch_idx -> (freq_hz, mag)
        self._sample_rate = 100000.0  # For spectral density calculation
        self._fft_size = 8192
        
        # Cursor setup
        self._cursors_enabled = False
        self._vline = pg.InfiniteLine(angle=90, movable=False, 
                                      pen=pg.mkPen('#FFAA00', width=1, style=Qt.PenStyle.DashLine))
        self._vline.setVisible(False)
        self.addItem(self._vline, ignoreBounds=True)
        
        # Floating cursor label
        self._cursor_label = pg.TextItem(color='#FFAA00', anchor=(0, 1))
        self._cursor_label.setFont(pg.QtGui.QFont('Consolas', 9))
        self._cursor_label.setVisible(False)
        self.addItem(self._cursor_label, ignoreBounds=True)
        
        # Mouse tracking
        self.setMouseTracking(True)
        self.scene().sigMouseMoved.connect(self._on_mouse_moved)
        
        # Averaging
        self._averaging_enabled = False
        self._avg_count = 1
        self._avg_buffers: Dict[int, List[np.ndarray]] = {}  # ch_idx -> list of spectra
        
    def set_sample_rate(self, rate: float):
        """Set sample rate for spectral density calculations"""
        self._sample_rate = rate
        
    def set_fft_size(self, size: int):
        """Set FFT size for spectral density calculations"""
        self._fft_size = size
        
    def set_x_log(self, log_scale: bool):
        """Set log scale VIEW for X axis (frequency)"""
        self._x_log = log_scale
        self.setLogMode(x=log_scale, y=self._y_log)
        self._update_all_displays()
        self.enableAutoRange()
        
    def set_y_log(self, log_scale: bool):
        """Set log scale VIEW for Y axis"""
        self._y_log = log_scale
        self.setLogMode(x=self._x_log, y=log_scale)
        self._update_all_displays()
        self.enableAutoRange()
        
    def set_use_db(self, use_db: bool):
        """Enable/disable dB display"""
        self._use_db = use_db
        self._update_y_label()
        self._update_all_displays()
        self.enableAutoRange()
        
    def set_display_mode(self, mode: FFTDisplayMode):
        """Set the spectral display mode (normal, spectral density, power, PSD)"""
        self._display_mode = mode
        self._update_y_label()
        self._update_all_displays()
        self.enableAutoRange()
    
    def set_power_correction(self, enabled: bool):
        """Set power correction mode (True=Power, False=Amplitude)"""
        self._power_correction = enabled
        self._update_y_label()
        self._update_all_displays()
        self.enableAutoRange()
        
    def set_freq_unit(self, unit: FrequencyUnit):
        """Set frequency axis units - rescales the frequency axis"""
        self._freq_unit = unit
        self.setLabel('bottom', 'Frequency', unit.label)
        self._update_all_displays()
        # Auto-range to show new data
        self.enableAutoRange()
        
    def _update_y_label(self):
        """Update Y axis label based on display mode and dB setting"""
        mode = self._display_mode
        
        if self._use_db:
            labels = {
                FFTDisplayMode.NORMAL: "dBV",
                FFTDisplayMode.SPECTRAL_DENSITY: "dBV/√Hz",
                FFTDisplayMode.POWER: "dBV²",
                FFTDisplayMode.PSD: "dBV²/Hz"
            }
        else:
            labels = {
                FFTDisplayMode.NORMAL: "V",
                FFTDisplayMode.SPECTRAL_DENSITY: "V/√Hz",
                FFTDisplayMode.POWER: "V²",
                FFTDisplayMode.PSD: "V²/Hz"
            }
        
        self.setLabel('left', 'Magnitude', labels.get(mode, "V"))
        
    def update_spectrum(self, ch_idx: int, freq: np.ndarray, mag: np.ndarray):
        """Update spectrum for a specific channel"""
        # Store raw data (voltage amplitude)
        self._channel_data[ch_idx] = (freq.copy(), mag.copy())
        self._update_channel_display(ch_idx)
        
    def _apply_averaging(self, ch_idx: int, data: np.ndarray) -> np.ndarray:
        """Apply averaging to FFT data"""
        if ch_idx not in self._avg_buffers:
            self._avg_buffers[ch_idx] = []
            
        buf = self._avg_buffers[ch_idx]
        
        if len(buf) == 0 or len(buf[-1]) != len(data):
            buf.clear()
            
        buf.append(data.copy())
        
        while len(buf) > self._avg_count:
            buf.pop(0)
            
        if len(buf) == 1:
            return data
        return np.mean(buf, axis=0).astype(np.float32)
        
    def set_averaging(self, enabled: bool, count: int = 1):
        """Enable/disable FFT averaging"""
        self._averaging_enabled = enabled
        self._avg_count = max(1, count)
        if not enabled:
            self._avg_buffers.clear()
            
    def clear_averaging(self):
        """Clear averaging buffers"""
        self._avg_buffers.clear()
        
    def _update_channel_display(self, ch_idx: int):
        """Update display for a single channel"""
        if ch_idx not in self._channel_data:
            return
            
        freq, mag = self._channel_data[ch_idx]
        
        if len(freq) == 0 or len(mag) == 0:
            return
        
        # Apply averaging if enabled
        if self._averaging_enabled and self._avg_count > 1:
            mag = self._apply_averaging(ch_idx, mag)
        
        # Convert based on display mode
        display_mag = self._convert_magnitude(mag.copy())
        
        # Apply dB if enabled
        if self._use_db:
            # Ensure positive values for log, use a reasonable floor
            display_mag = np.maximum(display_mag, 1e-12)
            display_mag = 20 * np.log10(display_mag)
            # Clip to reasonable dB range to avoid display issues
            display_mag = np.clip(display_mag, -200, 100)
        
        # Scale frequency by selected unit
        display_freq = freq / self._freq_unit.divisor
        
        # Filter for log scale (positive frequencies only)
        if self._x_log:
            mask = display_freq > 0
            display_freq = display_freq[mask]
            display_mag = display_mag[mask]
        
        # Create/update curve
        if ch_idx not in self._fft_curves:
            color = self.COLORS[ch_idx % len(self.COLORS)]
            self._fft_curves[ch_idx] = self.plot(pen=pg.mkPen(color, width=1.5))
            
        self._fft_curves[ch_idx].setData(display_freq, display_mag)
        
        # Update peaks
        if self._show_peaks:
            self._update_peaks(ch_idx, display_freq, display_mag)
            
    def _convert_magnitude(self, mag: np.ndarray) -> np.ndarray:
        """
        Convert magnitude based on display mode and correction type.
        
        Input 'mag' is amplitude-normalized (V) from DisplayManager._compute_fft()
        which already accounts for power_correction in computation.
        
        Conversion logic (matching UHFLI conventions):
        - Amplitude correction: amp = (2/W_sum)|X| is the reference
        - Power correction: psd = (2/(fs*W_sq))|X|² is the reference, amp derived as sqrt(psd*df)
        
        Display modes:
        - NORMAL: V (amplitude)
        - SPECTRAL_DENSITY: V/√Hz (ASD)
        - POWER: V²
        - PSD: V²/Hz
        """
        mode = self._display_mode
        df = self._sample_rate / self._fft_size
        
        if mode == FFTDisplayMode.NORMAL:
            # Amplitude (V) - already computed correctly
            return mag
            
        elif mode == FFTDisplayMode.SPECTRAL_DENSITY:
            # Amplitude Spectral Density (V/√Hz)
            if self._power_correction:
                # Power correction: ASD = sqrt(PSD) = amp / sqrt(df)
                # Since amp = sqrt(psd*df), then ASD = amp/sqrt(df)
                return mag / np.sqrt(df)
            else:
                # Amplitude correction: ASD = amp / sqrt(2*df)
                # Factor of 2 accounts for single-sided spectrum
                return mag / np.sqrt(2.0 * df)
                
        elif mode == FFTDisplayMode.POWER:
            # Power (V²)
            if self._power_correction:
                # Power = PSD * df = amp² (since amp = sqrt(psd*df))
                return mag ** 2
            else:
                # Amplitude correction: Power = amp²
                return mag ** 2
                
        elif mode == FFTDisplayMode.PSD:
            # Power Spectral Density (V²/Hz)
            if self._power_correction:
                # PSD = amp² / df (since amp = sqrt(psd*df), so amp² = psd*df)
                return (mag ** 2) / df
            else:
                # Amplitude correction: PSD = amp² / (2*df)
                return (mag ** 2) / (2.0 * df)
                
        return mag
        
    def _update_all_displays(self):
        """Update all channel displays"""
        for ch_idx in self._channel_data.keys():
            self._update_channel_display(ch_idx)
            
    def _update_peaks(self, ch_idx: int, freq: np.ndarray, display_mag: np.ndarray):
        """Update peak markers for a channel"""
        if len(freq) < 3:
            if ch_idx in self._peak_scatters:
                self._peak_scatters[ch_idx].setData([], [])
            return
            
        # Find peaks
        threshold = np.max(display_mag) - 40 if self._use_db else np.max(display_mag) * 0.01
        
        peaks = []
        for i in range(1, len(display_mag) - 1):
            if display_mag[i] > display_mag[i-1] and display_mag[i] > display_mag[i+1]:
                if display_mag[i] > threshold:
                    peaks.append((freq[i], display_mag[i]))
                    
        peaks.sort(key=lambda x: x[1], reverse=True)
        peaks = peaks[:5]
        
        # Create/update scatter
        if ch_idx not in self._peak_scatters:
            color = self.COLORS[ch_idx % len(self.COLORS)]
            self._peak_scatters[ch_idx] = pg.ScatterPlotItem(
                pen=pg.mkPen(color), brush=pg.mkBrush(color), size=8
            )
            self.addItem(self._peak_scatters[ch_idx])
            
        if peaks:
            pf = np.array([p[0] for p in peaks])
            pm = np.array([p[1] for p in peaks])
            self._peak_scatters[ch_idx].setData(pf, pm)
        else:
            self._peak_scatters[ch_idx].setData([], [])
        
    def toggle_peaks(self, show: bool):
        self._show_peaks = show
        for scatter in self._peak_scatters.values():
            scatter.setVisible(show)
        if show:
            self._update_all_displays()
            
    def clear_all(self):
        """Clear all channel data"""
        self._channel_data.clear()
        self._avg_buffers.clear()
        for curve in self._fft_curves.values():
            curve.setData([], [])
        for scatter in self._peak_scatters.values():
            scatter.setData([], [])
    
    def remove_channel(self, ch_idx: int):
        """Remove a channel's curve and data from the plot"""
        if ch_idx in self._fft_curves:
            self.removeItem(self._fft_curves[ch_idx])
            del self._fft_curves[ch_idx]
        if ch_idx in self._peak_scatters:
            self.removeItem(self._peak_scatters[ch_idx])
            del self._peak_scatters[ch_idx]
        if ch_idx in self._channel_data:
            del self._channel_data[ch_idx]
        if ch_idx in self._avg_buffers:
            del self._avg_buffers[ch_idx]
            
    # === Mouse Cursor Methods ===
    
    def set_cursors_enabled(self, enabled: bool):
        """Enable/disable crosshair cursor"""
        self._cursors_enabled = enabled
        self._vline.setVisible(enabled)
        self._cursor_label.setVisible(enabled)
        if not enabled:
            self.cursorMoved.emit({})
            
    def _on_mouse_moved(self, pos):
        """Handle mouse movement for crosshair"""
        if not self._cursors_enabled:
            return
            
        vb = self.getPlotItem().getViewBox()
        
        if self.sceneBoundingRect().contains(pos):
            mouse_point = vb.mapSceneToView(pos)
            x_view = mouse_point.x()  # View coordinate (may be log10 if log mode)
            y_view = mouse_point.y()  # View coordinate (may be log10 if log mode)
            
            # Convert view coordinate to actual display value
            # When log mode is on, mapSceneToView returns log10 values
            if self._x_log:
                x_display = 10 ** x_view  # Convert log10 back to linear
            else:
                x_display = x_view
            
            if self._y_log:
                y_display = 10 ** y_view
            else:
                y_display = y_view
            
            self._vline.setPos(x_view)  # Position uses view coords
            
            # Convert display frequency back to Hz
            freq_hz = x_display * self._freq_unit.divisor
            
            # Get actual magnitude at this frequency for all channels
            channel_mags = {}
            for ch_idx, (freq, mag) in self._channel_data.items():
                if len(freq) > 0:
                    idx = np.argmin(np.abs(freq - freq_hz))
                    channel_mags[ch_idx] = float(mag[idx])
            
            # Format floating label with frequency
            if freq_hz >= 1e6:
                f_str = f"F: {freq_hz/1e6:.3f} MHz"
            elif freq_hz >= 1e3:
                f_str = f"F: {freq_hz/1e3:.3f} kHz"
            else:
                f_str = f"F: {freq_hz:.1f} Hz"
            
            # Format magnitude based on current display mode and units
            lines = [f_str]
            for ch_idx, mag in sorted(channel_mags.items()):
                # Show in appropriate units based on settings
                if self._use_db:
                    # Convert to dB based on display mode
                    if self._display_mode == FFTDisplayMode.NORMAL:
                        mag_val = 20 * np.log10(mag + 1e-12)
                        lines.append(f"Ch{ch_idx+1}: {mag_val:.1f} dBV")
                    elif self._display_mode == FFTDisplayMode.SPECTRAL_DENSITY:
                        freq_res = self._sample_rate / self._fft_size
                        sd_val = mag / np.sqrt(freq_res)
                        mag_val = 20 * np.log10(sd_val + 1e-12)
                        lines.append(f"Ch{ch_idx+1}: {mag_val:.1f} dBV/√Hz")
                    elif self._display_mode == FFTDisplayMode.POWER:
                        mag_val = 10 * np.log10(mag**2 + 1e-24)
                        lines.append(f"Ch{ch_idx+1}: {mag_val:.1f} dBV²")
                    elif self._display_mode == FFTDisplayMode.PSD:
                        freq_res = self._sample_rate / self._fft_size
                        psd_val = (mag**2) / freq_res
                        mag_val = 10 * np.log10(psd_val + 1e-24)
                        lines.append(f"Ch{ch_idx+1}: {mag_val:.1f} dBV²/Hz")
                else:
                    # Linear units
                    if self._display_mode == FFTDisplayMode.NORMAL:
                        if mag < 1e-3:
                            lines.append(f"Ch{ch_idx+1}: {mag*1e6:.2f} µV")
                        elif mag < 1:
                            lines.append(f"Ch{ch_idx+1}: {mag*1e3:.2f} mV")
                        else:
                            lines.append(f"Ch{ch_idx+1}: {mag:.4f} V")
                    elif self._display_mode == FFTDisplayMode.SPECTRAL_DENSITY:
                        freq_res = self._sample_rate / self._fft_size
                        sd_val = mag / np.sqrt(freq_res)
                        lines.append(f"Ch{ch_idx+1}: {sd_val*1e6:.2f} µV/√Hz")
                    elif self._display_mode == FFTDisplayMode.POWER:
                        lines.append(f"Ch{ch_idx+1}: {mag**2:.2e} V²")
                    elif self._display_mode == FFTDisplayMode.PSD:
                        freq_res = self._sample_rate / self._fft_size
                        psd_val = (mag**2) / freq_res
                        lines.append(f"Ch{ch_idx+1}: {psd_val:.2e} V²/Hz")
            
            self._cursor_label.setText("\n".join(lines))
            
            # Position label near cursor (in view coords)
            view_range = self.viewRange()
            x_range = view_range[0][1] - view_range[0][0]
            y_range = view_range[1][1] - view_range[1][0]
            label_x = x_view + x_range * 0.02
            label_y = y_view + y_range * 0.02
            self._cursor_label.setPos(label_x, label_y)
            
            cursor_data = {
                "freq": freq_hz,
                "freq_display": x_display,
                "mag_display": y_display,
                "channel_mags": channel_mags,
                "unit": self._freq_unit.label
            }
            self.cursorMoved.emit(cursor_data)
            
    def get_cursor_values(self) -> dict:
        return {}


# =============================================================================
# CHANNEL CONFIG WIDGET
# =============================================================================

class ChannelConfigWidget(QWidget):
    """Widget for single channel configuration"""
    
    configChanged = Signal(int, object)
    
    def __init__(self, ch_idx: int, device_info: Optional[DAQDeviceInfo] = None, parent=None):
        super().__init__(parent)
        self.ch_idx = ch_idx
        self.device_info = device_info
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QGridLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        
        colors = ['#00BFFF', '#FF6B6B', '#4ECDC4', '#FFE66D']
        color = colors[self.ch_idx % len(colors)]
        
        # Enable
        self.enable_cb = QCheckBox(f"Ch {self.ch_idx + 1}")
        self.enable_cb.setChecked(True)
        self.enable_cb.setStyleSheet(f"color: {color};")
        self.enable_cb.toggled.connect(self._emit_config)
        layout.addWidget(self.enable_cb, 0, 0)
        
        # Physical channel
        self.phys_combo = QComboBox()
        if self.device_info and self.device_info.ai_channels:
            self.phys_combo.addItems([ch.split('/')[-1] for ch in self.device_info.ai_channels])
        else:
            self.phys_combo.addItems([f"ai{i}" for i in range(16)])
        if self.ch_idx < self.phys_combo.count():
            self.phys_combo.setCurrentIndex(self.ch_idx)
        self.phys_combo.currentTextChanged.connect(self._emit_config)
        layout.addWidget(self.phys_combo, 0, 1)
        
        # Terminal config
        layout.addWidget(QLabel("Term:"), 1, 0)
        self.term_combo = QComboBox()
        if self.device_info and self.device_info.terminal_configs:
            self.term_combo.addItems(self.device_info.terminal_configs)
        else:
            self.term_combo.addItems(["RSE", "NRSE", "DIFF"])
        self.term_combo.currentTextChanged.connect(self._emit_config)
        layout.addWidget(self.term_combo, 1, 1)
        
        # Voltage range
        layout.addWidget(QLabel("Range:"), 2, 0)
        range_layout = QHBoxLayout()
        
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(-100, 100)
        self.min_spin.setValue(-10)
        self.min_spin.setSuffix(" V")
        self.min_spin.valueChanged.connect(self._emit_config)
        range_layout.addWidget(self.min_spin)
        
        self.max_spin = QDoubleSpinBox()
        self.max_spin.setRange(-100, 100)
        self.max_spin.setValue(10)
        self.max_spin.setSuffix(" V")
        self.max_spin.valueChanged.connect(self._emit_config)
        range_layout.addWidget(self.max_spin)
        
        layout.addLayout(range_layout, 2, 1)
        
    def _emit_config(self):
        colors = ['#00BFFF', '#FF6B6B', '#4ECDC4', '#FFE66D']
        config = ChannelConfig(
            physical_channel=self.phys_combo.currentText(),
            enabled=self.enable_cb.isChecked(),
            min_voltage=self.min_spin.value(),
            max_voltage=self.max_spin.value(),
            terminal_config=self.term_combo.currentText(),
            color=colors[self.ch_idx % len(colors)],
            name=f"Channel {self.ch_idx + 1}"
        )
        self.configChanged.emit(self.ch_idx, config)
        
    def get_config(self) -> ChannelConfig:
        colors = ['#00BFFF', '#FF6B6B', '#4ECDC4', '#FFE66D']
        return ChannelConfig(
            physical_channel=self.phys_combo.currentText(),
            enabled=self.enable_cb.isChecked(),
            min_voltage=self.min_spin.value(),
            max_voltage=self.max_spin.value(),
            terminal_config=self.term_combo.currentText(),
            color=colors[self.ch_idx % len(colors)],
            name=f"Channel {self.ch_idx + 1}"
        )
        
    def update_device_info(self, device_info: DAQDeviceInfo):
        """Update with new device info"""
        self.device_info = device_info
        
        current = self.phys_combo.currentText()
        self.phys_combo.clear()
        if device_info and device_info.ai_channels:
            self.phys_combo.addItems([ch.split('/')[-1] for ch in device_info.ai_channels])
        else:
            self.phys_combo.addItems([f"ai{i}" for i in range(16)])
            
        idx = self.phys_combo.findText(current)
        if idx >= 0:
            self.phys_combo.setCurrentIndex(idx)
        elif self.ch_idx < self.phys_combo.count():
            self.phys_combo.setCurrentIndex(self.ch_idx)
            
        current_term = self.term_combo.currentText()
        self.term_combo.clear()
        if device_info and device_info.terminal_configs:
            self.term_combo.addItems(device_info.terminal_configs)
        else:
            self.term_combo.addItems(["RSE", "NRSE", "DIFF"])
            
        idx = self.term_combo.findText(current_term)
        if idx >= 0:
            self.term_combo.setCurrentIndex(idx)
