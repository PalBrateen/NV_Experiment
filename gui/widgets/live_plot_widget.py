"""
Live Plot Widget for NV Experiment GUI

Real-time plotting with pyqtgraph for voltage and intensity signals.
Supports both timeseries and FFT visualization modes.

Features:
- Timeseries mode: Rolling buffer display
- FFT mode: Frequency spectrum with peak detection
- Statistics display (mean, std, min, max)
- Auto-scaling or manual axes
- Export to PNG
- Pause/resume updates
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from collections import deque
from typing import Optional, Tuple

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QSpinBox, QRadioButton, QButtonGroup,
    QGroupBox, QFileDialog
)
from PySide6.QtCore import Qt, Signal, QTimer

import pyqtgraph as pg
from pyqtgraph import PlotWidget, PlotItem

from gui.utils.fft_processor import FFTProcessor, WindowFunction
from gui.styles.colors import *


class LivePlotWidget(QWidget):
    """
    Live plotting widget with timeseries and FFT modes.

    Signals:
        mode_changed: Emitted when switching between timeseries/FFT
    """

    mode_changed = Signal(str)  # 'timeseries' or 'fft'

    def __init__(
        self,
        title: str = "Live Plot",
        y_label: str = "Signal",
        y_unit: str = "V",
        buffer_size: int = 1000,
        parent=None
    ):
        """
        Initialize live plot widget.

        Args:
            title: Plot title
            y_label: Y-axis label
            y_unit: Y-axis unit
            buffer_size: Number of points to keep in rolling buffer
            parent: Parent widget
        """
        super().__init__(parent)

        self.title = title
        self.y_label = y_label
        self.y_unit = y_unit
        self.buffer_size = buffer_size

        # Data buffers
        self.time_buffer = deque(maxlen=buffer_size)
        self.data_buffer = deque(maxlen=buffer_size)

        # FFT processor
        self.fft_processor = FFTProcessor(use_gpu=False)

        # State
        self.mode = 'timeseries'  # 'timeseries' or 'fft'
        self.paused = False
        self.sampling_rate = 100.0  # Hz (default)
        self.update_rate = 10  # Hz

        # Statistics
        self.stats = {
            'current': 0.0,
            'mean': 0.0,
            'std': 0.0,
            'min': 0.0,
            'max': 0.0
        }

        # Setup UI
        self._create_ui()
        self._configure_plot()

        # Update timer
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self._update_plot)
        # Timer will be started when data arrives

    def _create_ui(self):
        """Create widget UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)

        # Header with mode selection
        header_layout = QHBoxLayout()

        title_label = QLabel(self.title)
        title_label.setProperty("class", "subheader")
        header_layout.addWidget(title_label)

        header_layout.addStretch()

        # Mode selection
        mode_label = QLabel("Display:")
        header_layout.addWidget(mode_label)

        self.mode_group = QButtonGroup(self)

        self.timeseries_radio = QRadioButton("Timeseries")
        self.timeseries_radio.setChecked(True)
        self.timeseries_radio.toggled.connect(self._on_mode_changed)
        self.mode_group.addButton(self.timeseries_radio)
        header_layout.addWidget(self.timeseries_radio)

        self.fft_radio = QRadioButton("FFT")
        self.fft_radio.toggled.connect(self._on_mode_changed)
        self.mode_group.addButton(self.fft_radio)
        header_layout.addWidget(self.fft_radio)

        layout.addLayout(header_layout)

        # Plot widget
        self.plot_widget = PlotWidget()
        layout.addWidget(self.plot_widget)

        # Statistics and controls
        controls_layout = QHBoxLayout()

        # Statistics display
        self.stats_label = QLabel("Current: 0.00 V | Mean: 0.00 V | Std: 0.00 V")
        self.stats_label.setProperty("class", "value")
        controls_layout.addWidget(self.stats_label)

        controls_layout.addStretch()

        # Buffer size control
        buffer_label = QLabel("Buffer:")
        controls_layout.addWidget(buffer_label)

        self.buffer_spinbox = QSpinBox()
        self.buffer_spinbox.setRange(100, 10000)
        self.buffer_spinbox.setValue(self.buffer_size)
        self.buffer_spinbox.setSuffix(" pts")
        self.buffer_spinbox.valueChanged.connect(self._on_buffer_size_changed)
        controls_layout.addWidget(self.buffer_spinbox)

        # Update rate control
        rate_label = QLabel("Update:")
        controls_layout.addWidget(rate_label)

        self.rate_spinbox = QSpinBox()
        self.rate_spinbox.setRange(1, 60)
        self.rate_spinbox.setValue(self.update_rate)
        self.rate_spinbox.setSuffix(" Hz")
        self.rate_spinbox.valueChanged.connect(self._on_update_rate_changed)
        controls_layout.addWidget(self.rate_spinbox)

        layout.addLayout(controls_layout)

        # Action buttons
        button_layout = QHBoxLayout()

        self.pause_btn = QPushButton("Pause")
        self.pause_btn.clicked.connect(self.toggle_pause)
        button_layout.addWidget(self.pause_btn)

        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self.clear_data)
        button_layout.addWidget(self.clear_btn)

        self.export_btn = QPushButton("Export PNG")
        self.export_btn.clicked.connect(self.export_plot)
        button_layout.addWidget(self.export_btn)

        button_layout.addStretch()

        layout.addLayout(button_layout)

        self.setLayout(layout)

    def _configure_plot(self):
        """Configure pyqtgraph plot appearance."""
        # Set background color
        self.plot_widget.setBackground(BG_DARK)

        # Get plot item
        plot_item = self.plot_widget.getPlotItem()

        # Set labels
        plot_item.setLabel('bottom', 'Time', units='s')
        plot_item.setLabel('left', self.y_label, units=self.y_unit)

        # Grid
        plot_item.showGrid(x=True, y=True, alpha=0.3)

        # Style axes
        axis_pen = pg.mkPen(color=PLOT_AXIS, width=1)
        plot_item.getAxis('bottom').setPen(axis_pen)
        plot_item.getAxis('left').setPen(axis_pen)

        # Style text
        text_color = TEXT_PRIMARY
        for axis in ['bottom', 'left']:
            plot_item.getAxis(axis).setTextPen(text_color)

        # Create plot curve
        self.curve = plot_item.plot(
            pen=pg.mkPen(color=PLOT_LINE1, width=2),
            name=self.title
        )

        # Auto-range
        plot_item.enableAutoRange()

    def add_data_point(self, timestamp: float, value: float):
        """
        Add a single data point to the plot.

        Args:
            timestamp: Time in seconds
            value: Signal value
        """
        if self.paused:
            return

        # Add to buffers
        self.time_buffer.append(timestamp)
        self.data_buffer.append(value)

        # Update statistics
        self._update_statistics()

        # Start update timer if not running
        if not self.update_timer.isActive():
            interval_ms = int(1000 / self.update_rate)
            self.update_timer.start(interval_ms)

    def add_data_batch(self, timestamps: np.ndarray, values: np.ndarray):
        """
        Add multiple data points at once.

        Args:
            timestamps: Array of timestamps
            values: Array of signal values
        """
        if self.paused:
            return

        # Add to buffers
        for t, v in zip(timestamps, values):
            self.time_buffer.append(t)
            self.data_buffer.append(v)

        # Update statistics
        self._update_statistics()

    def _update_plot(self):
        """Update plot with current data (called by timer)."""
        if len(self.data_buffer) == 0:
            return

        if self.mode == 'timeseries':
            self._plot_timeseries()
        elif self.mode == 'fft':
            self._plot_fft()

        # Update statistics label
        self._update_stats_label()

    def _plot_timeseries(self):
        """Plot timeseries data."""
        if len(self.time_buffer) == 0:
            return

        times = np.array(self.time_buffer)
        values = np.array(self.data_buffer)

        # Relative time (start from 0)
        times = times - times[0]

        self.curve.setData(times, values)

    def _plot_fft(self):
        """Plot FFT spectrum."""
        if len(self.data_buffer) < 10:
            return  # Not enough data for FFT

        data = np.array(self.data_buffer)

        # Compute FFT
        freqs, mags = self.fft_processor.compute_fft(
            data,
            self.sampling_rate,
            window=WindowFunction.HANNING,
            real_only=True,
            normalize=True
        )

        # Plot
        self.curve.setData(freqs, mags)

    def _update_statistics(self):
        """Calculate statistics from current buffer."""
        if len(self.data_buffer) == 0:
            return

        data = np.array(self.data_buffer)

        self.stats['current'] = data[-1]
        self.stats['mean'] = np.mean(data)
        self.stats['std'] = np.std(data)
        self.stats['min'] = np.min(data)
        self.stats['max'] = np.max(data)

    def _update_stats_label(self):
        """Update statistics label text."""
        text = (
            f"Current: {self.stats['current']:.4f} {self.y_unit} | "
            f"Mean: {self.stats['mean']:.4f} {self.y_unit} | "
            f"Std: {self.stats['std']:.4f} {self.y_unit} | "
            f"Min: {self.stats['min']:.4f} {self.y_unit} | "
            f"Max: {self.stats['max']:.4f} {self.y_unit}"
        )
        self.stats_label.setText(text)

    def _on_mode_changed(self):
        """Handle mode change (timeseries <-> FFT)."""
        if self.timeseries_radio.isChecked():
            self.mode = 'timeseries'
            # Update axis labels
            plot_item = self.plot_widget.getPlotItem()
            plot_item.setLabel('bottom', 'Time', units='s')
            plot_item.setLabel('left', self.y_label, units=self.y_unit)
        else:
            self.mode = 'fft'
            # Update axis labels
            plot_item = self.plot_widget.getPlotItem()
            plot_item.setLabel('bottom', 'Frequency', units='Hz')
            plot_item.setLabel('left', 'Magnitude', units='')

        self.mode_changed.emit(self.mode)
        self._update_plot()

    def _on_buffer_size_changed(self, new_size: int):
        """Handle buffer size change."""
        self.buffer_size = new_size
        # Recreate buffers with new size
        old_time = list(self.time_buffer)
        old_data = list(self.data_buffer)

        self.time_buffer = deque(old_time, maxlen=new_size)
        self.data_buffer = deque(old_data, maxlen=new_size)

    def _on_update_rate_changed(self, new_rate: int):
        """Handle update rate change."""
        self.update_rate = new_rate
        if self.update_timer.isActive():
            interval_ms = int(1000 / self.update_rate)
            self.update_timer.setInterval(interval_ms)

    def toggle_pause(self):
        """Toggle pause state."""
        self.paused = not self.paused
        if self.paused:
            self.pause_btn.setText("Resume")
            self.update_timer.stop()
        else:
            self.pause_btn.setText("Pause")
            interval_ms = int(1000 / self.update_rate)
            self.update_timer.start(interval_ms)

    def clear_data(self):
        """Clear all data from buffers."""
        self.time_buffer.clear()
        self.data_buffer.clear()
        self.curve.setData([], [])
        self._update_statistics()
        self._update_stats_label()

    def export_plot(self):
        """Export current plot to PNG."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Plot",
            f"{self.title.replace(' ', '_')}.png",
            "PNG Images (*.png)"
        )

        if file_path:
            exporter = pg.exporters.ImageExporter(self.plot_widget.plotItem)
            exporter.export(file_path)

    def set_sampling_rate(self, rate: float):
        """
        Set sampling rate for FFT calculation.

        Args:
            rate: Sampling rate in Hz
        """
        self.sampling_rate = rate

    def get_current_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get current buffered data.

        Returns:
            (times, values) as numpy arrays
        """
        return np.array(self.time_buffer), np.array(self.data_buffer)


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    # Create widget
    widget = LivePlotWidget(title="Test Signal", y_label="Voltage", y_unit="V")

    # Simulate data
    def add_test_data():
        import time
        t = time.time()
        value = np.sin(2 * np.pi * 1 * t) + 0.1 * np.random.randn()
        widget.add_data_point(t, value)

    timer = QTimer()
    timer.timeout.connect(add_test_data)
    timer.start(50)  # 20 Hz

    widget.show()
    sys.exit(app.exec())
