"""
Monitor Tab Widget for NV Experiment GUI

Real-time monitoring with single-point and continuous acquisition modes.
Allows quick field control and live signal observation.

Features:
- Single-point test acquisition
- Continuous streaming mode
- Magnetic field control (Bx, By, Bz sliders)
- Live signal display
- Connection to main window live plots
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Optional, Dict, Any
import numpy as np

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QSlider, QDoubleSpinBox, QGroupBox,
    QComboBox, QFrame
)
from PySide6.QtCore import Qt, Signal, QTimer

from gui.styles.colors import *


class MonitorTab(QWidget):
    """
    Monitor tab for single-point and continuous acquisition.

    Signals:
        single_point_requested: Emitted when single-point acquisition requested
        continuous_started: Emitted when continuous mode started
        continuous_stopped: Emitted when continuous mode stopped
        field_changed: Emitted when magnetic field changes (Bx, By, Bz)
        data_point_acquired: Emitted when data point is acquired (timestamp, voltage, intensity)
    """

    single_point_requested = Signal()
    continuous_started = Signal()
    continuous_stopped = Signal()
    field_changed = Signal(float, float, float)  # Bx, By, Bz in Gauss
    data_point_acquired = Signal(float, float, float)  # timestamp, voltage, intensity

    def __init__(self, parent=None):
        """Initialize monitor tab."""
        super().__init__(parent)

        # State
        self.continuous_mode = False
        self.acquiring = False

        # Field values (Gauss)
        self.Bx = 0.0
        self.By = 0.0
        self.Bz = 0.0

        # Latest values
        self.latest_voltage = 0.0
        self.latest_intensity = 0.0

        # Create UI
        self._create_ui()

        # Continuous acquisition timer (simulated for now)
        self.acquisition_timer = QTimer()
        self.acquisition_timer.timeout.connect(self._simulate_acquisition)

    def _create_ui(self):
        """Create monitor tab UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(15)

        # Header
        header = QLabel("Real-Time Monitoring")
        header.setProperty("class", "header")
        layout.addWidget(header)

        # Acquisition mode group
        mode_group = self._create_mode_group()
        layout.addWidget(mode_group)

        # Magnetic field control group
        field_group = self._create_field_control_group()
        layout.addWidget(field_group)

        # Current readings group
        readings_group = self._create_readings_group()
        layout.addWidget(readings_group)

        layout.addStretch()

        self.setLayout(layout)

    def _create_mode_group(self) -> QGroupBox:
        """Create acquisition mode control group."""
        group = QGroupBox("Acquisition Mode")
        layout = QVBoxLayout()

        # Mode description
        desc_label = QLabel(
            "Single-Point: Acquire one data point for quick testing\n"
            "Continuous: Stream data continuously to live plots"
        )
        desc_label.setStyleSheet("color: #888888; font-size: 9pt;")
        layout.addWidget(desc_label)

        # Buttons
        button_layout = QHBoxLayout()

        self.single_point_btn = QPushButton("Single Point Test")
        self.single_point_btn.setProperty("class", "primary")
        self.single_point_btn.clicked.connect(self._on_single_point_clicked)
        button_layout.addWidget(self.single_point_btn)

        self.continuous_btn = QPushButton("Start Continuous")
        self.continuous_btn.setProperty("class", "success")
        self.continuous_btn.clicked.connect(self._on_continuous_clicked)
        button_layout.addWidget(self.continuous_btn)

        layout.addLayout(button_layout)

        # Acquisition settings
        settings_layout = QHBoxLayout()

        settings_layout.addWidget(QLabel("Rate:"))

        self.rate_spinbox = QDoubleSpinBox()
        self.rate_spinbox.setRange(0.1, 100.0)
        self.rate_spinbox.setValue(10.0)
        self.rate_spinbox.setSuffix(" Hz")
        self.rate_spinbox.setDecimals(1)
        settings_layout.addWidget(self.rate_spinbox)

        settings_layout.addStretch()

        layout.addLayout(settings_layout)

        group.setLayout(layout)
        return group

    def _create_field_control_group(self) -> QGroupBox:
        """Create magnetic field control group."""
        group = QGroupBox("Magnetic Field Control")
        layout = QGridLayout()
        layout.setSpacing(10)

        # Bx control
        layout.addWidget(QLabel("Bx:"), 0, 0)

        self.bx_slider = QSlider(Qt.Horizontal)
        self.bx_slider.setRange(-1000, 1000)  # -100.0 to +100.0 Gauss (x10)
        self.bx_slider.setValue(0)
        self.bx_slider.setTickPosition(QSlider.TicksBelow)
        self.bx_slider.setTickInterval(200)
        self.bx_slider.valueChanged.connect(self._on_bx_changed)
        layout.addWidget(self.bx_slider, 0, 1)

        self.bx_spinbox = QDoubleSpinBox()
        self.bx_spinbox.setRange(-100.0, 100.0)
        self.bx_spinbox.setValue(0.0)
        self.bx_spinbox.setSuffix(" G")
        self.bx_spinbox.setDecimals(1)
        self.bx_spinbox.valueChanged.connect(self._on_bx_spinbox_changed)
        layout.addWidget(self.bx_spinbox, 0, 2)

        # By control
        layout.addWidget(QLabel("By:"), 1, 0)

        self.by_slider = QSlider(Qt.Horizontal)
        self.by_slider.setRange(-1000, 1000)
        self.by_slider.setValue(0)
        self.by_slider.setTickPosition(QSlider.TicksBelow)
        self.by_slider.setTickInterval(200)
        self.by_slider.valueChanged.connect(self._on_by_changed)
        layout.addWidget(self.by_slider, 1, 1)

        self.by_spinbox = QDoubleSpinBox()
        self.by_spinbox.setRange(-100.0, 100.0)
        self.by_spinbox.setValue(0.0)
        self.by_spinbox.setSuffix(" G")
        self.by_spinbox.setDecimals(1)
        self.by_spinbox.valueChanged.connect(self._on_by_spinbox_changed)
        layout.addWidget(self.by_spinbox, 1, 2)

        # Bz control
        layout.addWidget(QLabel("Bz:"), 2, 0)

        self.bz_slider = QSlider(Qt.Horizontal)
        self.bz_slider.setRange(-1000, 1000)
        self.bz_slider.setValue(0)
        self.bz_slider.setTickPosition(QSlider.TicksBelow)
        self.bz_slider.setTickInterval(200)
        self.bz_slider.valueChanged.connect(self._on_bz_changed)
        layout.addWidget(self.bz_slider, 2, 1)

        self.bz_spinbox = QDoubleSpinBox()
        self.bz_spinbox.setRange(-100.0, 100.0)
        self.bz_spinbox.setValue(0.0)
        self.bz_spinbox.setSuffix(" G")
        self.bz_spinbox.setDecimals(1)
        self.bz_spinbox.valueChanged.connect(self._on_bz_spinbox_changed)
        layout.addWidget(self.bz_spinbox, 2, 2)

        # Apply button
        apply_btn = QPushButton("Apply Field")
        apply_btn.setProperty("class", "primary")
        apply_btn.clicked.connect(self._on_apply_field)
        layout.addWidget(apply_btn, 3, 1, 1, 2)

        group.setLayout(layout)
        return group

    def _create_readings_group(self) -> QGroupBox:
        """Create current readings display group."""
        group = QGroupBox("Current Readings")
        layout = QGridLayout()

        # Voltage reading
        layout.addWidget(QLabel("Voltage:"), 0, 0)
        self.voltage_label = QLabel("0.0000 V")
        self.voltage_label.setProperty("class", "value")
        self.voltage_label.setStyleSheet(f"color: {PLOT_LINE1}; font-size: 14pt; font-weight: bold;")
        layout.addWidget(self.voltage_label, 0, 1)

        # Intensity reading
        layout.addWidget(QLabel("Intensity:"), 1, 0)
        self.intensity_label = QLabel("0 counts")
        self.intensity_label.setProperty("class", "value")
        self.intensity_label.setStyleSheet(f"color: {PLOT_LINE2}; font-size: 14pt; font-weight: bold;")
        layout.addWidget(self.intensity_label, 1, 1)

        # Timestamp
        layout.addWidget(QLabel("Last Update:"), 2, 0)
        self.timestamp_label = QLabel("N/A")
        self.timestamp_label.setStyleSheet("color: #888888;")
        layout.addWidget(self.timestamp_label, 2, 1)

        group.setLayout(layout)
        return group

    def _on_single_point_clicked(self):
        """Handle single-point acquisition button click."""
        if self.continuous_mode:
            return  # Don't allow if continuous is running

        self.single_point_btn.setEnabled(False)
        self.single_point_requested.emit()

        # Re-enable after short delay (simulated)
        QTimer.singleShot(500, lambda: self.single_point_btn.setEnabled(True))

    def _on_continuous_clicked(self):
        """Handle continuous acquisition button click."""
        if not self.continuous_mode:
            # Start continuous
            self.continuous_mode = True
            self.continuous_btn.setText("Stop Continuous")
            self.continuous_btn.setProperty("class", "danger")
            self.continuous_btn.setStyle(self.continuous_btn.style())  # Force refresh
            self.single_point_btn.setEnabled(False)

            # Start acquisition timer
            interval_ms = int(1000.0 / self.rate_spinbox.value())
            self.acquisition_timer.start(interval_ms)

            self.continuous_started.emit()
        else:
            # Stop continuous
            self.continuous_mode = False
            self.continuous_btn.setText("Start Continuous")
            self.continuous_btn.setProperty("class", "success")
            self.continuous_btn.setStyle(self.continuous_btn.style())  # Force refresh
            self.single_point_btn.setEnabled(True)

            # Stop acquisition timer
            self.acquisition_timer.stop()

            self.continuous_stopped.emit()

    def _on_bx_changed(self, value: int):
        """Handle Bx slider change."""
        self.Bx = value / 10.0  # Convert back to Gauss
        self.bx_spinbox.blockSignals(True)
        self.bx_spinbox.setValue(self.Bx)
        self.bx_spinbox.blockSignals(False)

    def _on_bx_spinbox_changed(self, value: float):
        """Handle Bx spinbox change."""
        self.Bx = value
        self.bx_slider.blockSignals(True)
        self.bx_slider.setValue(int(value * 10))
        self.bx_slider.blockSignals(False)

    def _on_by_changed(self, value: int):
        """Handle By slider change."""
        self.By = value / 10.0
        self.by_spinbox.blockSignals(True)
        self.by_spinbox.setValue(self.By)
        self.by_spinbox.blockSignals(False)

    def _on_by_spinbox_changed(self, value: float):
        """Handle By spinbox change."""
        self.By = value
        self.by_slider.blockSignals(True)
        self.by_slider.setValue(int(value * 10))
        self.by_slider.blockSignals(False)

    def _on_bz_changed(self, value: int):
        """Handle Bz slider change."""
        self.Bz = value / 10.0
        self.bz_spinbox.blockSignals(True)
        self.bz_spinbox.setValue(self.Bz)
        self.bz_spinbox.blockSignals(False)

    def _on_bz_spinbox_changed(self, value: float):
        """Handle Bz spinbox change."""
        self.Bz = value
        self.bz_slider.blockSignals(True)
        self.bz_slider.setValue(int(value * 10))
        self.bz_slider.blockSignals(False)

    def _on_apply_field(self):
        """Apply magnetic field values."""
        self.field_changed.emit(self.Bx, self.By, self.Bz)
        print(f"Applied field: Bx={self.Bx:.1f} G, By={self.By:.1f} G, Bz={self.Bz:.1f} G")

    def _simulate_acquisition(self):
        """
        Simulate data acquisition (for testing).
        In real implementation, this would trigger actual hardware acquisition.
        """
        import time

        timestamp = time.time()

        # Simulate voltage and intensity readings
        voltage = 2.5 + 0.3 * np.random.randn()
        intensity = 1000 + 50 * np.random.randn()

        self.update_readings(timestamp, voltage, intensity)
        self.data_point_acquired.emit(timestamp, voltage, intensity)

    def update_readings(self, timestamp: float, voltage: float, intensity: float):
        """
        Update current readings display.

        Args:
            timestamp: Acquisition timestamp
            voltage: Voltage value
            intensity: Intensity value (counts)
        """
        self.latest_voltage = voltage
        self.latest_intensity = intensity

        self.voltage_label.setText(f"{voltage:.4f} V")
        self.intensity_label.setText(f"{int(intensity)} counts")

        import time
        time_str = time.strftime("%H:%M:%S", time.localtime(timestamp))
        self.timestamp_label.setText(time_str)

    def set_acquiring(self, acquiring: bool):
        """
        Set acquisition state (external control).

        Args:
            acquiring: True if acquiring, False otherwise
        """
        self.acquiring = acquiring

        # Disable controls during acquisition
        self.single_point_btn.setEnabled(not acquiring and not self.continuous_mode)
        self.continuous_btn.setEnabled(not acquiring)


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    widget = MonitorTab()

    # Connect signals for testing
    widget.single_point_requested.connect(lambda: print("Single point requested"))
    widget.continuous_started.connect(lambda: print("Continuous started"))
    widget.continuous_stopped.connect(lambda: print("Continuous stopped"))
    widget.field_changed.connect(lambda bx, by, bz: print(f"Field: {bx}, {by}, {bz}"))
    widget.data_point_acquired.connect(
        lambda t, v, i: print(f"Data: t={t:.2f}, v={v:.4f}, i={i:.0f}")
    )

    widget.show()

    sys.exit(app.exec())
