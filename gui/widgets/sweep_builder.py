"""
Sweep Builder Widget for NV Experiment GUI

Visual parameter sweep configuration with multiple modes and real-time preview.

Features:
- Linear sweep mode
- Logarithmic sweep mode
- Manual array input
- Drag-drop parameter reordering (future)
- Real-time point count and time estimation
- Preview of sweep array
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import List, Tuple, Optional
import numpy as np

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QTextEdit, QRadioButton, QButtonGroup
)
from PySide6.QtCore import Qt, Signal

from gui.styles.colors import *


class SweepBuilderWidget(QWidget):
    """
    Sweep parameter builder with multiple modes.

    Signals:
        sweep_changed: Emitted when sweep configuration changes (sweep_array)
    """

    sweep_changed = Signal(np.ndarray)  # Sweep array

    def __init__(self, parent=None):
        """Initialize sweep builder."""
        super().__init__(parent)

        # State
        self.sweep_mode = 'linear'  # 'linear', 'log', 'manual'
        self.sweep_array = np.array([])

        # Create UI
        self._create_ui()
        self._update_sweep()

    def _create_ui(self):
        """Create sweep builder UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(10)

        # Mode selection
        mode_group = QGroupBox("Sweep Mode")
        mode_layout = QHBoxLayout()

        self.mode_button_group = QButtonGroup(self)

        self.linear_radio = QRadioButton("Linear")
        self.linear_radio.setChecked(True)
        self.linear_radio.toggled.connect(self._on_mode_changed)
        self.mode_button_group.addButton(self.linear_radio)
        mode_layout.addWidget(self.linear_radio)

        self.log_radio = QRadioButton("Logarithmic")
        self.log_radio.toggled.connect(self._on_mode_changed)
        self.mode_button_group.addButton(self.log_radio)
        mode_layout.addWidget(self.log_radio)

        self.manual_radio = QRadioButton("Manual")
        self.manual_radio.toggled.connect(self._on_mode_changed)
        self.mode_button_group.addButton(self.manual_radio)
        mode_layout.addWidget(self.manual_radio)

        mode_layout.addStretch()

        mode_group.setLayout(mode_layout)
        layout.addWidget(mode_group)

        # Linear/Log parameters
        self.params_group = QGroupBox("Parameters")
        params_layout = QGridLayout()

        # Start
        params_layout.addWidget(QLabel("Start:"), 0, 0)
        self.start_spinbox = QDoubleSpinBox()
        self.start_spinbox.setRange(-1e12, 1e12)
        self.start_spinbox.setValue(2.85e9)
        self.start_spinbox.setDecimals(6)
        self.start_spinbox.valueChanged.connect(self._update_sweep)
        params_layout.addWidget(self.start_spinbox, 0, 1)

        # Stop
        params_layout.addWidget(QLabel("Stop:"), 1, 0)
        self.stop_spinbox = QDoubleSpinBox()
        self.stop_spinbox.setRange(-1e12, 1e12)
        self.stop_spinbox.setValue(2.89e9)
        self.stop_spinbox.setDecimals(6)
        self.stop_spinbox.valueChanged.connect(self._update_sweep)
        params_layout.addWidget(self.stop_spinbox, 1, 1)

        # Points
        params_layout.addWidget(QLabel("Points:"), 2, 0)
        self.points_spinbox = QSpinBox()
        self.points_spinbox.setRange(2, 10000)
        self.points_spinbox.setValue(51)
        self.points_spinbox.valueChanged.connect(self._update_sweep)
        params_layout.addWidget(self.points_spinbox, 2, 1)

        # Step size (calculated)
        params_layout.addWidget(QLabel("Step:"), 3, 0)
        self.step_label = QLabel("0.0")
        self.step_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        params_layout.addWidget(self.step_label, 3, 1)

        self.params_group.setLayout(params_layout)
        layout.addWidget(self.params_group)

        # Manual input
        self.manual_group = QGroupBox("Manual Array")
        manual_layout = QVBoxLayout()

        manual_label = QLabel("Enter values (comma or space separated):")
        manual_label.setStyleSheet("font-size: 9pt;")
        manual_layout.addWidget(manual_label)

        self.manual_text = QTextEdit()
        self.manual_text.setPlaceholderText("e.g., 2.85e9, 2.86e9, 2.87e9, 2.88e9, 2.89e9")
        self.manual_text.setMaximumHeight(100)
        self.manual_text.textChanged.connect(self._update_sweep)
        manual_layout.addWidget(self.manual_text)

        self.manual_group.setLayout(manual_layout)
        self.manual_group.setVisible(False)
        layout.addWidget(self.manual_group)

        # Preview
        preview_group = QGroupBox("Preview")
        preview_layout = QVBoxLayout()

        self.preview_label = QLabel()
        self.preview_label.setStyleSheet(f"color: {TEXT_SECONDARY}; font-family: 'Consolas', monospace;")
        self.preview_label.setWordWrap(True)
        preview_layout.addWidget(self.preview_label)

        preview_group.setLayout(preview_layout)
        layout.addWidget(preview_group)

        # Statistics
        stats_group = QGroupBox("Statistics")
        stats_layout = QGridLayout()

        stats_layout.addWidget(QLabel("Total Points:"), 0, 0)
        self.total_points_label = QLabel("0")
        self.total_points_label.setStyleSheet(f"color: {STATUS_INFO}; font-weight: bold;")
        stats_layout.addWidget(self.total_points_label, 0, 1)

        stats_layout.addWidget(QLabel("Est. Time (1s/pt):"), 1, 0)
        self.est_time_label = QLabel("0 s")
        self.est_time_label.setStyleSheet(f"color: {STATUS_INFO};")
        stats_layout.addWidget(self.est_time_label, 1, 1)

        stats_group.setLayout(stats_layout)
        layout.addWidget(stats_group)

        layout.addStretch()

        self.setLayout(layout)

    def _on_mode_changed(self):
        """Handle sweep mode change."""
        if self.linear_radio.isChecked():
            self.sweep_mode = 'linear'
            self.params_group.setVisible(True)
            self.manual_group.setVisible(False)
        elif self.log_radio.isChecked():
            self.sweep_mode = 'log'
            self.params_group.setVisible(True)
            self.manual_group.setVisible(False)
        elif self.manual_radio.isChecked():
            self.sweep_mode = 'manual'
            self.params_group.setVisible(False)
            self.manual_group.setVisible(True)

        self._update_sweep()

    def _update_sweep(self):
        """Update sweep array based on current settings."""
        try:
            if self.sweep_mode == 'linear':
                start = self.start_spinbox.value()
                stop = self.stop_spinbox.value()
                points = self.points_spinbox.value()

                self.sweep_array = np.linspace(start, stop, points)

                # Calculate step
                if points > 1:
                    step = (stop - start) / (points - 1)
                    self.step_label.setText(f"{step:.6e}")
                else:
                    self.step_label.setText("N/A")

            elif self.sweep_mode == 'log':
                start = self.start_spinbox.value()
                stop = self.stop_spinbox.value()
                points = self.points_spinbox.value()

                if start <= 0 or stop <= 0:
                    self.sweep_array = np.array([])
                    self.step_label.setText("Error: Start and Stop must be > 0 for log")
                else:
                    self.sweep_array = np.logspace(
                        np.log10(start),
                        np.log10(stop),
                        points
                    )
                    self.step_label.setText("Logarithmic")

            elif self.sweep_mode == 'manual':
                text = self.manual_text.toPlainText()
                if text.strip():
                    # Parse comma or space separated values
                    text = text.replace(',', ' ')
                    values = [float(x) for x in text.split() if x.strip()]
                    self.sweep_array = np.array(values)
                else:
                    self.sweep_array = np.array([])

            # Update preview
            self._update_preview()

            # Update statistics
            self._update_statistics()

            # Emit signal
            self.sweep_changed.emit(self.sweep_array)

        except Exception as e:
            self.preview_label.setText(f"Error: {e}")
            self.sweep_array = np.array([])

    def _update_preview(self):
        """Update preview text."""
        if len(self.sweep_array) == 0:
            self.preview_label.setText("No sweep defined")
            return

        n = len(self.sweep_array)

        if n <= 10:
            # Show all values
            preview_text = ", ".join([f"{x:.6e}" for x in self.sweep_array])
        else:
            # Show first 5 and last 5
            first_5 = ", ".join([f"{x:.6e}" for x in self.sweep_array[:5]])
            last_5 = ", ".join([f"{x:.6e}" for x in self.sweep_array[-5:]])
            preview_text = f"{first_5} ... {last_5}"

        self.preview_label.setText(preview_text)

    def _update_statistics(self):
        """Update statistics display."""
        n_points = len(self.sweep_array)
        self.total_points_label.setText(str(n_points))

        # Estimate time (assuming 1 second per point)
        est_seconds = n_points * 1.0

        if est_seconds < 60:
            time_str = f"{est_seconds:.0f} s"
        elif est_seconds < 3600:
            time_str = f"{est_seconds / 60:.1f} min"
        else:
            time_str = f"{est_seconds / 3600:.2f} hr"

        self.est_time_label.setText(time_str)

    def get_sweep_array(self) -> np.ndarray:
        """
        Get current sweep array.

        Returns:
            Numpy array of sweep values
        """
        return self.sweep_array.copy()

    def set_sweep_array(self, array: np.ndarray):
        """
        Set sweep array (switches to manual mode).

        Args:
            array: Sweep values
        """
        self.manual_radio.setChecked(True)
        text = ", ".join([f"{x:.6e}" for x in array])
        self.manual_text.setPlainText(text)

    def set_linear_sweep(self, start: float, stop: float, points: int):
        """
        Configure linear sweep.

        Args:
            start: Start value
            stop: Stop value
            points: Number of points
        """
        self.linear_radio.setChecked(True)
        self.start_spinbox.setValue(start)
        self.stop_spinbox.setValue(stop)
        self.points_spinbox.setValue(points)


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    widget = SweepBuilderWidget()

    # Connect signal
    def on_sweep_changed(array):
        print(f"Sweep changed: {len(array)} points")
        if len(array) > 0:
            print(f"  Range: {array[0]:.3e} to {array[-1]:.3e}")

    widget.sweep_changed.connect(on_sweep_changed)

    widget.show()

    sys.exit(app.exec())
