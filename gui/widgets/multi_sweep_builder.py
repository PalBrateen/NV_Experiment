"""
Multi-Parameter Sweep Builder Widget for NV Experiment GUI

Build N-dimensional parameter sweeps with visual preview.

Features:
- Add/remove multiple sweep parameters
- Order parameters (inner loop = first)
- Visual sweep grid preview for 2D sweeps
- Total points calculation
- Estimated time display
- Compatible with ParameterSweep system

Supports:
- 1D sweeps (standard ESR, Rabi, etc.)
- 2D sweeps (frequency × power, tau × frequency, etc.)
- 3D+ sweeps (any combination)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import List, Tuple, Optional, Dict, Any
import numpy as np

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QListWidget, QListWidgetItem, QFrame
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from gui.styles.colors import *


class SweepParameter:
    """Represents one parameter in a multi-parameter sweep."""

    def __init__(self, name: str, start: float, stop: float, points: int, mode: str = 'linear'):
        """
        Initialize sweep parameter.

        Args:
            name: Parameter name (e.g., 'frequency', 'power', 'tau')
            start: Start value
            stop: Stop value
            points: Number of points
            mode: 'linear' or 'log'
        """
        self.name = name
        self.start = start
        self.stop = stop
        self.points = points
        self.mode = mode

    def get_array(self) -> np.ndarray:
        """Generate numpy array for this sweep."""
        if self.mode == 'linear':
            return np.linspace(self.start, self.stop, self.points)
        elif self.mode == 'log':
            if self.start <= 0 or self.stop <= 0:
                raise ValueError("Log sweep requires start and stop > 0")
            return np.logspace(np.log10(self.start), np.log10(self.stop), self.points)
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

    def __repr__(self):
        return f"SweepParameter({self.name}, {self.start:.3e} to {self.stop:.3e}, {self.points} pts, {self.mode})"


class MultiSweepBuilderWidget(QWidget):
    """
    Multi-parameter sweep builder.

    Signals:
        sweep_changed: Emitted when sweep configuration changes (List[SweepParameter])
    """

    sweep_changed = Signal(list)  # List of SweepParameter objects

    # Available parameters for different experiment types
    # AUTOMATE
    PARAMETER_OPTIONS = {
        'ESR': ['frequency', 'power', 'tau'],
        'Rabi': ['tau', 'frequency', 'power'],
        'T1': ['delay', 'frequency', 'power'],
        'T2': ['tau', 'frequency', 'power'],
        'Ramsey': ['tau', 'frequency', 'power', 'detuning']
    }

    def __init__(self, experiment_type: str = 'ESR', parent=None):
        """
        Initialize multi-sweep builder.

        Args:
            experiment_type: Type of experiment (determines available parameters)
            parent: Parent widget
        """
        super().__init__(parent)

        self.experiment_type = experiment_type
        self.sweep_parameters: List[SweepParameter] = []

        self._create_ui()

    def _create_ui(self):
        """Create multi-sweep builder UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(10)

        # Header
        header_layout = QHBoxLayout()
        header_label = QLabel("Multi-Parameter Sweep")
        header_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        header_layout.addWidget(header_label)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        # Add parameter controls
        add_group = QGroupBox("Add Sweep Parameter")
        add_layout = QGridLayout()

        # Parameter selector
        add_layout.addWidget(QLabel("Parameter:"), 0, 0)
        self.param_combo = QComboBox()
        self.param_combo.addItems(self.PARAMETER_OPTIONS.get(self.experiment_type, ['frequency', 'power']))
        self.param_combo.currentTextChanged.connect(self._on_param_selection_changed)
        add_layout.addWidget(self.param_combo, 0, 1)

        # Start value
        add_layout.addWidget(QLabel("Start:"), 1, 0)
        self.start_spinbox = QDoubleSpinBox()
        self.start_spinbox.setRange(-1e12, 1e12)
        self.start_spinbox.setValue(2.85e9)
        self.start_spinbox.setDecimals(6)
        add_layout.addWidget(self.start_spinbox, 1, 1)

        # Stop value
        add_layout.addWidget(QLabel("Stop:"), 2, 0)
        self.stop_spinbox = QDoubleSpinBox()
        self.stop_spinbox.setRange(-1e12, 1e12)
        self.stop_spinbox.setValue(2.89e9)
        self.stop_spinbox.setDecimals(6)
        add_layout.addWidget(self.stop_spinbox, 2, 1)

        # Points
        add_layout.addWidget(QLabel("Points:"), 3, 0)
        self.points_spinbox = QSpinBox()
        self.points_spinbox.setRange(2, 10000)
        self.points_spinbox.setValue(51)
        add_layout.addWidget(self.points_spinbox, 3, 1)

        # Mode
        add_layout.addWidget(QLabel("Mode:"), 4, 0)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Linear", "Logarithmic"])
        add_layout.addWidget(self.mode_combo, 4, 1)

        # Add button
        self.add_btn = QPushButton("Add Parameter")
        self.add_btn.setProperty("class", "primary")
        self.add_btn.clicked.connect(self._add_parameter)
        add_layout.addWidget(self.add_btn, 5, 0, 1, 2)

        add_group.setLayout(add_layout)
        layout.addWidget(add_group)

        # Current parameters list
        params_group = QGroupBox("Sweep Parameters (Inner Loop = First)")
        params_layout = QVBoxLayout()

        self.param_list = QListWidget()
        self.param_list.setMaximumHeight(150)
        self.param_list.setAlternatingRowColors(True)
        params_layout.addWidget(self.param_list)

        # Control buttons
        control_layout = QHBoxLayout()

        self.move_up_btn = QPushButton("↑ Move Up")
        self.move_up_btn.clicked.connect(self._move_up)
        control_layout.addWidget(self.move_up_btn)

        self.move_down_btn = QPushButton("↓ Move Down")
        self.move_down_btn.clicked.connect(self._move_down)
        control_layout.addWidget(self.move_down_btn)

        self.remove_btn = QPushButton("Remove")
        self.remove_btn.setProperty("class", "danger")
        self.remove_btn.clicked.connect(self._remove_parameter)
        control_layout.addWidget(self.remove_btn)

        params_layout.addLayout(control_layout)
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)

        # Statistics
        stats_group = QGroupBox("Sweep Statistics")
        stats_layout = QGridLayout()

        stats_layout.addWidget(QLabel("Dimensions:"), 0, 0)
        self.dim_label = QLabel("0D (Single Point)")
        self.dim_label.setStyleSheet(f"color: {STATUS_INFO}; font-weight: bold;")
        stats_layout.addWidget(self.dim_label, 0, 1)

        stats_layout.addWidget(QLabel("Total Points:"), 1, 0)
        self.total_points_label = QLabel("1")
        self.total_points_label.setStyleSheet(f"color: {STATUS_INFO}; font-weight: bold;")
        stats_layout.addWidget(self.total_points_label, 1, 1)

        stats_layout.addWidget(QLabel("Est. Time (1s/pt):"), 2, 0)
        self.est_time_label = QLabel("1 s")
        self.est_time_label.setStyleSheet(f"color: {STATUS_INFO};")
        stats_layout.addWidget(self.est_time_label, 2, 1)

        stats_group.setLayout(stats_layout)
        layout.addWidget(stats_group)

        # 2D sweep visualization (shown for 2D only)
        self.viz_group = QGroupBox("2D Sweep Visualization")
        viz_layout = QVBoxLayout()

        self.viz_label = QLabel("Grid: X × Y points")
        self.viz_label.setAlignment(Qt.AlignCenter)
        self.viz_label.setStyleSheet("font-family: monospace; font-size: 9pt;")
        viz_layout.addWidget(self.viz_label)

        self.viz_group.setLayout(viz_layout)
        self.viz_group.setVisible(False)
        layout.addWidget(self.viz_group)

        layout.addStretch()

        self.setLayout(layout)

        # Update parameter suggestions based on experiment type
        self._on_param_selection_changed(self.param_combo.currentText())

    def _on_param_selection_changed(self, param_name: str):
        """Update spinbox defaults based on selected parameter."""
        # Set sensible defaults for different parameters
        defaults = {
            'frequency': (2.85e9, 2.89e9, 51),
            'power': (0, 20, 11),  # dBm
            'tau': (0, 5000e-9, 51),  # 0-5000 ns
            'delay': (0, 100e-6, 51),  # 0-100 μs
            'detuning': (-5e6, 5e6, 21)  # ±5 MHz
        }

        if param_name in defaults:
            start, stop, points = defaults[param_name]
            self.start_spinbox.setValue(start)
            self.stop_spinbox.setValue(stop)
            self.points_spinbox.setValue(points)

    def _add_parameter(self):
        """Add parameter to sweep."""
        param_name = self.param_combo.currentText()
        start = self.start_spinbox.value()
        stop = self.stop_spinbox.value()
        points = self.points_spinbox.value()
        mode = self.mode_combo.currentText().lower()

        # Check if parameter already exists
        if any(p.name == param_name for p in self.sweep_parameters):
            return  # Already added

        # Create sweep parameter
        param = SweepParameter(param_name, start, stop, points, mode)
        self.sweep_parameters.append(param)

        # Update list
        self._update_parameter_list()
        self._update_statistics()

        # Emit signal
        self.sweep_changed.emit(self.sweep_parameters.copy())

    def _remove_parameter(self):
        """Remove selected parameter from sweep."""
        current_item = self.param_list.currentItem()
        if current_item is None:
            return

        index = self.param_list.row(current_item)
        if 0 <= index < len(self.sweep_parameters):
            del self.sweep_parameters[index]
            self._update_parameter_list()
            self._update_statistics()
            self.sweep_changed.emit(self.sweep_parameters.copy())

    def _move_up(self):
        """Move selected parameter up (towards inner loop)."""
        current_item = self.param_list.currentItem()
        if current_item is None:
            return

        index = self.param_list.row(current_item)
        if index > 0:
            self.sweep_parameters[index], self.sweep_parameters[index-1] = \
                self.sweep_parameters[index-1], self.sweep_parameters[index]
            self._update_parameter_list()
            self.param_list.setCurrentRow(index - 1)
            self.sweep_changed.emit(self.sweep_parameters.copy())

    def _move_down(self):
        """Move selected parameter down (towards outer loop)."""
        current_item = self.param_list.currentItem()
        if current_item is None:
            return

        index = self.param_list.row(current_item)
        if index < len(self.sweep_parameters) - 1:
            self.sweep_parameters[index], self.sweep_parameters[index+1] = \
                self.sweep_parameters[index+1], self.sweep_parameters[index]
            self._update_parameter_list()
            self.param_list.setCurrentRow(index + 1)
            self.sweep_changed.emit(self.sweep_parameters.copy())

    def _update_parameter_list(self):
        """Update parameter list display."""
        self.param_list.clear()

        for i, param in enumerate(self.sweep_parameters):
            level = "INNER LOOP" if i == 0 else f"Level {i+1}"
            text = f"{param.name}: {param.start:.3e} to {param.stop:.3e} ({param.points} pts, {param.mode}) [{level}]"
            item = QListWidgetItem(text)

            # Color code by level
            if i == 0:
                item.setForeground(QColor(STATUS_OK))  # Green for inner loop
            elif i == 1:
                item.setForeground(QColor(STATUS_WARNING))  # Yellow for 2nd level

            self.param_list.addItem(item)

    def _update_statistics(self):
        """Update statistics display."""
        n_params = len(self.sweep_parameters)

        # Dimension
        if n_params == 0:
            self.dim_label.setText("0D (Single Point)")
            total_points = 1
        elif n_params == 1:
            self.dim_label.setText("1D Sweep")
            total_points = self.sweep_parameters[0].points
        elif n_params == 2:
            self.dim_label.setText("2D Sweep")
            total_points = self.sweep_parameters[0].points * self.sweep_parameters[1].points
            self._update_2d_visualization()
        else:
            self.dim_label.setText(f"{n_params}D Sweep")
            total_points = 1
            for param in self.sweep_parameters:
                total_points *= param.points

        # Total points
        self.total_points_label.setText(f"{total_points:,}")

        # Estimated time
        est_seconds = total_points * 1.0
        if est_seconds < 60:
            time_str = f"{est_seconds:.0f} s"
        elif est_seconds < 3600:
            time_str = f"{est_seconds / 60:.1f} min"
        else:
            time_str = f"{est_seconds / 3600:.2f} hr"
        self.est_time_label.setText(time_str)

        # Show/hide 2D visualization
        self.viz_group.setVisible(n_params == 2)

    # REMOVE ??
    def _update_2d_visualization(self):
        """Update 2D sweep grid visualization."""
        if len(self.sweep_parameters) != 2:
            return

        inner = self.sweep_parameters[0]
        outer = self.sweep_parameters[1]

        viz_text = f"Grid: {inner.name} ({inner.points} pts) × {outer.name} ({outer.points} pts)\n"
        viz_text += f"Total: {inner.points * outer.points} points\n\n"
        viz_text += f"Inner loop (fast): {inner.name} sweeps {inner.points}× for each {outer.name}\n"
        viz_text += f"Outer loop (slow): {outer.name} changes {outer.points}×"

        self.viz_label.setText(viz_text)

    def get_sweep_parameters(self) -> List[SweepParameter]:
        """
        Get current sweep parameters.

        Returns:
            List of SweepParameter objects (first = inner loop)
        """
        return self.sweep_parameters.copy()

    def set_experiment_type(self, exp_type: str):
        """
        Set experiment type (updates available parameters).

        Args:
            exp_type: Experiment type (ESR, Rabi, T1, T2, Ramsey)
        """
        self.experiment_type = exp_type
        self.param_combo.clear()
        self.param_combo.addItems(self.PARAMETER_OPTIONS.get(exp_type, ['frequency', 'power']))

    def clear_parameters(self):
        """Clear all sweep parameters."""
        self.sweep_parameters.clear()
        self._update_parameter_list()
        self._update_statistics()
        self.sweep_changed.emit(self.sweep_parameters.copy())

    def set_single_parameter(self, param_name: str, start: float, stop: float, points: int, mode: str = 'linear'):
        """
        Set as single-parameter sweep (convenience method).

        Args:
            param_name: Parameter name
            start: Start value
            stop: Stop value
            points: Number of points
            mode: 'linear' or 'log'
        """
        self.clear_parameters()
        param = SweepParameter(param_name, start, stop, points, mode)
        self.sweep_parameters.append(param)
        self._update_parameter_list()
        self._update_statistics()
        self.sweep_changed.emit(self.sweep_parameters.copy())


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    widget = MultiSweepBuilderWidget(experiment_type='ESR')

    def on_sweep_changed(params):
        print(f"\nSweep changed: {len(params)} parameters")
        for i, param in enumerate(params):
            print(f"  [{i}] {param}")

    widget.sweep_changed.connect(on_sweep_changed)

    widget.show()

    sys.exit(app.exec())
