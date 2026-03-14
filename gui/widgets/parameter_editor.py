"""
Parameter Editor Widget for NV Experiment GUI

Unit-aware parameter editor with validation, templates, and save/load functionality.

Features:
- Unit-aware spinboxes (frequency, time, power, voltage)
- Automatic unit conversion
- Range validation
- Parameter templates (ESR, Rabi, T1, T2)
- Load/save configuration files
- Real-time validation feedback
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Dict, Any, Optional, Tuple, List
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QScrollArea, QFrame, QLineEdit, QFileDialog,
    QMessageBox
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QPalette

from gui.utils.unit_converter import UnitConverter, UnitType
from gui.styles.colors import *
import yaml
import json


class UnitAwareSpinBox(QWidget):
    """
    Spinbox with unit selection dropdown.

    Automatically converts values between units while maintaining
    the underlying base unit value.

    Signals:
        value_changed: Emitted when value changes (in base units)
    """

    value_changed = Signal(float)  # Value in base units

    def __init__(
        self,
        unit_type: UnitType,
        initial_value: float = 0.0,
        initial_unit: str = None,
        min_value: float = 0.0,
        max_value: float = 1e12,
        decimals: int = 6,
        parent=None
    ):
        """
        Initialize unit-aware spinbox.

        Args:
            unit_type: Type of unit (frequency, time, etc.)
            initial_value: Initial value in base units
            initial_unit: Initial display unit (auto-selected if None)
            min_value: Minimum value in base units
            max_value: Maximum value in base units
            decimals: Number of decimal places
            parent: Parent widget
        """
        super().__init__(parent)

        self.unit_type = unit_type
        self.converter = UnitConverter()
        self.min_value = min_value
        self.max_value = max_value
        self.decimals = decimals

        # Current value in base units
        self._base_value = initial_value

        # Get available units
        self.available_units = self.converter.get_available_units(unit_type)

        # Auto-select initial unit if not provided
        if initial_unit is None:
            _, initial_unit = self.converter.get_preferred_unit(initial_value, unit_type)

        self.current_unit = initial_unit

        # Create UI
        self._create_ui()
        self._update_spinbox_range()
        self.set_value(initial_value)

    def _create_ui(self):
        """Create spinbox and unit selector UI."""
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(5)

        # Spinbox for value
        self.spinbox = QDoubleSpinBox()
        self.spinbox.setDecimals(self.decimals)
        self.spinbox.setMinimumWidth(120)
        self.spinbox.valueChanged.connect(self._on_spinbox_changed)
        layout.addWidget(self.spinbox)

        # Unit selector
        self.unit_combo = QComboBox()
        self.unit_combo.addItems(self.available_units)
        self.unit_combo.setCurrentText(self.current_unit)
        self.unit_combo.currentTextChanged.connect(self._on_unit_changed)
        self.unit_combo.setMinimumWidth(70)
        layout.addWidget(self.unit_combo)

        self.setLayout(layout)

    def _update_spinbox_range(self):
        """Update spinbox min/max based on current unit."""
        min_display = self.converter.from_base(self.min_value, self.current_unit, self.unit_type)
        max_display = self.converter.from_base(self.max_value, self.current_unit, self.unit_type)

        self.spinbox.blockSignals(True)
        self.spinbox.setRange(min_display, max_display)
        self.spinbox.blockSignals(False)

    def _on_spinbox_changed(self, display_value: float):
        """Handle spinbox value change."""
        # Convert to base units
        self._base_value = self.converter.to_base(display_value, self.current_unit, self.unit_type)

        # Emit signal
        self.value_changed.emit(self._base_value)

    def _on_unit_changed(self, new_unit: str):
        """Handle unit change."""
        # Convert current base value to new unit
        old_unit = self.current_unit
        self.current_unit = new_unit

        # Update range
        self._update_spinbox_range()

        # Update displayed value (maintain base value)
        display_value = self.converter.from_base(self._base_value, new_unit, self.unit_type)

        self.spinbox.blockSignals(True)
        self.spinbox.setValue(display_value)
        self.spinbox.blockSignals(False)

    def set_value(self, value: float):
        """
        Set value in base units.

        Args:
            value: Value in base units
        """
        self._base_value = value

        # Convert to current display unit
        display_value = self.converter.from_base(value, self.current_unit, self.unit_type)

        self.spinbox.blockSignals(True)
        self.spinbox.setValue(display_value)
        self.spinbox.blockSignals(False)

    def get_value(self) -> float:
        """
        Get value in base units.

        Returns:
            Value in base units
        """
        return self._base_value

    def get_display_value(self) -> Tuple[float, str]:
        """
        Get current display value and unit.

        Returns:
            (value, unit)
        """
        return self.spinbox.value(), self.current_unit


class ParameterEditorWidget(QWidget):
    """
    Complete parameter editor for NV experiments.

    Features:
    - Unit-aware parameter inputs
    - Experiment type selection (ESR, Rabi, T1, T2)
    - Parameter validation
    - Template loading
    - Save/load configuration
    """

    parameters_changed = Signal(dict)  # Emitted when parameters change
    open_window_requested = Signal(dict)  # Emitted when Open Window clicked

    def __init__(self, parent=None):
        """Initialize parameter editor."""
        super().__init__(parent)

        # Parameter storage (base units)
        self.parameters = {}

        # Parameter widgets
        self.param_widgets = {}

        # Experiment type
        self.experiment_type = "ESR"

        # Create UI
        self._create_ui()
        self._load_default_parameters()

    def _create_ui(self):
        """Create parameter editor UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Header
        header_layout = QHBoxLayout()

        header_label = QLabel("Experiment Configuration")
        header_label.setProperty("class", "header")
        header_layout.addWidget(header_label)

        header_layout.addStretch()

        # Load/Save buttons
        self.load_btn = QPushButton("Load Config")
        self.load_btn.clicked.connect(self.load_configuration)
        header_layout.addWidget(self.load_btn)

        self.save_btn = QPushButton("Save Config")
        self.save_btn.setProperty("class", "success")
        self.save_btn.clicked.connect(self.save_configuration)
        header_layout.addWidget(self.save_btn)

        layout.addLayout(header_layout)

        # Experiment type selector
        type_group = QGroupBox("Experiment Type")
        type_layout = QHBoxLayout()

        type_label = QLabel("Type:")
        type_layout.addWidget(type_label)

        self.type_combo = QComboBox()
        self.type_combo.addItems(["ESR", "Rabi", "T1", "T2", "Ramsey"])
        self.type_combo.currentTextChanged.connect(self._on_experiment_type_changed)
        type_layout.addWidget(self.type_combo)

        self.load_template_btn = QPushButton("Load Template")
        self.load_template_btn.clicked.connect(self._load_template)
        type_layout.addWidget(self.load_template_btn)

        type_layout.addStretch()

        type_group.setLayout(type_layout)
        layout.addWidget(type_group)

        # Scrollable parameter area
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)

        # Parameter container
        param_widget = QWidget()
        self.param_layout = QVBoxLayout()
        self.param_layout.setSpacing(15)

        # Create parameter groups
        self._create_microwave_group()
        self._create_timing_group()
        self._create_acquisition_group()
        self._create_sweep_group()

        self.param_layout.addStretch()
        param_widget.setLayout(self.param_layout)
        scroll_area.setWidget(param_widget)

        layout.addWidget(scroll_area)

        # Apply and Open Window buttons
        apply_layout = QHBoxLayout()
        apply_layout.addStretch()

        self.apply_btn = QPushButton("Apply Parameters")
        self.apply_btn.setProperty("class", "success")
        self.apply_btn.clicked.connect(self._apply_parameters)
        apply_layout.addWidget(self.apply_btn)

        self.open_window_btn = QPushButton("Open Window")
        self.open_window_btn.setProperty("class", "primary")
        self.open_window_btn.setMinimumWidth(150)
        self.open_window_btn.clicked.connect(self._open_window)
        apply_layout.addWidget(self.open_window_btn)

        layout.addLayout(apply_layout)

        self.setLayout(layout)

    def _create_microwave_group(self):
        """Create microwave parameters group."""
        group = QGroupBox("Microwave Parameters")
        form_layout = QFormLayout()

        # Frequency
        self.freq_spinbox = UnitAwareSpinBox(
            UnitType.FREQUENCY,
            initial_value=2.87e9,
            initial_unit='GHz',
            min_value=0.0,
            max_value=20e9,
            decimals=6
        )
        self.param_widgets['frequency'] = self.freq_spinbox
        form_layout.addRow("Frequency:", self.freq_spinbox)

        # Power
        self.power_spinbox = QDoubleSpinBox()
        self.power_spinbox.setRange(-20.0, 20.0)
        self.power_spinbox.setValue(8.0)
        self.power_spinbox.setSuffix(" dBm")
        self.power_spinbox.setDecimals(2)
        self.param_widgets['power'] = self.power_spinbox
        form_layout.addRow("Power:", self.power_spinbox)

        # Modulation
        self.modulation_combo = QComboBox()
        self.modulation_combo.addItems(["Pulse", "CW", "FM", "AM"])
        self.param_widgets['modulation'] = self.modulation_combo
        form_layout.addRow("Modulation:", self.modulation_combo)

        group.setLayout(form_layout)
        self.param_layout.addWidget(group)

    def _create_timing_group(self):
        """Create timing parameters group."""
        group = QGroupBox("Timing Parameters")
        form_layout = QFormLayout()

        # t_AOM (AOM delay)
        self.t_aom_spinbox = UnitAwareSpinBox(
            UnitType.TIME,
            initial_value=8000e-9,  # 8000 ns = 8 us
            initial_unit='ns',
            min_value=0.0,
            max_value=1.0,  # 1 second
            decimals=3
        )
        self.param_widgets['t_aom'] = self.t_aom_spinbox
        form_layout.addRow("t_AOM:", self.t_aom_spinbox)

        # Read delay
        self.ro_delay_spinbox = UnitAwareSpinBox(
            UnitType.TIME,
            initial_value=100e-9,  # 100 ns
            initial_unit='ns',
            min_value=0.0,
            max_value=1e-3,  # 1 ms
            decimals=3
        )
        self.param_widgets['ro_delay'] = self.ro_delay_spinbox
        form_layout.addRow("Readout Delay:", self.ro_delay_spinbox)

        # Tau (for T1, T2, Rabi)
        self.tau_spinbox = UnitAwareSpinBox(
            UnitType.TIME,
            initial_value=1000e-9,  # 1000 ns = 1 us
            initial_unit='ns',
            min_value=0.0,
            max_value=1e-3,  # 1 ms
            decimals=3
        )
        self.param_widgets['tau'] = self.tau_spinbox
        form_layout.addRow("Tau (π-pulse):", self.tau_spinbox)

        group.setLayout(form_layout)
        self.param_layout.addWidget(group)

    def _create_acquisition_group(self):
        """Create acquisition parameters group."""
        group = QGroupBox("Acquisition Settings")
        form_layout = QFormLayout()

        # Hardware or simulated
        self.hw_mode_combo = QComboBox()
        self.hw_mode_combo.addItems(["Simulated", "Real Hardware"])
        self.param_widgets['hw_mode'] = self.hw_mode_combo
        form_layout.addRow("Hardware Mode:", self.hw_mode_combo)

        # Acquisition mode (Diode/Camera)
        self.acq_mode_combo = QComboBox()
        self.acq_mode_combo.addItems(["Diode", "Camera"])
        self.param_widgets['acquisition_mode'] = self.acq_mode_combo
        form_layout.addRow("Detector Mode:", self.acq_mode_combo)

        # Sampling rate
        self.sampling_rate_spinbox = UnitAwareSpinBox(
            UnitType.FREQUENCY,
            initial_value=100e3,  # 100 kS/s
            initial_unit='kHz',
            min_value=1e3,
            max_value=10e6,
            decimals=2
        )
        self.param_widgets['sampling_rate'] = self.sampling_rate_spinbox
        form_layout.addRow("Sampling Rate:", self.sampling_rate_spinbox)

        # Number of averages
        self.averages_spinbox = QSpinBox()
        self.averages_spinbox.setRange(1, 100000)
        self.averages_spinbox.setValue(1000)
        self.param_widgets['num_averages'] = self.averages_spinbox
        form_layout.addRow("Averages:", self.averages_spinbox)

        group.setLayout(form_layout)
        self.param_layout.addWidget(group)

    def _create_sweep_group(self):
        """Create sweep parameters group."""
        group = QGroupBox("Sweep Parameters")
        form_layout = QFormLayout()

        # Sweep parameter
        self.sweep_param_combo = QComboBox()
        self.sweep_param_combo.addItems(["Frequency", "Tau", "Power", "Magnetic Field"])
        self.param_widgets['sweep_parameter'] = self.sweep_param_combo
        form_layout.addRow("Sweep:", self.sweep_param_combo)

        # Start value
        self.sweep_start_spinbox = QDoubleSpinBox()
        self.sweep_start_spinbox.setRange(-1e12, 1e12)
        self.sweep_start_spinbox.setValue(2.85e9)
        self.sweep_start_spinbox.setDecimals(6)
        self.param_widgets['sweep_start'] = self.sweep_start_spinbox
        form_layout.addRow("Start:", self.sweep_start_spinbox)

        # Stop value
        self.sweep_stop_spinbox = QDoubleSpinBox()
        self.sweep_stop_spinbox.setRange(-1e12, 1e12)
        self.sweep_stop_spinbox.setValue(2.89e9)
        self.sweep_stop_spinbox.setDecimals(6)
        self.param_widgets['sweep_stop'] = self.sweep_stop_spinbox
        form_layout.addRow("Stop:", self.sweep_stop_spinbox)

        # Number of points
        self.sweep_points_spinbox = QSpinBox()
        self.sweep_points_spinbox.setRange(2, 10000)
        self.sweep_points_spinbox.setValue(51)
        self.param_widgets['sweep_points'] = self.sweep_points_spinbox
        form_layout.addRow("Points:", self.sweep_points_spinbox)

        # Sweep mode
        self.sweep_mode_combo = QComboBox()
        self.sweep_mode_combo.addItems(["Linear", "Log", "Manual"])
        self.param_widgets['sweep_mode'] = self.sweep_mode_combo
        form_layout.addRow("Mode:", self.sweep_mode_combo)

        group.setLayout(form_layout)
        self.param_layout.addWidget(group)

    def _load_default_parameters(self):
        """Load default parameter values."""
        self.parameters = {
            'experiment_type': 'ESR',
            'frequency': 2.87e9,
            'power': 8.0,
            'modulation': 'Pulse',
            't_aom': 8000e-9,
            'ro_delay': 100e-9,
            'tau': 1000e-9,
            'acquisition_mode': 'Diode',
            'sampling_rate': 100e3,
            'num_averages': 1000,
            'sweep_parameter': 'Frequency',
            'sweep_start': 2.85e9,
            'sweep_stop': 2.89e9,
            'sweep_points': 51,
            'sweep_mode': 'Linear'
        }

    def _on_experiment_type_changed(self, exp_type: str):
        """Handle experiment type change."""
        self.experiment_type = exp_type
        # Could adjust visible parameters based on type

    def _load_template(self):
        """Load parameter template for current experiment type."""
        templates = {
            'ESR': {
                'frequency': 2.87e9,
                'power': 8.0,
                't_aom': 8000e-9,
                'sweep_start': 2.85e9,
                'sweep_stop': 2.89e9,
                'sweep_points': 51,
                'num_averages': 1000
            },
            'Rabi': {
                'frequency': 2.87e9,
                'power': 8.0,
                't_aom': 8000e-9,
                'sweep_start': 0.0,
                'sweep_stop': 5000e-9,
                'sweep_points': 51,
                'num_averages': 2000
            },
            'T1': {
                'frequency': 2.87e9,
                'power': 8.0,
                't_aom': 8000e-9,
                'tau': 1000e-9,
                'sweep_start': 0.0,
                'sweep_stop': 100e-6,
                'sweep_points': 51,
                'num_averages': 2000
            },
            'T2': {
                'frequency': 2.87e9,
                'power': 8.0,
                't_aom': 8000e-9,
                'tau': 1000e-9,
                'sweep_start': 0.0,
                'sweep_stop': 10e-6,
                'sweep_points': 51,
                'num_averages': 3000
            }
        }

        if self.experiment_type in templates:
            template = templates[self.experiment_type]
            self.set_parameters(template)
            QMessageBox.information(self, "Template Loaded",
                                   f"Loaded default parameters for {self.experiment_type}")

    def _apply_parameters(self):
        """Apply current parameters."""
        params = self.get_parameters()
        self.parameters_changed.emit(params)
        QMessageBox.information(self, "Parameters Applied",
                               "Parameters have been applied successfully.")

    def _open_window(self):
        """Open experiment window with current parameters."""
        params = self.get_parameters()
        self.open_window_requested.emit(params)

    def get_parameters(self) -> Dict[str, Any]:
        """
        Get current parameter values (in base units).

        Returns:
            Dictionary of parameters
        """
        # Determine if using real hardware
        hw_mode = self.hw_mode_combo.currentText()
        acq_mode = 'real' if hw_mode == "Real Hardware" else 'simulated'

        params = {
            'experiment_type': self.type_combo.currentText(),
            'frequency': self.freq_spinbox.get_value(),
            'power': self.power_spinbox.value(),
            'modulation': self.modulation_combo.currentText(),
            't_aom': self.t_aom_spinbox.get_value(),
            'ro_delay': self.ro_delay_spinbox.get_value(),
            'tau': self.tau_spinbox.get_value(),
            'acquisition_mode': acq_mode,  # 'real' or 'simulated'
            'detector_mode': self.acq_mode_combo.currentText().lower(),  # 'diode' or 'camera'
            'sampling_rate': self.sampling_rate_spinbox.get_value(),
            'num_averages': self.averages_spinbox.value(),
            'sweep_parameter': self.sweep_param_combo.currentText(),
            'sweep_start': self.sweep_start_spinbox.value(),
            'sweep_stop': self.sweep_stop_spinbox.value(),
            'sweep_points': self.sweep_points_spinbox.value(),
            'sweep_mode': self.sweep_mode_combo.currentText(),
            'trial_run': ['y', 'y'] if acq_mode == 'simulated' else ['n', 'n']  # For hardware controller
        }
        return params

    def set_parameters(self, params: Dict[str, Any]):
        """
        Set parameter values from dictionary.

        Args:
            params: Dictionary of parameters (base units)
        """
        if 'experiment_type' in params:
            self.type_combo.setCurrentText(params['experiment_type'])

        if 'frequency' in params:
            self.freq_spinbox.set_value(params['frequency'])

        if 'power' in params:
            self.power_spinbox.setValue(params['power'])

        if 'modulation' in params:
            self.modulation_combo.setCurrentText(params['modulation'])

        if 't_aom' in params:
            self.t_aom_spinbox.set_value(params['t_aom'])

        if 'ro_delay' in params:
            self.ro_delay_spinbox.set_value(params['ro_delay'])

        if 'tau' in params:
            self.tau_spinbox.set_value(params['tau'])

        if 'acquisition_mode' in params:
            self.acq_mode_combo.setCurrentText(params['acquisition_mode'])

        if 'sampling_rate' in params:
            self.sampling_rate_spinbox.set_value(params['sampling_rate'])

        if 'num_averages' in params:
            self.averages_spinbox.setValue(params['num_averages'])

        if 'sweep_parameter' in params:
            self.sweep_param_combo.setCurrentText(params['sweep_parameter'])

        if 'sweep_start' in params:
            self.sweep_start_spinbox.setValue(params['sweep_start'])

        if 'sweep_stop' in params:
            self.sweep_stop_spinbox.setValue(params['sweep_stop'])

        if 'sweep_points' in params:
            self.sweep_points_spinbox.setValue(params['sweep_points'])

        if 'sweep_mode' in params:
            self.sweep_mode_combo.setCurrentText(params['sweep_mode'])

    def save_configuration(self):
        """Save current configuration to YAML file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Configuration",
            "",
            "YAML Files (*.yaml *.yml);;JSON Files (*.json)"
        )

        if file_path:
            params = self.get_parameters()

            try:
                if file_path.endswith('.json'):
                    with open(file_path, 'w') as f:
                        json.dump(params, f, indent=4)
                else:
                    with open(file_path, 'w') as f:
                        yaml.dump(params, f, default_flow_style=False)

                QMessageBox.information(self, "Success",
                                       f"Configuration saved to {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error",
                                    f"Failed to save configuration: {e}")

    def load_configuration(self):
        """Load configuration from YAML or JSON file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Configuration",
            "",
            "Config Files (*.yaml *.yml *.json)"
        )

        if file_path:
            try:
                if file_path.endswith('.json'):
                    with open(file_path, 'r') as f:
                        params = json.load(f)
                else:
                    with open(file_path, 'r') as f:
                        params = yaml.safe_load(f)

                self.set_parameters(params)
                QMessageBox.information(self, "Success",
                                       f"Configuration loaded from {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error",
                                    f"Failed to load configuration: {e}")


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    widget = ParameterEditorWidget()
    widget.show()

    sys.exit(app.exec())
