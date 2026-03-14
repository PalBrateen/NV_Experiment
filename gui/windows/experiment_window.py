"""
Generic Experiment Window for NV Experiment GUI

Single window class that adapts based on experiment type from parameters.
Supports ESR, Rabi, T1, T2, and other experiments.

Features:
- Adapts behavior based on experiment_type parameter
- Automatically determines required instruments
- Configures sweep parameter and axis labels
- Selects appropriate fit function
- All configuration from main window Config tab
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Dict, Any, Optional, List
import numpy as np
import time

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QCheckBox, QSplitter
)
from PySide6.QtCore import QTimer, Qt

from gui.windows.base_experiment_window import BaseExperimentWindow
from gui.widgets.sweep_builder import SweepBuilderWidget
from gui.widgets.quick_fit_panel import QuickFitPanel
from gui.widgets.scan_selector import ScanSelectorWidget
from gui.widgets.live_plot_widget import LivePlotWidget
from gui.gui_experiment_controller import GuiExperimentController


class ExperimentWindow(BaseExperimentWindow):
    """
    Generic experiment window that adapts based on parameters.

    Supported experiment types:
    - ESR: Frequency sweep, requires SG + PB + DAQ
    - Rabi: Tau sweep, requires SG + PB + DAQ
    - T1: Delay sweep, requires SG + PB + DAQ
    - T2: Tau sweep (Hahn echo), requires SG + PB + DAQ
    """

    # Experiment configurations
    EXPERIMENT_CONFIGS = {
        'ESR': {
            'required_instruments': ['signal_generator', 'pulseblaster', 'analog_input'],
            'sweep_parameter': 'frequency',
            'x_label': 'Frequency',
            'x_unit': 'GHz',
            'x_scale': 1e-9,  # Convert Hz to GHz for display
            'fit_type': 'lorentzian',
            'default_sweep': (2.85e9, 2.89e9, 51)  # start, stop, points
        },
        'Rabi': {
            'required_instruments': ['signal_generator', 'pulseblaster', 'analog_input'],
            'sweep_parameter': 'tau',
            'x_label': 'Pulse Width (τ)',
            'x_unit': 'ns',
            'x_scale': 1e9,  # Convert s to ns for display
            'fit_type': 'rabi',
            'default_sweep': (0.0, 5000e-9, 51)  # 0 to 5000 ns
        },
        'T1': {
            'required_instruments': ['signal_generator', 'pulseblaster', 'analog_input'],
            'sweep_parameter': 'delay',
            'x_label': 'Delay',
            'x_unit': 'μs',
            'x_scale': 1e6,  # Convert s to μs for display
            'fit_type': 'exponential',
            'default_sweep': (0.0, 100e-6, 51)  # 0 to 100 μs
        },
        'T2': {
            'required_instruments': ['signal_generator', 'pulseblaster', 'analog_input'],
            'sweep_parameter': 'tau',
            'x_label': 'Tau (2τ)',
            'x_unit': 'μs',
            'x_scale': 1e6,  # Convert s to μs for display
            'fit_type': 'exponential',
            'default_sweep': (0.0, 10e-6, 51)  # 0 to 10 μs
        },
        'Ramsey': {
            'required_instruments': ['signal_generator', 'pulseblaster', 'analog_input'],
            'sweep_parameter': 'tau',
            'x_label': 'Free Evolution Time',
            'x_unit': 'μs',
            'x_scale': 1e6,
            'fit_type': 'ramsey',
            'default_sweep': (0.0, 5e-6, 51)
        }
    }

    def __init__(self, window_id: str, parameters: Dict[str, Any], parent=None):
        """
        Initialize generic experiment window.

        Args:
            window_id: Unique identifier
            parameters: Parameter dictionary from main window Config tab
            parent: Parent widget
        """
        # Extract experiment type from parameters
        self.exp_type = parameters.get('experiment_type', 'ESR')

        # Get configuration for this experiment type
        if self.exp_type not in self.EXPERIMENT_CONFIGS:
            self.exp_type = 'ESR'  # Fallback

        self.config = self.EXPERIMENT_CONFIGS[self.exp_type]

        # Set required instruments before calling super().__init__
        self.required_instruments = self.config['required_instruments']

        super().__init__(
            window_id=window_id,
            experiment_type=self.exp_type,
            parent=parent
        )

        # Store all parameters
        self.set_parameters(parameters)

        # Acquisition mode: 'simulated' or 'real'
        self.acquisition_mode = parameters.get('acquisition_mode', 'simulated')

        # Simulated acquisition
        self.acquisition_timer = QTimer()
        self.acquisition_timer.timeout.connect(self._simulate_acquisition_step)
        self.current_sweep_index = 0

        # Real acquisition
        self.gui_controller: Optional[GuiExperimentController] = None
        self.use_real_hardware = (self.acquisition_mode == 'real')

        if self.use_real_hardware:
            self._initialize_hardware_controller()

        # Fit curve storage
        self.fit_curve_item = None  # pyqtgraph plot item for fit overlay

        # Configure initial sweep from parameters
        self._configure_sweep_from_parameters()

    def _create_controls_panel(self) -> QWidget:
        """Override to add experiment-specific controls."""
        # Get base controls
        widget = super()._create_controls_panel()

        # Display current experiment parameters
        self._add_parameter_display()

        # Sweep builder
        sweep_group = QGroupBox("Sweep Configuration")
        sweep_layout = QVBoxLayout()

        self.sweep_builder = SweepBuilderWidget()
        self.sweep_builder.sweep_changed.connect(self._on_sweep_changed)
        sweep_layout.addWidget(self.sweep_builder)

        sweep_group.setLayout(sweep_layout)
        self.parameters_layout.addWidget(sweep_group)

        # Auto-fit checkbox
        autofit_layout = QHBoxLayout()
        self.autofit_checkbox = QCheckBox(f"Auto-fit ({self.config['fit_type'].capitalize()})")
        self.autofit_checkbox.setChecked(True)
        autofit_layout.addWidget(self.autofit_checkbox)
        autofit_layout.addStretch()
        self.parameters_layout.addLayout(autofit_layout)

        # Scan selector
        scan_group = QGroupBox("Scan Management")
        scan_layout = QVBoxLayout()

        self.scan_selector = ScanSelectorWidget()
        self.scan_selector.selection_changed.connect(self._on_scan_selection_changed)
        self.scan_selector.color_changed.connect(self._on_scan_color_changed)
        scan_layout.addWidget(self.scan_selector)

        scan_group.setLayout(scan_layout)
        self.parameters_layout.addWidget(scan_group)

        return widget

    def _add_parameter_display(self):
        """Add parameter display labels to controls."""
        # Experiment type
        exp_layout = QHBoxLayout()
        exp_layout.addWidget(QLabel("Experiment:"))
        self.exp_type_label = QLabel(self.exp_type)
        self.exp_type_label.setStyleSheet("font-weight: bold; color: #4ec9b0;")
        exp_layout.addWidget(self.exp_type_label)
        exp_layout.addStretch()
        self.parameters_layout.addLayout(exp_layout)

        # Sweep parameter
        param_layout = QHBoxLayout()
        param_layout.addWidget(QLabel("Sweep Parameter:"))
        self.sweep_param_label = QLabel(self.config['sweep_parameter'].capitalize())
        self.sweep_param_label.setStyleSheet("font-weight: bold;")
        param_layout.addWidget(self.sweep_param_label)
        param_layout.addStretch()
        self.parameters_layout.addLayout(param_layout)

        # Frequency (for ESR)
        if self.exp_type == 'ESR':
            freq_layout = QHBoxLayout()
            freq_layout.addWidget(QLabel("Center Freq:"))
            self.center_freq_label = QLabel("N/A")
            self.center_freq_label.setStyleSheet("font-weight: bold;")
            freq_layout.addWidget(self.center_freq_label)
            freq_layout.addStretch()
            self.parameters_layout.addLayout(freq_layout)

        # Number of averages
        avg_layout = QHBoxLayout()
        avg_layout.addWidget(QLabel("Averages:"))
        self.averages_label = QLabel(str(self.parameters.get('num_averages', 1000)))
        self.averages_label.setStyleSheet("font-weight: bold;")
        avg_layout.addWidget(self.averages_label)
        avg_layout.addStretch()
        self.parameters_layout.addLayout(avg_layout)

    def _create_plot_panel(self) -> QWidget:
        """Override to add fit panel alongside plot."""
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        # Splitter for plot and fit panel
        splitter = QSplitter(Qt.Vertical)

        # Live plot widget
        self.live_plot = LivePlotWidget(
            title=f"{self.exp_type} Spectrum",
            y_label="Signal",
            y_unit="V",
            buffer_size=1000
        )

        # Override x-axis label based on experiment type
        plot_item = self.live_plot.plot_widget.getPlotItem()
        plot_item.setLabel('bottom', self.config['x_label'], units=self.config['x_unit'])

        splitter.addWidget(self.live_plot)

        # Quick fit panel
        self.fit_panel = QuickFitPanel()
        self.fit_panel.set_fit_type(self.config['fit_type'])
        self.fit_panel.fit_completed.connect(self._on_fit_completed)
        self.fit_panel.fit_failed.connect(self._on_fit_failed)
        splitter.addWidget(self.fit_panel)

        # Set splitter sizes (70% plot, 30% fit panel)
        splitter.setSizes([700, 300])

        layout.addWidget(splitter)

        widget.setLayout(layout)
        return widget

    def _configure_sweep_from_parameters(self):
        """Configure sweep builder from parameters."""
        # Try to get sweep configuration from parameters
        start = self.parameters.get('sweep_start')
        stop = self.parameters.get('sweep_stop')
        points = self.parameters.get('sweep_points')

        # If not in parameters, use defaults from config
        if start is None or stop is None or points is None:
            start, stop, points = self.config['default_sweep']

        # Set sweep
        self.sweep_builder.set_linear_sweep(start, stop, points)

    def _initialize_hardware_controller(self):
        """Initialize hardware controller for real acquisition."""
        try:
            # Create controller
            self.gui_controller = GuiExperimentController(
                owner_id=self.window_id,
                acquisition_mode='diode'  # Default to diode mode
            )

            # Connect signals
            self.gui_controller.progress_updated.connect(self._on_hardware_progress)
            self.gui_controller.data_point_acquired.connect(self._on_hardware_data_point)
            self.gui_controller.experiment_finished.connect(self._on_hardware_finished)
            self.gui_controller.error_occurred.connect(self._on_hardware_error)

            print(f"Hardware controller initialized for {self.window_id}")

        except Exception as e:
            print(f"Error initializing hardware controller: {e}")
            self.use_real_hardware = False
            self.gui_controller = None

    def _on_sweep_changed(self, sweep_array: np.ndarray):
        """Handle sweep configuration change."""
        if len(sweep_array) > 0 and self.exp_type == 'ESR':
            # Update center frequency display for ESR
            center_freq = (sweep_array[0] + sweep_array[-1]) / 2
            self.center_freq_label.setText(f"{center_freq / 1e9:.6f} GHz")

    def _start_acquisition(self):
        """Override to implement generic acquisition."""
        super()._start_acquisition()

        # Get sweep array
        sweep_array = self.sweep_builder.get_sweep_array()

        if len(sweep_array) == 0:
            self.status_label.setText("Status: Error - No sweep defined")
            self._stop_acquisition()
            return

        # Initialize acquisition
        self.current_sweep_index = 0
        self.sweep_params = []
        self.sweep_signals = []

        # Clear fit curve
        if self.fit_curve_item is not None:
            self.live_plot.plot_widget.removeItem(self.fit_curve_item)
            self.fit_curve_item = None

        # Update plot axes based on experiment type
        plot_item = self.live_plot.plot_widget.getPlotItem()
        plot_item.setLabel('bottom', self.config['x_label'], units=self.config['x_unit'])
        plot_item.setLabel('left', 'Signal', units='V')

        # Start acquisition (real or simulated)
        if self.use_real_hardware and self.gui_controller is not None:
            self._start_hardware_acquisition()
        else:
            self._start_simulated_acquisition()

    def _start_simulated_acquisition(self):
        """Start simulated acquisition."""
        self.acquisition_timer.start(50)  # 50ms per point (simulated)

    def _start_hardware_acquisition(self):
        """Start real hardware acquisition."""
        try:
            # Build params_dict from parameters and sweep
            sweep_array = self.sweep_builder.get_sweep_array()

            # Create params_dict in expected format
            params_dict = self._build_params_dict(sweep_array)

            # Configure controller
            trial_run = self.parameters.get('trial_run', ['n', 'n'])
            success, error = self.gui_controller.configure(params_dict, trial_run)

            if not success:
                self.status_label.setText(f"Status: Configuration failed - {error}")
                self._stop_acquisition()
                return

            # Initialize instruments
            success, error = self.gui_controller.initialize_instruments()
            if not success:
                self.status_label.setText(f"Status: Initialization failed - {error}")
                self._stop_acquisition()
                return

            # Start acquisition
            success, error = self.gui_controller.start_acquisition()
            if not success:
                self.status_label.setText(f"Status: Start failed - {error}")
                self._stop_acquisition()
                return

            self.status_label.setText("Status: Running (Hardware)")

        except Exception as e:
            self.status_label.setText(f"Status: Error - {str(e)}")
            self._stop_acquisition()

    def _build_params_dict(self, sweep_array: np.ndarray) -> Dict[str, Any]:
        """Build params_dict for hardware controller."""
        # Map experiment type to sequence name
        sequence_map = {
            'ESR': 'esr',
            'Rabi': 'rabi',
            'T1': 't1',
            'T2': 't2_hahn',
            'Ramsey': 'ramsey'
        }

        params_dict = {
            'seq': {
                'sequence': sequence_map.get(self.exp_type, 'esr')
            },
            'scan': {
                'Nruns': self.parameters.get('num_averages', 1000)
            },
            'mw': {
                'freq': sweep_array if self.exp_type == 'ESR' else self.parameters.get('frequency', 2.87e9),
                'power': self.parameters.get('power', 8.0)
            },
            'daq': {
                'sampling_rate': self.parameters.get('sampling_rate', 100000),
                'Nsamples': 1000,
                'use_ao': False
            }
        }

        # Add sweep configuration
        if self.exp_type == 'ESR':
            params_dict['sweep'] = {'frequency': sweep_array}
        elif self.exp_type in ['Rabi', 'T2', 'Ramsey']:
            params_dict['sweep'] = {'tau': sweep_array}
        elif self.exp_type == 'T1':
            params_dict['sweep'] = {'delay': sweep_array}

        return params_dict

    def _simulate_acquisition_step(self):
        """Simulate one acquisition step (for testing)."""
        if not self.acquisition_running or self.acquisition_paused_state:
            return

        sweep_array = self.sweep_builder.get_sweep_array()

        if self.current_sweep_index >= len(sweep_array):
            # Acquisition complete
            self.acquisition_timer.stop()
            self._finish_acquisition()
            return

        # Get current parameter value
        param_value = sweep_array[self.current_sweep_index]

        # Simulate signal based on experiment type
        signal = self._simulate_signal(param_value)

        # Store data
        self.sweep_params.append(param_value)
        self.sweep_signals.append(signal)

        # Update plot (apply scale for display)
        display_x = param_value * self.config['x_scale']
        self.live_plot.add_data_point(display_x, signal)

        # Update progress
        progress = int(100 * (self.current_sweep_index + 1) / len(sweep_array))
        self.progress_bar.setValue(progress)
        self.progress_label.setText(f"Point {self.current_sweep_index + 1} / {len(sweep_array)}")

        # Emit signals
        self.data_point_acquired.emit(self.current_sweep_index, param_value, signal)
        self.progress_updated.emit(self.current_sweep_index + 1, len(sweep_array))

        # Move to next point
        self.current_sweep_index += 1

    def _simulate_signal(self, param_value: float) -> float:
        """
        Simulate signal based on experiment type.

        Args:
            param_value: Current sweep parameter value

        Returns:
            Simulated signal value
        """
        baseline = 2.5
        noise_level = 0.05

        if self.exp_type == 'ESR':
            # Lorentzian dip centered at 2.87 GHz
            center = 2.87e9
            width = 5e6  # 5 MHz
            amplitude = 0.5
            signal = baseline - amplitude / (1 + ((param_value - center) / width)**2)

        elif self.exp_type == 'Rabi':
            # Rabi oscillation
            period = 2000e-9  # 2000 ns period
            amplitude = 0.5
            signal = baseline - amplitude * (1 - np.cos(2 * np.pi * param_value / period)) / 2

        elif self.exp_type == 'T1':
            # Exponential decay
            T1 = 20e-6  # 20 μs
            amplitude = 0.5
            signal = baseline - amplitude * np.exp(-param_value / T1)

        elif self.exp_type == 'T2':
            # Exponential decay (Hahn echo)
            T2 = 5e-6  # 5 μs
            amplitude = 0.5
            signal = baseline - amplitude * np.exp(-param_value / T2)

        elif self.exp_type == 'Ramsey':
            # Ramsey fringes (oscillation + decay)
            T2_star = 2e-6  # 2 μs
            detuning = 1e6  # 1 MHz detuning
            amplitude = 0.5
            signal = baseline - amplitude * np.exp(-param_value / T2_star) * np.cos(2 * np.pi * detuning * param_value)

        else:
            # Default
            signal = baseline

        # Add noise
        signal += noise_level * np.random.randn()

        return signal

    def _finish_acquisition(self):
        """Override to add scan to selector and trigger fitting."""
        # Call base implementation first (stores scan, updates UI)
        super()._finish_acquisition()

        # Add scan to selector
        if len(self.sweep_params) > 0:
            metadata = {
                'experiment_type': self.exp_type,
                'num_points': len(self.sweep_params),
                'num_averages': self.parameters.get('num_averages', 1000)
            }
            self.scan_selector.add_scan(
                (np.array(self.sweep_params), np.array(self.sweep_signals)),
                metadata
            )

            # Auto-fit if enabled
            if self.autofit_checkbox.isChecked():
                self._perform_fit()

    def _perform_fit(self):
        """Perform fitting on current data."""
        if len(self.sweep_params) == 0:
            return

        # Set data in fit panel
        self.fit_panel.set_data(
            np.array(self.sweep_params),
            np.array(self.sweep_signals)
        )

        # Trigger fit
        self.fit_panel._perform_fit()

    def _on_fit_completed(self, params: dict, uncertainties: dict, curve: tuple):
        """Handle fit completion."""
        # Remove previous fit curve if exists
        if self.fit_curve_item is not None:
            self.live_plot.plot_widget.removeItem(self.fit_curve_item)

        # Add new fit curve
        x_fit, y_fit = curve

        # Apply display scaling
        x_fit_display = x_fit * self.config['x_scale']

        # Plot fit curve (red dashed line)
        self.fit_curve_item = self.live_plot.plot_widget.plot(
            x_fit_display,
            y_fit,
            pen={'color': '#ff0000', 'width': 2, 'style': Qt.DashLine},
            name='Fit'
        )

    def _on_fit_failed(self, error: str):
        """Handle fit failure."""
        self.status_label.setText(f"Status: Fit failed - {error}")

    def _on_scan_selection_changed(self, selected_indices: list):
        """Handle scan selection change from selector."""
        # Clear plot
        self.live_plot.clear_data()

        # Clear fit curve
        if self.fit_curve_item is not None:
            self.live_plot.plot_widget.removeItem(self.fit_curve_item)
            self.fit_curve_item = None

        # Plot selected scans
        for scan_index in selected_indices:
            scan_info = self.scan_selector.scans[scan_index]
            x_data, y_data = scan_info['data']
            color = self.scan_selector.get_scan_color(scan_index)

            # Apply display scaling
            x_display = x_data * self.config['x_scale']

            # Plot scan with its color
            self.live_plot.plot_widget.plot(
                x_display,
                y_data,
                pen={'color': color.name(), 'width': 2},
                symbol='o',
                symbolSize=5,
                symbolBrush=color,
                name=f"Scan #{scan_index + 1}"
            )

    def _on_scan_color_changed(self, scan_index: int, color):
        """Handle scan color change."""
        # Refresh plot to show new color
        selected_indices = self.scan_selector.get_selected_scan_indices()
        self._on_scan_selection_changed(selected_indices)

    def _on_hardware_progress(self, current: int, total: int):
        """Handle progress update from hardware acquisition."""
        progress = int(100 * current / total) if total > 0 else 0
        self.progress_bar.setValue(progress)
        self.progress_label.setText(f"Point {current} / {total}")
        self.progress_updated.emit(current, total)

    def _on_hardware_data_point(self, point_idx: int, param_values: dict, signal: float, metadata: dict):
        """Handle data point from hardware acquisition."""
        # Extract the sweep parameter value
        sweep_param = self.config['sweep_parameter']
        param_value = param_values.get(sweep_param, 0.0)

        # Store data
        self.sweep_params.append(param_value)
        self.sweep_signals.append(signal)

        # Update plot (apply scale for display)
        display_x = param_value * self.config['x_scale']
        self.live_plot.add_data_point(display_x, signal)

        # Emit signal
        self.data_point_acquired.emit(point_idx, param_value, signal)

    def _on_hardware_finished(self):
        """Handle hardware acquisition completion."""
        self._finish_acquisition()

    def _on_hardware_error(self, error_msg: str):
        """Handle hardware acquisition error."""
        self.status_label.setText(f"Status: Error - {error_msg}")
        self._stop_acquisition()

    def _stop_acquisition(self):
        """Override to stop acquisition timer."""
        self.acquisition_timer.stop()

        # Stop hardware controller if running
        if self.use_real_hardware and self.gui_controller is not None:
            self.gui_controller.stop_acquisition()

        super()._stop_acquisition()

    def _pause_acquisition(self):
        """Override to handle hardware pause."""
        super()._pause_acquisition()

        if self.use_real_hardware and self.gui_controller is not None:
            if self.acquisition_paused_state:
                self.gui_controller.pause_acquisition()
            else:
                self.gui_controller.resume_acquisition()


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    # Create experiment window with ESR parameters
    params = {
        'experiment_type': 'ESR',
        'frequency': 2.87e9,
        'power': 8.0,
        'num_averages': 1000,
        'sweep_start': 2.85e9,
        'sweep_stop': 2.89e9,
        'sweep_points': 51
    }

    window = ExperimentWindow(window_id="ESR #1", parameters=params)
    window.show()

    sys.exit(app.exec())
