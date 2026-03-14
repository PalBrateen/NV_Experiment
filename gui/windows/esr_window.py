"""
ESR Experiment Window

Electron Spin Resonance experiment with frequency sweep and Lorentzian fitting.

Features:
- Frequency sweep configuration
- Live spectrum display
- On-the-go Lorentzian fitting
- Multiple scan accumulation
- Scan selector for comparison
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Dict, Any, Optional
import numpy as np
import time

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QCheckBox
)
from PySide6.QtCore import QTimer

from gui.windows.base_experiment_window import BaseExperimentWindow
from gui.widgets.sweep_builder import SweepBuilderWidget


class ESRWindow(BaseExperimentWindow):
    """
    ESR experiment window with frequency sweep.

    Typical parameters:
    - Frequency: 2.85 - 2.89 GHz
    - Power: 8 dBm
    - Points: 51
    - Averages: 1000
    """

    def __init__(self, window_id: str, parent=None):
        """
        Initialize ESR window.

        Args:
            window_id: Unique identifier
            parent: Parent widget
        """
        # Set required instruments before calling super().__init__
        self.required_instruments = ['signal_generator', 'pulseblaster', 'analog_input']

        super().__init__(
            window_id=window_id,
            experiment_type="ESR",
            parent=parent
        )

        # ESR-specific parameters
        self.center_frequency = 2.87e9  # Hz
        self.span_frequency = 40e6  # Hz (40 MHz span)
        self.num_points = 51
        self.num_averages = 1000

        # Acquisition simulation timer
        self.acquisition_timer = QTimer()
        self.acquisition_timer.timeout.connect(self._simulate_acquisition_step)
        self.current_sweep_index = 0

    def _create_controls_panel(self) -> QWidget:
        """Override to add ESR-specific controls."""
        # Get base controls
        widget = super()._create_controls_panel()

        # Add ESR-specific parameters to parameters_layout
        # Frequency center
        freq_layout = QHBoxLayout()
        freq_layout.addWidget(QLabel("Center Freq:"))
        self.center_freq_label = QLabel("2.870 GHz")
        self.center_freq_label.setStyleSheet("font-weight: bold;")
        freq_layout.addWidget(self.center_freq_label)
        freq_layout.addStretch()
        self.parameters_layout.addLayout(freq_layout)

        # Frequency span
        span_layout = QHBoxLayout()
        span_layout.addWidget(QLabel("Span:"))
        self.span_label = QLabel("40 MHz")
        self.span_label.setStyleSheet("font-weight: bold;")
        span_layout.addWidget(self.span_label)
        span_layout.addStretch()
        self.parameters_layout.addLayout(span_layout)

        # Number of points
        points_layout = QHBoxLayout()
        points_layout.addWidget(QLabel("Points:"))
        self.points_label = QLabel("51")
        self.points_label.setStyleSheet("font-weight: bold;")
        points_layout.addWidget(self.points_label)
        points_layout.addStretch()
        self.parameters_layout.addLayout(points_layout)

        # Averages
        avg_layout = QHBoxLayout()
        avg_layout.addWidget(QLabel("Averages:"))
        self.averages_label = QLabel("1000")
        self.averages_label.setStyleSheet("font-weight: bold;")
        avg_layout.addWidget(self.averages_label)
        avg_layout.addStretch()
        self.parameters_layout.addLayout(avg_layout)

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
        self.autofit_checkbox = QCheckBox("Auto-fit (Lorentzian)")
        self.autofit_checkbox.setChecked(True)
        autofit_layout.addWidget(self.autofit_checkbox)
        autofit_layout.addStretch()
        self.parameters_layout.addLayout(autofit_layout)

        return widget

    def _on_sweep_changed(self, sweep_array: np.ndarray):
        """Handle sweep configuration change."""
        if len(sweep_array) > 0:
            self.num_points = len(sweep_array)
            self.points_label.setText(str(self.num_points))

            # Update center and span
            self.center_frequency = (sweep_array[0] + sweep_array[-1]) / 2
            self.span_frequency = sweep_array[-1] - sweep_array[0]

            self.center_freq_label.setText(f"{self.center_frequency / 1e9:.6f} GHz")
            self.span_label.setText(f"{self.span_frequency / 1e6:.2f} MHz")

    def set_parameters(self, params: Dict[str, Any]):
        """
        Override to extract ESR-specific parameters.

        Args:
            params: Parameter dictionary from main window
        """
        super().set_parameters(params)

        # Extract relevant parameters
        if 'frequency' in params:
            self.center_frequency = params['frequency']

        if 'num_averages' in params:
            self.num_averages = params['num_averages']
            self.averages_label.setText(str(self.num_averages))

        if 'sweep_start' in params and 'sweep_stop' in params and 'sweep_points' in params:
            self.sweep_builder.set_linear_sweep(
                params['sweep_start'],
                params['sweep_stop'],
                params['sweep_points']
            )

        # Update displays
        self.center_freq_label.setText(f"{self.center_frequency / 1e9:.6f} GHz")

    def _start_acquisition(self):
        """Override to implement ESR-specific acquisition."""
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

        # Update plot axes
        plot_item = self.live_plot.plot_widget.getPlotItem()
        plot_item.setLabel('bottom', 'Frequency', units='GHz')
        plot_item.setLabel('left', 'Signal', units='V')

        # Start simulated acquisition (replace with real acquisition later)
        self.acquisition_timer.start(50)  # 50ms per point (simulated)

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

        # Get current frequency
        freq = sweep_array[self.current_sweep_index]

        # Simulate ESR signal (Lorentzian dip centered at 2.87 GHz)
        center = 2.87e9
        width = 5e6  # 5 MHz
        amplitude = 0.5
        baseline = 2.5
        noise_level = 0.05

        signal = baseline - amplitude / (1 + ((freq - center) / width)**2) + noise_level * np.random.randn()

        # Store data
        self.sweep_params.append(freq)
        self.sweep_signals.append(signal)

        # Update plot (convert frequency to GHz for display)
        self.live_plot.add_data_point(freq / 1e9, signal)

        # Update progress
        progress = int(100 * (self.current_sweep_index + 1) / len(sweep_array))
        self.progress_bar.setValue(progress)
        self.progress_label.setText(f"Point {self.current_sweep_index + 1} / {len(sweep_array)}")

        # Emit signal
        self.data_point_acquired.emit(self.current_sweep_index, freq, signal)
        self.progress_updated.emit(self.current_sweep_index + 1, len(sweep_array))

        # Move to next point
        self.current_sweep_index += 1

    def _stop_acquisition(self):
        """Override to stop acquisition timer."""
        self.acquisition_timer.stop()
        super()._stop_acquisition()


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    # Create ESR window
    window = ESRWindow(window_id="ESR #1")

    # Set some parameters
    params = {
        'frequency': 2.87e9,
        'power': 8.0,
        'num_averages': 1000,
        'sweep_start': 2.85e9,
        'sweep_stop': 2.89e9,
        'sweep_points': 51
    }
    window.set_parameters(params)

    window.show()

    sys.exit(app.exec())
