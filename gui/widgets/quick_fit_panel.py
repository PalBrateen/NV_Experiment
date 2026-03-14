"""
Quick Fit Panel Widget for NV Experiment GUI

On-the-go fitting for experiment data with multiple fit types.

Supported fits:
- Lorentzian: ESR experiments (dip/peak)
- Exponential Decay: T1, T2 experiments
- Rabi: Damped cosine oscillation
- Ramsey: Decaying oscillation with detuning

Features:
- Automatic fit on data completion (optional)
- Manual fit button
- Display fit parameters with uncertainties
- Overlay fit curve on plot
- Export fit parameters
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Dict, Any, Optional, Tuple, Callable
import numpy as np
from scipy.optimize import curve_fit

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QPushButton, QTableWidget, QTableWidgetItem,
    QCheckBox, QHeaderView
)
from PySide6.QtCore import Qt, Signal

from gui.styles.colors import *


class QuickFitPanel(QWidget):
    """
    Quick fitting panel for experiment data.

    Signals:
        fit_completed: Emitted when fit completes (fit_params, uncertainties, fit_curve)
        fit_failed: Emitted when fit fails (error_message)
    """

    fit_completed = Signal(dict, dict, tuple)  # params, uncertainties, (x_fit, y_fit)
    fit_failed = Signal(str)

    def __init__(self, parent=None):
        """Initialize quick fit panel."""
        super().__init__(parent)

        # State
        self.current_fit_type = 'lorentzian'
        self.data_x = None
        self.data_y = None
        self.fit_params = {}
        self.fit_uncertainties = {}
        self.fit_curve = None

        # Create UI
        self._create_ui()

    def _create_ui(self):
        """Create fit panel UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(10)

        # Fit controls
        control_layout = QHBoxLayout()

        self.fit_btn = QPushButton("Fit Data")
        self.fit_btn.setProperty("class", "primary")
        self.fit_btn.setEnabled(False)
        self.fit_btn.clicked.connect(self._perform_fit)
        control_layout.addWidget(self.fit_btn)

        self.clear_fit_btn = QPushButton("Clear Fit")
        self.clear_fit_btn.setEnabled(False)
        self.clear_fit_btn.clicked.connect(self._clear_fit)
        control_layout.addWidget(self.clear_fit_btn)

        control_layout.addStretch()

        layout.addLayout(control_layout)

        # Fit parameters table
        param_group = QGroupBox("Fit Parameters")
        param_layout = QVBoxLayout()

        self.param_table = QTableWidget()
        self.param_table.setColumnCount(3)
        self.param_table.setHorizontalHeaderLabels(["Parameter", "Value", "Uncertainty"])
        self.param_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.param_table.setMaximumHeight(200)
        self.param_table.setAlternatingRowColors(True)
        param_layout.addWidget(self.param_table)

        param_group.setLayout(param_layout)
        layout.addWidget(param_group)

        # Fit quality
        quality_layout = QHBoxLayout()
        quality_layout.addWidget(QLabel("R²:"))
        self.r2_label = QLabel("N/A")
        self.r2_label.setStyleSheet(f"color: {STATUS_INFO}; font-weight: bold;")
        quality_layout.addWidget(self.r2_label)
        quality_layout.addStretch()
        layout.addLayout(quality_layout)

        layout.addStretch()

        self.setLayout(layout)

    def set_fit_type(self, fit_type: str):
        """
        Set the fit function type.

        Args:
            fit_type: One of 'lorentzian', 'exponential', 'rabi', 'ramsey'
        """
        self.current_fit_type = fit_type.lower()

    def set_data(self, x: np.ndarray, y: np.ndarray):
        """
        Set data for fitting.

        Args:
            x: X data (parameter values)
            y: Y data (signal values)
        """
        self.data_x = np.array(x)
        self.data_y = np.array(y)

        # Enable fit button if we have enough data
        if len(self.data_x) >= 5:  # Need at least 5 points
            self.fit_btn.setEnabled(True)
        else:
            self.fit_btn.setEnabled(False)

    def _perform_fit(self):
        """Perform fitting based on current fit type."""
        if self.data_x is None or self.data_y is None:
            self.fit_failed.emit("No data available for fitting")
            return

        if len(self.data_x) < 5:
            self.fit_failed.emit("Not enough data points (need at least 5)")
            return

        try:
            # Select fit function and initial parameters
            if self.current_fit_type == 'lorentzian':
                fit_func, p0, param_names = self._lorentzian_setup()
            elif self.current_fit_type == 'exponential':
                fit_func, p0, param_names = self._exponential_setup()
            elif self.current_fit_type == 'rabi':
                fit_func, p0, param_names = self._rabi_setup()
            elif self.current_fit_type == 'ramsey':
                fit_func, p0, param_names = self._ramsey_setup()
            else:
                self.fit_failed.emit(f"Unknown fit type: {self.current_fit_type}")
                return

            # Perform fit
            popt, pcov = curve_fit(fit_func, self.data_x, self.data_y, p0=p0, maxfev=10000)

            # Calculate uncertainties
            perr = np.sqrt(np.diag(pcov))

            # Store results
            self.fit_params = {name: val for name, val in zip(param_names, popt)}
            self.fit_uncertainties = {name: val for name, val in zip(param_names, perr)}

            # Generate fit curve (use more points for smooth curve)
            x_fit = np.linspace(self.data_x.min(), self.data_x.max(), 500)
            y_fit = fit_func(x_fit, *popt)
            self.fit_curve = (x_fit, y_fit)

            # Calculate R²
            y_pred = fit_func(self.data_x, *popt)
            ss_res = np.sum((self.data_y - y_pred)**2)
            ss_tot = np.sum((self.data_y - np.mean(self.data_y))**2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

            # Update UI
            self._update_parameter_table()
            self.r2_label.setText(f"{r2:.6f}")
            self.clear_fit_btn.setEnabled(True)

            # Emit signal
            self.fit_completed.emit(self.fit_params, self.fit_uncertainties, self.fit_curve)

        except Exception as e:
            self.fit_failed.emit(f"Fit failed: {str(e)}")

    def _lorentzian_setup(self) -> Tuple[Callable, list, list]:
        """Setup for Lorentzian fit (ESR)."""
        # Lorentzian: y = baseline - amplitude / (1 + ((x - center) / width)^2)
        def lorentzian(x, baseline, amplitude, center, width):
            return baseline - amplitude / (1 + ((x - center) / width)**2)

        # Initial guess
        baseline_guess = np.max(self.data_y)
        amplitude_guess = baseline_guess - np.min(self.data_y)
        center_guess = self.data_x[np.argmin(self.data_y)]
        width_guess = (self.data_x[-1] - self.data_x[0]) / 10

        p0 = [baseline_guess, amplitude_guess, center_guess, width_guess]
        param_names = ['baseline', 'amplitude', 'center', 'width']

        return lorentzian, p0, param_names

    def _exponential_setup(self) -> Tuple[Callable, list, list]:
        """Setup for exponential decay fit (T1, T2)."""
        # Exponential: y = baseline - amplitude * exp(-x / tau)
        def exponential(x, baseline, amplitude, tau):
            return baseline - amplitude * np.exp(-x / tau)

        # Initial guess
        baseline_guess = self.data_y[-1]  # Final value
        amplitude_guess = self.data_y[0] - baseline_guess
        tau_guess = self.data_x[-1] / 3  # Roughly 1/3 of range

        p0 = [baseline_guess, amplitude_guess, tau_guess]
        param_names = ['baseline', 'amplitude', 'tau']

        return exponential, p0, param_names

    def _rabi_setup(self) -> Tuple[Callable, list, list]:
        """Setup for Rabi oscillation fit."""
        # Rabi: y = baseline - amplitude * (1 - cos(2*pi*x / period)) / 2
        def rabi(x, baseline, amplitude, period, phase):
            return baseline - amplitude * (1 - np.cos(2*np.pi*x / period + phase)) / 2

        # Initial guess (need to estimate period from data)
        baseline_guess = np.mean(self.data_y)
        amplitude_guess = np.max(self.data_y) - np.min(self.data_y)

        # Estimate period from zero crossings or FFT
        # Simple approach: use range / 2 as initial guess
        period_guess = (self.data_x[-1] - self.data_x[0]) / 2
        phase_guess = 0.0

        p0 = [baseline_guess, amplitude_guess, period_guess, phase_guess]
        param_names = ['baseline', 'amplitude', 'period', 'phase']

        return rabi, p0, param_names

    def _ramsey_setup(self) -> Tuple[Callable, list, list]:
        """Setup for Ramsey fringes fit."""
        # Ramsey: y = baseline - amplitude * exp(-x / T2_star) * cos(2*pi*detuning*x + phase)
        def ramsey(x, baseline, amplitude, T2_star, detuning, phase):
            return baseline - amplitude * np.exp(-x / T2_star) * np.cos(2*np.pi*detuning*x + phase)

        # Initial guess
        baseline_guess = np.mean(self.data_y)
        amplitude_guess = np.max(self.data_y) - np.min(self.data_y)
        T2_star_guess = self.data_x[-1] / 3
        detuning_guess = 1e6  # 1 MHz default
        phase_guess = 0.0

        p0 = [baseline_guess, amplitude_guess, T2_star_guess, detuning_guess, phase_guess]
        param_names = ['baseline', 'amplitude', 'T2_star', 'detuning', 'phase']

        return ramsey, p0, param_names

    def _update_parameter_table(self):
        """Update parameter table with fit results."""
        self.param_table.setRowCount(len(self.fit_params))

        for i, (name, value) in enumerate(self.fit_params.items()):
            # Parameter name
            name_item = QTableWidgetItem(name.replace('_', ' ').title())
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self.param_table.setItem(i, 0, name_item)

            # Value (scientific notation for small/large values)
            if abs(value) < 0.01 or abs(value) > 1000:
                value_str = f"{value:.6e}"
            else:
                value_str = f"{value:.6f}"
            value_item = QTableWidgetItem(value_str)
            value_item.setFlags(value_item.flags() & ~Qt.ItemIsEditable)
            self.param_table.setItem(i, 1, value_item)

            # Uncertainty
            uncert = self.fit_uncertainties.get(name, 0)
            if abs(uncert) < 0.01 or abs(uncert) > 1000:
                uncert_str = f"{uncert:.6e}"
            else:
                uncert_str = f"{uncert:.6f}"
            uncert_item = QTableWidgetItem(f"± {uncert_str}")
            uncert_item.setFlags(uncert_item.flags() & ~Qt.ItemIsEditable)
            self.param_table.setItem(i, 2, uncert_item)

    def _clear_fit(self):
        """Clear current fit results."""
        self.fit_params = {}
        self.fit_uncertainties = {}
        self.fit_curve = None
        self.param_table.setRowCount(0)
        self.r2_label.setText("N/A")
        self.clear_fit_btn.setEnabled(False)

    def get_fit_results(self) -> Dict[str, Any]:
        """
        Get current fit results.

        Returns:
            Dictionary with fit parameters, uncertainties, and curve
        """
        return {
            'parameters': self.fit_params.copy(),
            'uncertainties': self.fit_uncertainties.copy(),
            'curve': self.fit_curve
        }


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    panel = QuickFitPanel()

    # Generate test ESR data
    x_data = np.linspace(2.85e9, 2.89e9, 51)
    center = 2.87e9
    width = 5e6
    amplitude = 0.5
    baseline = 2.5
    y_data = baseline - amplitude / (1 + ((x_data - center) / width)**2) + 0.05 * np.random.randn(len(x_data))

    panel.set_fit_type('lorentzian')
    panel.set_data(x_data, y_data)

    def on_fit_completed(params, uncertainties, curve):
        print("Fit completed!")
        print("Parameters:", params)
        print("Uncertainties:", uncertainties)

    def on_fit_failed(error):
        print(f"Fit failed: {error}")

    panel.fit_completed.connect(on_fit_completed)
    panel.fit_failed.connect(on_fit_failed)

    panel.show()

    sys.exit(app.exec())
