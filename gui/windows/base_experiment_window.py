"""
Base Experiment Window for NV Experiment GUI

Base class for all experiment-specific windows (ESR, Rabi, T1, T2, etc.).
Handles resource management, parameter storage, acquisition control, and plotting.

Features:
- Resource claim/release on window open/close
- Independent parameter storage (copied from main window)
- Live plotting with FFT
- Progress tracking
- Pause/resume/stop controls
- Integration with resource manager
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Dict, Any, Optional, List
import numpy as np
from datetime import datetime

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QProgressBar, QGroupBox,
    QSplitter, QMessageBox
)
from PySide6.QtCore import Qt, Signal, QTimer

from gui.resource_manager import get_instrument_manager, InstrumentState
from gui.widgets.live_plot_widget import LivePlotWidget
from gui.styles.colors import *


class BaseExperimentWindow(QMainWindow):
    """
    Base class for experiment windows.

    Signals:
        acquisition_started: Emitted when acquisition starts
        acquisition_paused: Emitted when acquisition pauses
        acquisition_resumed: Emitted when acquisition resumes
        acquisition_stopped: Emitted when acquisition stops
        acquisition_finished: Emitted when acquisition completes
        progress_updated: Emitted with (current, total) progress
        data_point_acquired: Emitted with (index, param_value, signal_value)
    """

    acquisition_started = Signal()
    acquisition_paused = Signal()
    acquisition_resumed = Signal()
    acquisition_stopped = Signal()
    acquisition_finished = Signal()
    progress_updated = Signal(int, int)  # current, total
    data_point_acquired = Signal(int, float, float)  # index, param, signal

    def __init__(
        self,
        window_id: str,
        experiment_type: str,
        parent=None
    ):
        """
        Initialize base experiment window.

        Args:
            window_id: Unique identifier for this window
            experiment_type: Type of experiment (ESR, Rabi, T1, T2, etc.)
            parent: Parent widget (main window)
        """
        super().__init__(parent)

        self.window_id = window_id
        self.experiment_type = experiment_type

        # Window title
        self.setWindowTitle(f"{experiment_type} Experiment - {window_id}")
        self.setGeometry(150, 150, 1200, 800)

        # Resource manager
        self.resource_manager = get_instrument_manager()

        # Parameters (copied from main window)
        self.parameters = {}

        # Required instruments (to be set by subclasses)
        self.required_instruments = []

        # State
        self.instruments_claimed = False
        self.acquisition_running = False
        self.acquisition_paused_state = False

        # Data storage
        self.sweep_params = []
        self.sweep_signals = []
        self.current_scan_number = 0
        self.all_scans = []  # List of (params, signals) tuples

        # Create UI
        self._create_ui()

        # Connect to resource manager signals
        self.resource_manager.conflict_detected.connect(self._on_resource_conflict)

    def _create_ui(self):
        """Create base UI layout."""
        # Central widget
        central_widget = QWidget()
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        # Header with experiment type
        header_layout = QHBoxLayout()

        header_label = QLabel(f"{self.experiment_type} Experiment")
        header_label.setProperty("class", "header")
        header_layout.addWidget(header_label)

        header_layout.addStretch()

        # Window ID
        id_label = QLabel(f"ID: {self.window_id}")
        id_label.setStyleSheet("color: #888888;")
        header_layout.addWidget(id_label)

        main_layout.addLayout(header_layout)

        # Main content splitter (left: controls, right: plot)
        splitter = QSplitter(Qt.Horizontal)

        # Left panel: Controls
        self.controls_widget = self._create_controls_panel()
        splitter.addWidget(self.controls_widget)

        # Right panel: Plot
        self.plot_widget = self._create_plot_panel()
        splitter.addWidget(self.plot_widget)

        # Set splitter sizes (40% controls, 60% plot)
        splitter.setSizes([400, 720])

        main_layout.addWidget(splitter)

        # Status bar
        self.status_bar_widget = self._create_status_bar()
        main_layout.addWidget(self.status_bar_widget)

        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

    def _create_controls_panel(self) -> QWidget:
        """Create left panel with controls (to be extended by subclasses)."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Resource status group
        resource_group = QGroupBox("Resource Status")
        resource_layout = QVBoxLayout()

        self.resource_status_label = QLabel("Instruments: Not claimed")
        self.resource_status_label.setStyleSheet(f"color: {TEXT_DISABLED};")
        resource_layout.addWidget(self.resource_status_label)

        self.claim_btn = QPushButton("Claim Instruments")
        self.claim_btn.setProperty("class", "primary")
        self.claim_btn.clicked.connect(self._claim_instruments)
        resource_layout.addWidget(self.claim_btn)

        self.release_btn = QPushButton("Release Instruments")
        self.release_btn.setEnabled(False)
        self.release_btn.clicked.connect(self._release_instruments)
        resource_layout.addWidget(self.release_btn)

        resource_group.setLayout(resource_layout)
        layout.addWidget(resource_group)

        # Parameters group (placeholder - to be filled by subclasses)
        self.parameters_group = QGroupBox("Parameters")
        self.parameters_layout = QVBoxLayout()
        self.parameters_group.setLayout(self.parameters_layout)
        layout.addWidget(self.parameters_group)

        # Acquisition control group
        control_group = QGroupBox("Acquisition Control")
        control_layout = QVBoxLayout()

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        control_layout.addWidget(self.progress_bar)

        # Progress label
        self.progress_label = QLabel("Ready")
        control_layout.addWidget(self.progress_label)

        # Control buttons
        button_layout = QHBoxLayout()

        self.start_btn = QPushButton("▶ START")
        self.start_btn.setProperty("class", "primary")
        self.start_btn.setMinimumHeight(40)
        self.start_btn.clicked.connect(self._start_acquisition)
        button_layout.addWidget(self.start_btn)

        self.pause_btn = QPushButton("⏸ PAUSE")
        self.pause_btn.setProperty("class", "warning")
        self.pause_btn.setMinimumHeight(40)
        self.pause_btn.setEnabled(False)
        self.pause_btn.clicked.connect(self._pause_acquisition)
        button_layout.addWidget(self.pause_btn)

        self.stop_btn = QPushButton("⏹ STOP")
        self.stop_btn.setProperty("class", "danger")
        self.stop_btn.setMinimumHeight(40)
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_acquisition)
        button_layout.addWidget(self.stop_btn)

        control_layout.addLayout(button_layout)

        control_group.setLayout(control_layout)
        layout.addWidget(control_group)

        layout.addStretch()

        widget.setLayout(layout)
        return widget

    def _create_plot_panel(self) -> QWidget:
        """Create right panel with plot."""
        widget = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        # Live plot widget
        self.live_plot = LivePlotWidget(
            title=f"{self.experiment_type} Spectrum",
            y_label="Signal",
            y_unit="V",
            buffer_size=1000
        )

        # For experiment windows, we typically plot parameter vs signal
        # Override x-axis label
        plot_item = self.live_plot.plot_widget.getPlotItem()
        plot_item.setLabel('bottom', 'Parameter', units='')

        layout.addWidget(self.live_plot)

        widget.setLayout(layout)
        return widget

    def _create_status_bar(self) -> QWidget:
        """Create status bar at bottom."""
        widget = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)

        self.status_label = QLabel("Status: Idle")
        layout.addWidget(self.status_label)

        layout.addStretch()

        # Scan counter
        self.scan_label = QLabel("Scan: 0")
        layout.addWidget(self.scan_label)

        widget.setLayout(layout)
        return widget

    def _claim_instruments(self):
        """Claim required instruments."""
        if not self.required_instruments:
            QMessageBox.warning(
                self,
                "No Instruments Required",
                "This experiment does not require any instruments."
            )
            return

        success, error = self.resource_manager.claim_instruments(
            requester=self.window_id,
            instrument_names=self.required_instruments
        )

        if success:
            self.instruments_claimed = True
            self.resource_status_label.setText(
                f"Instruments: Claimed ({', '.join(self.required_instruments)})"
            )
            self.resource_status_label.setStyleSheet(f"color: {STATUS_OK};")
            self.claim_btn.setEnabled(False)
            self.release_btn.setEnabled(True)
            self.start_btn.setEnabled(True)
            self.status_label.setText("Status: Ready to acquire")
        else:
            QMessageBox.critical(
                self,
                "Resource Conflict",
                f"Failed to claim instruments:\n{error}"
            )

    def _release_instruments(self):
        """Release claimed instruments."""
        if not self.instruments_claimed:
            return

        self.resource_manager.release_instruments(
            requester=self.window_id,
            instrument_names=self.required_instruments
        )

        self.instruments_claimed = False
        self.resource_status_label.setText("Instruments: Released")
        self.resource_status_label.setStyleSheet(f"color: {TEXT_DISABLED};")
        self.claim_btn.setEnabled(True)
        self.release_btn.setEnabled(False)
        self.start_btn.setEnabled(False)
        self.status_label.setText("Status: Idle")

    def _on_resource_conflict(self, requester: str, instrument: str):
        """Handle resource conflict notification."""
        if self.instruments_claimed and instrument in self.required_instruments:
            QMessageBox.warning(
                self,
                "Resource Conflict",
                f"Another window ({requester}) is requesting {instrument}.\n"
                f"Your acquisition may be affected."
            )

    def _start_acquisition(self):
        """Start acquisition (to be implemented by subclasses)."""
        if not self.instruments_claimed:
            QMessageBox.warning(
                self,
                "Instruments Not Claimed",
                "Please claim instruments before starting acquisition."
            )
            return

        # Set UI state
        self.acquisition_running = True
        self.start_btn.setEnabled(False)
        self.pause_btn.setEnabled(True)
        self.stop_btn.setEnabled(True)
        self.status_label.setText("Status: Acquiring...")

        # Clear previous data
        self.sweep_params = []
        self.sweep_signals = []
        self.live_plot.clear_data()

        # Emit signal
        self.acquisition_started.emit()

        # Subclasses should override and implement actual acquisition

    def _pause_acquisition(self):
        """Pause/resume acquisition."""
        if self.acquisition_paused_state:
            # Resume
            self.acquisition_paused_state = False
            self.pause_btn.setText("⏸ PAUSE")
            self.status_label.setText("Status: Acquiring...")
            self.acquisition_resumed.emit()
        else:
            # Pause
            self.acquisition_paused_state = True
            self.pause_btn.setText("▶ RESUME")
            self.status_label.setText("Status: Paused")
            self.acquisition_paused.emit()

    def _stop_acquisition(self):
        """Stop acquisition."""
        self.acquisition_running = False
        self.acquisition_paused_state = False

        # Reset UI
        self.start_btn.setEnabled(True)
        self.pause_btn.setEnabled(False)
        self.pause_btn.setText("⏸ PAUSE")
        self.stop_btn.setEnabled(False)
        self.status_label.setText("Status: Stopped")

        # Emit signal
        self.acquisition_stopped.emit()

    def _finish_acquisition(self):
        """Called when acquisition completes successfully."""
        self.acquisition_running = False
        self.acquisition_paused_state = False

        # Store scan data
        if len(self.sweep_params) > 0:
            self.all_scans.append((
                np.array(self.sweep_params),
                np.array(self.sweep_signals)
            ))
            self.current_scan_number += 1
            self.scan_label.setText(f"Scan: {self.current_scan_number}")

        # Reset UI
        self.start_btn.setEnabled(True)
        self.pause_btn.setEnabled(False)
        self.pause_btn.setText("⏸ PAUSE")
        self.stop_btn.setEnabled(False)
        self.progress_bar.setValue(100)
        self.progress_label.setText("Complete!")
        self.status_label.setText("Status: Finished")

        # Emit signal
        self.acquisition_finished.emit()

    def set_parameters(self, params: Dict[str, Any]):
        """
        Set experiment parameters (copied from main window).

        Args:
            params: Parameter dictionary
        """
        self.parameters = params.copy()

    def get_parameters(self) -> Dict[str, Any]:
        """
        Get current parameters.

        Returns:
            Parameter dictionary
        """
        return self.parameters.copy()

    def closeEvent(self, event):
        """Handle window close event."""
        # Check if acquisition is running
        if self.acquisition_running:
            reply = QMessageBox.question(
                self,
                "Acquisition Running",
                "Acquisition is still running. Stop and close?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )

            if reply == QMessageBox.No:
                event.ignore()
                return

            self._stop_acquisition()

        # Release instruments
        if self.instruments_claimed:
            self._release_instruments()

        event.accept()


# Example usage for testing
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    # Create test window
    window = BaseExperimentWindow(
        window_id="Test Window #1",
        experiment_type="Test"
    )
    window.required_instruments = ['signal_generator', 'analog_input']
    window.show()

    sys.exit(app.exec())
