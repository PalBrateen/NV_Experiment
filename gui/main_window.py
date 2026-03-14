"""
NV Experiment GUI - Main Window

Main application window with:
- Tabbed left panel (Instruments, Config, Monitor, Data)
- Dockable IPython console (right side)
- Live plots (bottom, Voltage and Intensity)
- Status bar with connection indicators
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QPushButton, QStatusBar,
    QDockWidget, QSplitter, QFrame
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction

from gui.widgets.ipython_console_widget import create_nv_console
from gui.widgets.live_plot_widget import LivePlotWidget
from gui.widgets.instrument_status_card import (
    SignalGeneratorCard, PulseBlasterCard, DAQCard, CameraCard
)
from gui.widgets.parameter_editor import ParameterEditorWidget
from gui.widgets.monitor_tab import MonitorTab
from gui.widgets.multi_sweep_builder import MultiSweepBuilderWidget
from gui.resource_manager import get_instrument_manager, InstrumentState, InstrumentType


class NVExperimentGUI(QMainWindow):
    """
    Main GUI window for NV center experiment control.

    Layout:
    - Left panel (60%): Tabbed interface for instruments, config, monitor, data
    - Right panel (40%): IPython console (dockable)
    - Bottom: Live plots for voltage (left) and intensity (right) (dockable)
    - Status bar: Connection status, current state
    """

    def __init__(self):
        super().__init__()

        self.setWindowTitle("NV Experiment Control")
        self.setGeometry(100, 100, 1600, 1000)  # x, y, width, height

        # Resource manager
        self.resource_manager = get_instrument_manager()

        # Instrument instances (will be populated on connection)
        self.instruments = {}

        # Instrument status cards
        self.instrument_cards = {}

        # Live plot widgets
        self.voltage_plot = None
        self.intensity_plot = None

        # Experiment windows
        self.experiment_windows = []
        self.experiment_window_counter = 0

        # Setup UI
        self._create_menu_bar()
        self._create_main_layout()
        self._create_ipython_console_dock()
        self._create_live_plots_dock()
        self._create_status_bar()

        # Initialize instruments (simulation mode for now)
        self._initialize_simulated_instruments()

        # Start status update timer
        self.status_timer = QTimer()
        self.status_timer.timeout.connect(self._update_status)
        self.status_timer.start(1000)  # Update every 1 second

        # Test data generator (for demonstration)
        self.test_data_timer = QTimer()
        self.test_data_timer.timeout.connect(self._generate_test_data)
        self.test_data_enabled = False  # Will be enabled via menu
        self.test_time = 0.0

    def _create_menu_bar(self):
        """Create menu bar."""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("File")

        load_config_action = QAction("Load Configuration...", self)
        file_menu.addAction(load_config_action)

        save_config_action = QAction("Save Configuration...", self)
        file_menu.addAction(save_config_action)

        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Tools menu
        tools_menu = menubar.addMenu("Tools")

        # Experiment windows submenu
        exp_windows_menu = tools_menu.addMenu("Experiment Windows")

        open_esr_action = QAction("ESR Window", self)
        open_esr_action.triggered.connect(lambda: self._open_experiment_window('ESR'))
        exp_windows_menu.addAction(open_esr_action)

        open_rabi_action = QAction("Rabi Window", self)
        open_rabi_action.triggered.connect(lambda: self._open_experiment_window('Rabi'))
        exp_windows_menu.addAction(open_rabi_action)

        open_t1_action = QAction("T1 Window", self)
        open_t1_action.triggered.connect(lambda: self._open_experiment_window('T1'))
        exp_windows_menu.addAction(open_t1_action)

        open_t2_action = QAction("T2 Window", self)
        open_t2_action.triggered.connect(lambda: self._open_experiment_window('T2'))
        exp_windows_menu.addAction(open_t2_action)

        open_ramsey_action = QAction("Ramsey Window", self)
        open_ramsey_action.triggered.connect(lambda: self._open_experiment_window('Ramsey'))
        exp_windows_menu.addAction(open_ramsey_action)

        tools_menu.addSeparator()

        simulation_mode_action = QAction("Simulation Mode", self)
        simulation_mode_action.setCheckable(True)
        tools_menu.addAction(simulation_mode_action)

        test_data_action = QAction("Generate Test Data", self)
        test_data_action.setCheckable(True)
        test_data_action.toggled.connect(self._toggle_test_data)
        tools_menu.addAction(test_data_action)

        # View menu
        view_menu = menubar.addMenu("View")

        toggle_console_action = QAction("Toggle IPython Console", self)
        toggle_console_action.triggered.connect(self._toggle_console_dock)
        view_menu.addAction(toggle_console_action)

        toggle_plots_action = QAction("Toggle Live Plots", self)
        toggle_plots_action.triggered.connect(self._toggle_plots_dock)
        view_menu.addAction(toggle_plots_action)

        # Help menu
        help_menu = menubar.addMenu("Help")

        about_action = QAction("About", self)
        help_menu.addAction(about_action)

        keyboard_shortcuts_action = QAction("Keyboard Shortcuts", self)
        help_menu.addAction(keyboard_shortcuts_action)

    def _create_main_layout(self):
        """Create main central widget with tabbed interface."""
        central_widget = QWidget()
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Create tab widget for left panel
        self.tab_widget = QTabWidget()

        # Tab 1: Instruments
        self.instruments_tab = self._create_instruments_tab()
        self.tab_widget.addTab(self.instruments_tab, "Instruments")

        # Tab 2: Configuration
        self.config_tab = self._create_config_tab()
        self.tab_widget.addTab(self.config_tab, "Configuration")

        # Tab 3: Monitor (single-point/continuous)
        self.monitor_tab = self._create_monitor_tab()
        self.tab_widget.addTab(self.monitor_tab, "Monitor")

        # Tab 4: Data Management
        self.data_tab = self._create_data_tab()
        self.tab_widget.addTab(self.data_tab, "Data")

        main_layout.addWidget(self.tab_widget)
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

    def _create_instruments_tab(self) -> QWidget:
        """Create instruments status tab."""
        tab = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Header
        header = QLabel("Instrument Status")
        header.setProperty("class", "header")
        layout.addWidget(header)

        # Create status cards
        self.sg_card = SignalGeneratorCard()
        self.sg_card.reconnect_requested.connect(lambda: self._reconnect_instrument('signal_generator'))
        layout.addWidget(self.sg_card)
        self.instrument_cards['signal_generator'] = self.sg_card

        self.pb_card = PulseBlasterCard()
        self.pb_card.reconnect_requested.connect(lambda: self._reconnect_instrument('pulseblaster'))
        layout.addWidget(self.pb_card)
        self.instrument_cards['pulseblaster'] = self.pb_card

        self.daq_card = DAQCard()
        self.daq_card.reconnect_requested.connect(lambda: self._reconnect_instrument('analog_input'))
        layout.addWidget(self.daq_card)
        self.instrument_cards['analog_input'] = self.daq_card

        self.camera_card = CameraCard()
        self.camera_card.reconnect_requested.connect(lambda: self._reconnect_instrument('camera'))
        layout.addWidget(self.camera_card)
        self.instrument_cards['camera'] = self.camera_card

        layout.addStretch()
        tab.setLayout(layout)
        return tab

    def _create_config_tab(self) -> QWidget:
        """Create configuration tab."""
        tab = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        # Parameter editor widget
        self.parameter_editor = ParameterEditorWidget()
        self.parameter_editor.parameters_changed.connect(self._on_parameters_changed)
        self.parameter_editor.open_window_requested.connect(self._on_open_window_requested)
        layout.addWidget(self.parameter_editor)

        # Multi-parameter sweep builder (added below parameter editor)
        self.multi_sweep_builder = MultiSweepBuilderWidget(experiment_type='ESR')
        self.multi_sweep_builder.sweep_changed.connect(self._on_multi_sweep_changed)
        layout.addWidget(self.multi_sweep_builder)

        tab.setLayout(layout)
        return tab

    def _create_monitor_tab(self) -> QWidget:
        """Create monitor tab for single-point/continuous acquisition."""
        # Create monitor tab widget
        self.monitor_tab = MonitorTab()

        # Connect signals
        self.monitor_tab.single_point_requested.connect(self._on_single_point_acquisition)
        self.monitor_tab.continuous_started.connect(self._on_continuous_started)
        self.monitor_tab.continuous_stopped.connect(self._on_continuous_stopped)
        self.monitor_tab.field_changed.connect(self._on_field_changed)
        self.monitor_tab.data_point_acquired.connect(self._on_monitor_data_acquired)

        return self.monitor_tab

    def _create_data_tab(self) -> QWidget:
        """Create data management tab."""
        tab = QWidget()
        layout = QVBoxLayout()

        # Header
        header = QLabel("Data Management")
        header.setProperty("class", "header")
        layout.addWidget(header)

        # Placeholder
        placeholder = QLabel("Data management interface will appear here.\n\n"
                           "Features:\n"
                           "• Save directory configuration\n"
                           "• Filename prefix and timestamp\n"
                           "• Format selection (YAML, NPY, PNG)\n"
                           "• Auto-save settings\n"
                           "• Recent files browser")
        placeholder.setAlignment(Qt.AlignCenter)
        placeholder.setStyleSheet("color: #888888; font-size: 11pt; padding: 50px;")
        layout.addWidget(placeholder)

        layout.addStretch()
        tab.setLayout(layout)
        return tab

    def _create_ipython_console_dock(self):
        """Create dockable IPython console."""
        # Create dock widget
        self.console_dock = QDockWidget("IPython Console", self)
        self.console_dock.setAllowedAreas(Qt.RightDockWidgetArea | Qt.LeftDockWidgetArea)

        # Create console widget
        self.ipython_console = create_nv_console(
            gui_instance=self,
            instruments=self.instruments,
            parent=self.console_dock
        )

        self.console_dock.setWidget(self.ipython_console)

        # Add to main window (right side)
        self.addDockWidget(Qt.RightDockWidgetArea, self.console_dock)

    def _create_live_plots_dock(self):
        """Create dockable live plots (voltage and intensity)."""
        # Create dock widget
        self.plots_dock = QDockWidget("Live Plots", self)
        self.plots_dock.setAllowedAreas(Qt.BottomDockWidgetArea | Qt.TopDockWidgetArea)

        # Container widget
        plots_widget = QWidget()
        plots_layout = QHBoxLayout()
        plots_layout.setContentsMargins(5, 5, 5, 5)
        plots_layout.setSpacing(10)

        # Voltage plot (left)
        self.voltage_plot = LivePlotWidget(
            title="Live Voltage",
            y_label="Voltage",
            y_unit="V",
            buffer_size=1000
        )
        plots_layout.addWidget(self.voltage_plot)

        # Intensity plot (right)
        self.intensity_plot = LivePlotWidget(
            title="Live Intensity",
            y_label="Intensity",
            y_unit="counts",
            buffer_size=1000
        )
        plots_layout.addWidget(self.intensity_plot)

        plots_widget.setLayout(plots_layout)
        self.plots_dock.setWidget(plots_widget)

        # Add to main window (bottom)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.plots_dock)

        # Set initial size (1/3 of window height)
        self.plots_dock.setMinimumHeight(300)

    def _create_status_bar(self):
        """Create status bar with connection indicators."""
        self.status_bar = QStatusBar()

        # Status message (left side)
        self.status_label = QLabel("Status: Idle")
        self.status_bar.addWidget(self.status_label)

        # Connection indicators (right side)
        self.sg_indicator = QLabel("⚫ SG")
        self.sg_indicator.setToolTip("Signal Generator: Disconnected")
        self.status_bar.addPermanentWidget(self.sg_indicator)

        self.pb_indicator = QLabel("⚫ PB")
        self.pb_indicator.setToolTip("PulseBlaster: Disconnected")
        self.status_bar.addPermanentWidget(self.pb_indicator)

        self.daq_indicator = QLabel("⚫ DAQ")
        self.daq_indicator.setToolTip("DAQ: Disconnected")
        self.status_bar.addPermanentWidget(self.daq_indicator)

        self.camera_indicator = QLabel("⚫ CAM")
        self.camera_indicator.setToolTip("Camera: Disconnected")
        self.status_bar.addPermanentWidget(self.camera_indicator)

        self.setStatusBar(self.status_bar)

    def _toggle_console_dock(self):
        """Toggle IPython console visibility."""
        if self.console_dock.isVisible():
            self.console_dock.hide()
        else:
            self.console_dock.show()

    def _toggle_plots_dock(self):
        """Toggle live plots visibility."""
        if self.plots_dock.isVisible():
            self.plots_dock.hide()
        else:
            self.plots_dock.show()

    def _update_status(self):
        """Update status bar indicators (called every 1 second)."""
        # Get instrument statuses from resource manager
        instruments_status = self.resource_manager.get_all_instruments()

        # Update status bar indicators
        self._update_indicator(
            self.sg_indicator,
            'signal_generator',
            instruments_status,
            "Signal Generator"
        )
        self._update_indicator(
            self.pb_indicator,
            'pulseblaster',
            instruments_status,
            "PulseBlaster"
        )
        self._update_indicator(
            self.daq_indicator,
            'analog_input',
            instruments_status,
            "DAQ"
        )
        self._update_indicator(
            self.camera_indicator,
            'camera',
            instruments_status,
            "Camera"
        )

    def _update_indicator(self, label: QLabel, inst_name: str, statuses: dict, display_name: str):
        """Update a single status indicator."""
        if inst_name in statuses:
            state = statuses[inst_name]['state']
            owner = statuses[inst_name]['owner']

            if state == InstrumentState.FREE:
                label.setText(f"🟢 {display_name.split()[0]}")
                tooltip = f"{display_name}: Connected (FREE)"
            elif state == InstrumentState.CLAIMED:
                label.setText(f"🟡 {display_name.split()[0]}")
                tooltip = f"{display_name}: Claimed by {owner}"
            elif state == InstrumentState.IN_USE:
                label.setText(f"🟡 {display_name.split()[0]}")
                tooltip = f"{display_name}: In use by {owner}"
            elif state == InstrumentState.ERROR:
                label.setText(f"🔴 {display_name.split()[0]}")
                tooltip = f"{display_name}: Error"
            else:  # DISCONNECTED
                label.setText(f"⚫ {display_name.split()[0]}")
                tooltip = f"{display_name}: Disconnected"

            label.setToolTip(tooltip)
        else:
            label.setText(f"⚫ {display_name.split()[0]}")
            label.setToolTip(f"{display_name}: Not registered")

    def _initialize_simulated_instruments(self):
        """Initialize simulated instruments for testing."""
        # Register simulated instruments with resource manager
        # These will be replaced with real instruments when connected

        from SGcontrol import SignalGenerator_sim
        from PBcontrol import PulseBlaster

        # Signal Generator (simulation)
        try:
            sg_sim = SignalGenerator_sim(name="sg1")
            self.resource_manager.register_instrument(
                'signal_generator',
                sg_sim,
                InstrumentType.SIGNAL_GENERATOR,
                InstrumentState.DISCONNECTED  # Start as disconnected
            )
            self.instruments['sg'] = sg_sim
        except Exception as e:
            print(f"⚠️ Could not initialize SG simulation: {e}")

        # PulseBlaster, DAQ, Camera will be registered when connected
        # For now, just register placeholders
        self.resource_manager.register_instrument(
            'pulseblaster',
            None,
            InstrumentType.PULSEBLASTER,
            InstrumentState.DISCONNECTED
        )

        self.resource_manager.register_instrument(
            'analog_input',
            None,
            InstrumentType.ANALOG_INPUT,
            InstrumentState.DISCONNECTED
        )

        self.resource_manager.register_instrument(
            'camera',
            None,
            InstrumentType.CAMERA,
            InstrumentState.DISCONNECTED
        )

        # Update instrument cards
        self._update_instrument_cards()

        # Update console namespace
        self.update_console_namespace({'sg': self.instruments.get('sg')})

    def _reconnect_instrument(self, instrument_name: str):
        """
        Reconnect to an instrument.

        Args:
            instrument_name: Name of instrument to reconnect
        """
        print(f"Reconnecting to {instrument_name}...")
        # TODO: Implement actual reconnection logic
        # For now, just update state to FREE
        self.resource_manager.set_instrument_state(
            instrument_name,
            InstrumentState.FREE
        )
        self._update_instrument_cards()

    def _update_instrument_cards(self):
        """Update all instrument status cards."""
        instruments_status = self.resource_manager.get_all_instruments()

        for inst_name, card in self.instrument_cards.items():
            if inst_name in instruments_status:
                status = instruments_status[inst_name]
                card.update_state(status['state'], status['owner'])

                # Update specific parameters
                if inst_name == 'signal_generator' and 'sg' in self.instruments:
                    sg = self.instruments['sg']
                    if hasattr(sg, 'frequency'):
                        card.update_sg_parameters(
                            frequency=getattr(sg, 'frequency', 2.87e9),
                            power=getattr(sg, 'power', 8.0),
                            output_enabled=True,
                            modulation="Pulse"
                        )

    def _toggle_test_data(self, enabled: bool):
        """Toggle test data generation for live plots."""
        import numpy as np

        self.test_data_enabled = enabled
        if enabled:
            self.test_time = 0.0
            self.voltage_plot.set_sampling_rate(100.0)  # 100 Hz
            self.intensity_plot.set_sampling_rate(100.0)
            self.test_data_timer.start(50)  # 20 Hz update
            self.status_label.setText("Status: Generating test data...")
        else:
            self.test_data_timer.stop()
            self.status_label.setText("Status: Idle")

    def _generate_test_data(self):
        """Generate test data for live plots."""
        import numpy as np
        import time

        # Generate simulated voltage signal (10 Hz sine + noise)
        t = self.test_time
        voltage = 2.5 + 0.5 * np.sin(2 * np.pi * 1.0 * t) + 0.1 * np.random.randn()

        # Generate simulated intensity signal (5 Hz sine + noise)
        intensity = 1000 + 200 * np.sin(2 * np.pi * 0.5 * t + np.pi/4) + 20 * np.random.randn()

        # Add to plots
        self.voltage_plot.add_data_point(t, voltage)
        self.intensity_plot.add_data_point(t, intensity)

        # Increment time
        self.test_time += 0.05  # 50ms timestep

    def _on_parameters_changed(self, params: dict):
        """
        Handle parameter changes from editor.

        Args:
            params: Updated parameters dictionary
        """
        print(f"Parameters updated: {params.get('experiment_type', 'Unknown')}")
        # Store parameters for use by experiment windows
        self.current_parameters = params

        # Update multi-sweep builder experiment type
        exp_type = params.get('experiment_type', 'ESR')
        self.multi_sweep_builder.set_experiment_type(exp_type)

    def _on_multi_sweep_changed(self, sweep_params: list):
        """
        Handle multi-parameter sweep configuration changes.

        Args:
            sweep_params: List of SweepParameter objects
        """
        n_params = len(sweep_params)
        if n_params == 0:
            print("Multi-sweep: No parameters configured (single point)")
        elif n_params == 1:
            param = sweep_params[0]
            print(f"Multi-sweep: 1D sweep of {param.name} ({param.points} points)")
        elif n_params == 2:
            inner = sweep_params[0]
            outer = sweep_params[1]
            total = inner.points * outer.points
            print(f"Multi-sweep: 2D sweep - {inner.name} × {outer.name} ({total} points)")
        else:
            total = 1
            for p in sweep_params:
                total *= p.points
            print(f"Multi-sweep: {n_params}D sweep ({total} total points)")

        # Store for use when opening experiment windows
        if hasattr(self, 'current_parameters'):
            self.current_parameters['multi_sweep_params'] = sweep_params

    def _on_single_point_acquisition(self):
        """Handle single-point acquisition request from monitor tab."""
        print("Single-point acquisition requested")
        # TODO: Implement actual single-point acquisition
        # For now, just simulate
        import time
        import numpy as np
        timestamp = time.time()
        voltage = 2.5 + 0.3 * np.random.randn()
        intensity = 1000 + 50 * np.random.randn()

        self.monitor_tab.update_readings(timestamp, voltage, intensity)
        self.status_label.setText("Status: Single-point acquisition complete")

    def _on_continuous_started(self):
        """Handle continuous acquisition start from monitor tab."""
        print("Continuous acquisition started")
        self.status_label.setText("Status: Continuous acquisition running...")

    def _on_continuous_stopped(self):
        """Handle continuous acquisition stop from monitor tab."""
        print("Continuous acquisition stopped")
        self.status_label.setText("Status: Idle")

    def _on_field_changed(self, bx: float, by: float, bz: float):
        """
        Handle magnetic field change from monitor tab.

        Args:
            bx, by, bz: Magnetic field components in Gauss
        """
        print(f"Magnetic field changed: Bx={bx:.1f} G, By={by:.1f} G, Bz={bz:.1f} G")
        # TODO: Apply field to DAQ output
        self.status_label.setText(f"Status: Field applied (Bx={bx:.1f}, By={by:.1f}, Bz={bz:.1f} G)")

    def _on_monitor_data_acquired(self, timestamp: float, voltage: float, intensity: float):
        """
        Handle data point acquired from monitor tab.
        Send to live plots.

        Args:
            timestamp: Acquisition timestamp
            voltage: Voltage value
            intensity: Intensity value
        """
        # Add to live plots
        self.voltage_plot.add_data_point(timestamp, voltage)
        self.intensity_plot.add_data_point(timestamp, intensity)

    def _on_open_window_requested(self, params: dict):
        """
        Handle Open Window button click from Config tab.

        Args:
            params: Parameters from parameter editor
        """
        from gui.windows.experiment_window import ExperimentWindow

        # Generate unique window ID
        exp_type = params.get('experiment_type', 'ESR')
        self.experiment_window_counter += 1
        window_id = f"{exp_type} #{self.experiment_window_counter}"

        # Create generic experiment window with parameters
        window = ExperimentWindow(window_id=window_id, parameters=params, parent=self)

        # Track window
        self.experiment_windows.append(window)

        # Clean up when window closes
        def on_window_closed():
            if window in self.experiment_windows:
                self.experiment_windows.remove(window)

        window.destroyed.connect(on_window_closed)

        # Show window
        window.show()

        self.status_label.setText(f"Status: Opened {window_id}")

    def _open_experiment_window(self, exp_type: str):
        """
        Open a new experiment window (legacy method - kept for menu compatibility).

        Args:
            exp_type: Experiment type ('ESR', 'Rabi', 'T1', 'T2')
        """
        # Get current parameters from editor
        if hasattr(self, 'parameter_editor'):
            params = self.parameter_editor.get_parameters()
            # Override experiment type
            params['experiment_type'] = exp_type
            self._on_open_window_requested(params)
        else:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self,
                "No Parameters",
                "Please configure parameters in Config tab first."
            )

    def update_console_namespace(self, updates: dict):
        """
        Update IPython console namespace with new variables.

        Args:
            updates: Dictionary of variables to add/update
        """
        if hasattr(self, 'ipython_console'):
            self.ipython_console.update_namespace(updates)

    def closeEvent(self, event):
        """Handle window close event."""
        # Stop timers
        if hasattr(self, 'status_timer'):
            self.status_timer.stop()
        if hasattr(self, 'test_data_timer'):
            self.test_data_timer.stop()

        # Clean up IPython console
        if hasattr(self, 'ipython_console'):
            self.ipython_console.shutdown()

        # Clean up resource manager
        self.resource_manager.cleanup()

        event.accept()


# For testing standalone
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication
    from pathlib import Path

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    window = NVExperimentGUI()
    window.show()

    sys.exit(app.exec())
