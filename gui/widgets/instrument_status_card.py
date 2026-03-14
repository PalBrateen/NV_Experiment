"""
Instrument Status Card Widget

Visual card displaying instrument connection status, current parameters,
ownership information, and control buttons.

Features:
- Color-coded status indicators (green/yellow/red/gray)
- Real-time parameter display
- Ownership tracking (which window is using it)
- Reconnect/disconnect buttons
- Simulation mode toggle
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import Optional, Dict, Any, Tuple, List
from PySide6.QtWidgets import (
    QWidget, QFrame, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QCheckBox
)
from PySide6.QtCore import Qt, Signal

from gui.resource_manager import InstrumentState, InstrumentType
from gui.styles.colors import *


class InstrumentStatusCard(QFrame):
    """
    Status card for a single instrument.

    Signals:
        reconnect_requested: Emitted when user clicks reconnect
        disconnect_requested: Emitted when user clicks disconnect
        simulation_toggled: Emitted when simulation mode is toggled
        test_requested: Emitted when user clicks test connection
    """

    reconnect_requested = Signal()
    disconnect_requested = Signal()
    simulation_toggled = Signal(bool)  # True if simulation enabled
    test_requested = Signal()

    def __init__(
        self,
        instrument_name: str,
        instrument_type: InstrumentType,
        parent=None
    ):
        """
        Initialize instrument status card.

        Args:
            instrument_name: Display name (e.g., "Signal Generator")
            instrument_type: Type of instrument
            parent: Parent widget
        """
        super().__init__(parent)

        self.instrument_name = instrument_name
        self.instrument_type = instrument_type

        # Current state
        self.state = InstrumentState.DISCONNECTED
        self.owner = None
        self.parameters = {}

        # Setup UI
        self.setProperty("class", "card")
        self.setFrameStyle(QFrame.Box | QFrame.Raised)
        self._create_ui()

    def _create_ui(self):
        """Create card UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)

        # Header with name and status indicator
        header_layout = QHBoxLayout()

        # Status indicator (colored circle)
        self.status_indicator = QLabel("●")
        self.status_indicator.setProperty("class", "status_indicator")
        self.status_indicator.setStyleSheet(f"color: {TEXT_DISABLED};")  # Gray by default
        header_layout.addWidget(self.status_indicator)

        # Instrument name
        name_label = QLabel(self.instrument_name)
        name_label.setProperty("class", "subheader")
        header_layout.addWidget(name_label)

        header_layout.addStretch()

        # State label
        self.state_label = QLabel("Disconnected")
        self.state_label.setProperty("class", "value")
        header_layout.addWidget(self.state_label)

        layout.addLayout(header_layout)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        layout.addWidget(separator)

        # Parameters grid
        self.params_grid = QGridLayout()
        self.params_grid.setSpacing(5)

        # Owner row
        owner_label = QLabel("Owner:")
        owner_label.setStyleSheet("font-weight: bold;")
        self.owner_value = QLabel("None")
        self.owner_value.setProperty("class", "value")
        self.params_grid.addWidget(owner_label, 0, 0)
        self.params_grid.addWidget(self.owner_value, 0, 1)

        # State row
        state_label = QLabel("State:")
        state_label.setStyleSheet("font-weight: bold;")
        self.state_value = QLabel("FREE")
        self.state_value.setProperty("class", "value")
        self.params_grid.addWidget(state_label, 1, 0)
        self.params_grid.addWidget(self.state_value, 1, 1)

        # Parameter rows (will be added dynamically)
        self.param_labels = {}
        self.param_values = {}
        self.next_param_row = 2

        layout.addLayout(self.params_grid)

        layout.addStretch()

        # Control buttons
        button_layout = QHBoxLayout()

        self.reconnect_btn = QPushButton("Reconnect")
        self.reconnect_btn.clicked.connect(self.reconnect_requested.emit)
        button_layout.addWidget(self.reconnect_btn)

        self.simulate_checkbox = QCheckBox("Simulate")
        self.simulate_checkbox.toggled.connect(self.simulation_toggled.emit)
        button_layout.addWidget(self.simulate_checkbox)

        self.test_btn = QPushButton("Test")
        self.test_btn.clicked.connect(self.test_requested.emit)
        button_layout.addWidget(self.test_btn)

        layout.addLayout(button_layout)

        self.setLayout(layout)

    def update_state(self, state: InstrumentState, owner: Optional[str] = None):
        """
        Update instrument state.

        Args:
            state: New instrument state
            owner: Current owner (if claimed or in use)
        """
        self.state = state
        self.owner = owner

        # Update status indicator color
        if state == InstrumentState.FREE:
            color = STATUS_OK  # Green
            status_text = "Connected"
        elif state == InstrumentState.CLAIMED:
            color = STATUS_WARNING  # Yellow
            status_text = "Claimed"
        elif state == InstrumentState.IN_USE:
            color = STATUS_WARNING  # Yellow
            status_text = "In Use"
        elif state == InstrumentState.ERROR:
            color = STATUS_ERROR  # Red
            status_text = "Error"
        elif state == InstrumentState.DISCONNECTED:
            color = TEXT_DISABLED  # Gray
            status_text = "Disconnected"
        else:
            color = TEXT_DISABLED
            status_text = "Unknown"

        self.status_indicator.setStyleSheet(f"color: {color};")
        self.state_label.setText(status_text)
        self.state_value.setText(state.value.upper())

        # Update owner
        if owner:
            self.owner_value.setText(owner)
            self.owner_value.setStyleSheet(f"color: {STATUS_WARNING};")
        else:
            self.owner_value.setText("None")
            self.owner_value.setStyleSheet(f"color: {TEXT_SECONDARY};")

        # Update button states
        if state == InstrumentState.DISCONNECTED or state == InstrumentState.ERROR:
            self.reconnect_btn.setEnabled(True)
            self.test_btn.setEnabled(False)
        else:
            self.reconnect_btn.setEnabled(False)
            self.test_btn.setEnabled(True)

    def update_parameters(self, parameters: Dict[str, Any]):
        """
        Update instrument parameters display.

        Args:
            parameters: Dictionary of parameter name -> value
        """
        self.parameters = parameters

        # Add or update parameter rows
        for param_name, param_value in parameters.items():
            if param_name not in self.param_labels:
                # Create new row
                label = QLabel(f"{param_name}:")
                label.setStyleSheet("font-weight: bold;")

                value = QLabel(str(param_value))
                value.setProperty("class", "value")

                self.params_grid.addWidget(label, self.next_param_row, 0)
                self.params_grid.addWidget(value, self.next_param_row, 1)

                self.param_labels[param_name] = label
                self.param_values[param_name] = value

                self.next_param_row += 1
            else:
                # Update existing row
                self.param_values[param_name].setText(str(param_value))

    def set_simulation_mode(self, enabled: bool):
        """
        Set simulation mode checkbox state.

        Args:
            enabled: True to check, False to uncheck
        """
        self.simulate_checkbox.setChecked(enabled)


class SignalGeneratorCard(InstrumentStatusCard):
    """Status card specifically for Signal Generator."""

    def __init__(self, parent=None):
        super().__init__("Signal Generator (SRS SG384)", InstrumentType.SIGNAL_GENERATOR, parent)

    def update_sg_parameters(
        self,
        frequency: Optional[float] = None,
        power: Optional[float] = None,
        output_enabled: Optional[bool] = None,
        modulation: Optional[str] = None
    ):
        """
        Update SG-specific parameters.

        Args:
            frequency: Frequency in Hz
            power: Power in dBm
            output_enabled: Output state
            modulation: Modulation type
        """
        params = {}

        if frequency is not None:
            # Convert to GHz for display
            params['Frequency'] = f"{frequency / 1e9:.6f} GHz"

        if power is not None:
            params['Power'] = f"{power:.2f} dBm"

        if output_enabled is not None:
            params['Output'] = "ON" if output_enabled else "OFF"

        if modulation is not None:
            params['Modulation'] = modulation

        self.update_parameters(params)


class PulseBlasterCard(InstrumentStatusCard):
    """Status card specifically for PulseBlaster."""

    def __init__(self, parent=None):
        super().__init__("PulseBlaster (SpinCore)", InstrumentType.PULSEBLASTER, parent)

    def update_pb_parameters(
        self,
        sequence: Optional[str] = None,
        t_aom: Optional[float] = None,
        last_program_time: Optional[str] = None,
        status: Optional[str] = None
    ):
        """
        Update PB-specific parameters.

        Args:
            sequence: Current sequence name
            t_aom: AOM time in microseconds
            last_program_time: Last program timestamp
            status: PB status
        """
        params = {}

        if sequence is not None:
            params['Sequence'] = sequence

        if t_aom is not None:
            params['t_AOM'] = f"{t_aom:.2f} μs"

        if last_program_time is not None:
            params['Last Program'] = last_program_time

        if status is not None:
            params['Status'] = status

        self.update_parameters(params)


class DAQCard(InstrumentStatusCard):
    """Status card for DAQ (Analog Input/Output)."""

    def __init__(self, parent=None):
        super().__init__("DAQ (NI P6363)", InstrumentType.ANALOG_INPUT, parent)

    def update_daq_parameters(
        self,
        sampling_rate: Optional[float] = None,
        voltage_range: Optional[Tuple[float, float]] = None,
        trigger_source: Optional[str] = None,
        buffer_status: Optional[str] = None
    ):
        """
        Update DAQ-specific parameters.

        Args:
            sampling_rate: Sampling rate in S/s
            voltage_range: (min, max) voltage range
            trigger_source: Trigger source name
            buffer_status: Buffer status
        """
        params = {}

        if sampling_rate is not None:
            if sampling_rate >= 1e6:
                params['Sampling Rate'] = f"{sampling_rate / 1e6:.2f} MS/s"
            elif sampling_rate >= 1e3:
                params['Sampling Rate'] = f"{sampling_rate / 1e3:.2f} kS/s"
            else:
                params['Sampling Rate'] = f"{sampling_rate:.2f} S/s"

        if voltage_range is not None:
            params['Voltage Range'] = f"{voltage_range[0]:.1f} to {voltage_range[1]:.1f} V"

        if trigger_source is not None:
            params['Trigger'] = trigger_source

        if buffer_status is not None:
            params['Buffer'] = buffer_status

        self.update_parameters(params)


class CameraCard(InstrumentStatusCard):
    """Status card for Camera."""

    def __init__(self, parent=None):
        super().__init__("Camera (Hamamatsu)", InstrumentType.CAMERA, parent)

    def update_camera_parameters(
        self,
        exposure_time: Optional[float] = None,
        roi: Optional[Tuple[int, int, int, int]] = None,
        trigger_mode: Optional[str] = None,
        temperature: Optional[float] = None
    ):
        """
        Update camera-specific parameters.

        Args:
            exposure_time: Exposure time in ms
            roi: (x, y, width, height)
            trigger_mode: Trigger mode name
            temperature: Sensor temperature in °C
        """
        params = {}

        if exposure_time is not None:
            params['Exposure'] = f"{exposure_time:.2f} ms"

        if roi is not None:
            params['ROI'] = f"{roi[2]}×{roi[3]} @ ({roi[0]},{roi[1]})"

        if trigger_mode is not None:
            params['Trigger Mode'] = trigger_mode

        if temperature is not None:
            params['Temperature'] = f"{temperature:.1f} °C"

        self.update_parameters(params)


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication, QVBoxLayout, QWidget

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    # Create container
    container = QWidget()
    layout = QVBoxLayout()

    layout_horz1 = QHBoxLayout()
    layout_horz2 = QHBoxLayout()

    # Create cards
    sg_card = SignalGeneratorCard()
    sg_card.update_state(InstrumentState.FREE)
    sg_card.update_sg_parameters(
        frequency=2.87e9,
        power=8.0,
        output_enabled=True,
        modulation="Pulse"
    )
    layout_horz1.addWidget(sg_card)

    pb_card = PulseBlasterCard()
    pb_card.update_state(InstrumentState.CLAIMED, owner="ESR Window #1")
    pb_card.update_pb_parameters(
        sequence="esr_seq",
        t_aom=8000.0,
        last_program_time="14:32:45",
        status="Idle"
    )
    layout_horz1.addWidget(pb_card)

    daq_card = DAQCard()
    daq_card.update_state(InstrumentState.FREE)
    daq_card.update_daq_parameters(
        sampling_rate = 2e6,
        voltage_range = (-0.2, 0.2),
        trigger_source = None,
        buffer_status = None,
    )
    layout_horz2.addWidget(daq_card)

    layout.addLayout(layout_horz1)
    layout.addLayout(layout_horz2)
    
    layout.addStretch()
    container.setLayout(layout)
    container.show()

    sys.exit(app.exec())
