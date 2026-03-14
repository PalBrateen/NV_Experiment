#!/usr/bin/env python3
"""
Verdi V5 Laser Control GUI
===========================
PyQt5-based interface for Coherent Verdi V5 laser control via RS-232.

Features:
- Main tab with essential info (fast polling ~1 Hz)
- Secondary tabs with detailed info (lazy polling - only when visible)
- Mirrors the Verdi display menu structure

Author: Claude (Anthropic)
"""

import sys
from typing import Optional, Dict, Any, Callable
from dataclasses import dataclass
from enum import IntEnum

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QTabWidget, QLabel, QLineEdit, QPushButton,
    QGroupBox, QFrame, QStatusBar, QComboBox, QMessageBox,
    QSplitter, QSizePolicy
)
from PySide6.QtCore import QTimer, Qt, Signal, QThread
from PySide6.QtGui import QFont, QColor, QPalette, QDoubleValidator

# Try importing pyvisa - will use mock if not available
try:
    import pyvisa
    PYVISA_AVAILABLE = True
except ImportError:
    PYVISA_AVAILABLE = False
    print("⚠️  pyvisa not found - running in simulation mode")

# %gui qt

# =============================================================================
# Enums and Constants
# =============================================================================

class ServoState(IntEnum):
    """Servo status codes from Verdi documentation."""
    OPEN = 0
    LOCKED = 1
    SEEKING = 2
    FAULT = 3
    OPTIMIZING = 4
    CPEAKING = 5

    @classmethod
    def to_string(cls, value: int) -> str:
        names = {0: "open", 1: "lock", 2: "seek", 3: "FAULT", 4: "opt", 5: "cpeak"}
        return names.get(value, f"?{value}")


class LaserState(IntEnum):
    """Laser on/off status codes."""
    STANDBY = 0
    ON = 1
    FAULT = 2

    @classmethod
    def to_string(cls, value: int) -> str:
        names = {0: "STANDBY", 1: "ON", 2: "FAULT"}
        return names.get(value, f"?{value}")


class LightRegState(IntEnum):
    """Light regulation status codes."""
    OPEN = 0
    LOCKED = 1
    SEEKING = 2
    FAULT = 3

    @classmethod
    def to_string(cls, value: int) -> str:
        names = {0: "current", 1: "lock", 2: "seek", 3: "FAULT"}
        return names.get(value, f"?{value}")


# Fault code definitions from documentation
FAULT_CODES = {
    1: "Laser Head Interlock",
    2: "External Interlock",
    3: "PS Cover Interlock",
    4: "LBO Temperature",
    5: "LBO Not Locked at Set Temp",
    6: "Vanadate Temp.",
    7: "Etalon Temp.",
    8: "Diode 1 Temp.",
    9: "Diode 2 Temp.",
    10: "Baseplate Temp.",
    11: "Heatsink 1 Temp.",
    12: "Heatsink 2 Temp.",
    16: "Diode 1 Over Current",
    17: "Diode 2 Over Current",
    18: "Over Current",
    19: "Diode 1 Under Volt",
    20: "Diode 2 Under Volt",
    21: "Diode 1 Over Volt",
    22: "Diode 2 Over Volt",
    25: "Diode 1 EEPROM",
    26: "Diode 2 EEPROM",
    27: "Laser Head EEPROM",
    28: "Power Supply EEPROM",
    29: "PS-Head Mismatch",
    31: "Shutter State Mismatch",
    40: "Head-Diode Mismatch",
    47: "Vanadate2 Temp.",
}


# =============================================================================
# Communication Layer
# =============================================================================

class VerdiComm:
    """
    Communication interface for Verdi V5 laser via RS-232.
    
    Uses pyvisa for VISA resource management. Falls back to simulation
    mode if pyvisa is not available or no device is connected.
    """
    
    DEFAULT_TIMEOUT = 2000  # ms
    DEFAULT_BAUDRATE = 19200
    
    def __init__(self, resource_name: Optional[str] = None):
        self.resource_name = resource_name
        self.instrument = None
        self.connected = False
        self.simulation_mode = not PYVISA_AVAILABLE
        
        # Simulation state
        self._sim_state = {
            'power': 0.01,
            'set_power': 5.00,
            'current': 0.0,
            'd1_current': 0.0,
            'd2_current': 0.0,
            'baseplate_temp': 23.45,
            'laser_state': 0,
            'keyswitch': 0,
            'shutter': 0,
            'lbo_temp': 22.34,
            'lbo_set_temp': 147.90,
            'lbo_drive': 8191,
            'lbo_servo': 2,
            'lbo_heater': 1,
            'vanadate_temp': 23.44,
            'vanadate_set_temp': 30.00,
            'vanadate_drive': 0,
            'vanadate_servo': 0,
            'etalon_temp': 33.60,
            'etalon_set_temp': 33.60,
            'etalon_drive': 1050,
            'etalon_servo': 1,
            'diode1_temp': 23.82,
            'diode1_set_temp': 23.82,
            'diode1_drive': -199,
            'diode1_servo': 1,
            'diode1_heatsink': 26.11,
            'diode1_photocell': 0.0,
            'diode2_temp': 23.80,
            'diode2_set_temp': 23.80,
            'diode2_drive': -180,
            'diode2_servo': 1,
            'diode2_heatsink': 26.05,
            'diode2_photocell': 0.0,
            'light_reg_status': 2,
            'head_hours': 790.23,
            'diode1_hours': 1442.08,
            'diode2_hours': 1440.00,
            'ps_hours': 850.5,
            'software_version': '8.33',
            'faults': '',
            'fault_history': 'SYSTEM OK',
        }
    
    def connect(self, resource_name: str, baudrate: int = DEFAULT_BAUDRATE) -> bool:
        """Connect to the Verdi laser."""
        if self.simulation_mode:
            self.resource_name = resource_name
            self.connected = True
            return True
        
        try:
            rm = pyvisa.ResourceManager()
            self.instrument = rm.open_resource(resource_name)
            self.instrument.baud_rate = baudrate
            self.instrument.data_bits = 8
            self.instrument.stop_bits = pyvisa.constants.StopBits.one
            self.instrument.parity = pyvisa.constants.Parity.none
            self.instrument.timeout = self.DEFAULT_TIMEOUT
            self.instrument.read_termination = '\r\n'
            self.instrument.write_termination = '\r\n'
            
            self.resource_name = resource_name
            self.connected = True
            return True
        except Exception as e:
            print(f"❌ Connection failed: {e}")
            self.connected = False
            return False
    
    def disconnect(self):
        """Disconnect from the laser."""
        if self.instrument:
            try:
                self.instrument.close()
            except:
                pass
        self.instrument = None
        self.connected = False
    
    def query(self, command: str) -> str:
        """Send a query and return the response."""
        if self.simulation_mode:
            return self._simulate_query(command)
        
        if not self.connected:
            raise ConnectionError("Not connected to laser")
        
        try:
            response = self.instrument.query(command)
            return response.strip()
        except Exception as e:
            raise IOError(f"Query failed: {e}")
    
    def write(self, command: str) -> bool:
        """Send a command (no response expected)."""
        if self.simulation_mode:
            return self._simulate_write(command)
        
        if not self.connected:
            raise ConnectionError("Not connected to laser")
        
        try:
            self.instrument.write(command)
            # Read acknowledgment
            _ = self.instrument.read()
            return True
        except Exception as e:
            raise IOError(f"Write failed: {e}")
    
    def _simulate_query(self, command: str) -> str:
        """Simulate laser responses for testing."""
        cmd = command.upper().strip()
        
        # Map queries to simulated values
        query_map = {
            '?P': lambda: f"{self._sim_state['power']:.3f}",
            '?SP': lambda: f"{self._sim_state['set_power']:.4f}",
            '?C': lambda: f"{self._sim_state['current']:.1f}",
            '?D1C': lambda: f"{self._sim_state['d1_current']:.1f}",
            '?D2C': lambda: f"{self._sim_state['d2_current']:.1f}",
            '?BT': lambda: f"{self._sim_state['baseplate_temp']:.2f}",
            '?L': lambda: str(self._sim_state['laser_state']),
            '?K': lambda: str(self._sim_state['keyswitch']),
            '?S': lambda: str(self._sim_state['shutter']),
            '?LRS': lambda: str(self._sim_state['light_reg_status']),
            '?LBOT': lambda: f"{self._sim_state['lbo_temp']:.2f}",
            '?LBOST': lambda: f"{self._sim_state['lbo_set_temp']:.2f}",
            '?LBOD': lambda: str(self._sim_state['lbo_drive']),
            '?LBOSS': lambda: str(self._sim_state['lbo_servo']),
            '?LBOH': lambda: str(self._sim_state['lbo_heater']),
            '?VT': lambda: f"{self._sim_state['vanadate_temp']:.2f}",
            '?VST': lambda: f"{self._sim_state['vanadate_set_temp']:.2f}",
            '?VD': lambda: str(self._sim_state['vanadate_drive']),
            '?VSS': lambda: str(self._sim_state['vanadate_servo']),
            '?ET': lambda: f"{self._sim_state['etalon_temp']:.2f}",
            '?EST': lambda: f"{self._sim_state['etalon_set_temp']:.2f}",
            '?ED': lambda: str(self._sim_state['etalon_drive']),
            '?ESS': lambda: str(self._sim_state['etalon_servo']),
            '?D1T': lambda: f"{self._sim_state['diode1_temp']:.2f}",
            '?D1ST': lambda: f"{self._sim_state['diode1_set_temp']:.2f}",
            '?D1TD': lambda: str(self._sim_state['diode1_drive']),
            '?D1SS': lambda: str(self._sim_state['diode1_servo']),
            '?D1HST': lambda: f"{self._sim_state['diode1_heatsink']:.2f}",
            '?D1PC': lambda: f"{self._sim_state['diode1_photocell']:.2f}",
            '?D2T': lambda: f"{self._sim_state['diode2_temp']:.2f}",
            '?D2ST': lambda: f"{self._sim_state['diode2_set_temp']:.2f}",
            '?D2TD': lambda: str(self._sim_state['diode2_drive']),
            '?D2SS': lambda: str(self._sim_state['diode2_servo']),
            '?D2HST': lambda: f"{self._sim_state['diode2_heatsink']:.2f}",
            '?D2PC': lambda: f"{self._sim_state['diode2_photocell']:.2f}",
            '?HH': lambda: f"{self._sim_state['head_hours']:.2f}",
            '?D1H': lambda: f"{self._sim_state['diode1_hours']:.2f}",
            '?D2H': lambda: f"{self._sim_state['diode2_hours']:.2f}",
            '?PSH': lambda: f"{self._sim_state['ps_hours']:.1f}",
            '?SV': lambda: self._sim_state['software_version'],
            '?F': lambda: self._sim_state['faults'],
            '?FH': lambda: self._sim_state['fault_history'],
            '?M': lambda: '1',  # Light mode
            '?B': lambda: '19200',
        }
        
        if cmd in query_map:
            return query_map[cmd]()
        
        return "Command Error: " + command
    
    def _simulate_write(self, command: str) -> bool:
        """Simulate command execution."""
        cmd = command.upper().strip()
        
        if cmd.startswith('L=') or cmd.startswith('LASER='):
            val = int(cmd.split('=')[1].strip())
            self._sim_state['laser_state'] = val
        elif cmd.startswith('S=') or cmd.startswith('SHUTTER='):
            val = int(cmd.split('=')[1].strip())
            self._sim_state['shutter'] = val
        elif cmd.startswith('P=') or cmd.startswith('POWER='):
            val = float(cmd.split('=')[1].strip())
            self._sim_state['set_power'] = val
        elif cmd.startswith('LBOH='):
            val = int(cmd.split('=')[1].strip())
            self._sim_state['lbo_heater'] = val
        
        return True
    
    @staticmethod
    def list_resources() -> list:
        """List available VISA resources."""
        if not PYVISA_AVAILABLE:
            return ["SIM::Verdi_V5::INSTR"]
        try:
            rm = pyvisa.ResourceManager()
            return list(rm.list_resources())
        except:
            return []


# =============================================================================
# Custom Widgets
# =============================================================================

class StatusIndicator(QLabel):
    """A colored status indicator widget."""
    
    COLORS = {
        'ok': '#27ae60',       # Green
        'warning': '#e67e22',  # Orange
        'error': '#e74c3c',    # Red
        'inactive': '#555555', # Dark gray
        'seeking': '#3498db',  # Blue
    }
    
    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setAlignment(Qt.AlignCenter)
        self.setMinimumWidth(60)
        self.set_status('inactive')
    
    def set_status(self, status: str, text: Optional[str] = None):
        """Set the indicator status and optionally update text."""
        color = self.COLORS.get(status, self.COLORS['inactive'])
        self.setStyleSheet(f"""
            QLabel {{
                background-color: {color};
                color: white;
                font-weight: bold;
                padding: 4px 8px;
                border-radius: 4px;
                font-family: 'Consolas', 'Monaco', monospace;
            }}
        """)
        if text is not None:
            self.setText(text)


class ValueDisplay(QFrame):
    """A labeled value display widget."""
    
    def __init__(self, label: str, unit: str = "", parent=None):
        super().__init__(parent)
        self.unit = unit
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(2)
        
        self.label = QLabel(label)
        self.label.setStyleSheet("color: #888888; font-size: 11px;")
        
        self.value = QLabel("---")
        self.value.setStyleSheet("""
            font-size: 18px;
            font-weight: bold;
            font-family: 'Consolas', 'Monaco', monospace;
            color: #00dd88;
        """)
        
        layout.addWidget(self.label)
        layout.addWidget(self.value)
        
        self.setStyleSheet("""
            QFrame {
                background-color: #2a2a32;
                border-radius: 6px;
                border: 1px solid #444444;
            }
        """)
    
    def set_value(self, value: Any, fmt: str = "{}"):
        """Update the displayed value."""
        try:
            text = fmt.format(value)
            if self.unit:
                text += f" {self.unit}"
            self.value.setText(text)
        except:
            self.value.setText("---")


class ServoStatusRow(QWidget):
    """A row in the servo status table."""
    
    def __init__(self, name: str, parent=None):
        super().__init__(parent)
        
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 2, 0, 2)
        layout.setSpacing(8)
        
        self.name_label = QLabel(name)
        self.name_label.setMinimumWidth(80)
        self.name_label.setStyleSheet("font-weight: bold;")
        
        self.status = StatusIndicator()
        self.temp = QLabel("---")
        self.temp.setMinimumWidth(70)
        self.temp.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.temp.setStyleSheet("font-family: monospace; color: #dddddd;")
        
        self.drive = QLabel("---")
        self.drive.setMinimumWidth(70)
        self.drive.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.drive.setStyleSheet("font-family: monospace; color: #dddddd;")
        
        layout.addWidget(self.name_label)
        layout.addWidget(self.status)
        layout.addWidget(self.temp)
        layout.addWidget(self.drive)
    
    def update_values(self, state: int, temp: float, drive: float):
        """Update all values in the row."""
        state_str = ServoState.to_string(state)
        
        if state == ServoState.LOCKED:
            self.status.set_status('ok', state_str)
        elif state == ServoState.SEEKING:
            self.status.set_status('seeking', state_str)
        elif state == ServoState.FAULT:
            self.status.set_status('error', state_str)
        elif state == ServoState.OPEN:
            self.status.set_status('inactive', state_str)
        else:
            self.status.set_status('warning', state_str)
        
        self.temp.setText(f"{temp:.2f} °C")
        self.drive.setText(f"{drive:.0f}")


# =============================================================================
# Tab Widgets
# =============================================================================

class MainStatusTab(QWidget):
    """
    Main status tab - displays essential laser info.
    This tab is polled at high frequency (~1 Hz).
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        
        # === Top Status Bar ===
        status_frame = QFrame()
        status_frame.setStyleSheet("""
            QFrame {
                background-color: #1a1a22;
                border-radius: 8px;
                padding: 8px;
                border: 1px solid #444444;
            }
        """)
        status_layout = QHBoxLayout(status_frame)
        
        self.laser_status = StatusIndicator("STANDBY")
        self.laser_status.setMinimumWidth(100)
        
        self.shutter_status = StatusIndicator("CLOSED")
        self.light_reg_status = StatusIndicator("---")
        
        status_layout.addWidget(QLabel("Laser:"))
        status_layout.addWidget(self.laser_status)
        status_layout.addSpacing(20)
        status_layout.addWidget(QLabel("Shutter:"))
        status_layout.addWidget(self.shutter_status)
        status_layout.addSpacing(20)
        status_layout.addWidget(QLabel("Light Reg:"))
        status_layout.addWidget(self.light_reg_status)
        status_layout.addStretch()
        
        for label in status_frame.findChildren(QLabel):
            if not isinstance(label, StatusIndicator):
                label.setStyleSheet("color: white; font-weight: bold;")
        
        layout.addWidget(status_frame)
        
        # === Power Display ===
        power_group = QGroupBox("Power Control")
        power_layout = QGridLayout(power_group)
        
        self.power_display = ValueDisplay("Output Power", "W")
        self.set_power_display = ValueDisplay("Set Power", "W")
        
        self.power_entry = QLineEdit()
        self.power_entry.setPlaceholderText("Enter power (W)")
        self.power_entry.setValidator(QDoubleValidator(0.0, 10.0, 4))
        self.power_entry.setMaximumWidth(120)
        
        self.set_power_btn = QPushButton("Set Power")
        self.set_power_btn.clicked.connect(self._on_set_power)
        
        power_layout.addWidget(self.power_display, 0, 0)
        power_layout.addWidget(self.set_power_display, 0, 1)
        power_layout.addWidget(self.power_entry, 1, 0)
        power_layout.addWidget(self.set_power_btn, 1, 1)
        
        layout.addWidget(power_group)
        
        # === Current & Temperature ===
        info_group = QGroupBox("Diode & Temperature")
        info_layout = QHBoxLayout(info_group)
        
        self.current_display = ValueDisplay("Avg Current", "A")
        self.d1_current_display = ValueDisplay("Diode 1 Current", "A")
        self.d2_current_display = ValueDisplay("Diode 2 Current", "A")
        self.baseplate_display = ValueDisplay("Baseplate Temp", "°C")
        
        info_layout.addWidget(self.current_display)
        info_layout.addWidget(self.d1_current_display)
        info_layout.addWidget(self.d2_current_display)
        info_layout.addWidget(self.baseplate_display)
        
        layout.addWidget(info_group)
        
        # === Control Buttons ===
        ctrl_group = QGroupBox("Control")
        ctrl_layout = QHBoxLayout(ctrl_group)
        
        self.laser_on_btn = QPushButton("LASER ON")
        self.laser_on_btn.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-weight: bold;
                padding: 12px 24px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #2ecc71; }
            QPushButton:pressed { background-color: #1e8449; }
        """)
        self.laser_on_btn.clicked.connect(lambda: self._send_command("L=1"))
        
        self.laser_off_btn = QPushButton("STANDBY")
        self.laser_off_btn.setStyleSheet("""
            QPushButton {
                background-color: #e74c3c;
                color: white;
                font-weight: bold;
                padding: 12px 24px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #ec7063; }
            QPushButton:pressed { background-color: #c0392b; }
        """)
        self.laser_off_btn.clicked.connect(lambda: self._send_command("L=0"))
        
        self.shutter_btn = QPushButton("Toggle Shutter")
        self.shutter_btn.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                font-weight: bold;
                padding: 12px 24px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #5dade2; }
        """)
        self.shutter_btn.clicked.connect(self._toggle_shutter)
        
        self.flash_btn = QPushButton("Flash (Etalon)")
        self.flash_btn.setStyleSheet("""
            QPushButton {
                background-color: #9b59b6;
                color: white;
                padding: 12px 24px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #a569bd; }
        """)
        self.flash_btn.clicked.connect(lambda: self._send_command("FL=1"))
        
        ctrl_layout.addWidget(self.laser_on_btn)
        ctrl_layout.addWidget(self.laser_off_btn)
        ctrl_layout.addWidget(self.shutter_btn)
        ctrl_layout.addWidget(self.flash_btn)
        
        layout.addWidget(ctrl_group)
        layout.addStretch()
    
    def _on_set_power(self):
        """Handle set power button click."""
        try:
            power = float(self.power_entry.text())
            if 0 <= power <= 10:
                self._send_command(f"P={power:.4f}")
                self.power_entry.clear()
            else:
                QMessageBox.warning(self, "Invalid Power", "Power must be between 0 and 10 W")
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number")
    
    def _toggle_shutter(self):
        """Toggle the external shutter."""
        try:
            current = int(self.comm.query("?S"))
            new_state = 0 if current == 1 else 1
            self._send_command(f"S={new_state}")
        except Exception as e:
            print(f"❌ Shutter toggle failed: {e}")
    
    def _send_command(self, cmd: str):
        """Send a command to the laser."""
        try:
            self.comm.write(cmd)
        except Exception as e:
            QMessageBox.critical(self, "Command Error", f"Failed to send command: {e}")
    
    def update_data(self):
        """Poll and update all displayed values."""
        try:
            # Power
            power = float(self.comm.query("?P"))
            self.power_display.set_value(power, "{:.3f}")
            
            set_power = float(self.comm.query("?SP"))
            self.set_power_display.set_value(set_power, "{:.4f}")
            
            # Current
            current = float(self.comm.query("?C"))
            self.current_display.set_value(current, "{:.1f}")
            
            d1_current = float(self.comm.query("?D1C"))
            self.d1_current_display.set_value(d1_current, "{:.1f}")
            
            d2_current = float(self.comm.query("?D2C"))
            self.d2_current_display.set_value(d2_current, "{:.1f}")
            
            # Temperature
            baseplate = float(self.comm.query("?BT"))
            self.baseplate_display.set_value(baseplate, "{:.2f}")
            
            # Status indicators
            laser_state = int(self.comm.query("?L"))
            state_str = LaserState.to_string(laser_state)
            if laser_state == LaserState.ON:
                self.laser_status.set_status('ok', state_str)
            elif laser_state == LaserState.FAULT:
                self.laser_status.set_status('error', state_str)
            else:
                self.laser_status.set_status('inactive', state_str)
            
            shutter = int(self.comm.query("?S"))
            if shutter == 1:
                self.shutter_status.set_status('warning', "OPEN")
            else:
                self.shutter_status.set_status('inactive', "CLOSED")
            
            light_reg = int(self.comm.query("?LRS"))
            lrs_str = LightRegState.to_string(light_reg)
            if light_reg == LightRegState.LOCKED:
                self.light_reg_status.set_status('ok', lrs_str)
            elif light_reg == LightRegState.SEEKING:
                self.light_reg_status.set_status('seeking', lrs_str)
            elif light_reg == LightRegState.FAULT:
                self.light_reg_status.set_status('error', lrs_str)
            else:
                self.light_reg_status.set_status('inactive', lrs_str)
                
        except Exception as e:
            print(f"❌ Main tab update error: {e}")


class ServoStatusTab(QWidget):
    """
    Servo Status Screen tab - mirrors the Verdi servo status display.
    Shows state, temperature, and drive for all servos.
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Header
        header = QLabel("Servo Status Screen")
        header.setStyleSheet("""
            font-size: 18px;
            font-weight: bold;
            color: #dddddd;
            padding: 8px;
        """)
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # Column headers
        header_widget = QWidget()
        header_layout = QHBoxLayout(header_widget)
        header_layout.setContentsMargins(0, 0, 0, 0)
        
        headers = [("Servo", 80), ("State", 60), ("Temp (°C)", 70), ("Drive", 70)]
        for text, width in headers:
            lbl = QLabel(text)
            lbl.setMinimumWidth(width)
            lbl.setStyleSheet("font-weight: bold; color: #888888;")
            if text != "Servo":
                lbl.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            header_layout.addWidget(lbl)
        
        layout.addWidget(header_widget)
        
        # Servo rows
        self.servo_rows = {}
        
        servo_frame = QFrame()
        servo_frame.setStyleSheet("""
            QFrame {
                background-color: #2a2a32;
                border-radius: 8px;
                padding: 8px;
                border: 1px solid #444444;
            }
        """)
        servo_layout = QVBoxLayout(servo_frame)
        
        for name in ["Laser", "LBO", "Vanadate", "Etalon", "Diode 1", "Diode 2"]:
            row = ServoStatusRow(name)
            self.servo_rows[name] = row
            servo_layout.addWidget(row)
        
        layout.addWidget(servo_frame)
        layout.addStretch()
    
    def update_data(self):
        """Poll and update servo status."""
        try:
            # Laser (light regulation)
            lrs = int(self.comm.query("?LRS"))
            power = float(self.comm.query("?P"))
            self.servo_rows["Laser"].update_values(lrs, power, 0)
            
            # LBO
            lbo_ss = int(self.comm.query("?LBOSS"))
            lbo_t = float(self.comm.query("?LBOT"))
            lbo_d = float(self.comm.query("?LBOD"))
            self.servo_rows["LBO"].update_values(lbo_ss, lbo_t, lbo_d)
            
            # Vanadate
            v_ss = int(self.comm.query("?VSS"))
            v_t = float(self.comm.query("?VT"))
            v_d = float(self.comm.query("?VD"))
            self.servo_rows["Vanadate"].update_values(v_ss, v_t, v_d)
            
            # Etalon
            e_ss = int(self.comm.query("?ESS"))
            e_t = float(self.comm.query("?ET"))
            e_d = float(self.comm.query("?ED"))
            self.servo_rows["Etalon"].update_values(e_ss, e_t, e_d)
            
            # Diode 1
            d1_ss = int(self.comm.query("?D1SS"))
            d1_t = float(self.comm.query("?D1T"))
            d1_d = float(self.comm.query("?D1TD"))
            self.servo_rows["Diode 1"].update_values(d1_ss, d1_t, d1_d)
            
            # Diode 2
            d2_ss = int(self.comm.query("?D2SS"))
            d2_t = float(self.comm.query("?D2T"))
            d2_d = float(self.comm.query("?D2TD"))
            self.servo_rows["Diode 2"].update_values(d2_ss, d2_t, d2_d)
            
        except Exception as e:
            print(f"❌ Servo status update error: {e}")


class TemperatureTab(QWidget):
    """
    Temperature Set Points tab - shows detailed temp info for each servo.
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        header = QLabel("Temperature Set Points")
        header.setStyleSheet("font-size: 18px; font-weight: bold; padding: 8px; color: #dddddd;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # Create sections for each temperature-controlled component
        self.temp_displays = {}
        
        components = [
            ("LBO", "?LBOT", "?LBOST", "?LBOD", "?LBOSS"),
            ("Vanadate", "?VT", "?VST", "?VD", "?VSS"),
            ("Etalon", "?ET", "?EST", "?ED", "?ESS"),
            ("Diode 1", "?D1T", "?D1ST", "?D1TD", "?D1SS"),
            ("Diode 2", "?D2T", "?D2ST", "?D2TD", "?D2SS"),
        ]
        
        for name, temp_q, set_q, drive_q, status_q in components:
            group = QGroupBox(f"{name} Temperature")
            grid = QGridLayout(group)
            
            displays = {
                'temp': ValueDisplay("Read T", "°C"),
                'set': ValueDisplay("Set Pt", "°C"),
                'drive': ValueDisplay("Drive", ""),
                'status': StatusIndicator("---"),
            }
            
            grid.addWidget(displays['temp'], 0, 0)
            grid.addWidget(displays['set'], 0, 1)
            grid.addWidget(displays['drive'], 0, 2)
            grid.addWidget(QLabel("Status:"), 0, 3)
            grid.addWidget(displays['status'], 0, 4)
            
            self.temp_displays[name] = {
                'widgets': displays,
                'queries': (temp_q, set_q, drive_q, status_q)
            }
            
            layout.addWidget(group)
        
        layout.addStretch()
    
    def update_data(self):
        """Poll and update temperature data."""
        try:
            for name, data in self.temp_displays.items():
                widgets = data['widgets']
                temp_q, set_q, drive_q, status_q = data['queries']
                
                temp = float(self.comm.query(temp_q))
                set_pt = float(self.comm.query(set_q))
                drive = float(self.comm.query(drive_q))
                status = int(self.comm.query(status_q))
                
                widgets['temp'].set_value(temp, "{:.2f}")
                widgets['set'].set_value(set_pt, "{:.2f}")
                widgets['drive'].set_value(drive, "{:.1f}")
                
                status_str = ServoState.to_string(status)
                if status == ServoState.LOCKED:
                    widgets['status'].set_status('ok', status_str)
                elif status == ServoState.SEEKING:
                    widgets['status'].set_status('seeking', status_str)
                elif status == ServoState.FAULT:
                    widgets['status'].set_status('error', status_str)
                else:
                    widgets['status'].set_status('inactive', status_str)
                    
        except Exception as e:
            print(f"❌ Temperature tab update error: {e}")


class DiodeParametersTab(QWidget):
    """
    Diode Parameters Screen - voltage, current, photocell readings.
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        header = QLabel("Diode Parameters Screen")
        header.setStyleSheet("font-size: 18px; font-weight: bold; padding: 8px; color: #dddddd;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # Diode 1
        d1_group = QGroupBox("Diode 1")
        d1_layout = QHBoxLayout(d1_group)
        
        self.d1_current = ValueDisplay("Current", "A")
        self.d1_photocell = ValueDisplay("Photocell", "V")
        self.d1_heatsink = ValueDisplay("Heatsink T", "°C")
        self.d1_temp = ValueDisplay("Diode T", "°C")
        
        d1_layout.addWidget(self.d1_current)
        d1_layout.addWidget(self.d1_photocell)
        d1_layout.addWidget(self.d1_heatsink)
        d1_layout.addWidget(self.d1_temp)
        
        layout.addWidget(d1_group)
        
        # Diode 2
        d2_group = QGroupBox("Diode 2")
        d2_layout = QHBoxLayout(d2_group)
        
        self.d2_current = ValueDisplay("Current", "A")
        self.d2_photocell = ValueDisplay("Photocell", "V")
        self.d2_heatsink = ValueDisplay("Heatsink T", "°C")
        self.d2_temp = ValueDisplay("Diode T", "°C")
        
        d2_layout.addWidget(self.d2_current)
        d2_layout.addWidget(self.d2_photocell)
        d2_layout.addWidget(self.d2_heatsink)
        d2_layout.addWidget(self.d2_temp)
        
        layout.addWidget(d2_group)
        layout.addStretch()
    
    def update_data(self):
        """Poll and update diode parameters."""
        try:
            # Diode 1
            self.d1_current.set_value(float(self.comm.query("?D1C")), "{:.2f}")
            self.d1_photocell.set_value(float(self.comm.query("?D1PC")), "{:.2f}")
            self.d1_heatsink.set_value(float(self.comm.query("?D1HST")), "{:.2f}")
            self.d1_temp.set_value(float(self.comm.query("?D1T")), "{:.2f}")
            
            # Diode 2
            self.d2_current.set_value(float(self.comm.query("?D2C")), "{:.2f}")
            self.d2_photocell.set_value(float(self.comm.query("?D2PC")), "{:.2f}")
            self.d2_heatsink.set_value(float(self.comm.query("?D2HST")), "{:.2f}")
            self.d2_temp.set_value(float(self.comm.query("?D2T")), "{:.2f}")
            
        except Exception as e:
            print(f"❌ Diode parameters update error: {e}")


class LaserStatusTab(QWidget):
    """
    Laser Status Screen - hours, software version, etc.
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        header = QLabel("Laser Status Screen")
        header.setStyleSheet("font-size: 18px; font-weight: bold; padding: 8px; color: #dddddd;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        info_group = QGroupBox("System Information")
        info_layout = QGridLayout(info_group)
        
        self.sw_version = ValueDisplay("S/W Version", "")
        self.head_hours = ValueDisplay("HEAD Hours", "hrs")
        self.ps_hours = ValueDisplay("PS Hours", "hrs")
        self.d1_hours = ValueDisplay("Diode 1 Hours", "hrs")
        self.d2_hours = ValueDisplay("Diode 2 Hours", "hrs")
        self.baseplate = ValueDisplay("Baseplate T", "°C")
        self.heatsink1 = ValueDisplay("Heatsink 1 T", "°C")
        
        info_layout.addWidget(self.sw_version, 0, 0)
        info_layout.addWidget(self.head_hours, 0, 1)
        info_layout.addWidget(self.ps_hours, 0, 2)
        info_layout.addWidget(self.d1_hours, 1, 0)
        info_layout.addWidget(self.d2_hours, 1, 1)
        info_layout.addWidget(self.baseplate, 1, 2)
        info_layout.addWidget(self.heatsink1, 2, 0)
        
        layout.addWidget(info_group)
        layout.addStretch()
    
    def update_data(self):
        """Poll and update laser status."""
        try:
            self.sw_version.set_value(self.comm.query("?SV"))
            self.head_hours.set_value(float(self.comm.query("?HH")), "{:.2f}")
            self.ps_hours.set_value(float(self.comm.query("?PSH")), "{:.1f}")
            self.d1_hours.set_value(float(self.comm.query("?D1H")), "{:.2f}")
            self.d2_hours.set_value(float(self.comm.query("?D2H")), "{:.2f}")
            self.baseplate.set_value(float(self.comm.query("?BT")), "{:.2f}")
            self.heatsink1.set_value(float(self.comm.query("?D1HST")), "{:.2f}")
            
        except Exception as e:
            print(f"❌ Laser status update error: {e}")


class LBOSettingsTab(QWidget):
    """
    LBO Settings tab - temperature control and optimization.
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        header = QLabel("LBO Settings")
        header.setStyleSheet("font-size: 18px; font-weight: bold; padding: 8px; color: #dddddd;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # LBO Status
        status_group = QGroupBox("LBO Status")
        status_layout = QHBoxLayout(status_group)
        
        self.lbo_temp = ValueDisplay("Temperature", "°C")
        self.lbo_set = ValueDisplay("Set Point", "°C")
        self.lbo_drive = ValueDisplay("Drive", "")
        self.lbo_status = StatusIndicator("---")
        
        status_layout.addWidget(self.lbo_temp)
        status_layout.addWidget(self.lbo_set)
        status_layout.addWidget(self.lbo_drive)
        status_layout.addWidget(QLabel("Status:"))
        status_layout.addWidget(self.lbo_status)
        
        layout.addWidget(status_group)
        
        # Heater Control
        heater_group = QGroupBox("Heater Control")
        heater_layout = QHBoxLayout(heater_group)
        
        self.heater_status = StatusIndicator("---")
        self.heater_status.setMinimumWidth(120)
        
        self.heat_btn = QPushButton("Start HEATING")
        self.heat_btn.setStyleSheet("""
            QPushButton {
                background-color: #e74c3c;
                color: white;
                font-weight: bold;
                padding: 10px 20px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #ec7063; }
        """)
        self.heat_btn.clicked.connect(lambda: self._send_command("LBOH=1"))
        
        self.cool_btn = QPushButton("Start COOLING")
        self.cool_btn.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                font-weight: bold;
                padding: 10px 20px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #5dade2; }
        """)
        self.cool_btn.clicked.connect(lambda: self._send_command("LBOH=0"))
        
        self.optimize_btn = QPushButton("LBO Optimize")
        self.optimize_btn.setStyleSheet("""
            QPushButton {
                background-color: #9b59b6;
                color: white;
                font-weight: bold;
                padding: 10px 20px;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #a569bd; }
        """)
        self.optimize_btn.clicked.connect(lambda: self._send_command("LBOOPT=1"))
        
        heater_layout.addWidget(QLabel("Heater:"))
        heater_layout.addWidget(self.heater_status)
        heater_layout.addStretch()
        heater_layout.addWidget(self.heat_btn)
        heater_layout.addWidget(self.cool_btn)
        heater_layout.addWidget(self.optimize_btn)
        
        layout.addWidget(heater_group)
        layout.addStretch()
    
    def _send_command(self, cmd: str):
        """Send a command to the laser."""
        try:
            self.comm.write(cmd)
        except Exception as e:
            QMessageBox.critical(self, "Command Error", f"Failed: {e}")
    
    def update_data(self):
        """Poll and update LBO settings."""
        try:
            self.lbo_temp.set_value(float(self.comm.query("?LBOT")), "{:.2f}")
            self.lbo_set.set_value(float(self.comm.query("?LBOST")), "{:.2f}")
            self.lbo_drive.set_value(float(self.comm.query("?LBOD")), "{:.0f}")
            
            status = int(self.comm.query("?LBOSS"))
            status_str = ServoState.to_string(status)
            if status == ServoState.LOCKED:
                self.lbo_status.set_status('ok', status_str)
            elif status == ServoState.SEEKING:
                self.lbo_status.set_status('seeking', status_str)
            elif status == ServoState.FAULT:
                self.lbo_status.set_status('error', status_str)
            else:
                self.lbo_status.set_status('inactive', status_str)
            
            heater = int(self.comm.query("?LBOH"))
            if heater == 1:
                self.heater_status.set_status('warning', "HEATING")
            else:
                self.heater_status.set_status('inactive', "COOLING")
                
        except Exception as e:
            print(f"❌ LBO settings update error: {e}")


class FaultsTab(QWidget):
    """
    Faults tab - displays current faults and fault history.
    """
    
    def __init__(self, comm: VerdiComm, parent=None):
        super().__init__(parent)
        self.comm = comm
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        header = QLabel("Faults")
        header.setStyleSheet("font-size: 18px; font-weight: bold; padding: 8px; color: #dddddd;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # Current Faults
        current_group = QGroupBox("Current Faults")
        current_layout = QVBoxLayout(current_group)
        
        self.current_faults_label = QLabel("SYSTEM OK")
        self.current_faults_label.setStyleSheet("""
            font-size: 14px;
            padding: 16px;
            background-color: #1a3d2e;
            border-radius: 6px;
            color: #27ae60;
            border: 1px solid #27ae60;
        """)
        self.current_faults_label.setWordWrap(True)
        current_layout.addWidget(self.current_faults_label)
        
        layout.addWidget(current_group)
        
        # Fault History
        history_group = QGroupBox("Fault History (since last LASER ON)")
        history_layout = QVBoxLayout(history_group)
        
        self.fault_history_label = QLabel("SYSTEM OK")
        self.fault_history_label.setStyleSheet("""
            font-size: 14px;
            padding: 16px;
            background-color: #2a2a32;
            border-radius: 6px;
            color: #aaaaaa;
            border: 1px solid #444444;
        """)
        self.fault_history_label.setWordWrap(True)
        history_layout.addWidget(self.fault_history_label)
        
        layout.addWidget(history_group)
        layout.addStretch()
    
    def update_data(self):
        """Poll and update fault information."""
        try:
            # Current faults
            faults = self.comm.query("?F")
            if faults and faults != "SYSTEM OK":
                fault_list = faults.split("&")
                fault_texts = []
                for f in fault_list:
                    try:
                        code = int(f)
                        fault_texts.append(f"⚠️ {FAULT_CODES.get(code, f'Unknown ({code})')}")
                    except:
                        fault_texts.append(f"⚠️ {f}")
                self.current_faults_label.setText("\n".join(fault_texts))
                self.current_faults_label.setStyleSheet("""
                    font-size: 14px;
                    padding: 16px;
                    background-color: #3d1a1a;
                    border-radius: 6px;
                    color: #e74c3c;
                    border: 1px solid #e74c3c;
                """)
            else:
                self.current_faults_label.setText("✅ SYSTEM OK")
                self.current_faults_label.setStyleSheet("""
                    font-size: 14px;
                    padding: 16px;
                    background-color: #1a3d2e;
                    border-radius: 6px;
                    color: #27ae60;
                    border: 1px solid #27ae60;
                """)
            
            # Fault history
            history = self.comm.query("?FH")
            if history and history != "SYSTEM OK":
                fault_list = history.split("&")
                fault_texts = []
                for f in fault_list:
                    try:
                        code = int(f)
                        fault_texts.append(f"• {FAULT_CODES.get(code, f'Unknown ({code})')}")
                    except:
                        fault_texts.append(f"• {f}")
                self.fault_history_label.setText("\n".join(fault_texts))
            else:
                self.fault_history_label.setText("✅ SYSTEM OK")
                
        except Exception as e:
            print(f"❌ Faults update error: {e}")


# =============================================================================
# Main Window
# =============================================================================

class VerdiV5MainWindow(QMainWindow):
    """
    Main application window for Verdi V5 control.
    """
    
    # Polling intervals (ms)
    FAST_POLL_INTERVAL = 1000    # 1 Hz for main tab
    SLOW_POLL_INTERVAL = 5000    # 0.2 Hz for secondary tabs when visible
    
    def __init__(self):
        super().__init__()
        
        self.comm = VerdiComm()
        self.current_tab_index = 0
        
        self._setup_ui()
        self._setup_timers()
        self._setup_connections()
        
        # Auto-connect in simulation mode
        if self.comm.simulation_mode:
            self._on_connect()
    
    def _setup_ui(self):
        self.setWindowTitle("Verdi V5 Laser Control")
        self.setMinimumSize(800, 600)
        
        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        
        # Connection bar
        conn_frame = QFrame()
        conn_frame.setStyleSheet("""
            QFrame {
                background-color: #1a1a22;
                padding: 8px;
                border-bottom: 1px solid #444444;
            }
        """)
        conn_layout = QHBoxLayout(conn_frame)
        
        port_label = QLabel("Port:")
        port_label.setStyleSheet("color: #dddddd; font-weight: bold;")
        
        self.port_combo = QComboBox()
        self.port_combo.setMinimumWidth(200)
        self.port_combo.setEditable(True)
        self._refresh_ports()
        
        self.refresh_btn = QPushButton("↻")
        self.refresh_btn.setMaximumWidth(40)
        self.refresh_btn.clicked.connect(self._refresh_ports)
        
        self.connect_btn = QPushButton("Connect")
        self.connect_btn.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-weight: bold;
                padding: 6px 16px;
                border-radius: 4px;
            }
            QPushButton:hover { background-color: #2ecc71; }
        """)
        self.connect_btn.clicked.connect(self._on_connect)
        
        self.disconnect_btn = QPushButton("Disconnect")
        self.disconnect_btn.setStyleSheet("""
            QPushButton {
                background-color: #e74c3c;
                color: white;
                font-weight: bold;
                padding: 6px 16px;
                border-radius: 4px;
            }
            QPushButton:hover { background-color: #ec7063; }
        """)
        self.disconnect_btn.clicked.connect(self._on_disconnect)
        self.disconnect_btn.setEnabled(False)
        
        self.conn_status = StatusIndicator("Disconnected")
        
        conn_layout.addWidget(port_label)
        conn_layout.addWidget(self.port_combo)
        conn_layout.addWidget(self.refresh_btn)
        conn_layout.addWidget(self.connect_btn)
        conn_layout.addWidget(self.disconnect_btn)
        conn_layout.addStretch()
        conn_layout.addWidget(self.conn_status)
        
        layout.addWidget(conn_frame)
        
        # Tab widget
        self.tabs = QTabWidget()
        
        # Create tabs
        self.main_tab = MainStatusTab(self.comm)
        self.servo_tab = ServoStatusTab(self.comm)
        self.temp_tab = TemperatureTab(self.comm)
        self.diode_tab = DiodeParametersTab(self.comm)
        self.laser_status_tab = LaserStatusTab(self.comm)
        self.lbo_tab = LBOSettingsTab(self.comm)
        self.faults_tab = FaultsTab(self.comm)
        
        self.tabs.addTab(self.main_tab, "📊 Main Status")
        self.tabs.addTab(self.servo_tab, "🔄 Servo Status")
        self.tabs.addTab(self.temp_tab, "🌡️ Temperatures")
        self.tabs.addTab(self.diode_tab, "💡 Diode Parameters")
        self.tabs.addTab(self.laser_status_tab, "ℹ️ Laser Info")
        self.tabs.addTab(self.lbo_tab, "⚙️ LBO Settings")
        self.tabs.addTab(self.faults_tab, "⚠️ Faults")
        
        layout.addWidget(self.tabs)
        
        # Status bar
        self.statusBar().showMessage("Ready - Connect to laser to begin")
    
    def _setup_timers(self):
        """Setup polling timers."""
        # Fast timer for main tab
        self.fast_timer = QTimer(self)
        self.fast_timer.timeout.connect(self._fast_poll)
        
        # Slow timer for secondary tabs
        self.slow_timer = QTimer(self)
        self.slow_timer.timeout.connect(self._slow_poll)
    
    def _setup_connections(self):
        """Setup signal connections."""
        self.tabs.currentChanged.connect(self._on_tab_changed)
    
    def _refresh_ports(self):
        """Refresh the list of available ports."""
        self.port_combo.clear()
        ports = VerdiComm.list_resources()
        self.port_combo.addItems(ports)
        
        if not ports:
            self.port_combo.addItem("No ports found")
    
    def _on_connect(self):
        """Handle connect button click."""
        port = self.port_combo.currentText()
        
        if self.comm.connect(port):
            self.conn_status.set_status('ok', "Connected")
            self.connect_btn.setEnabled(False)
            self.disconnect_btn.setEnabled(True)
            self.port_combo.setEnabled(False)
            
            # Start polling
            self.fast_timer.start(self.FAST_POLL_INTERVAL)
            self.slow_timer.start(self.SLOW_POLL_INTERVAL)
            
            # Initial update
            self._fast_poll()
            self._slow_poll()
            
            mode = "(simulation)" if self.comm.simulation_mode else ""
            self.statusBar().showMessage(f"Connected to {port} {mode}")
        else:
            self.conn_status.set_status('error', "Failed")
            QMessageBox.critical(self, "Connection Error", 
                                f"Failed to connect to {port}")
    
    def _on_disconnect(self):
        """Handle disconnect button click."""
        self.fast_timer.stop()
        self.slow_timer.stop()
        
        self.comm.disconnect()
        
        self.conn_status.set_status('inactive', "Disconnected")
        self.connect_btn.setEnabled(True)
        self.disconnect_btn.setEnabled(False)
        self.port_combo.setEnabled(True)
        
        self.statusBar().showMessage("Disconnected")
    
    def _on_tab_changed(self, index: int):
        """Handle tab change - update the newly visible tab immediately."""
        self.current_tab_index = index
        
        if self.comm.connected:
            # Immediately update the newly visible tab
            self._slow_poll()
    
    def _fast_poll(self):
        """Fast polling for main tab (always runs)."""
        if not self.comm.connected:
            return
        
        try:
            self.main_tab.update_data()
        except Exception as e:
            print(f"❌ Fast poll error: {e}")
    
    def _slow_poll(self):
        """Slow polling for secondary tabs (only updates visible tab)."""
        if not self.comm.connected:
            return
        
        index = self.current_tab_index
        
        try:
            if index == 1:
                self.servo_tab.update_data()
            elif index == 2:
                self.temp_tab.update_data()
            elif index == 3:
                self.diode_tab.update_data()
            elif index == 4:
                self.laser_status_tab.update_data()
            elif index == 5:
                self.lbo_tab.update_data()
            elif index == 6:
                self.faults_tab.update_data()
        except Exception as e:
            print(f"❌ Slow poll error: {e}")
    
    def closeEvent(self, event):
        """Handle window close."""
        self.fast_timer.stop()
        self.slow_timer.stop()
        self.comm.disconnect()
        event.accept()


# =============================================================================
# Entry Point
# =============================================================================

def main():
    """Main entry point - handles both Spyder and standalone execution."""
    
    # Check if QApplication already exists (e.g., in Spyder's IPython console)
    app = QApplication.instance()
    new_app = app is None
    
    if new_app:
        app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    # Dark mode palette
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.Window, QColor(35, 35, 40))
    dark_palette.setColor(QPalette.WindowText, QColor(220, 220, 220))
    dark_palette.setColor(QPalette.Base, QColor(25, 25, 30))
    dark_palette.setColor(QPalette.AlternateBase, QColor(45, 45, 50))
    dark_palette.setColor(QPalette.ToolTipBase, QColor(220, 220, 220))
    dark_palette.setColor(QPalette.ToolTipText, QColor(220, 220, 220))
    dark_palette.setColor(QPalette.Text, QColor(220, 220, 220))
    dark_palette.setColor(QPalette.Button, QColor(45, 45, 50))
    dark_palette.setColor(QPalette.ButtonText, QColor(220, 220, 220))
    dark_palette.setColor(QPalette.BrightText, QColor(255, 50, 50))
    dark_palette.setColor(QPalette.Link, QColor(90, 150, 220))
    dark_palette.setColor(QPalette.Highlight, QColor(80, 120, 180))
    dark_palette.setColor(QPalette.HighlightedText, QColor(240, 240, 240))
    dark_palette.setColor(QPalette.Disabled, QPalette.Text, QColor(127, 127, 127))
    dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(127, 127, 127))
    
    app.setPalette(dark_palette)
    
    # Dark mode stylesheet for finer control
    app.setStyleSheet("""
        QToolTip { 
            color: #dcdcdc; 
            background-color: #2d2d30; 
            border: 1px solid #555555; 
        }
        QGroupBox {
            border: 1px solid #555555;
            border-radius: 6px;
            margin-top: 8px;
            padding-top: 8px;
            color: #dcdcdc;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 5px;
        }
        QTabWidget::pane {
            border: 1px solid #555555;
            border-radius: 4px;
            background-color: #23232a;
        }
        QTabBar::tab {
            padding: 8px 16px;
            margin-right: 2px;
            background-color: #2d2d35;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
            color: #aaaaaa;
        }
        QTabBar::tab:selected {
            background-color: #23232a;
            color: #ffffff;
            font-weight: bold;
        }
        QLineEdit {
            background-color: #2d2d35;
            border: 1px solid #555555;
            border-radius: 4px;
            padding: 4px;
            color: #dcdcdc;
        }
        QComboBox {
            background-color: #2d2d35;
            border: 1px solid #555555;
            border-radius: 4px;
            padding: 4px;
            color: #dcdcdc;
        }
        QComboBox::drop-down {
            border: none;
        }
        QComboBox QAbstractItemView {
            background-color: #2d2d35;
            color: #dcdcdc;
            selection-background-color: #4a6a9a;
        }
    """)
    
    # Close any existing Verdi window before creating new one
    for widget in app.topLevelWidgets():
        if isinstance(widget, VerdiV5MainWindow):
            # Stop timers first
            widget.fast_timer.stop()
            widget.slow_timer.stop()
            # Disconnect from hardware
            widget.comm.disconnect()
            # Close and schedule for deletion
            widget.close()
            widget.deleteLater()
    
    # Process events to ensure cleanup completes
    app.processEvents()
    
    # Create and show main window
    window = VerdiV5MainWindow()
    window.setAttribute(Qt.WA_DeleteOnClose, False)  # Don't auto-delete, we manage it
    window.show()
    
    # Store reference on app to prevent garbage collection
    app._verdi_window = window
    
    # Only start event loop if we created the app (standalone mode)
    if new_app:
        sys.exit(app.exec())
    
    return window  # Return for Spyder's benefit


if __name__ == "__main__":
    main()
