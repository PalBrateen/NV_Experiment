"""
SG384 Signal Generator Control Panel - PySide6 GUI

A comprehensive GUI matching the Stanford Research Systems SG384 hardware panel.
Features:
- Dark mode theme for light-sensitive experiments
- N-Type and BNC output sections with collapsible panels
- RF Doubler and Clock output sections
- Unit-aware input fields with suffix parsing (e.g., "2.87G", "100M", "1k")
- Real-time hardware synchronization
- Status bar with connection info

Author: Lab Control System
"""
# TODO: interface should show the latest change in the status bar temporary error display area
# TODO: when run, the file runs from a different folder while the session.py runs correctly... why? Same with pyPBLV.py. Also the python environment selected by default is not NV_py3.12...

import sys, re, ctypes, json, os
from pathlib import Path
from typing import Optional, Union, Tuple, Dict, Any
from enum import IntEnum
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QLineEdit, QComboBox,
    QGroupBox, QFrame, QStatusBar, QSizePolicy, QSpacerItem,
    QToolButton, QStyle, QScrollArea, QCheckBox
)
from PySide6.QtCore import Qt, Signal, Slot, QTimer, QSize, QUrl
from PySide6.QtGui import QPalette, QColor, QFont, QIcon, QPainter, QFontMetrics, QDesktopServices

expt_dir = os.path.abspath(r'D:\Brateen\NV_Experiment')    #absolute path to the directory containing the module
sys.path.append(expt_dir)      # Add the directory to sys.path
from SGcontrol import SignalGenerator, ModulationFunction, ModulationType, SGDISPLAY, ErrorCodes    # Import the module

# =============================================================================
# CONFIGURATION
# =============================================================================
WORKING_DIRECTORY = r'D:\Brateen\Saved_Data\SavedStates\SGStates'  # Default working directory for state files
DEFAULT_STATE_FILE = 'last_state.json'

# =============================================================================
# UNIT PARSING UTILITIES
# =============================================================================

class UnitParser:
    """
    Parse user input with SI prefixes and units.
    
    Supports:
    - SI prefixes: p, n, u/μ, m, k, M, G, T
    - Scientific notation: 1e9, 2.87E9
    - Combined: "2.87GHz", "100MHz", "1kHz"
    """
    
    SI_PREFIXES = {
        'T': 1e12, 'G': 1e9, 'M': 1e6, 'k': 1e3, 'K': 1e3,
        '': 1, 'm': 1e-3, 'u': 1e-6, 'μ': 1e-6, 'n': 1e-9, 'p': 1e-12
    }
    
    # Unit aliases (case-insensitive matching)
    UNIT_ALIASES = {
        'hz': 'Hz', 'khz': 'Hz', 'mhz': 'Hz', 'ghz': 'Hz',
        'dbm': 'dBm', 'db': 'dBm', 'v': 'V', 'vpp': 'Vpp',
        's': 's', 'ms': 's', 'us': 's', 'ns': 's',
        '%': '%', 'deg': 'deg', '°': 'deg'
    }
    
    @classmethod
    def parse(cls, text: str, default_unit: str = '') -> Tuple[float, str]:
        """
        Parse a value string with optional SI prefix and unit.
        
        Args:
            text: Input string like "2.87GHz", "100M", "1e9", "-10dBm"
            default_unit: Unit to assume if none provided
            
        Returns:
            Tuple of (value_in_base_units, unit_string)
        """
        text = text.strip()
        if not text:
            raise ValueError("Empty input")
        
        # Try scientific notation first (handles negative exponents)
        sci_match = re.match(r'^([+-]?\d*\.?\d+)[eE]([+-]?\d+)(.*)$', text)
        if sci_match:
            mantissa = float(sci_match.group(1))
            exponent = int(sci_match.group(2))
            remainder = sci_match.group(3).strip()
            value = mantissa * (10 ** exponent)
            unit = cls._extract_unit(remainder) or default_unit
            return value, unit
        
        # Match: number + optional SI prefix + optional unit
        pattern = r'^([+-]?\d*\.?\d+)\s*([TGMkKmuμnp]?)([A-Za-z%°]*)$'
        match = re.match(pattern, text)
        
        if match:
            number = float(match.group(1))
            prefix = match.group(2)
            unit_str = match.group(3)
            
            multiplier = cls.SI_PREFIXES.get(prefix, 1)
            value = number * multiplier
            unit = cls._extract_unit(unit_str) or default_unit
            
            return value, unit
        
        # Last resort: try to extract just a number
        try:
            return float(text), default_unit
        except ValueError:
            raise ValueError(f"Cannot parse: {text}")
    
    @classmethod
    def _extract_unit(cls, unit_str: str) -> str:
        """Normalize unit string."""
        if not unit_str:
            return ''
        unit_lower = unit_str.lower()
        return cls.UNIT_ALIASES.get(unit_lower, unit_str)
    
    @classmethod
    def format_value(cls, value: float, unit: str = '', precision: int = 6) -> str:
        """
        Format a value with appropriate SI prefix.
        
        Args:
            value: Value in base units
            unit: Unit string (Hz, dBm, etc.)
            precision: Number of significant figures
            
        Returns:
            Formatted string like "2.87 GHz"
        """
        if value == 0:
            return f"0 {unit}".strip()
        
        abs_val = abs(value)
        
        # Find appropriate prefix
        for prefix, multiplier in [('T', 1e12), ('G', 1e9), ('M', 1e6), 
                                    ('k', 1e3), ('', 1), ('m', 1e-3),
                                    ('μ', 1e-6), ('n', 1e-9), ('p', 1e-12)]:
            if abs_val >= multiplier * 0.999 or multiplier == 1e-12:
                scaled = value / multiplier
                # Format with appropriate precision
                if abs(scaled) >= 100:
                    formatted = f"{scaled:.{max(0, precision-3)}f}"
                elif abs(scaled) >= 10:
                    formatted = f"{scaled:.{max(0, precision-2)}f}"
                elif abs(scaled) >= 1:
                    formatted = f"{scaled:.{max(0, precision-1)}f}"
                else:
                    formatted = f"{scaled:.{precision}f}"
                # Remove trailing zeros after decimal
                if '.' in formatted:
                    formatted = formatted.rstrip('0').rstrip('.')
                return f"{formatted} {prefix}{unit}".strip()
        
        return f"{value} {unit}".strip()

class StatusBarLongLabel(QLabel):
    def __init__(self, parent=None):
        super().__init__(parent)
        # Set a minimum width so it doesn't disappear, 
        # but no maximum so it can take up available space.
        self.setMinimumWidth(100)

    def paintEvent(self, event):
        """Custom paint event to draw elided text."""
        painter = QPainter(self)
        metrics = QFontMetrics(self.font())
        
        # Calculate the elided text based on CURRENT label width
        # Qt.ElideRight puts the '...' at the end.
        elided_text = metrics.elidedText(self.text(), Qt.TextElideMode.ElideRight, self.width())
        
        # Draw the text within the label's rectangle
        painter.drawText(self.rect(), self.alignment(), elided_text)
        painter.end()

    # def set_error(self, message):
    #     """Update text and set the full message as a tooltip."""
    #     self.setText(message)
    #     self.setStyleSheet("color: #ff6b6b;") # Light red for dark mode errors
    #     self.setToolTip(message) # The OS handles the popup on hover

    def set_msg(self, message, error=False):
        """Update text and set the full message as a tooltip."""
        self.setText(message)
        if error:
            self.setStyleSheet("color: #ff6b6b;") # Light red for dark mode errors
        else:
            self.setStyleSheet("color: #4ae34a;") # Light green for messages
        self.setToolTip(message) # The OS handles the popup on hover

# =============================================================================
# CUSTOM WIDGETS
# =============================================================================

class UnitLineEdit(QLineEdit):
    """
    Line edit with unit parsing and display.
    
    Features:
    - Accepts input with SI prefixes (e.g., "2.87G", "100M")
    - Displays formatted value with appropriate prefix
    - Validates input on focus out
    """
    
    valueChanged = Signal(float)  # Emitted with value in base units
    
    def __init__(self, 
                 default_value: float = 0.,
                 unit: str = '',
                 min_value: float = 0.,
                 max_value: float = 0.,
                 parent=None):
        super().__init__(parent)
        
        self.unit = unit
        self.min_value = min_value
        self.max_value = max_value
        self._value = default_value
        self._is_editing = False
        
        # Display initial value
        self._update_display()
        
        # Connect signals
        self.editingFinished.connect(self._on_editing_finished)
        self.textChanged.connect(self._on_text_changed)
    
    def focusInEvent(self, event):
        """Select all text when gaining focus."""
        super().focusInEvent(event)
        self._is_editing = True
        QTimer.singleShot(0, self.selectAll)
    
    def focusOutEvent(self, event):
        """Validate and format on focus out."""
        self._is_editing = False
        super().focusOutEvent(event)
    
    def _on_text_changed(self, text):
        """Track that user is editing."""
        pass
    
    def _on_editing_finished(self):
        """Parse input and update value."""
        try:
            value, _ = UnitParser.parse(self.text(), self.unit)
            
            # Clamp to range if specified
            if self.min_value is not None:
                value = max(self.min_value, value)
            if self.max_value is not None:
                value = min(self.max_value, value)
            
            if value != self._value:
                self._value = value
                self.valueChanged.emit(value)
            
            self._update_display()
            
        except ValueError:
            # Revert to previous value
            self._update_display()
    
    def _update_display(self):
        """Update display with formatted value."""
        if not self._is_editing:
            self.setText(UnitParser.format_value(self._value, self.unit))
    
    def value(self) -> float:
        """Get current value in base units."""
        return self._value
    
    def setValue(self, value: float):
        """Set value (in base units)."""
        if self.min_value is not None:
            value = max(self.min_value, value)
        if self.max_value is not None:
            value = min(self.max_value, value)
        
        self._value = value
        self._update_display()


class ToggleButton(QPushButton):
    """
    Toggle button that changes color based on state.
    
    Green when ON, Red (or grey) when OFF.
    """
    
    toggled_state = Signal(bool)
    
    def __init__(self, text_on: str = "ON", text_off: str = "OFF", parent=None):
        super().__init__(parent)
        self.text_on = text_on
        self.text_off = text_off
        self._is_on = False
        
        self.setCheckable(True)
        self.clicked.connect(self._on_clicked)
        self._update_appearance()
    
    def _on_clicked(self):
        """Handle click."""
        self._is_on = self.isChecked()
        self._update_appearance()
        self.toggled_state.emit(self._is_on)
    
    def _update_appearance(self):
        """Update button appearance based on state."""
        if self._is_on:
            self.setText(self.text_on)
            self.setStyleSheet("""
                QPushButton {
                    background-color: #009de0;
                    color: white;
                    border: none;
                    border-radius: 4px;
                    padding: 2px 2px;
                    /*min-width: 40px;*/
                }

                QPushButton:hover {
                    background-color: #02b0fa;
                }
            """)
        else:
            self.setText(self.text_off)
            self.setStyleSheet("""
                QPushButton {
                    background-color: #4a4a4a;
                    border: 2px solid #666666;
                    color: #cccccc;
                    font-weight: bold;
                    padding: 8px 16px;
                    border-radius: 4px;
                }
                QPushButton:hover {
                    background-color: #5a5a5a;
                }
            """)
    
    def isOn(self) -> bool:
        """Get current state."""
        return self._is_on
    
    def setOn(self, state: bool):
        """Set state programmatically."""
        self._is_on = state
        self.setChecked(state)
        self._update_appearance()


class CollapsibleSection(QWidget):
    """
    A collapsible section with header and content.
    
    Click the header to expand/collapse the content.
    """
    
    def __init__(self, title: str, parent=None, collapsed: bool = False):
        super().__init__(parent)
        
        self._is_collapsed = collapsed
        
        # Main layout
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # Header button
        self.header = QPushButton(title)
        self.header.setCheckable(True)
        self.header.setChecked(not collapsed)
        self.header.clicked.connect(self._toggle)
        self.header.setStyleSheet("""
            QPushButton {
                background-color: #3a3a3f;
                border: 1px solid #555555;
                border-radius: 4px;
                color: #e0e0e0;
                font-weight: bold;
                font-size: 13px;
                padding: 8px 12px;
                text-align: left;
            }
            QPushButton:hover {
                background-color: #454550;
            }
            QPushButton:checked {
                background-color: #454550;
                border-bottom-left-radius: 0px;
                border-bottom-right-radius: 0px;
            }
        """)
        main_layout.addWidget(self.header)
        
        # Content widget
        self.content = QWidget()
        self.content_layout = QVBoxLayout(self.content)
        self.content_layout.setContentsMargins(8, 8, 8, 8)
        self.content.setStyleSheet("""
            QWidget {
                background-color: #2d2d32;
                border: 1px solid #555555;
                
                border-bottom-left-radius: 4px;
                border-bottom-right-radius: 4px;
            }
        """)
        main_layout.addWidget(self.content)
        
        # Set initial state
        self.content.setVisible(not collapsed)
        self._update_arrow()
    
    def _toggle(self):
        """Toggle collapsed state."""
        self._is_collapsed = not self.header.isChecked()
        self.content.setVisible(not self._is_collapsed)
        self._update_arrow()
    
    def _update_arrow(self):
        """Update header arrow indicator."""
        arrow = "▼" if not self._is_collapsed else "▶"
        # Get text without existing arrow
        text = self.header.text()
        if text.startswith("▼ ") or text.startswith("▶ "):
            text = text[2:]
        self.header.setText(f"{arrow} {text}")
    
    def addWidget(self, widget):
        """Add widget to content area."""
        self.content_layout.addWidget(widget)
    
    def addLayout(self, layout):
        """Add layout to content area."""
        self.content_layout.addLayout(layout)


# =============================================================================
# OUTPUT SECTION WIDGET
# =============================================================================

class DisplayToggleButton(QPushButton):
    """
    Small toggle button to switch the instrument display to show this parameter.
    
    Only one can be active at a time (managed by parent).
    """
    
    displayRequested = Signal(int)  # Emits SGDISPLAY value
    
    def __init__(self, display_mode: int, tooltip: str = "", parent=None):
        super().__init__("", parent)
        self.display_mode = display_mode
        self.setCheckable(True)
        self.setFixedSize(24, 24)
        self.setToolTip(tooltip or f"Show on display (DISP {display_mode})")
        
        self.clicked.connect(self._on_clicked)
        self._update_style()
    
    def _on_clicked(self):
        """Handle click - emit display request."""
        if self.isChecked():
            self.displayRequested.emit(self.display_mode)
        self._update_style()
    
    def _update_style(self):
        """Update button appearance based on state."""
        if self.isChecked():
            self.setStyleSheet("""
                QPushButton {
                    background-color: #009de0;
                    color: white;
                    border: none;
                    border-radius: 4px;
                    padding: 2px 2px;
                    /*min-width: 40px;*/
                }

                QPushButton:hover {
                    background-color: #02b0fa;
                }
            """)
        else:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #3a3a3f;
                    border: 1px solid #555555;
                    border-radius: 4px;
                    color: #909090;
                    font-weight: bold;
                    font-size: 10px;
                }
                QPushButton:hover {
                    background-color: #009de0;
                    color: #ffffff;
                }
            """)
    
    def setActive(self, active: bool):
        """Set the active state (called by parent to manage exclusivity)."""
        self.setChecked(active)
        self._update_style()


class OutputSection(QWidget):
    """
    Widget for a single output section (N-Type or BNC).
    
    Contains:
    - Output enable toggle
    - Frequency input with display toggle
    - Amplitude input with display toggle
    - Modulation settings (type, function, deviation, rate) with display toggles
    """
    
    # Signals
    outputToggled = Signal(bool)
    frequencyChanged = Signal(float)
    amplitudeChanged = Signal(float)
    modulationToggled = Signal(bool)
    modulationTypeChanged = Signal(str)
    modulationFunctionChanged = Signal(str)
    modulationRateChanged = Signal(float)
    modulationDeviationChanged = Signal(float)
    displayChangeRequested = Signal(int)  # SGDISPLAY value
    
    def __init__(self, 
                 title: str,
                 amp_unit: str = 'dBm',
                 amp_range: Tuple[float, float] = (-110, 16.5),
                 freq_range: Tuple[float, float] = (0, 4.05e9),
                 is_ntype: bool = True,  # True for N-Type, False for BNC
                 parent=None):
        super().__init__(parent)
        
        self.title = title
        self.amp_unit = amp_unit
        self.amp_range = amp_range
        self.freq_range = freq_range
        self.is_ntype = is_ntype
        
        # Track all display buttons for exclusivity
        self._display_buttons: list = []
        
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the UI layout."""

        pre_layout = QHBoxLayout(self)
        pre_layout.setSpacing(5)
        pre_layout.setContentsMargins(5, 5, 5, 5)

        output_switch_layout = QVBoxLayout()
        output_switch_layout.addStretch()
        # Output toggle
        lbl = QLabel("Output")
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        output_switch_layout.addWidget(lbl)

        self.output_btn = ToggleButton("ON", "OFF")
        self.output_btn.setFixedSize(70, 36)
        self.output_btn.toggled_state.connect(self.outputToggled.emit)
        output_switch_layout.addWidget(self.output_btn, Qt.AlignmentFlag.AlignCenter)
        output_switch_layout.addStretch()
        pre_layout.addLayout(output_switch_layout)
        # pre_layout.addStretch()

        layout1 = QGridLayout()
        # layout1.setSpacing(10)
        layout1.setContentsMargins(8, 5, 0, 5)
        
        col = 0
        # dummy
        lbl = QLabel("")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("height: 4px; color: #ffffff; font-size: 4px; font-weight: bold; border: none;")
        layout1.addWidget(lbl, 0, col)

        lbl = QLabel("Frequency")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        layout1.addWidget(lbl, 1, col)

        # Frequency input + display button
        self.freq_input = UnitLineEdit(
            default_value=2.87e9,
            unit='Hz',
            min_value=self.freq_range[0],
            max_value=self.freq_range[1]
        )
        self.freq_input.setFixedWidth(120)
        self.freq_input.valueChanged.connect(self.frequencyChanged.emit)
        self.freq_input.setStyleSheet("background-color: #484848; color: #ffffff;")
        layout1.addWidget(self.freq_input, 2, col)
        
        self.freq_disp_btn = DisplayToggleButton(
            SGDISPLAY.FREQUENCY, "Show FREQUENCY on display"
        )
        self.freq_disp_btn.displayRequested.connect(self._on_display_requested)
        self._display_buttons.append(self.freq_disp_btn)
        layout1.addWidget(self.freq_disp_btn, 2, col+1)
        
        lbl = QLabel("Amplitude")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        layout1.addWidget(lbl, 3, col)

        # Amplitude input + display button
        self.amp_input = UnitLineEdit(
            default_value=0.0,
            unit=self.amp_unit,
            min_value=self.amp_range[0],
            max_value=self.amp_range[1]
        )
        self.amp_input.setFixedWidth(120)
        self.amp_input.valueChanged.connect(self.amplitudeChanged.emit)
        self.amp_input.setStyleSheet("background-color: #484848; color: #ffffff;")
        layout1.addWidget(self.amp_input, 4, col)
        
        amp_display = SGDISPLAY.AMPLITUDE_NTYPE if self.is_ntype else SGDISPLAY.AMPLITUDE_BNC
        self.amp_disp_btn = DisplayToggleButton(
            amp_display, f"Show {'N-TYPE' if self.is_ntype else 'BNC'} AMPLITUDE on display"
        )
        self.amp_disp_btn.displayRequested.connect(self._on_display_requested)
        self._display_buttons.append(self.amp_disp_btn)
        layout1.addWidget(self.amp_disp_btn, 4, col+1)
        
        # dummy
        lbl = QLabel("")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("height: 4px; color: #ffffff; font-size: 4px; font-weight: bold; border: none;")
        layout1.addWidget(lbl, 5, col)

        pre_layout.addLayout(layout1)
        pre_layout.addSpacing(5)

        mod_group = QGroupBox("Modulation")
        mod_layout = QHBoxLayout(mod_group)

        layout3 = QVBoxLayout()
        layout3.addStretch()
        
        # mod_label = QLabel("Modulation")
        # mod_label.setFixedSize(80, 18)
        # mod_label.setStyleSheet("color: #ffffff; border: None; font-size: 13px; font-weight: bold;")
        # layout3.addWidget(mod_label, Qt.AlignmentFlag.AlignRight)
        
        self.mod_enable_btn = ToggleButton("ON", "OFF")
        self.mod_enable_btn.setFixedSize(70, 36)
        self.mod_enable_btn.toggled_state.connect(self.modulationToggled.emit)
        layout3.addWidget(self.mod_enable_btn, Qt.AlignmentFlag.AlignRight)

        combo_layout = QGridLayout()
        combo_layout.setSpacing(5)
        # combo_layout.setContentsMargins(5, 5, 5, 5)

        col = 0
        lbl = QLabel("Type")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        combo_layout.addWidget(lbl, 0, col)#, Qt.AlignmentFlag.AlignCenter)

        # # Modulation type combo + display button
        self.mod_type_combo = QComboBox()
        self.mod_type_combo.addItems(ModulationType._member_names_)
        self.mod_type_combo.currentTextChanged.connect(self._on_mod_type_changed)
        self.mod_type_combo.setFixedWidth(100)
        self.mod_type_combo.setStyleSheet("background-color: #484848; color: #ffffff;")
        combo_layout.addWidget(self.mod_type_combo, 1, col)#, Qt.AlignmentFlag.AlignCenter)

        self.mod_type_disp_btn = DisplayToggleButton(
            SGDISPLAY.MODULATION_TYPE, "Show MODULATION TYPE on display"
        )
        self.mod_type_disp_btn.displayRequested.connect(self._on_display_requested)
        self._display_buttons.append(self.mod_type_disp_btn)
        combo_layout.addWidget(self.mod_type_disp_btn, 1, col+1)

        layout3.addLayout(combo_layout)
        layout3.addStretch()
        mod_layout.addLayout(layout3)
        # mod_group
        
        # pre_layout.addLayout(layout3)
        # pre_layout.addStretch()
        
        layout2 = QGridLayout()
        # layout2.setSpacing(5)
        # layout2.setContentsMargins(5, 5, 0, 5)

        col = 0
        lbl = QLabel("Function")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        layout2.addWidget(lbl, 2, col)

        # Modulation function combo + display button
        self.mod_func_combo = QComboBox()
        self.mod_func_combo.addItems(ModulationFunction._member_names_)
        self.mod_func_combo.currentTextChanged.connect(self.modulationFunctionChanged.emit)
        self.mod_func_combo.setFixedWidth(100)
        self.mod_func_combo.setStyleSheet("background-color: #484848; color: #ffffff;")
        layout2.addWidget(self.mod_func_combo, 3, col)
        
        self.mod_func_disp_btn = DisplayToggleButton(
            SGDISPLAY.MODULATION_FUNCTION, "Show MODULATION FUNCTION on display"
        )
        self.mod_func_disp_btn.displayRequested.connect(self._on_display_requested)
        self._display_buttons.append(self.mod_func_disp_btn)
        layout2.addWidget(self.mod_func_disp_btn, 3, col+1)
        # col += 2
        
        lbl = QLabel("Deviation")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        layout2.addWidget(lbl, 0, col)
        
        # Modulation Deviation input + display button
        self.mod_dev_input = UnitLineEdit(
            default_value=200e3,
            unit='Hz',
            min_value=1.,
            max_value=50e3
        )
        self.mod_dev_input.setFixedWidth(100)
        self.mod_dev_input.valueChanged.connect(self.modulationDeviationChanged.emit)
        self.mod_dev_input.setStyleSheet("background-color: #484848; color: #ffffff;")
        layout2.addWidget(self.mod_dev_input, 1, col)
        
        # Deviation display toggle
        self.mod_dev_disp_btn = DisplayToggleButton(
            SGDISPLAY.MODULATION_DEVIATION, "Show MODULATION DEVIATION on display"
        )
        self.mod_dev_disp_btn.displayRequested.connect(self._on_display_requested)
        self._display_buttons.append(self.mod_dev_disp_btn)
        layout2.addWidget(self.mod_dev_disp_btn, 1, col+1)

        # Rate label
        lbl = QLabel("Rate")
        lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lbl.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold; border: none;")
        layout2.addWidget(lbl, 4, col)

        # Modulation rate input + display button
        self.mod_rate_input = UnitLineEdit(
            default_value=1e3,
            unit='Hz',
            min_value=0.001,
            max_value=50e3
        )
        self.mod_rate_input.setFixedWidth(100)
        self.mod_rate_input.valueChanged.connect(self.modulationRateChanged.emit)
        self.mod_rate_input.setStyleSheet("background-color: #484848; color: #ffffff;")
        layout2.addWidget(self.mod_rate_input, 5, col)
        
        self.mod_rate_disp_btn = DisplayToggleButton(
            SGDISPLAY.MODULATION_RATE, "Show MODULATION RATE on display"
        )
        self.mod_rate_disp_btn.displayRequested.connect(self._on_display_requested)
        self._display_buttons.append(self.mod_rate_disp_btn)
        layout2.addWidget(self.mod_rate_disp_btn, 5, col+1)
        
        mod_layout.addLayout(layout2)
        pre_layout.addWidget(mod_group)
    
    def _on_display_requested(self, display_mode: int):
        """Handle display button click - ensure only one is active."""
        # Uncheck all other display buttons
        for btn in self._display_buttons:
            if btn.display_mode != display_mode:
                btn.setActive(False)
        
        # Emit signal to change instrument display
        self.displayChangeRequested.emit(display_mode)
    
    def _on_mod_type_changed(self, text: str):
        """Handle modulation type change."""
        # Map display text to enum name
        # type_map = {
        #     "NONE": "NONE", "AM": "AMPLITUDE", "FM": "FREQUENCY",
        #     "PM": "PHASE", "SWEEP": "SWEEP", "PULSE": "PULSE",
        #     "BLANK": "BLANK", "IQ": "IQ"
        # }
        # self.modulationTypeChanged.emit(type_map.get(text, text))
        self.modulationTypeChanged.emit(ModulationType[text.upper()].name)

    
    # Getters for current values
    def getOutputState(self) -> bool:
        return self.output_btn.isOn()
    
    def getFrequency(self) -> float:
        return self.freq_input.value()
    
    def getAmplitude(self) -> float:
        return self.amp_input.value()
    
    def getModulationType(self) -> str:
        return self.mod_type_combo.currentText()
    
    def getModulationFunction(self) -> str:
        return self.mod_func_combo.currentText()
    
    def getModulationDeviation(self) -> float:
        return self.mod_dev_input.value()

    def getModulationRate(self) -> float:
        return self.mod_rate_input.value()
    
    def getModulationState(self) -> bool:
        return self.mod_enable_btn.isOn()
    
    # Setters for values from hardware
    def setOutputState(self, state: bool):
        self.output_btn.setOn(state)
    
    def setFrequency(self, value: float):
        self.freq_input.setValue(value)
    
    def setAmplitude(self, value: float):
        self.amp_input.setValue(value)
    
    def setModulationType(self, type_str: str):
        # Map enum name to display text
        # type_map = {
        #     "NONE": "NONE", "AMPLITUDE": "AM", "FREQUENCY": "FM",
        #     "PHASE": "PM", "SWEEP": "SWEEP", "PULSE": "PULSE",
        #     "BLANK": "BLANK", "IQ": "IQ"
        # }
        # display_text = ModulationType.get(type_str.upper(), type_str)
        display_text = ModulationType[type_str.upper()].name
        idx = self.mod_type_combo.findText(display_text)
        if idx >= 0:
            self.mod_type_combo.setCurrentIndex(idx)
    
    def setModulationFunction(self, func_str: str):
        idx = self.mod_func_combo.findText(func_str.upper())
        if idx >= 0:
            self.mod_func_combo.setCurrentIndex(idx)
    
    def setModulationDeviation(self, value: float):
        self.mod_dev_input.setValue(value)
    
    def setModulationRate(self, value: float):
        self.mod_rate_input.setValue(value)
    
    def setModulationState(self, state: bool):
        self.mod_enable_btn.setOn(state)


# =============================================================================
# AUXILIARY OUTPUT SECTION
# =============================================================================

class AuxOutputSection(QWidget):
    """
    Compact widget for auxiliary outputs (RF Doubler, Clock).
    
    Contains:
    - Amplitude input
    - Offset input (if applicable)
    """
    
    amplitudeChanged = Signal(float)
    offsetChanged = Signal(float)
    
    def __init__(self, 
                 title: str,
                 has_offset: bool = True,
                 amp_unit: str = 'dBm',
                 offset_unit: str = 'V',
                 parent=None):
        super().__init__(parent)
        
        self.has_offset = has_offset
        
        layout = QHBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(5)
        
        # Amplitude
        amp_layout = QVBoxLayout()
        amp_label = QLabel("Amplitude")
        amp_label.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold;")
        amp_layout.addWidget(amp_label)
        
        self.amp_input = UnitLineEdit(default_value=0.0, unit=amp_unit)
        self.amp_input.setMinimumWidth(100)
        self.amp_input.valueChanged.connect(self.amplitudeChanged.emit)
        amp_layout.addWidget(self.amp_input)
        
        amp_hint = QLabel(f"[{amp_unit}]")
        amp_hint.setStyleSheet("color: #707070; font-size: 10px;")
        amp_layout.addWidget(amp_hint)
        
        layout.addLayout(amp_layout)
        
        # Offset (if applicable)
        if has_offset:
            offset_layout = QVBoxLayout()
            offset_label = QLabel("Offset")
            offset_label.setStyleSheet("color: #ffffff; font-size: 13px; font-weight: bold;")
            offset_layout.addWidget(offset_label)
            
            self.offset_input = UnitLineEdit(default_value=0.0, unit=offset_unit)
            self.offset_input.setMinimumWidth(100)
            self.offset_input.valueChanged.connect(self.offsetChanged.emit)
            offset_layout.addWidget(self.offset_input)
            
            offset_hint = QLabel(f"[{offset_unit}]")
            offset_hint.setStyleSheet("color: #707070; font-size: 10px;")
            offset_layout.addWidget(offset_hint)
            
            layout.addLayout(offset_layout)
        
        layout.addStretch()

# =============================================================================
# MAIN WINDOW
# =============================================================================

class SG384ControlPanel(QMainWindow):
    """
    Main control panel for the SG384 Signal Generator.
    
    Replicates the hardware front panel layout with:
    - N-Type output section (main, expanded)
    - BNC output section (collapsed by default)
    - RF Doubler output section (collapsed)
    - Clock output section (collapsed)
    - Status bar with connection info
    """
    
    def __init__(self, sg_instance, auto_connect_hardware=False, parent=None):
        super().__init__(parent)
        
        self.sg: SignalGenerator = sg_instance  # SignalGenerator instance (optional)
        self._connected = False
        self._address = "Not connected"
        self._signals_connected = False  # Track if signals are connected
        
        # Working directory for state files
        self.working_directory = Path(WORKING_DIRECTORY)
        self.working_directory.mkdir(parents=True, exist_ok=True)
        
        self.setWindowTitle("SG384 Signal Generator Control")
        self.setMinimumSize(575, 500)
        
        self._apply_dark_theme()
        self._setup_ui()
        self._setup_status_bar()
        
        # Initialize button states
        self._update_connect_button()
        
        # Refresh state files list
        self._refresh_state_files()
        self._load_state()
        
        # If we have a signal generator instance and auto_connect is True, connect to it
        if self.sg is not None and auto_connect_hardware:
            self._connect()
    
    def _apply_dark_theme(self):
        """Apply dark mode theme to the application."""
        app = QApplication.instance()
        app.setStyle('Fusion')
        
        dark_palette = QPalette()
        
        # Base colors
        dark_palette.setColor(QPalette.ColorRole.Window, "#232328")
        dark_palette.setColor(QPalette.ColorRole.WindowText, "#dcdcdc")
        dark_palette.setColor(QPalette.ColorRole.Base, "#19191e")
        dark_palette.setColor(QPalette.ColorRole.AlternateBase, "#2d2d32")
        dark_palette.setColor(QPalette.ColorRole.ToolTipBase, "#2d2d32")
        dark_palette.setColor(QPalette.ColorRole.ToolTipText, "#dcdcdc")
        dark_palette.setColor(QPalette.ColorRole.Text, "#dcdcdc")
        dark_palette.setColor(QPalette.ColorRole.Button, "#323237")
        dark_palette.setColor(QPalette.ColorRole.ButtonText, "#dcdcdc")
        dark_palette.setColor(QPalette.ColorRole.BrightText, "#ff3232")
        dark_palette.setColor(QPalette.ColorRole.Link, "#5a96dc")
        dark_palette.setColor(QPalette.ColorRole.Highlight, "#5078b4")
        dark_palette.setColor(QPalette.ColorRole.HighlightedText, "#f0f0f0")
        
        # Disabled colors
        dark_palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, "#7f7f7f")
        dark_palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, "#7f7f7f")
        
        app.setPalette(dark_palette)
        
        # Global stylesheet
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background-color: #242424;
                color: #ffffff;
                font-family: 'Segoe UI', Arial, sans-serif;
                font-size: 12px;
                /*border: 3px dash #3f3f46*/
            }
            QLineEdit, QSpinBox, QComboBox {
                background-color: #484848;
                color: #ffffff;
                border: 1px solid #3f3f46;
                border-radius: 3px;
                padding: 4px;
                selection-background-color: #264f78;
                font-size: 12px;
                font-weight: bold;
            }

            QLineEdit:focus, QSpinBox:focus, QComboBox:focus {
                border: 1px solid #009de0;
            }

            QLineEdit:disabled, QSpinBox:disabled {
                background-color: #252526;
                color: #6d6d6d;
            }
            /*QLineEdit {
                background-color: #2a2a30;
                border: 1px solid #555555;
                border-radius: 4px;
                padding: 6px 10px;
                color: #00ff88;
                font-family: 'Consolas', 'Courier New', monospace;
                font-size: 13px;
                font-weight: bold;
            }
            QLineEdit:focus {
                border: 2px solid #4a90d9;
                background-color: #35353a;
            }*/
            QComboBox {
                padding-right: 20px;
            }

            QComboBox::drop-down {
                border: none;
                width: 20px;
            }

            QComboBox::down-arrow {
                width: 0px; height: 0px;
                border-left: 4px solid #484848;
                border-right: 4px solid #484848;
                border-top: 5px solid #ffffff;
                margin-right: 5px;
            }

            QComboBox QAbstractItemView {
                background-color: #2d2d30;
                color: #ffffff;
                selection-background-color: #094771;
                border: 1px solid #3f3f46;
            }
            QLabel {
                color: #c0c0c0;
            }
            QStatusBar {
                background-color: #1e1e22;
                color: #909090;
                border-top: 2px solid #878787;
            }
            QGroupBox {
                font-weight: bold;
                border: 1px solid #404040;
                border-radius: 6px;
                margin-top: 12px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #ffffff;
            }
                           QGroupBox {
    border: 3px solid #3f3f46;
    border-radius: 5px;
    margin-top: 8px;
    padding-top: 0px;
    color: #ffffff;
}

QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 0px;
}
        """)
    
    def _setup_ui(self):
        """Setup the main UI layout."""
        # Create scroll area for the entire content
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setStyleSheet("""
            QScrollArea {
                border: none;
                background-color: #232328;
            }
            QScrollBar:vertical {
                background-color: #2a2a30;
                width: 12px;
                border-radius: 6px;
                margin: 2px;
            }
            QScrollBar::handle:vertical {
                background-color: #555555;
                border-radius: 5px;
                min-height: 30px;
            }
            QScrollBar::handle:vertical:hover {
                background-color: #666666;
            }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                height: 0px;
            }
            QScrollBar:horizontal {
                background-color: #2a2a30;
                height: 12px;
                border-radius: 6px;
                margin: 2px;
            }
            QScrollBar::handle:horizontal {
                background-color: #555555;
                border-radius: 5px;
                min-width: 30px;
            }
            QScrollBar::handle:horizontal:hover {
                background-color: #666666;
            }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                width: 0px;
            }
        """)
        
        # Content widget inside scroll area
        content_widget = QWidget()
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(content_widget)
        main_layout.setSpacing(5)
        # main_layout.setContentsMargins(16, 16, 16, 16)
        
        # State management and connection options panel
        options_layout = QHBoxLayout()
        options_layout.setSpacing(15)
        
        # State file management group
        state_group = QGroupBox("State Management")
        state_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #3f3f46;
                border-radius: 5px;
                margin-top: 8px;
                padding-top: 5px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #b0b0b0;
            }
        """)
        state_layout = QVBoxLayout(state_group)
        state_layout.setSpacing(5)
        state_layout.setContentsMargins(8, 12, 8, 8)
        
        # File name row with editable combo
        filename_row = QHBoxLayout()
        filename_row.addWidget(QLabel("File:"))
        
        self.state_combo = QComboBox()
        self.state_combo.setEditable(True)
        self.state_combo.setMinimumWidth(150)
        self.state_combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.state_combo.lineEdit().setPlaceholderText("Enter name or select...")
        filename_row.addWidget(self.state_combo)
        
        open_folder_btn = QPushButton("...")
        open_folder_btn.setFixedSize(25, 25)
        open_folder_btn.setStyleSheet("""
            QPushButton {
                background-color: #484848;
                border: 1px solid #555555;
                border-radius: 3px;
                color: #ffffff;
            }
            QPushButton:hover {
                background-color: #5a5a5a;
            }
        """)
        open_folder_btn.clicked.connect(lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(WORKING_DIRECTORY)))
        filename_row.addWidget(open_folder_btn)

        # Refresh file list button
        self.refresh_files_btn = QPushButton("↻")
        self.refresh_files_btn.setFixedSize(25, 25)
        self.refresh_files_btn.setToolTip("Refresh file list")
        self.refresh_files_btn.clicked.connect(self._refresh_state_files)
        self.refresh_files_btn.setStyleSheet("""
            QPushButton {
                background-color: #484848;
                border: 1px solid #555555;
                border-radius: 3px;
                color: #ffffff;
            }
            QPushButton:hover {
                background-color: #5a5a5a;
            }
        """)
        filename_row.addWidget(self.refresh_files_btn)
        
        state_layout.addLayout(filename_row)
        
        # Load/Save buttons row
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        
        self.load_state_btn = QPushButton("Load")
        self.load_state_btn.setFixedSize(60, 25)
        self.load_state_btn.clicked.connect(self._load_state)
        self.load_state_btn.setStyleSheet("""
            QPushButton {
                background-color: #484848;
                border: 1px solid #555555;
                border-radius: 3px;
                color: #ffffff;
            }
            QPushButton:hover {
                background-color: #009de0;
            }
        """)
        btn_row.addWidget(self.load_state_btn)
        
        self.save_state_btn = QPushButton("Save")
        self.save_state_btn.setFixedSize(60, 25)
        self.save_state_btn.clicked.connect(self._save_state)
        self.save_state_btn.setStyleSheet("""
            QPushButton {
                background-color: #d3dbde;
                border: 1px solid #555555;
                border-radius: 3px;
                color: #484848;
            }
            QPushButton:hover {
                background-color: #009de0;
                color: #ffffff;
            }
        """)
        btn_row.addWidget(self.save_state_btn)
        
        state_layout.addLayout(btn_row)
        options_layout.addWidget(state_group)
        
        # Connection options group
        conn_group = QGroupBox("Connection Options")
        conn_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #3f3f46;
                border-radius: 5px;
                margin-top: 8px;
                padding-top: 5px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #b0b0b0;
            }
        """)
        conn_layout = QVBoxLayout(conn_group)
        conn_layout.setSpacing(8)
        conn_layout.setContentsMargins(8, 12, 8, 8)
        
        # Checkbox: Initialize SG with GUI parameters (overwrite SG)
        self.init_with_gui_cb = QCheckBox("Init SG with GUI params")
        self.init_with_gui_cb.setToolTip(
            "When checked: On connect, send GUI values TO the signal generator.\n"
            "When unchecked: On connect, read values FROM the signal generator to GUI."
        )
        self.init_with_gui_cb.setChecked(False)
        self.init_with_gui_cb.setStyleSheet("""
            QCheckBox {
                color: #b0b0b0;
                spacing: 5px;
            }
            QCheckBox::indicator {
                width: 16px;
                height: 16px;
                border: 1px solid #555555;
                border-radius: 3px;
                background-color: #2a2a30;
            }
            QCheckBox::indicator:checked {
                background-color: #4a9a4a;
                border-color: #6aba6a;
            }
            QCheckBox::indicator:hover {
                border-color: #777777;
            }
        """)
        conn_layout.addWidget(self.init_with_gui_cb)
        
        # Info label
        info_label = QLabel("↑ Overwrite SG on connect")
        info_label.setStyleSheet("color: #707070; font-size: 10px;")
        conn_layout.addWidget(info_label)
        
        conn_layout.addStretch()
        options_layout.addWidget(conn_group)
        
        options_layout.addStretch()
        # main_layout.addLayout(options_layout)

        # Title bar with model info
        top_layout = QVBoxLayout()
        
        # SRS Logo placeholder
        # self.logo_label = QLabel("SRS SG384")
        # self.logo_label.setStyleSheet("""
        #     font-size: 18px;
        #     font-weight: bold;
        #     color: #ffffff;
        #     padding: 4px 8px;
        #     border: 2px solid #ffffff;
        #     border-radius: 4px;
        # """)
        # top_layout.addWidget(self.logo_label)
        # top_layout.addStretch()
        
        # Connect/Disconnect button
        self.connect_btn = QPushButton("Connect")
        # self.connect_btn.setStyleSheet("""
        #     QPushButton {
        #         background-color: #009de0;
        #         border: 1px solid #4a9a4a;
        #         border-radius: 4px;
        #         color: #ffffff;
        #         padding: 6px 16px;
        #         font-weight: bold;
        #         min-width: 100px;
        #     }
        #     QPushButton:hover {
        #         background-color: #009de0;
        #     }
        # """)
        self.connect_btn.setStyleSheet("""
            QPushButton {
            font-size: 18px;
            font-weight: bold;
            color: #ffffff;
            padding: 4px 8px;
            border: 2px solid #ffffff;
            border-radius: 4px;
            }
        """)
        self.connect_btn.clicked.connect(self._toggle_connection)
        top_layout.addWidget(self.connect_btn)
        
        # Refresh button
        self.refresh_btn = QPushButton("🔄 Refresh")
        self.refresh_btn.setStyleSheet("""
            QPushButton {
                background-color: #3a5a8a;
                border: 1px solid #5a7aaa;
                border-radius: 4px;
                color: #e0e0e0;
                padding: 6px 12px;
            }
            QPushButton:hover {
                background-color: #4a6a9a;
            }
            QPushButton:disabled {
                background-color: #2a2a30;
                border: 1px solid #404040;
                color: #606060;
            }
        """)
        self.refresh_btn.setFixedSize(130, 30)
        self.refresh_btn.clicked.connect(self._sync_from_hardware)
        self.refresh_btn.setEnabled(False)  # Disabled until connected
        top_layout.addWidget(self.refresh_btn)
        
        # main_layout.addLayout(top_layout)
        options_layout.addLayout(top_layout)
        main_layout.addLayout(options_layout)
        
        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("background-color: #404040;")
        main_layout.addWidget(sep)
        
        # N-Type Output Section (expanded by default)
        self.ntype_section = CollapsibleSection("N-TYPE OUTPUT (950 kHz to 4.05 GHz)", collapsed=False)
        self.ntype_output = OutputSection(
            "N-Type",
            amp_unit='dBm',
            amp_range=(-110, 16.5),
            freq_range=(950e3, 4.05e9),
            is_ntype=True
        )
        self.ntype_section.addWidget(self.ntype_output)
        main_layout.addWidget(self.ntype_section)
        
        # BNC Output Section (collapsed by default)
        self.bnc_section = CollapsibleSection("BNC OUTPUT (DC to 62.5 MHz)", collapsed=True)
        self.bnc_output = OutputSection(
            "BNC",
            amp_unit='Vpp',
            amp_range=(-47, 13),  # dBm range, but can also be Vpp
            freq_range=(0, 62.5e6),
            is_ntype=False
        )
        self.bnc_section.addWidget(self.bnc_output)
        main_layout.addWidget(self.bnc_section)
        
        # Auxiliary outputs in a horizontal layout
        aux_layout = QHBoxLayout()
        
        # RF Doubler Section (collapsed)
        self.doubler_section = CollapsibleSection("RF DOUBLER (4 to 8 GHz)", collapsed=True)
        self.doubler_output = AuxOutputSection(
            "RF Doubler",
            has_offset=False,
            amp_unit='dBm'
        )
        self.doubler_section.addWidget(self.doubler_output)
        aux_layout.addWidget(self.doubler_section)
        
        # Clock Output Section (collapsed)
        self.clock_section = CollapsibleSection("CLOCK OUTPUT", collapsed=True)
        self.clock_output = AuxOutputSection(
            "Clock",
            has_offset=True,
            amp_unit='dBm',
            offset_unit='V'
        )
        self.clock_section.addWidget(self.clock_output)
        aux_layout.addWidget(self.clock_section)
        
        main_layout.addLayout(aux_layout)
        
        # Spacer to push content up
        main_layout.addStretch()
        
        # QWidget combining
        # Set scroll area as central widget
        self.setCentralWidget(scroll_area)
    
    def _setup_status_bar(self):
        """Setup the status bar with connection info and error display."""
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        
        # Address label
        self.addr_label = StatusBarLongLabel()
        self.addr_label.setStyleSheet("color: #aaaaaa;")
        self.addr_label.setFixedWidth(180)
        self.status_bar.addWidget(self.addr_label, 1)

        # Connection status label
        self.conn_status_label = QLabel()
        self.conn_status_label.setFixedWidth(100)
        self.conn_status_label.setStyleSheet("padding-right: 5px;")
        self.status_bar.addWidget(self.conn_status_label, 1)

        # # Error display label (left side, temporary messages)
        # self.error_label = StatusBarLongLabel()
        # self.error_label.setStyleSheet("""
        #     color: #ff6b6b;
        #     font-weight: bold;
        #     padding: 2px 10px;
        # """)
        # self.status_bar.addWidget(self.error_label, 1)  # stretch=1 to take available space
        
        # Last error code display
        self.last_error_label = StatusBarLongLabel("Last Error: None")#QLabel("Last Error: None")
        self.last_error_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.last_error_label.setStyleSheet("""
            color: #909090;
            padding-right: 5px;
        """)
        self.status_bar.addWidget(self.last_error_label, 1)

        self._update_status_bar()
        
        # Timer to clear temporary error messages
        self._error_clear_timer = QTimer(self)
        self._error_clear_timer.setSingleShot(True)
        self._error_clear_timer.timeout.connect(self._clear_error_message)
    
    def _update_status_bar(self):
        """Update status bar with current connection info."""
        if self._connected:
            self.conn_status_label.setText("● Connected")
            self.conn_status_label.setStyleSheet("""
                color: #4ae34a;
                font-weight: bold;
                padding-right: 5px;
            """)
        else:
            self.conn_status_label.setText("○ Disconnected")
            self.conn_status_label.setStyleSheet("""
                color: #e34a4a;
                font-weight: bold;
                padding-right: 5px;
            """)
        
        self.addr_label.set_msg(f"Address: {self._address}")
        self.addr_label.setStyleSheet("color: #aaaaaa;")
    
    def _show_error(self, error_code: int, context: str = ""):
        """Display an error in the status bar."""
        if error_code == 0:
            return  # No error
        
        try:
            error_name = ErrorCodes(error_code).name
        except ValueError:
            error_name = f"UNKNOWN_{error_code}"
        
        description = ErrorCodes.get_description(error_code)
        
        # Update last error display
        # self.last_error_label.set_msg(f"Last Error: {error_code} ({error_name})", error=True)
        # self.last_error_label.setStyleSheet("""
        #     color: #ff6b6b;
        #     font-weight: bold;
        #     padding-right: 15px;
        # """)
        
        # Show temporary error message
        context_str = f" [{context}]" if context else ""
        self.last_error_label.set_msg(f"⚠ Error {error_code}: {description}{context_str}", error=True)
        self.last_error_label.setStyleSheet("""
            color: #ff6b6b;
            font-weight: bold;
            padding-right: 15px;
        """)
        
        # Clear message after 5 seconds
        self._error_clear_timer.start(5000)
        
        # print(f"⚠ SG384 Error {error_code} ({error_name}): {description}{context_str}")
    
    def _clear_error_message(self):
        """Clear the temporary error message."""
        # self.error_label.set_error("")
        pass
    
    def _clear_last_error(self):
        """Clear the last error display (call after successful operation)."""
        self.last_error_label.set_msg("Last Error: None")
        self.last_error_label.setStyleSheet("""
            color: #4ae34a;
            padding-right: 15px;
        """)
    
    def _check_and_report_error(self, context: str = "") -> int:
        """
        Check for errors from the signal generator and report them.
        
        Returns:
            Error code (0 if no error)
        """
        if self.sg is None or not self._connected:
            return 0
        
        try:
            # Query last error
            error_code = int(self.sg._instr.query('LERR?'))
            
            if error_code != 0:
                self._show_error(error_code, context)
            else:
                # Successful operation - could optionally clear last error
                pass
            
            return error_code
            
        except Exception as e:
            # print(f"Error checking SG status: {e}")
            self.last_error_label.set_msg(f"Error checking SG status: {e}", error=True)
            return -1
    
    def _safe_execute(self, func, *args, context: str = "", **kwargs) -> bool:
        """
        Safely execute a signal generator command with error checking.
        
        Args:
            func: The function to execute
            *args: Arguments for the function
            context: Description of the operation for error messages
            **kwargs: Keyword arguments for the function
            
        Returns:
            True if successful, False if error occurred
        """
        if self.sg is None or not self._connected:
            self._show_error(-1, "Not connected")
            return False
        
        try:
            # Execute the command
            func(*args, **kwargs)
            
            # Check for errors
            error_code = self._check_and_report_error(context)
            if error_code == 0:
                self.last_error_label.set_msg(f"{context}: {args}")
                return True
            else:
                return False
            # return error_code == 0
            
        except Exception as e:
            self.last_error_label.set_msg(f"⚠ Exception: {str(e)}", error=True)
            # self._error_clear_timer.start(5000)
            # print(f"Exception during {context}: {e}")
            return False
    
    def _update_connect_button(self):
        """Update connect button appearance based on connection state."""
        if self._connected:
            self.connect_btn.setText("Connected")
            # self.connect_btn.setStyleSheet("""
            #     QPushButton {
            #         background-color: #009de0;
            #         color: #ffffff;
            #         border: 1px solid #009de0;
            #         border-radius: 4px;
            #         padding: 6px 16px;
            #         font-weight: bold;
            #         min-width: 100px;
            #     }
            #     QPushButton:hover {
            #         background-color: #484848;
            #     }
            # """)
            self.connect_btn.setStyleSheet("""
                QPushButton {
                    font-size: 18px;
                    font-weight: bold;
                    color: #36d21d;
                    padding: 4px 8px;
                    border: 2px solid #36d21d;
                    border-radius: 4px;
                }
                QPushButton:hover {
                    color: #db8972;
                    border: 2px solid #db8972;
                }
            """)
            # self.logo_label.setStyleSheet("""
            #     QLabel {
            #         font-size: 18px;
            #         font-weight: bold;
            #         color: #36d21d;
            #         padding: 4px 8px;
            #         border: 2px solid #36d21d;
            #         border-radius: 4px;
            #                 }
            # """)
            self.refresh_btn.setEnabled(True)
        else:
            self.connect_btn.setText("Disconnected")
            # self.connect_btn.setStyleSheet("""
            #     QPushButton {
            #         background-color: #484848;
            #         border: 1px solid #484848;
            #         border-radius: 4px;
            #         color: #ffffff;
            #         padding: 6px 16px;
            #         font-weight: bold;
            #         min-width: 100px;
            #     }
            #     
            # """)
            self.connect_btn.setStyleSheet("""
                QPushButton {
                    font-size: 18px;
                    font-weight: bold;
                    color: #ffffff;
                    padding: 4px 8px;
                    border: 2px solid #ffffff;
                    border-radius: 4px;
                }
                QPushButton:hover {
                    color: #7ede6f;
                    border: 2px solid #7ede6f;
                }
            """)
            # self.logo_label.setStyleSheet("""
            #     QLabel {
            #         font-size: 18px;
            #         font-weight: bold;
            #         color: #ffffff;
            #         padding: 4px 8px;
            #         border: 2px solid #ffffff;
            #         border-radius: 4px;
            #     }
            # """)
            self.refresh_btn.setEnabled(False)
    
    def _toggle_connection(self):
        """Toggle connection to the signal generator."""
        if self._connected:
            self._disconnect()
        else:
            self._connect()
    
    def _connect(self):
        """Establish connection to the signal generator."""
        if self.sg is not None and self._connected:
            # Already connected
            return
        
        init_with_gui = self.init_with_gui_cb.isChecked()
        
        if self.sg is not None:
            # Have an existing instance (passed to __init__)
            try:
                # Initialize hardware connection if not already done
                if not hasattr(self.sg, '_instr') or self.sg._instr is None:
                    from connectionConfig import sg_addr
                    self.sg.init(sg_addr)
                
                if init_with_gui:
                    # Push GUI values TO the signal generator
                    self._connected = True
                    self._address = getattr(self.sg, 'addr', 'Unknown')
                    self._connect_signals()
                    self._push_gui_to_hardware()
                    self.last_error_label.set_msg("✓ Connected to SG384 (initialized with GUI params)")
                    self.last_error_label.styleSheet()
                else:
                    # Read values FROM the signal generator to GUI
                    self._sync_from_hardware()
                    if self._connected:
                        self._connect_signals()
                        self.last_error_label.set_msg("✓ Connected to SG384 (synced from hardware)")
                        self.last_error_label.styleSheet()
                    
            except Exception as e:
                self.last_error_label.set_msg(f"✗ Connection failed: {e}", error=True)
                self._connected = False
        else:
            # Try to create a new instance
            try:
                
                if init_with_gui:
                    # Create SG without auto-init, then push GUI values
                    self.sg = SignalGenerator(auto_init_hardware=False)
                    from connectionConfig import sg_addr
                    self.sg.init(sg_addr)
                    self._connected = True
                    self._address = getattr(self.sg, 'addr', 'Unknown')
                    self._connect_signals()
                    self._push_gui_to_hardware()
                    self.last_error_label.set_msg("✓ Connected to SG384 (initialized with GUI params)")
                    self.last_error_label.styleSheet()
                else:
                    # Normal init - will auto-initialize and we sync from hardware
                    self.sg = SignalGenerator(auto_init_hardware=True)
                    self._sync_from_hardware()
                    self._connect_signals()
                    self.last_error_label.set_msg("✓ Connected to SG384 (synced from hardware)")
                    self.last_error_label.styleSheet()
                    
            except ImportError:
                self.last_error_label.set_msg("✗ Could not import SGcontrol module", error=True)
                self._connected = False
            except Exception as e:
                self.last_error_label.set_msg(f"✗ Connection failed: {e}", error=True)
                self._connected = False
                self.sg = None
        
        self._update_connect_button()
        self._update_status_bar()
    
    def _disconnect(self):
        """Disconnect from the signal generator."""
        if self.sg is not None:
            try:
                self.sg.uninit()
                self.last_error_label.set_msg("✓ Disconnected from SG384")
                self.last_error_label.styleSheet()
            except Exception as e:
                self.last_error_label.set_msg(f"⚠ Error during disconnect: {e}", error=True)
            finally:
                self.sg = None
        
        self._connected = False
        self._address = "Not connected"
        self._update_connect_button()
        self._update_status_bar()
    
    def _sync_from_hardware(self):
        """Sync GUI values from hardware."""
        if self.sg is None:
            self.last_error_label.set_msg("No signal generator instance connected", error=True)
            return
        
        try:
            # Query all parameters
            self.sg.query_all()
            
            # Check for any errors during query
            error_code = self._check_and_report_error("Sync from hardware")
            
            # Update connection info
            self._connected = True
            self._address = getattr(self.sg, 'addr', 'Unknown')
            
            # Update N-Type section
            self.ntype_output.setOutputState(bool(self.sg.status_ntype))
            self.ntype_output.setFrequency(self.sg.freq)
            self.ntype_output.setAmplitude(self.sg.amp_rf)
            self.ntype_output.setModulationState(bool(self.sg.mod_status))
            if self.sg.mod_status:
                self.ntype_output.setModulationType(self.sg.mod_type)
                self.ntype_output.setModulationFunction(self.sg.mod_func)
                self.ntype_output.setModulationDeviation(self.sg.mod_dev)
                self.ntype_output.setModulationRate(self.sg.mod_rate)

            # Update BNC section
            self.bnc_output.setOutputState(bool(self.sg.status_bnc))
            self.bnc_output.setAmplitude(self.sg.amp_bnc)
            
            self._update_status_bar()
            
            if error_code == 0:
                self._clear_last_error()
                self.last_error_label.set_msg("✓ Synced from hardware")
            else:
                self.last_error_label.set_msg(f"⚠ Synced with errors (code {error_code})", error=True)
            
        except Exception as e:
            self.last_error_label.set_msg(f"⚠ Sync failed: {str(e)}", error=True)
            # self._error_clear_timer.start(5000)
            self._connected = False
            self._update_status_bar()
    
    def _connect_signals(self):
        """Connect GUI signals to hardware control methods with error checking."""
        if self.sg is None:
            return
        
        # Prevent duplicate connections
        if self._signals_connected:
            return
        self._signals_connected = True
        
        # N-Type connections with error handling
        self.ntype_output.outputToggled.connect(
            lambda state: self._safe_execute(
                self.sg.enable_ntype, state, context="N-Type Output"
            )
        )
        self.ntype_output.frequencyChanged.connect(
            lambda freq: self._safe_execute(
                self.sg.set_freq, freq, context="Set Frequency"
            )
        )
        self.ntype_output.amplitudeChanged.connect(
            lambda amp: self._safe_execute(
                self.sg.set_amp_rf, amp, context="Set N-Type Amplitude"
            )
        )
        self.ntype_output.modulationToggled.connect(
            lambda state: self._safe_execute(
                self.sg.enable_modulation, state, context="Modulation Enable"
            )
        )
        self.ntype_output.modulationTypeChanged.connect(
            lambda type: self._safe_execute(
                self.sg.set_mod_type, type, context="Modulation Type"
            )
        )
        self.ntype_output.modulationFunctionChanged.connect(
            lambda func: self._safe_execute(
                self.sg.set_mod_func, func, context="Modulation Function"
            )
        )
        self.ntype_output.modulationDeviationChanged.connect(
            lambda dev: self._safe_execute(
                self.sg.set_mod_dev, dev, context="Modulation Deviation"
            )
        )
        self.ntype_output.modulationRateChanged.connect(
            lambda rate: self._safe_execute(
                self.sg.set_mod_rate, rate, context="Modulation Rate"
            )
        )
        
        # BNC connections with error handling
        self.bnc_output.outputToggled.connect(
            lambda state: self._safe_execute(
                self.sg.enable_bnc, state, context="BNC Output"
            )
        )
        self.bnc_output.amplitudeChanged.connect(
            lambda amp: self._safe_execute(
                self.sg.set_amp_bnc, amp, context="Set BNC Amplitude"
            )
        )
        
        # Display change connections (both sections share the same display)
        self.ntype_output.displayChangeRequested.connect(self._on_display_change)
        self.bnc_output.displayChangeRequested.connect(self._on_display_change)
    
    def _on_display_change(self, display_mode: int):
        """Handle display change request from any output section."""
        # Uncheck display buttons in the OTHER section to maintain global exclusivity
        sender = self.sender()
        if sender == self.ntype_output:
            for btn in self.bnc_output._display_buttons:
                btn.setActive(False)
        else:
            for btn in self.ntype_output._display_buttons:
                btn.setActive(False)
        
        # Send command to instrument
        self._safe_execute(
            self.sg.set_display, display_mode, context="Set Display"
        )

        self.last_error_label.set_msg(f"Display: {self.sg.display}")
    
    # =========================================================================
    # STATE MANAGEMENT
    # =========================================================================
    
    def _refresh_state_files(self):
        """Refresh the list of available state files."""
        self.state_combo.clear()
        
        try:
            if self.working_directory.exists():
                json_files = sorted(self.working_directory.glob("*.json"))
                for f in json_files:
                    self.state_combo.addItem(f.stem, str(f))
        except Exception as e:
            self.last_error_label.set_msg(f"Error refreshing state files: {e}", error=True)
    
    def _get_gui_state(self) -> dict:
        """Get current GUI state as a dictionary."""
        state = {
            'ntype': {
                'output_enabled': self.ntype_output.getOutputState(),
                'frequency': self.ntype_output.getFrequency(),
                'amplitude': self.ntype_output.getAmplitude(),
                'mod_enabled': self.ntype_output.getModulationState(),
                'mod_type': self.ntype_output.getModulationType(),
                'mod_func': self.ntype_output.getModulationFunction(),
                'mod_dev': self.ntype_output.getModulationDeviation(),
                'mod_rate': self.ntype_output.getModulationRate(),
            },
            'bnc': {
                'output_enabled': self.bnc_output.getOutputState(),
                'frequency': self.bnc_output.getFrequency(),
                'amplitude': self.bnc_output.getAmplitude(),
                'mod_enabled': self.bnc_output.getModulationState(),
                'mod_type': self.bnc_output.getModulationType(),
                'mod_func': self.bnc_output.getModulationFunction(),
                'mod_dev': self.ntype_output.getModulationDeviation(),
                'mod_rate': self.bnc_output.getModulationRate(),
            },
            'window': {
                'width': self.width(),
                'height': self.height(),
            },
            'options': {
                'init_with_gui': self.init_with_gui_cb.isChecked(),
            }
        }
        return state
    
    def _set_gui_state(self, state: dict):
        """Set GUI state from a dictionary."""
        try:
            # N-Type settings
            if 'ntype' in state:
                ntype = state['ntype']
                self.ntype_output.setOutputState(ntype.get('output_enabled', False))
                self.ntype_output.setFrequency(ntype.get('frequency', 2.87e9))
                self.ntype_output.setAmplitude(ntype.get('amplitude', 0.0))
                self.ntype_output.setModulationState(ntype.get('mod_enabled', False))
                self.ntype_output.setModulationType(ntype.get('mod_type'))
                self.ntype_output.setModulationFunction(ntype.get('mod_func', 'EXTERNAL'))
                self.ntype_output.setModulationDeviation(ntype.get('mod_dev', 100e3))
                self.ntype_output.setModulationRate(ntype.get('mod_rate', 1e3))
            
            # BNC settings
            if 'bnc' in state:
                bnc = state['bnc']
                self.bnc_output.setOutputState(bnc.get('output_enabled', False))
                self.bnc_output.setFrequency(bnc.get('frequency', 10e6))
                self.bnc_output.setAmplitude(bnc.get('amplitude', 0.0))
                self.bnc_output.setModulationState(bnc.get('mod_enabled', False))
                self.bnc_output.setModulationType(bnc.get('mod_type'))
                self.bnc_output.setModulationFunction(bnc.get('mod_func', 'EXTERNAL'))
                self.bnc_output.setModulationDeviation(bnc.get('mod_dev', 100e3))
                self.bnc_output.setModulationRate(bnc.get('mod_rate', 1e3))
            
            # Window settings
            if 'window' in state:
                w = state['window'].get('width', 550)
                h = state['window'].get('height', 500)
                self.resize(w, h)
            
            # Options
            if 'options' in state:
                self.init_with_gui_cb.setChecked(state['options'].get('init_with_gui', False))
                
        except Exception as e:
            self.last_error_label.set_msg(f"Error setting GUI state: {e}", error=True)
    
    def _save_state(self):
        """Save current GUI state to file."""
        try:
            filename = Path(self.state_combo.currentText().strip())
            
            if not filename or str(filename) == '.':
                self._show_error(-1, "No filename specified")
                return
            
            filepath = self.working_directory / filename.with_suffix(".json")
            
            state = self._get_gui_state()
            
            with open(filepath, 'w') as f:
                json.dump(state, f, indent=2)
            
            self.last_error_label.set_msg(f"✓ Saved: {filename.stem}")
            # self._error_clear_timer.start(3000)
            self.last_error_label.set_msg(f"✓ State saved to {filepath}")
            
            # Refresh file list and select saved file
            self._refresh_state_files()
            idx = self.state_combo.findText(filename.stem)
            if idx >= 0:
                self.state_combo.setCurrentIndex(idx)
                
        except Exception as e:
            self.last_error_label.set_msg(f"Error saving state: {e}", error=True)
            self._show_error(-1, f"Save failed: {str(e)}")
    
    def _load_state(self):
        """Load GUI state from file."""
        try:
            filename = Path(self.state_combo.currentText().strip())
            
            if not filename or str(filename) == '.':
                self._show_error(-1, "No filename specified")
                return
            
            filepath = self.working_directory / filename.with_suffix(".json")
            
            if not filepath.exists():
                self._show_error(-1, f"File not found: {filename}")
                return
            
            with open(filepath, 'r') as f:
                state = json.load(f)
            
            self._set_gui_state(state)
            
            self.last_error_label.set_msg(f"✓ Loaded: {filename.stem}")
            # self._error_clear_timer.start(3000)
            self.last_error_label.set_msg(f"✓ State loaded from {filepath}")
            
        except Exception as e:
            self.last_error_label.set_msg(f"Error loading state: {e}", error=True)
            self._show_error(-1, f"Load failed: {str(e)}")
    
    def _push_gui_to_hardware(self):
        """Push all GUI values to the signal generator."""
        if self.sg is None or not self._connected:
            return
        
        self.last_error_label.set_msg("Pushing GUI values to hardware...", error=True)
        
        # N-Type settings
        self._safe_execute(self.sg.enable_ntype, self.ntype_output.getOutputState(), context="N-Type Output")
        self._safe_execute(self.sg.set_freq, self.ntype_output.getFrequency(), context="Set Frequency")
        self._safe_execute(self.sg.set_amp_rf, self.ntype_output.getAmplitude(), context="Set N-Type Amplitude")
        
        # Modulation settings
        mod_type = self.ntype_output.getModulationType()
        if mod_type != "NONE":
            # type_map = {
            #     "NONE": "NONE", "AM": "AMPLITUDE", "FM": "FREQUENCY",
            #     "PM": "PHASE", "SWEEP": "SWEEP", "PULSE": "PULSE",
            #     "BLANK": "BLANK", "IQ": "IQ"
            # }
            self._safe_execute(self.sg.set_mod_type, ModulationType[mod_type.upper()].name, context="Modulation Type")
            self._safe_execute(self.sg.set_mod_func, self.ntype_output.getModulationFunction(), context="Modulation Function")
            self._safe_execute(self.sg.set_mod_dev, self.ntype_output.getModulationDeviation(), context="Modulation Deviation")
            self._safe_execute(self.sg.set_mod_rate, self.ntype_output.getModulationRate(), context="Modulation Rate")
        
        self._safe_execute(self.sg.enable_modulation, self.ntype_output.getModulationState(), context="Modulation Enable")
        
        # BNC settings
        self._safe_execute(self.sg.enable_bnc, self.bnc_output.getOutputState(), context="BNC Output")
        self._safe_execute(self.sg.set_amp_bnc, self.bnc_output.getAmplitude(), context="Set BNC Amplitude")
        
        self.last_error_label.set_msg("✓ GUI values pushed to hardware")
    
    def closeEvent(self, event):
        """Handle window close - auto-save last state."""
        try:
            # Save last state
            filepath = self.working_directory / DEFAULT_STATE_FILE
            state = self._get_gui_state()
            with open(filepath, 'w') as f:
                json.dump(state, f, indent=2)
            print(f"✓ Auto-saved state to {filepath}")
        except Exception as e:
            print(f"Warning: Could not auto-save state: {e}")
        
        event.accept()
    
    def setSignalGenerator(self, sg):
        """Set the signal generator instance and sync."""
        self.sg = sg
        self._sync_from_hardware()
        self._connect_signals()


# =============================================================================
# STANDALONE MODE
# =============================================================================

def main():
    """Run the GUI in standalone mode (without hardware)."""
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Create window (no hardware connection)
    window = SG384ControlPanel(sg_instance=None)
    window.show()
    
    # Set some demo values
    window.ntype_output.setFrequency(2.87e9)
    window.ntype_output.setAmplitude(-10)
    window.bnc_output.setFrequency(10e6)
    
    return app.exec()


def main_with_hardware():
    """Run the GUI with actual hardware connection."""
    # Must be called before creating QApplication
    myappid = 'aglab.sg'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Import and create signal generator
    try:
        from SGcontrol import SignalGenerator
        sg = SignalGenerator(auto_init_hardware=False)
    except Exception as e:
        print(f"Could not connect to hardware: {e}")
        sg = None
    
    window = SG384ControlPanel(sg_instance=sg, auto_connect_hardware=False)

    app_icon = QIcon(expt_dir + r"\gui\sg\gui_icon.png") # .ico is preferred for Windows
    window.setWindowIcon(app_icon)
    app.setWindowIcon(app_icon) # Sets it for the whole application

    window.show()
    return app.exec()

if __name__ == '__main__':
    sys.exit(main_with_hardware())
