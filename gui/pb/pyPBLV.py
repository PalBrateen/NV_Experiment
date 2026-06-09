# -*- coding: utf-8 -*-
"""
PulseBlaster SpinAPI Controller - Enhanced Version
Based on pyPBLV_qbuttongroup.py with improvements:
- Dark theme for light-sensitive experiments
- QSpinBox for # instructions and # channels
- Status bar for errors/status and hardware info
- Channel naming system with file persistence
- State file management via dropdown with working directory
- Filter checkbox to show only named channels
- Green highlighting for active instruction buttons

@author: brate
"""
# TODO: check if the input can be reflected in the status bar for the user..

import sys, os, logging, json, ctypes
from pathlib import Path
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QLineEdit,
    QVBoxLayout, QHBoxLayout, QWidget, QGridLayout,
    QLabel, QScrollArea, QComboBox, QButtonGroup, QSpinBox,
    QCheckBox, QStatusBar, QFrame, QSizePolicy, QGroupBox
)
from PySide6.QtCore import Qt, Slot, QTimer
from PySide6.QtGui import QColor, QPalette, QPainter, QFontMetrics, QIcon
from functools import reduce
# from dark_mode_and_extended_error import DARK_STYLESHEET, StatusBarLongLabel

# Get the absolute path to the directory containing the module
expt_dir = os.path.abspath(r'D:\Brateen\NV_Experiment')
# Add the directory to sys.path
sys.path.append(expt_dir)

# Try to import spinapi, fall back to simulation mode if not available
try:
    import spinapi as sp
    SPINAPI_AVAILABLE = True
except ImportError:
    SPINAPI_AVAILABLE = False
    print("Warning: spinapi not available, running in simulation mode")

# ============================================================================
# Configuration
# ============================================================================
CLOCK_FREQ = 500
SPINAPI_DLL_PATH = r'C:\SpinCore\SpinAPI\lib\spinapi64.dll'
STATEFILE_DIRECTORY = r'D:\Brateen\Saved_Data\SavedStates\PBStates'  # Default working directory for state files
CHANNEL_CONFIG_FILE = 'pb_channels.json'
DEFAULT_STATE_FILE = 'last_state.json'
MAX_CHANNELS = 21  # PulseBlaster has 21 channels (1-21)

# Default channel names (commonly used)
DEFAULT_CHANNEL_NAMES = {
    1: "Laser",
    2: "MW",
    5: "SampleCLK",
    8: "Camera",
}

DARK_STYLESHEET = """
QMainWindow, QWidget {
    background-color: #242424;
    color: #ffffff;
    font-family: 'Segoe UI', Arial, sans-serif;
    font-size: 10pt;
}
QLabel {
    color: #ffffff;
    padding: 2px;
}
QLineEdit, QSpinBox, QComboBox {
    background-color: #484848;
    color: #ffffff;
    border: 1px solid #3f3f46;
    border-radius: 3px;
    padding: 4px;
    selection-background-color: #264f78;
}
QLineEdit:focus, QSpinBox:focus, QComboBox:focus {
    border: 1px solid #009de0;
}
QLineEdit:disabled, QSpinBox:disabled {
    background-color: #252526;
    color: #6d6d6d;
}
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
QPushButton:pressed {
    background-color: #094771;
}
QPushButton:disabled {
    background-color: #3f3f46;
    color: #6d6d6d;
}
/* Checkable channel buttons - unchecked state */
QPushButton[checkable="true"] {
    background-color: #3c3c3c;
    color: #ffffff;
    border: 1px solid #555555;
    padding: 4px 8px;
    min-width: 30px;
}
/* Checkable channel buttons - checked (active) state - GREEN */
QPushButton[checkable="true"]:checked {
    background-color: #2e7d32;
    color: white;
    border: 1px solid #4caf50;
}
QPushButton[checkable="true"]:hover {
    background-color: #4a4a4a;
}
QPushButton[checkable="true"]:checked:hover {
    background-color: #388e3c;
}
/* The Main Box */
QSpinBox {
    background-color: #2d2d2d;
    color: #ffffff;
    border: 1px solid #555555;
    border-radius: 4px;
    padding-right: 5px; /* Leave space for buttons */
    selection-background-color: #444444;
}
/* The Buttons Container */
QSpinBox::up-button, QSpinBox::down-button {
    background-color: #3d3d3d;
    border-left: 1px solid #555555;
    width: 20px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {
    background-color: #4d4d4d;
}
/* The Arrows (Triangle Hack) */
QSpinBox::up-arrow {
    width: 0px; height: 0px;
    border-left: 4px solid #3d3d3d;
    border-right: 4px solid #3d3d3d;
    border-bottom: 5px solid #ffffff;
}
QSpinBox::down-arrow {
    width: 0px; height: 0px;
    border-left: 4px solid #3d3d3d;
    border-right: 4px solid #3d3d3d;
    border-top: 5px solid #ffffff;
}
/* Disabled State - Critical for Logic */
QSpinBox:disabled {
    background-color: #1e1e1e;
    color: #777777;
}
QSpinBox::up-arrow:disabled, QSpinBox::down-arrow:disabled {
    border-bottom-color: #555555;
    border-top-color: #555555;
}
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
    background-color: #2d2d30; color: #ffffff;
    selection-background-color: #094771;
    border: 1px solid #3f3f46;
}
QScrollArea {
    border: 3px solid #3f3f46;
    border-radius: 5px;
    background-color: #252526;
}
QScrollBar:vertical {
    background-color: #1e1e1e;
    width: 14px;
    margin: 0;
}
QScrollBar::handle:vertical {
    background-color: #5a5a5a;
    min-height: 30px;
    border-radius: 7px;
    margin: 2px;
}
QScrollBar::handle:vertical:hover {
    background-color: #787878;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}
QScrollBar:horizontal {
    background-color: #1e1e1e;
    height: 14px;
    margin: 0;
}
QScrollBar::handle:horizontal {
    background-color: #5a5a5a;
    min-width: 30px;
    border-radius: 7px;
    margin: 2px;
}
QScrollBar::handle:horizontal:hover {
    background-color: #787878;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}
QGroupBox {
    border: 3px solid #3f3f46;
    border-radius: 5px;
    margin-top: 8px;
    padding-top: 8px;
    color: #ffffff;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 5px;
}
QCheckBox {
    color: #ffffff;
    spacing: 8px;
}
QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 1px solid #555555;
    border-radius: 3px;
    background-color: #2d2d30;
}
QCheckBox::indicator:checked {
    background-color: #0e639c;
    border-color: #859199;
}
QStatusBar {
    background-color: #303030;
    color: white;
    border-top: 2px solid #878787;
}
QStatusBar QLabel {
    color: white;
    padding: 2px 8px;
}
QFrame#separator {
    background-color: #3f3f46;
}
"""

class StatusBarLongLabel(QLabel):
    def __init__(self, parent=None):
        super().__init__(parent)
        # Set a minimum width so it doesn't disappear, 
        # but no maximum so it can take up available space.
        self.setMinimumWidth(100)
        self.setStyleSheet("color: #ff6b6b;") # Light red for dark mode errors

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

    def set_error(self, message):
        """Update text and set the full message as a tooltip."""
        self.setText(message)
        self.setToolTip(message) # The OS handles the popup on hover

# ============================================================================
# Logging Setup
# ============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def setup_logging():
    logger = logging.getLogger('pyPBLV')
    logger.setLevel(logging.INFO)
    
    if logger.hasHandlers():
        logger.handlers.clear()
    
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)s: %(message)s')
    console.setFormatter(formatter)
    logger.addHandler(console)
    
    return logger


# ============================================================================
# Simulation classes for spinapi when not available
# ============================================================================
if not SPINAPI_AVAILABLE:
    class Inst:
        CONTINUE = 0
        STOP = 1
        LOOP = 2
        END_LOOP = 3
        JSR = 4
        RTS = 5
        BRANCH = 6
        LONG_DELAY = 7
        WAIT = 8
    
    ns = 1
    us = 1000
    ms = 1000000
    s = 1000000000
    PULSE_PROGRAM = 0
    
    def pb_count_boards(): return 1
    def pb_get_error(): return "Simulation mode"
    def pb_close(): return 0
    def pb_select_board(n): return 0
    def pb_init(): return 0
    def pb_get_version(): return "Simulation"
    def pb_core_clock(freq): return 0
    def pb_start_programming(mode): return 0
    def pb_stop_programming(): return 0
    def pb_inst_pbonly(flags, inst, data, length): return 0
    def pb_start(): return 0
    def pb_stop(): return 0
else:
    from spinapi import (
        Inst, ns, us, ms, s, PULSE_PROGRAM,
        pb_count_boards, pb_get_error, pb_close, pb_select_board,
        pb_init, pb_get_version, pb_core_clock, pb_start_programming,
        pb_stop_programming, pb_inst_pbonly, pb_start, pb_stop
    )


# ============================================================================
# Main GUI Class
# ============================================================================
class SpinAPIGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        # self.setStyleSheet("QMainWindow { border: 50px solid #ffffff; }")
        self.spinapi_version = ''
        self.logger = logging.getLogger('pyPBLV')
        
        # Working directory for state files
        self.statefile_directory = Path(STATEFILE_DIRECTORY)
        self.statefile_directory.mkdir(parents=True, exist_ok=True)
        
        # Channel configuration
        self.channel_names = {}
        self.load_channel_names()
        
        # GUI state
        self.n_channels = 8  # Default visible channels
        self.button_states = {}
        self.button_groups = {}
        self.column_data = {}
        self.channel_name_widgets = {}  # Store QLineEdit widgets for channel names
        self.d_struct = []
        
        self.window_width: int = 500
        self.window_height: int = 700
        # Program state
        self.prg_order = 0

        # self.status_bar.setStyleSheet("QStatusBar{background-color: yellow; color: black;}")
        self.initUI()
        self.load_state()

    def load_channel_names(self):
        """Load channel names from config file in working directory."""
        config_path = self.statefile_directory / CHANNEL_CONFIG_FILE
        try:
            if config_path.exists():
                with open(config_path, 'r') as f:
                    loaded = json.load(f)
                    # Convert string keys to int
                    self.channel_names = {int(k): v for k, v in loaded.items()}
                self.logger.info(f'Loaded channel names from {config_path}')
            else:
                self.channel_names = DEFAULT_CHANNEL_NAMES.copy()
                self.save_channel_names()
        except Exception as e:
            self.logger.error(f'Error loading channel names: {e}')
            self.channel_names = DEFAULT_CHANNEL_NAMES.copy()

    def save_channel_names(self):
        """Save channel names to config file in working directory."""
        config_path = self.statefile_directory / CHANNEL_CONFIG_FILE
        try:
            # Convert int keys to string for JSON
            to_save = {str(k): v for k, v in self.channel_names.items() if v.strip()}
            with open(config_path, 'w') as f:
                json.dump(to_save, f, indent=2)
            self.logger.info(f'Saved channel names to {config_path}')
        except Exception as e:
            self.logger.error(f'Error saving channel names: {e}')

    def initUI(self):
        self.setWindowTitle('SpinAPI Controller')
        self.setGeometry(100, 100, self.window_width, self.window_height)
        self.setStyleSheet(DARK_STYLESHEET)
        
        # Create central widget and main layout
        central_widget = QWidget()
        # QWidget combining
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(10)
        main_layout.setContentsMargins(10, 10, 10, 10)
        
        # Top section with controls
        top_layout = QHBoxLayout()
        
        # Left side: Control buttons
        self._create_control_buttons(top_layout)
        
        # Right side: State file management and settings
        self._create_settings_panel(top_layout)
        
        main_layout.addLayout(top_layout)
        
        # Separator
        separator = QFrame()
        separator.setObjectName("separator")
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFixedHeight(2)
        main_layout.addWidget(separator)
        
        # Instruction grid with scroll area
        self._create_instruction_grid(main_layout)
        
        # QWidget combining
        # Status bar
        self._create_status_bar()
        
        # Initialize grid
        self.update_grid()
        
        # Initialize board connection
        QTimer.singleShot(100, self.board_num_Callback)
        
        self.logger.info('UI Initialized')

    def _create_control_buttons(self, parent_layout):
        """Create the main control buttons (Load/Start/Stop/Change Board)."""
        button_group = QGroupBox("Board Control")
        button_layout = QVBoxLayout(button_group)
        
        # Row 1: Load and Start/Restart
        # row1 = QHBoxLayout()
        self.load_board_btn = QPushButton('🔁')     # LOAD
        self.load_board_btn.setFixedSize(100, 30)
        self.load_board_btn.setStyleSheet("QPushButton { font-size: 20px; }")
        self.load_board_btn.clicked.connect(self.load_pushbutton_Callback)
        button_layout.addWidget(self.load_board_btn)
        
        self.start_restart_btn = QPushButton('▶')   # START
        self.start_restart_btn.setFixedSize(100, 30)
        self.start_restart_btn.setStyleSheet("QPushButton { font-size: 30px; }")
        self.start_restart_btn.clicked.connect(self.start_pushbutton_Callback)
        button_layout.addWidget(self.start_restart_btn)
        # button_layout.addLayout(row1)
        
        # Row 2: Change Board and Stop
        # row2 = QHBoxLayout()
        # self.change_board_btn = QPushButton('CHANGE BOARD')
        # self.change_board_btn.setFixedSize(100, 30)
        # self.change_board_btn.clicked.connect(self.change_board)
        # row2.addWidget(self.change_board_btn)
        # |🔄🔃🔀
        self.stop_btn = QPushButton('🛑') # STOP
        self.stop_btn.setFixedSize(100, 30)
        self.stop_btn.setStyleSheet("QPushButton { font-size: 30px; }")
        self.stop_btn.setStyleSheet("""
            QPushButton {
                background-color: #c62828;
            }
            QPushButton:hover {
                background-color: #d32f2f;
            }
            QPushButton:pressed {
                background-color: #b71c1c;
            }
        """)
        self.stop_btn.clicked.connect(self.stop_pushbutton_Callback)
        button_layout.addWidget(self.stop_btn)
        # button_layout.addLayout(row2)
        
        parent_layout.addWidget(button_group)

    def _create_settings_panel(self, parent_layout):
        """Create the settings panel with state management and display options."""
        settings_group = QGroupBox("Settings")
        settings_layout = QVBoxLayout(settings_group)
        
        # State file management section (LabOne style - editable combo)
        state_layout = QVBoxLayout()
        
        # File Name row with editable combo
        filename_row = QHBoxLayout()
        filename_row.addWidget(QLabel("File Name:"))
        
        self.state_combo = QComboBox()
        self.state_combo.setEditable(True)  # Allow typing new names
        self.state_combo.setMinimumWidth(180)
        self.state_combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)  # Don't auto-add typed text
        self.state_combo.lineEdit().setPlaceholderText("Enter name or select...")
        self.refresh_state_files()
        filename_row.addWidget(self.state_combo)
        
        state_layout.addLayout(filename_row)
        
        # Load and Save buttons row
        btn_row = QHBoxLayout()
        btn_row.addStretch(20)
        
        self.load_state_btn = QPushButton('Load')
        self.load_state_btn.setStyleSheet("""
            QPushButton {
            background-color: #484848;
            color: #ffffff;
            }
            QPushButton:hover {
                background-color: #009de0;
            }
        """)
        self.load_state_btn.setFixedSize(60, 25)
        self.load_state_btn.clicked.connect(self.load_state)
        btn_row.addWidget(self.load_state_btn)

        self.save_state_btn = QPushButton('Save')
        self.save_state_btn.setStyleSheet("""
            QPushButton {
            background-color: #d3dbde;
            color: #484848;
            }
            QPushButton:hover {
                background-color: #009de0;
                color: #ffffff;
            }
        """)
        self.save_state_btn.setFixedSize(60, 25)
        self.save_state_btn.clicked.connect(self.save_state)
        btn_row.addWidget(self.save_state_btn)

        # Refresh button inline with combo
        self.refresh_btn = QPushButton('↻')
        self.refresh_btn.setStyleSheet("""
            QPushButton {
            background-color: #484848;
            color: #ffffff;
            }
            QPushButton:hover {
                background-color: #009de0;
                color: #ffffff;
            }
        """)
        self.refresh_btn.setFixedSize(25, 25)
        self.refresh_btn.setToolTip("Refresh file list")
        self.refresh_btn.clicked.connect(self.refresh_state_files)
        btn_row.addWidget(self.refresh_btn)
        
        state_layout.addLayout(btn_row)
        settings_layout.addLayout(state_layout)
        
        # Grid settings
        grid_settings = QHBoxLayout()
        
        # Number of instructions
        inst_col = QVBoxLayout()
        inst_col.addWidget(QLabel("# Instructions:"))
        self.num_instructions = QSpinBox()
        self.num_instructions.setFixedSize(70, 25)
        self.num_instructions.setRange(1, 50)
        self.num_instructions.setValue(6)
        self.num_instructions.valueChanged.connect(self.update_grid)
        inst_col.addWidget(self.num_instructions)
        grid_settings.addLayout(inst_col)
        
        # Number of channels
        chan_col = QVBoxLayout()
        chan_col.addWidget(QLabel("# Channels:"))
        self.num_channels_spin = QSpinBox()
        self.num_channels_spin.setFixedSize(70, 25)
        self.num_channels_spin.setRange(0, MAX_CHANNELS-1)
        self.num_channels_spin.setValue(8)
        self.num_channels_spin.valueChanged.connect(self.on_channels_changed)
        chan_col.addWidget(self.num_channels_spin)
        grid_settings.addLayout(chan_col)

        grid_settings.addStretch()
        settings_layout.addLayout(grid_settings)
        
        # Filter checkbox
        filter_layout = QHBoxLayout()
        self.filter_named_checkbox = QCheckBox("Show only named channels")
        self.filter_named_checkbox.stateChanged.connect(self.update_grid)
        filter_layout.addWidget(self.filter_named_checkbox)
        filter_layout.addStretch()
        settings_layout.addLayout(filter_layout)
        
        parent_layout.addWidget(settings_group)
        parent_layout.addStretch()

    def _create_instruction_grid(self, parent_layout):
        """Create the scrollable instruction grid."""
        self.grid_layout = QGridLayout()
        self.grid_layout.setSpacing(4)
        self.grid_layout.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)  # Align to top-left
        
        self.grid_widget = QWidget()
        self.grid_widget.setLayout(self.grid_layout)
        
        # Wrap in another widget to ensure top alignment in scroll area
        scroll_content = QWidget()
        scroll_content_layout = QVBoxLayout(scroll_content)
        scroll_content_layout.setContentsMargins(0, 0, 0, 0)
        scroll_content_layout.addWidget(self.grid_widget, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        scroll_content_layout.addStretch(1)  # Push content to top
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(scroll_content)
        scroll_area.setMinimumHeight(400)
        
        parent_layout.addWidget(scroll_area)

    def _create_status_bar(self):
        """Create the status bar with hardware info and error display."""
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        
        # DLL path label
        self.dll_label = StatusBarLongLabel()
        self.dll_label.setStyleSheet("color: #aaaaaa;")
        self.status_bar.addWidget(self.dll_label,1)
        self.dll_label.set_error(f"DLL: {SPINAPI_DLL_PATH}")
        
        # # Separator
        # sep = QLabel("|")
        # sep.setStyleSheet("color: #666666;")
        # self.status_bar.addWidget(sep)
        
        # Clock frequency label
        self.clock_label = QLabel(f"Clock: {CLOCK_FREQ} MHz")
        self.status_bar.addWidget(self.clock_label)
        
        # # Separator
        # sep2 = QLabel("|")
        # sep2.setStyleSheet("color: #666666;")
        # self.status_bar.addWidget(sep2)
        
        # Board number label
        self.board_label = QLabel("Board: 0")
        self.status_bar.addWidget(self.board_label)
        
        # # Stretch to push status to right
        # self.status_bar.addWidget(QWidget(), 1)
        
        # Status message (right side)
        self.status_message = StatusBarLongLabel() #QLabel("OK")
        self.status_message.setStyleSheet("color: #4caf50; font-weight: bold;")
        self.status_message.setStyleSheet("padding-right: 20px;")
        # self.status_bar.addPermanentWidget(self.status_message)
        self.status_bar.addWidget(self.status_message, 1)

    def set_status(self, message, is_error=False):
        """Set status bar message with appropriate color."""
        self.status_message.set_error(message)
        if is_error:
            self.status_message.setStyleSheet("color: #f44336; font-weight: bold;")
        else:
            self.status_message.setStyleSheet("color: #4caf50; font-weight: bold;")

    def refresh_state_files(self):
        """Refresh the list of state files (.json) in the state save directory."""
        # Save current text if user typed something
        current_text = self.state_combo.currentText()
        
        self.state_combo.clear()
        
        try:
            self.state_combo.addItem(Path(DEFAULT_STATE_FILE).stem, DEFAULT_STATE_FILE)
            if self.statefile_directory.exists():
                json_files = sorted(self.statefile_directory.glob("*.json"))
                json_files.remove(self.statefile_directory / DEFAULT_STATE_FILE)
                for f in json_files:
                    if f.name != CHANNEL_CONFIG_FILE:
                        self.state_combo.addItem(f.stem, str(f))
        except Exception as e:
            self.logger.error(f'Error refreshing state files: {e}')
        
        # Restore typed text or set placeholder
        if current_text and current_text not in [self.state_combo.itemText(i) for i in range(self.state_combo.count())]:
            self.state_combo.setCurrentText(current_text)
        elif self.state_combo.count() > 0:
            self.state_combo.setCurrentIndex(0)

    def on_channels_changed(self, value):
        """Handle change in number of channels."""
        self.n_channels = value
        self.update_grid()

    def _safe_update_grid(self):
        """Wrapper for update_grid with exception handling."""
        try:
            self.update_grid()
        except Exception as e:
            self.logger.error(f"Grid update error: {e}")

    def get_visible_channels(self):
        """Get list of channel indices to display based on filter settings."""
        if self.filter_named_checkbox.isChecked():
            # Return only channels that have names
            named = [ch for ch in range(0, MAX_CHANNELS) if ch in self.channel_names and self.channel_names[ch].strip()]
            return named if named else list(range(0, self.n_channels + 1))
        else:
            return list(range(0, self.n_channels + 1))

    @Slot()
    def update_grid(self):
        """Update the instruction grid based on current settings."""
        # Prevent re-entrant calls during rapid changes
        if hasattr(self, '_updating_grid') and self._updating_grid:
            return
        self._updating_grid = True
        
        try:
            # Block signals from spinboxes during update to prevent cascade
            self.num_instructions.blockSignals(True)
            self.num_channels_spin.blockSignals(True)
            self.filter_named_checkbox.blockSignals(True)
            
            # Clear existing grid - disconnect signals first, then delete
            for i in reversed(range(self.grid_layout.count())):
                item = self.grid_layout.itemAt(i)
                if item:
                    widget = item.widget()
                    if widget:
                        # Disconnect all signals to prevent callbacks during destruction
                        try:
                            widget.blockSignals(True)
                        except:
                            pass
                        widget.setParent(None)
                        widget.deleteLater()
            
            self.button_groups.clear()
            self.channel_name_widgets.clear()
            
            # Process pending deletions before creating new widgets
            QApplication.processEvents()
            
            try:
                num_cols = self.num_instructions.value()
            except:
                num_cols = 6
            
            visible_channels = self.get_visible_channels()
            
            # Initialize column data for new columns
            for col in range(num_cols):
                if col not in self.column_data:
                    self.column_data[col] = {
                        'time': '100m',
                        'buttons': set(),
                        'opcode': 'CONTINUE',
                        'instruct_data': '0'
                    }
            
            # Remove data for columns beyond current count
            cols_to_remove = [c for c in list(self.column_data.keys()) if c >= num_cols]
            for c in cols_to_remove:
                del self.column_data[c]
            
            # Calculate total rows for separators
            total_rows = len(visible_channels) + 4  # header + time + channels + opcode + inst_data
            
            # Helper function to get actual grid column (accounting for separators)
            # Layout: [Channel names col 0] [sep] [Inst1 col 2] [sep col 3] [Inst2 col 4] ...
            def get_grid_col(inst_col):
                return 2 + inst_col * 2  # Each instruction takes 2 columns (content + separator)
            
            # Row 0: Column headers (Name column + instruction columns)
            channel_header = QLabel("Channel")
            channel_header.setStyleSheet("font-weight: bold;")
            self.grid_layout.addWidget(channel_header, 0, 0)
            
            for col in range(num_cols):
                header = QLabel(f"Inst {col + 1}")
                header.setAlignment(Qt.AlignmentFlag.AlignCenter)
                header.setStyleSheet("font-weight: bold;")
                self.grid_layout.addWidget(header, 0, get_grid_col(col))
            
            # Row 1: Time inputs with [s] label inline
            self.grid_layout.addWidget(QLabel(""), 1, 0)  # Empty for name column
            for col in range(num_cols):
                time_widget = QWidget()
                time_layout = QHBoxLayout(time_widget)
                time_layout.setContentsMargins(0, 0, 0, 0)
                time_layout.setSpacing(2)
                
                time_input = QLineEdit(self.column_data[col]['time'])
                time_input.setFixedWidth(60)
                time_input.setProperty("column", col)
                time_input.textChanged.connect(self._on_time_changed)
                time_layout.addWidget(time_input)
                
                unit_label = QLabel("[s]")
                unit_label.setStyleSheet("color: #888888; font-size: 9pt;")
                time_layout.addWidget(unit_label)
                time_layout.addStretch()
                
                self.grid_layout.addWidget(time_widget, 1, get_grid_col(col))
            
            # Channel rows with name inputs
            for row_idx, channel in enumerate(visible_channels):
                grid_row = row_idx + 2
                
                # Channel name input (column 0)
                name_widget = QWidget()
                name_layout = QHBoxLayout(name_widget)
                name_layout.setContentsMargins(0, 0, 0, 0)
                name_layout.setSpacing(4)
                
                # Channel number label
                ch_label = QLabel(f"{channel}:")
                ch_label.setFixedWidth(25)
                ch_label.setStyleSheet("color: #888888;")
                name_layout.addWidget(ch_label)
                
                # Editable name field
                name_input = QLineEdit(self.channel_names.get(channel, ""))
                name_input.setFixedWidth(80)
                name_input.setPlaceholderText("name")
                name_input.setProperty("channel", channel)
                name_input.textChanged.connect(self._on_channel_name_changed)
                name_layout.addWidget(name_input)
                
                self.channel_name_widgets[channel] = name_input
                self.grid_layout.addWidget(name_widget, grid_row, 0)
                
                # Instruction buttons for this channel
                for col in range(num_cols):
                    btn = QPushButton(str(channel))
                    btn.setFixedSize(50, 28)
                    btn.setCheckable(True)
                    btn.setProperty("checkable", True)  # For stylesheet
                    
                    # Convert channel to 0-based for internal storage
                    channel_idx = channel - 1
                    btn.setProperty("channel_idx", channel_idx)
                    btn.setProperty("column", col)
                    btn.setChecked(channel_idx in self.column_data[col].get('buttons', set()))
                    btn.clicked.connect(self._on_button_clicked)
                    
                    self.grid_layout.addWidget(btn, grid_row, get_grid_col(col))
            
            # Opcode row
            opcode_row = len(visible_channels) + 2
            opcode_label = QLabel("Opcode")
            opcode_label.setStyleSheet("color: #888888;")
            self.grid_layout.addWidget(opcode_label, opcode_row, 0)
            for col in range(num_cols):
                opcode_combo = QComboBox()
                opcode_combo.addItems(['CONTINUE', 'STOP', 'LOOP', 'END_LOOP',
                                       'JSR', 'RTS', 'BRANCH', 'LONG_DELAY', 'WAIT'])
                opcode_combo.setCurrentText(self.column_data[col].get('opcode', 'CONTINUE'))
                opcode_combo.setProperty("column", col)
                opcode_combo.currentTextChanged.connect(self._on_opcode_changed)
                opcode_combo.setFixedWidth(110)
                self.grid_layout.addWidget(opcode_combo, opcode_row, get_grid_col(col))
            
            # Instruction data row
            inst_data_row = opcode_row + 1
            inst_data_label = QLabel("Inst Data")
            inst_data_label.setStyleSheet("color: #888888;")
            self.grid_layout.addWidget(inst_data_label, inst_data_row, 0)
            for col in range(num_cols):
                inst_spin = QSpinBox()
                inst_spin.setRange(0, 1000000)
                inst_spin.setValue(int(self.column_data[col].get('instruct_data', 0)))
                inst_spin.setProperty("column", col)
                inst_spin.valueChanged.connect(self._on_inst_data_changed)
                inst_spin.setFixedWidth(80)
                self.grid_layout.addWidget(inst_spin, inst_data_row, get_grid_col(col))
            
            # Add vertical separators between instruction columns
            # Separator after channel names column (column 1)
            for row in range(total_rows):
                sep = QFrame()
                sep.setFrameShape(QFrame.Shape.VLine)
                sep.setStyleSheet("background-color: #3f3f46;")
                sep.setFixedWidth(2)
                self.grid_layout.addWidget(sep, row, 1)
            
            # Separators between instruction columns
            for col in range(num_cols - 1):  # No separator after last column
                sep_col = get_grid_col(col) + 1  # Separator column after each instruction
                for row in range(total_rows):
                    sep = QFrame()
                    sep.setFrameShape(QFrame.Shape.VLine)
                    sep.setStyleSheet("background-color: #3f3f46;")
                    sep.setFixedWidth(2)
                    self.grid_layout.addWidget(sep, row, sep_col)
            
            self.grid_widget.adjustSize()
            
        finally:
            # Always unblock signals and clear flag
            self.num_instructions.blockSignals(False)
            self.num_channels_spin.blockSignals(False)
            self.filter_named_checkbox.blockSignals(False)
            self._updating_grid = False
    
    # Slot methods that use sender() to get widget properties instead of lambdas
    @Slot(str)
    def _on_time_changed(self, text):
        """Handle time input change."""
        sender = self.sender()
        if sender:
            col = sender.property("column")
            if col is not None and col in self.column_data:
                self.update_column_data(col, 'time', text)
    
    @Slot(str)
    def _on_channel_name_changed(self, text):
        """Handle channel name change."""
        sender = self.sender()
        if sender:
            channel = sender.property("channel")
            if channel is not None:
                self.update_channel_name(channel, text)
    
    @Slot()
    def _on_button_clicked(self):
        """Handle instruction button click."""
        sender = self.sender()
        if sender:
            channel_idx = sender.property("channel_idx")
            col = sender.property("column")
            checked = sender.isChecked()
            if channel_idx is not None and col is not None and col in self.column_data:
                self.update_button_state(channel_idx, col, checked)
    
    @Slot(str)
    def _on_opcode_changed(self, text):
        """Handle opcode combo change."""
        sender = self.sender()
        if sender:
            col = sender.property("column")
            if col is not None and col in self.column_data:
                self.update_column_data(col, 'opcode', text)
    
    @Slot(int)
    def _on_inst_data_changed(self, value):
        """Handle instruction data spinbox change."""
        sender = self.sender()
        if sender:
            col = sender.property("column")
            if col is not None and col in self.column_data:
                self.update_column_data(col, 'instruct_data', str(value))

    def update_channel_name(self, channel, name):
        """Update channel name and mark for saving."""
        if name.strip():
            self.channel_names[channel] = name.strip()
        elif channel in self.channel_names:
            del self.channel_names[channel]

    def update_button_state(self, channel_idx, col, checked):
        """Update button state in column data."""
        if checked:
            self.column_data[col]['buttons'].add(channel_idx)
        else:
            self.column_data[col]['buttons'].discard(channel_idx)
        self.logger.debug(f"Column {col} buttons: {self.column_data[col]['buttons']}")

    def update_column_data(self, col, key, value):
        """Update column data for a specific key."""
        self.column_data[col][key] = value
        self.logger.debug(f"Column {col} {key}: {value}")

    def save_state(self):
        """Save current state to a JSON file."""
        try:
            # Get filename from combo (typed or selected)
            filename = Path(self.state_combo.currentText().strip())# if filename is None else Path(filename)
            
            if not filename:
                # Create default filename if empty
                from datetime import datetime
                filename = Path(f"pb_state_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            
            # Remove .json extension if user typed it
            filepath = self.statefile_directory / filename.with_suffix('.json')
            
            state_data = {
                'window_width': self.width(),
                'window_height': self.height(),
                'num_instructions': self.num_instructions.value(),
                'num_channels': self.num_channels_spin.value(),
                'clock_freq': CLOCK_FREQ,
                'columns': {}
            }
            
            for col, data in self.column_data.items():
                state_data['columns'][str(col)] = {
                    'time': data['time'],
                    'buttons': list(data['buttons']),
                    'opcode': data['opcode'],
                    'instruct_data': data['instruct_data']
                }
            
            with open(filepath, 'w') as f:
                json.dump(state_data, f, indent=2)
            
            self.set_status(f"Saved: {filename}")
            self.logger.info(f'State saved to {filepath}')
            self.refresh_state_files()
            
            # Select the saved file in dropdown
            idx = self.state_combo.findText(filename.stem)
            if idx >= 0:
                self.state_combo.setCurrentIndex(idx)
            
        except Exception as e:
            self.logger.error(f'Error saving state: {e}')
            self.set_status(f"Save failed: {str(e)}", is_error=True)

    def load_state(self):
        """Load state from selected or typed JSON file."""
        try:
            filename = Path(self.state_combo.currentText().strip())# if filename is None else Path(filename)
            
            if not filename:
                self.set_status("No file specified", is_error=True)
                return
            
            # # Try to get filepath from combo data first (for selected items)
            # filepath = Path(self.state_combo.currentData())
            
            # If no data (user typed), construct path
            # if not filepath:
            print(f"filename: {filename}")
            filepath = self.statefile_directory / filename.with_suffix(".json")
            # else:
            #     filepath = Path(filepath)
            
            if not filepath.exists():
                self.set_status(f"File not found: {filename}", is_error=True)
                return
            
            with open(filepath, 'r') as f:
                state_data = json.load(f)
            
            # Set window
            self.window_width = state_data.get('window_width',500)
            self.window_height = state_data.get('window_height',700)
            self.resize(self.window_width, self.window_height)

            # Set settings
            self.num_instructions.setValue(state_data.get('num_instructions', 6))
            self.num_channels_spin.setValue(state_data.get('num_channels', 8))
            
            # Clear and load column data
            self.column_data.clear()
            for col_str, data in state_data.get('columns', {}).items():
                col = int(col_str)
                self.column_data[col] = {
                    'time': data['time'],
                    'buttons': set(data['buttons']),
                    'opcode': data['opcode'],
                    'instruct_data': data['instruct_data']
                }
            
            self.update_grid()
            self.set_status(f"Loaded: {filename}")
            self.logger.info(f'State loaded from {filepath}')
            
        except Exception as e:
            self.logger.error(f'Error loading state: {e}')
            self.set_status(f"Load failed: {str(e)}", is_error=True)

    def generate_d_struct(self):
        """Generate the instruction data structure for PulseBlaster programming."""
        d_mult = {'n': ns, 'u': us, 'm': ms, '': s}
        
        self.d_struct = []
        
        def calc_flag_val(flag_list):
            return sum(2**(bit+1) for bit in flag_list) if flag_list else 0
        
        for col, data in self.column_data.items():
            d = {'inst_num': col}
            
            time_str = data['time']
            try:
                time_unit = int(time_str[-1])
                time_unit = ''
                time_val = float(time_str)
            except ValueError:
                time_unit = time_str[-1]
                time_val = float(time_str[:-1])
            
            d['dur'] = time_val
            d['s_mult'] = time_unit
            d['d_mult'] = d_mult.get(time_unit, s)
            # print(f"Buttons = {list(data['buttons'])}")
            d['flags'] = calc_flag_val(list(data['buttons']))
            
            try:
                if SPINAPI_AVAILABLE:
                    d['inst'] = Inst[data['opcode']]#.value
                else:
                    d['inst'] = getattr(Inst, data['opcode'])
            except:
                d['inst'] = 0
            
            d['instruct_data'] = int(data['instruct_data'])
            
            self.d_struct.append(d)
            self.logger.debug(f"Instruction {col}: {d}")

    def board_num_Callback(self):
        """Handle board number callback - initialize board connection."""
        try:
            count = pb_count_boards()
            self.logger.debug(f'Board count: {count}')
            
            if count < 0:
                self.set_status(f"Board error: {pb_get_error()}", is_error=True)
                self.prg_order = -1
            elif count == 0:
                self.set_status("NO board", is_error=True)
                self.prg_order = -1
            else:
                pb_close()
                if pb_select_board(0) < 0:
                    self.set_status(f"Select board failed: {pb_get_error()}", is_error=True)
                    return
                pb_init()
                self.prg_order = 0
                self.spinapi_version = pb_get_version()
                self.set_status(f"v{self.spinapi_version} Ready")
                self.board_label.setText(f"Board: 0")
                pb_close()
                
        except Exception as e:
            self.logger.error(f'Board init error: {e}')
            self.set_status(f"Init error: {str(e)}", is_error=True)
            self.prg_order = -1

    # def clock_freq_callback(self):
    #     """Handle clock frequency change."""
    #     global CLOCK_FREQ
    #     try:
    #         CLOCK_FREQ = int(self.clock_freq.text())
    #         self.clock_label.setText(f"Clock: {CLOCK_FREQ} MHz")
    #     except:
    #         pass

    def change_board(self):
        """Handle change board button."""
        self.set_status("Change board: not implemented")

    def start_pushbutton_Callback(self):
        """Start or restart the PulseBlaster program."""
        pb_init()
        if self.prg_order == 0:
            self.set_status("Must Load Board First", is_error=True)
        elif self.prg_order == 1:
            if pb_start() < 0:
                self.set_status(f"Start failed: {pb_get_error()}", is_error=True)
            else:
                self.prg_order = 2
                self.set_status("Running")
        elif self.prg_order == 2:
            self.set_status("Running already", is_error=True)
        elif self.prg_order < 0:
            self.set_status("NO board!", is_error=True)
        pb_close()

    def stop_pushbutton_Callback(self):
        """Stop the PulseBlaster program."""
        pb_init()
        if self.prg_order == 0:
            self.set_status("Must Load Board First", is_error=True)
        elif self.prg_order == 1:
            self.set_status("Stopped already")
        elif self.prg_order == 2:
            if pb_stop() < 0:
                self.set_status(f"Stop failed: {pb_get_error()}", is_error=True)
            else:
                self.prg_order = 1
                self.set_status("Stopped")
        elif self.prg_order < 0:
            self.set_status("NO board!", is_error=True)
        pb_close()

    def load_pushbutton_Callback(self):
        """Load the program into the PulseBlaster."""
        global CLOCK_FREQ
        self.generate_d_struct()
        
        self.logger.debug(f'd_struct: {self.d_struct}')
        
        if self.prg_order >= 0:
            try:
                pb_init()
                pb_core_clock(CLOCK_FREQ)
                pb_start_programming(PULSE_PROGRAM)
                
                # print(self.d_struct)
                for s_inst in self.d_struct:
                    self.logger.debug(f"Programming: flags={s_inst['flags']}, inst={s_inst['inst']}, "
                                     f"data={s_inst['instruct_data']}, dur={s_inst['dur'] * s_inst['d_mult']}")
                    pb_inst_pbonly(s_inst['flags'], s_inst['inst'], 
                                  s_inst['instruct_data'], s_inst['dur'] * s_inst['d_mult'])
                    # print(s_inst['flags'], s_inst['inst'], s_inst['instruct_data'], s_inst['dur'], s_inst['d_mult'])
                    # print(type(s_inst['flags']), type(s_inst['inst']), type(s_inst['instruct_data']), type(s_inst['dur']), type(s_inst['d_mult']))
                    err = pb_get_error()
                    # print(f"error = {err}")
                    if err and 'ok' not in err.lower():
                        self.logger.warning(f'Instruction warning: {err}')
                
                pb_stop_programming()
                self.prg_order = 1
                self.set_status("Program loaded")
                pb_close()
                
            except Exception as e:
                # print(e)
                self.logger.error(f'Load error: {e}')
                self.set_status(f"Load failed: {str(e)}", is_error=True)
                pb_close()
        else:
            self.set_status("NO board!", is_error=True)

    def closeEvent(self, event):
        """Handle window close - save channel names and cleanup."""
        # Save channel names before closing
        self.save_channel_names()
        
        # # Stop board if running
        # if self.prg_order == 2:
        #     self.stop_pushbutton_Callback()
        
        # Save last state before exiting
        self.save_state()

        event.accept()
        QApplication.quit()

def set_dark_theme(app):
    """Set application-wide dark theme"""
    app.setStyle("Fusion")
    
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ColorRole.WindowText, QColor(255, 255, 255))
    dark_palette.setColor(QPalette.ColorRole.Base, QColor(25, 25, 25))
    dark_palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(30, 30, 30))  # Dark tooltip
    dark_palette.setColor(QPalette.ColorRole.ToolTipText, QColor(212, 212, 212))  # Light text
    dark_palette.setColor(QPalette.ColorRole.Text, QColor(255, 255, 255))
    dark_palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ColorRole.ButtonText, QColor(255, 255, 255))
    dark_palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 0, 0))
    dark_palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.ColorRole.HighlightedText, QColor(0, 0, 0))
    # Disabled colors
    dark_palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, "#7f7f7f")
    dark_palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, "#7f7f7f")
    app.setPalette(dark_palette)
# ============================================================================
# Main Entry Point
# ============================================================================
if __name__ == '__main__':
    
    # Must be called before creating QApplication
    myappid = 'aglab.pb'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    setup_logging()
    app = QApplication(sys.argv)
    set_dark_theme(app)

    pbgui = SpinAPIGUI()

    app_icon = QIcon(expt_dir + r"\gui\pb\gui_icon.png") # .ico is preferred for Windows
    pbgui.setWindowIcon(app_icon)
    app.setWindowIcon(app_icon) # Sets it for the whole application

    pbgui.show()
    
    sys.exit(app.exec())

