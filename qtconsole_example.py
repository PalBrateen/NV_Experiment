import sys, os

# Fix for Windows DPI issue - MUST BE BEFORE Qt imports
# if sys.platform == 'win32':
#     # os.environ['QT_ENABLE_HIGHDPI_SCALING'] = '0'
#     # Alternative options (try one at a time if first doesn't work):
#     # os.environ['QT_AUTO_SCREEN_SCALE_FACTOR'] = '1'
#     os.environ['QT_SCALE_FACTOR'] = '1'

from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, 
                               QVBoxLayout, QHBoxLayout, QPushButton, 
                               QLabel, QSpinBox)
from PySide6.QtCore import QTimer
from PySide6.QtGui import QPalette, QColor
from qtconsole.rich_ipython_widget import RichIPythonWidget
from qtconsole.inprocess import QtInProcessKernelManager

class SimpleDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Qt Console Demo - Dark Mode")
        self.setGeometry(100, 100, 900, 600)
        
        # Some "experimental" variables you might want to access
        self.counter = 0
        self.measurement_value = 42.0
        self.is_running = False
        
        self.setup_ui()
        
        # Timer to simulate ongoing measurement/display update
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_display)
        self.timer.start(1000)  # Update every second
        
    def setup_ui(self):
        central_widget = QWidget()
        main_layout = QHBoxLayout()
        
        # Left side: Your "GUI controls"
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        
        # Display with dark mode styling
        self.counter_label = QLabel(f"Counter: {self.counter}")
        self.counter_label.setStyleSheet("""
            font-size: 24px; 
            padding: 20px;
            color: #00ff00;
            background-color: #1e1e1e;
            border-radius: 5px;
        """)
        left_layout.addWidget(self.counter_label)
        
        self.value_label = QLabel(f"Measurement: {self.measurement_value:.2f}")
        self.value_label.setStyleSheet("""
            font-size: 18px; 
            padding: 10px;
            color: #4ec9b0;
            background-color: #1e1e1e;
            border-radius: 5px;
        """)
        left_layout.addWidget(self.value_label)
        
        # Controls with dark styling
        btn_start = QPushButton("Start Counter")
        btn_start.clicked.connect(self.start_counter)
        btn_start.setStyleSheet("""
            QPushButton {
                background-color: #0e639c;
                color: white;
                padding: 10px;
                border: none;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #1177bb;
            }
            QPushButton:pressed {
                background-color: #0d5589;
            }
        """)
        left_layout.addWidget(btn_start)
        
        btn_stop = QPushButton("Stop Counter")
        btn_stop.clicked.connect(self.stop_counter)
        btn_stop.setStyleSheet("""
            QPushButton {
                background-color: #c5000b;
                color: white;
                padding: 10px;
                border: none;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #e81123;
            }
            QPushButton:pressed {
                background-color: #a80000;
            }
        """)
        left_layout.addWidget(btn_stop)
        
        btn_reset = QPushButton("Reset Counter")
        btn_reset.clicked.connect(self.reset_counter)
        btn_reset.setStyleSheet("""
            QPushButton {
                background-color: #3a3a3a;
                color: white;
                padding: 10px;
                border: none;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #505050;
            }
            QPushButton:pressed {
                background-color: #2a2a2a;
            }
        """)
        left_layout.addWidget(btn_reset)
        
        # SpinBox for measurement value
        spinbox_layout = QHBoxLayout()
        spinbox_label = QLabel("Set Measurement:")
        spinbox_label.setStyleSheet("color: #cccccc; font-size: 14px;")
        spinbox_layout.addWidget(spinbox_label)
        
        self.spinbox = QSpinBox()
        self.spinbox.setRange(0, 100)
        self.spinbox.setValue(int(self.measurement_value))
        self.spinbox.valueChanged.connect(self.set_measurement_value)
        self.spinbox.setStyleSheet("""
            QSpinBox {
                background-color: #3c3c3c;
                color: white;
                border: 1px solid #555555;
                padding: 5px;
                border-radius: 3px;
            }
            QSpinBox::up-button, QSpinBox::down-button {
                background-color: #505050;
                border: 1px solid #555555;
            }
            QSpinBox::up-button:hover, QSpinBox::down-button:hover {
                background-color: #606060;
            }
        """)
        spinbox_layout.addWidget(self.spinbox)
        left_layout.addLayout(spinbox_layout)
        
        left_layout.addStretch()
        
        # Style the left panel
        left_panel.setLayout(left_layout)
        left_panel.setStyleSheet("background-color: #252526;")
        
        # Right side: IPython Console with colors
        console = self.create_console()
        
        # Add both to main layout
        main_layout.addWidget(left_panel, stretch=1)
        main_layout.addWidget(console, stretch=1)
        
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)
        
        # Set main window background
        central_widget.setStyleSheet("background-color: #1e1e1e;")
    
    def create_console(self):
        """Create the embedded IPython console with proper dark colors"""
        # Create kernel
        kernel_manager = QtInProcessKernelManager()
        kernel_manager.start_kernel()
        kernel_client = kernel_manager.client()
        kernel_client.start_channels()
        
        # Create console
        console_widget = RichIPythonWidget()
        console_widget.kernel_manager = kernel_manager
        console_widget.kernel_client = kernel_client
        
        # MANUAL COLOR CONFIGURATION - This should work!
        # Set ANSI colors manually (this is what controls syntax highlighting)
        console_widget.ansi_codes = True  # Enable ANSI color codes
        
        # Configure the color palette manually
        from PySide6.QtGui import QColor
        
        # Set individual colors for syntax elements
        # These correspond to pygments token types
        console_widget._ansi_color_names = [
            'black', 'darkred', 'darkgreen', 'brown',
            'darkblue', 'darkviolet', 'steelblue', 'grey',
            'lightgrey', 'red', 'green', 'yellow',
            'blue', 'violet', 'lightblue', 'white'
        ]
        
        # Set the style sheet with comprehensive styling
        # background-color: ; #004052; #757170
        style_sheet = """
            QPlainTextEdit, QTextEdit {
                background-color: #1e1e1e;
                color: #d4d4d4;
                selection-background-color: #264f78;
                selection-color: #ffffff;
                font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
                font-size: 11pt;
            }
            QToolTip {
                background-color: #1e1e1e;
                color: #d4d4d4;
                border: 1px solid #454545;
            }
        """
        console_widget.setStyleSheet(style_sheet)
        
        # Set console to use dark background
        console_widget.kind = 'rich'  # Use rich text formatting
        
        # Configure pygments style after widget is created
        console_widget.syntax_style = 'native'  # Dark theme
        
        # Manually set color mappings for different text types
        # This is the key part that makes colors work!
        console_widget._ansi_foreground_colors = {
            0: QColor(0, 0, 0),           # black
            1: QColor(205, 0, 0),         # red
            2: QColor(0, 205, 0),         # green  
            3: QColor(205, 205, 0),       # yellow
            4: QColor(92, 92, 255),       # blue
            5: QColor(205, 0, 205),       # magenta
            6: QColor(0, 205, 205),       # cyan
            7: QColor(229, 229, 229),     # white
            8: QColor(127, 127, 127),     # bright black (grey)
            9: QColor(255, 0, 0),         # bright red
            10: QColor(0, 255, 0),        # bright green
            11: QColor(255, 255, 0),      # bright yellow
            12: QColor(92, 92, 255),      # bright blue
            13: QColor(255, 0, 255),      # bright magenta
            14: QColor(0, 255, 255),      # bright cyan
            15: QColor(255, 255, 255),    # bright white
        }
        
        console_widget._ansi_background_colors = {
            0: QColor(0, 0, 0),
            1: QColor(205, 0, 0),
            2: QColor(0, 205, 0),
            3: QColor(205, 205, 0),
            4: QColor(0, 0, 238),
            5: QColor(205, 0, 205),
            6: QColor(0, 205, 205),
            7: QColor(229, 229, 229),
        }
        
        # Push variables
        kernel_manager.kernel.shell.push({
            'gui': self,
        })
        
        # Print colorful welcome message using ANSI codes
        console_widget.execute("""
from IPython.core.getipython import get_ipython
ipython = get_ipython()

# Configure IPython to use colors
ipython.colors = 'Linux'  # or 'Neutral', 'NoColor', 'LightBG'

print("\\033[92m" + "="*50 + "\\033[0m")
print("\\033[96mConsole Ready!\\033[0m")
print("\\033[92m" + "="*50 + "\\033[0m")
print("\\033[93mTry these commands:\\033[0m")
print("  gui.counter")
print("  gui.start_counter()")
print("  gui.stop_counter()")
print()

# Test syntax highlighting
import numpy as np
x = 42
s = "hello world"
print("\\033[92mSyntax highlighting should now work!\\033[0m")
        """)
        
        return console_widget
    
    def start_counter(self):
        self.is_running = True
        
    def stop_counter(self):
        self.is_running = False
        
    def reset_counter(self):
        self.counter = 0
        self.update_display()
        
    def set_measurement_value(self, value):
        self.measurement_value = float(value)
        self.update_display()
        
    def update_display(self):
        if self.is_running:
            self.counter += 1
        
        self.counter_label.setText(f"Counter: {self.counter}")
        self.value_label.setText(f"Measurement: {self.measurement_value:.2f}")

def set_dark_theme(app):
    """Set application-wide dark theme"""
    app.setStyle("Fusion")
    
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.WindowText, QColor(255, 255, 255))
    dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
    dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ToolTipBase, QColor(30, 30, 30))  # Dark tooltip
    dark_palette.setColor(QPalette.ToolTipText, QColor(212, 212, 212))  # Light text
    dark_palette.setColor(QPalette.Text, QColor(255, 255, 255))
    dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ButtonText, QColor(255, 255, 255))
    dark_palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
    dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.HighlightedText, QColor(0, 0, 0))
    
    app.setPalette(dark_palette)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Apply dark theme to entire application
    set_dark_theme(app)
    
    window = SimpleDemo()
    window.show()
    sys.exit(app.exec())