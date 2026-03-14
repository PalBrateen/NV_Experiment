"""
IPython Console Widget for NV Experiment GUI

Embedded IPython console with dark mode styling for live troubleshooting
and direct access to instruments and data during experiments.

Based on qtconsole_example.py with enhanced dark theme.
"""

from PySide6.QtWidgets import QWidget
from PySide6.QtGui import QColor
from PySide6.QtCore import Qt
from qtconsole.rich_ipython_widget import RichIPythonWidget
from qtconsole.inprocess import QtInProcessKernelManager


class IPythonConsoleWidget(RichIPythonWidget):
    """
    Embedded IPython console with dark mode and syntax highlighting.

    Features:
    - Dark background matching GUI theme
    - Syntax highlighting for Python code
    - Auto-complete and history
    - Direct access to GUI, instruments, and data
    - Thread-safe kernel communication
    """

    def __init__(self, namespace: dict = None, parent=None):
        """
        Initialize IPython console widget.

        Args:
            namespace: Dictionary of variables to expose in console
                      (e.g., {'gui': main_window, 'sg': signal_generator, ...})
            parent: Parent widget
        """
        super().__init__(parent)

        # Create in-process kernel
        self.kernel_manager = QtInProcessKernelManager()
        self.kernel_manager.start_kernel()

        # Create kernel client
        self.kernel_client = self.kernel_manager.client()
        self.kernel_client.start_channels()

        # Configure console styling
        self._configure_dark_theme()

        # Push namespace variables to kernel
        if namespace:
            self.push_namespace(namespace)

        # Print welcome message
        self._print_welcome()

    def _configure_dark_theme(self):
        """Configure dark mode colors and styling."""
        # Enable ANSI color codes
        self.ansi_codes = True

        # Set kind to rich (supports formatting)
        self.kind = 'rich'

        # Set pygments syntax highlighting style to dark theme
        self.syntax_style = 'native'  # Dark theme for code

        # Configure stylesheet
        style_sheet = """
            QPlainTextEdit, QTextEdit {
                background-color: #1e1e1e;
                color: #d4d4d4;
                selection-background-color: #264f78;
                selection-color: #ffffff;
                font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
                font-size: 10pt;
                border: 1px solid #454545;
            }
            QToolTip {
                background-color: #1e1e1e;
                color: #d4d4d4;
                border: 1px solid #454545;
            }
        """
        self.setStyleSheet(style_sheet)

        # Set ANSI foreground colors (for colored output)
        self._ansi_foreground_colors = {
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

        # Set ANSI background colors
        self._ansi_background_colors = {
            0: QColor(0, 0, 0),
            1: QColor(205, 0, 0),
            2: QColor(0, 205, 0),
            3: QColor(205, 205, 0),
            4: QColor(0, 0, 238),
            5: QColor(205, 0, 205),
            6: QColor(0, 205, 205),
            7: QColor(229, 229, 229),
        }

    def push_namespace(self, namespace: dict):
        """
        Push variables into the IPython kernel namespace.

        Args:
            namespace: Dictionary of {name: object} to expose
        """
        if self.kernel_manager and self.kernel_manager.kernel:
            self.kernel_manager.kernel.shell.push(namespace)

    def execute_command(self, command: str, hidden: bool = False):
        """
        Execute a command in the IPython kernel.

        Args:
            command: Python code to execute
            hidden: If True, don't show command in console
        """
        if hidden:
            self.kernel_client.execute(command, silent=True)
        else:
            self.execute(command)

    def _print_welcome(self):
        """Print welcome message to console."""
        welcome_code = """
from IPython.core.getipython import get_ipython
ipython = get_ipython()

# Configure IPython to use colors
ipython.colors = 'Linux'  # Dark terminal color scheme

print("\\033[92m" + "="*60 + "\\033[0m")
print("\\033[96m  NV Experiment Control - IPython Console\\033[0m")
print("\\033[92m" + "="*60 + "\\033[0m")
print("\\033[93mAvailable variables:\\033[0m")
print("  gui       - Main window instance")
print("  sg        - Signal Generator (if connected)")
print("  pb        - PulseBlaster (if connected)")
print("  daq       - DAQ tasks (if connected)")
print("  camera    - Camera (if connected)")
print("  data      - Last acquired dataset")
print("  np        - NumPy module")
print("  plt       - Matplotlib pyplot")
print()
print("\\033[93mUseful commands:\\033[0m")
print("  %whos               - List all variables")
print("  %history            - Show command history")
print("  gui.<tab>           - Autocomplete GUI methods")
print("  sg.frequency        - Access instrument parameters")
print()
print("\\033[92mReady for troubleshooting!\\033[0m")
print("\\033[92m" + "="*60 + "\\033[0m")
print()

# Import common modules into namespace
import numpy as np
import matplotlib.pyplot as plt
        """
        self.execute_command(welcome_code, hidden=False)

    def update_namespace(self, updates: dict):
        """
        Update namespace with new or changed variables.

        Args:
            updates: Dictionary of variables to update
        """
        self.push_namespace(updates)

    def shutdown(self):
        """Shutdown the kernel and clean up."""
        if self.kernel_client:
            self.kernel_client.stop_channels()

        if self.kernel_manager:
            self.kernel_manager.shutdown_kernel()

    def clear_console(self):
        """Clear the console output."""
        self.clear()


# Convenience function for creating console with common namespace
def create_nv_console(gui_instance, instruments: dict = {}, parent=None) -> IPythonConsoleWidget:
    """
    Create IPython console with NV experiment namespace.

    Args:
        gui_instance: Main GUI window instance
        instruments: Dictionary of instrument instances (sg, pb, daq, camera, etc.)
        parent: Parent widget

    Returns:
        Configured IPythonConsoleWidget
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Build namespace
    namespace = {
        'gui': gui_instance,
        'np': np,
        'plt': plt,
    }

    # Add instruments if provided
    if instruments:
        namespace.update(instruments)

    # Add placeholder for data
    namespace['data'] = None

    return IPythonConsoleWidget(namespace=namespace, parent=parent)
