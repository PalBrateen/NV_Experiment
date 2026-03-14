"""
Color constants for NV Experiment GUI dark theme.
Based on verdi_v5_gui.py and qtconsole_example.py.
"""

from PySide6.QtGui import QColor

# Main background colors
BG_DARK = "#1e1e1e"           # Main background
BG_MEDIUM = "#252526"         # Panel background
BG_LIGHT = "#2d2d30"          # Widget background
BG_LIGHTER = "#3c3c3c"        # Input background

# Border colors
BORDER_DARK = "#454545"       # Subtle borders
BORDER_MEDIUM = "#555555"     # Input borders
BORDER_LIGHT = "#707070"      # Hover borders

# Text colors
TEXT_PRIMARY = "#d4d4d4"      # Main text
TEXT_SECONDARY = "#cccccc"    # Secondary text
TEXT_DISABLED = "#858585"     # Disabled text
TEXT_BRIGHT = "#ffffff"       # Emphasized text

# Status indicator colors
STATUS_OK = "#00ff00"         # Green - Connected/Success
STATUS_WARNING = "#ffff00"    # Yellow - Warning/Claimed
STATUS_ERROR = "#ff0000"      # Red - Error/Disconnected
STATUS_INFO = "#4ec9b0"       # Cyan - Info

# Syntax highlighting (for IPython console)
SYNTAX_KEYWORD = "#569cd6"    # Blue - Keywords
SYNTAX_STRING = "#ce9178"     # Orange - Strings
SYNTAX_NUMBER = "#b5cea8"     # Green - Numbers
SYNTAX_COMMENT = "#6a9955"    # Green - Comments
SYNTAX_FUNCTION = "#dcdcaa"   # Yellow - Functions
SYNTAX_CLASS = "#4ec9b0"      # Cyan - Classes

# Button colors
BTN_PRIMARY_BG = "#0e639c"    # Blue button background
BTN_PRIMARY_HOVER = "#1177bb" # Blue button hover
BTN_PRIMARY_PRESSED = "#0d5589" # Blue button pressed

BTN_DANGER_BG = "#c5000b"     # Red button background
BTN_DANGER_HOVER = "#e81123"  # Red button hover
BTN_DANGER_PRESSED = "#a80000" # Red button pressed

BTN_NEUTRAL_BG = "#3a3a3a"    # Gray button background
BTN_NEUTRAL_HOVER = "#505050" # Gray button hover
BTN_NEUTRAL_PRESSED = "#2a2a2a" # Gray button pressed

# Selection colors
SELECTION_BG = "#264f78"      # Selection background
SELECTION_TEXT = "#ffffff"    # Selection text

# Plot colors (for pyqtgraph)
PLOT_BG = "#1e1e1e"          # Plot background
PLOT_GRID = "#3c3c3c"        # Grid lines
PLOT_AXIS = "#cccccc"        # Axis lines/text
PLOT_LINE1 = "#00ff00"       # First line color (green)
PLOT_LINE2 = "#ff00ff"       # Second line color (magenta)
PLOT_LINE3 = "#00ffff"       # Third line color (cyan)
PLOT_LINE4 = "#ffff00"       # Fourth line color (yellow)
PLOT_MARKER = "#ff0000"      # Current point marker (red)

# Progress bar colors
PROGRESS_BG = "#3c3c3c"      # Progress bar background
PROGRESS_CHUNK = "#0e639c"   # Progress bar fill

# Tooltip colors
TOOLTIP_BG = "#1e1e1e"       # Tooltip background
TOOLTIP_TEXT = "#d4d4d4"     # Tooltip text
TOOLTIP_BORDER = "#454545"   # Tooltip border


def get_qcolor(hex_color: str) -> QColor:
    """Convert hex color string to QColor object."""
    return QColor(hex_color)


def get_rgb(hex_color: str) -> tuple[int, int, int]:
    """Convert hex color string to RGB tuple."""
    color = QColor(hex_color)
    return color.red(), color.green(), color.blue()
