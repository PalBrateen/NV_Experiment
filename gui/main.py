"""
NV Experiment GUI - Main Entry Point

Launch the NV center experiment control GUI with dark theme.

Usage:
    python gui/main.py
"""

import sys
import os

# Add parent directory to path to import experiment modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QPalette, QColor
from PySide6.QtCore import Qt
from pathlib import Path

# Import main window (will be created next)
from gui.main_window import NVExperimentGUI


def set_dark_palette(app: QApplication):
    """
    Set application-wide dark color palette.

    Based on verdi_v5_gui.py and qtconsole_example.py.
    """
    app.setStyle("Fusion")

    dark_palette = QPalette()

    # Window colors
    dark_palette.setColor(QPalette.Window, QColor(30, 30, 30))  # #1e1e1e
    dark_palette.setColor(QPalette.WindowText, QColor(212, 212, 212))  # #d4d4d4

    # Base colors (input fields, etc.)
    dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
    dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))

    # Tooltip colors
    dark_palette.setColor(QPalette.ToolTipBase, QColor(30, 30, 30))
    dark_palette.setColor(QPalette.ToolTipText, QColor(212, 212, 212))

    # Text colors
    dark_palette.setColor(QPalette.Text, QColor(212, 212, 212))
    dark_palette.setColor(QPalette.PlaceholderText, QColor(133, 133, 133))

    # Button colors
    dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ButtonText, QColor(212, 212, 212))

    # Bright text (emphasis)
    dark_palette.setColor(QPalette.BrightText, QColor(255, 255, 255))

    # Link colors
    dark_palette.setColor(QPalette.Link, QColor(14, 99, 156))  # #0e639c
    dark_palette.setColor(QPalette.LinkVisited, QColor(78, 201, 176))  # #4ec9b0

    # Selection colors
    dark_palette.setColor(QPalette.Highlight, QColor(38, 79, 120))  # #264f78
    dark_palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))

    # Disabled colors
    dark_palette.setColor(QPalette.Disabled, QPalette.Text, QColor(133, 133, 133))
    dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(133, 133, 133))
    dark_palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(133, 133, 133))

    app.setPalette(dark_palette)


def load_stylesheet(app: QApplication):
    """Load QSS stylesheet."""
    qss_path = Path(__file__).parent / "styles" / "dark_theme.qss"

    if qss_path.exists():
        with open(qss_path, 'r') as f:
            stylesheet = f.read()
            app.setStyleSheet(stylesheet)
    else:
        print(f"⚠️  Warning: Stylesheet not found at {qss_path}")


def main():
    """Main application entry point."""
    # Create application
    app = QApplication(sys.argv)

    # Set application metadata
    app.setApplicationName("NV Experiment Control")
    app.setApplicationVersion("1.0.0")
    app.setOrganizationName("NV Lab")

    # Apply dark theme
    set_dark_palette(app)
    load_stylesheet(app)

    # Create and show main window
    window = NVExperimentGUI()
    window.show()

    # Start event loop
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
