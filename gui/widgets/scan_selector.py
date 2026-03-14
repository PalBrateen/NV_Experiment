"""
Scan Selector Widget for NV Experiment GUI

Select and compare multiple scans from an experiment.

Features:
- Checkbox list of completed scans
- Individual scan selection/deselection
- Overlay multiple scans on plot
- Color assignment for each scan
- Scan metadata display (timestamp, parameters)
- Export selected scans
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from typing import List, Dict, Any, Tuple
import numpy as np
from datetime import datetime

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QPushButton, QListWidget, QListWidgetItem,
    QCheckBox, QColorDialog
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from gui.styles.colors import *


class ScanSelectorWidget(QWidget):
    """
    Scan selector for comparing multiple experiment scans.

    Signals:
        selection_changed: Emitted when scan selection changes (selected_indices)
        color_changed: Emitted when scan color changes (scan_index, color)
    """

    selection_changed = Signal(list)  # List of selected scan indices
    color_changed = Signal(int, QColor)  # scan_index, color

    # Default color palette for scans
    DEFAULT_COLORS = [
        "#00ff00",  # Green
        "#ff0000",  # Red
        "#0000ff",  # Blue
        "#ffff00",  # Yellow
        "#ff00ff",  # Magenta
        "#00ffff",  # Cyan
        "#ffa500",  # Orange
        "#ff1493",  # Deep pink
        "#00ff7f",  # Spring green
        "#ff4500",  # Orange red
    ]

    def __init__(self, parent=None):
        """Initialize scan selector."""
        super().__init__(parent)

        # State
        self.scans = []  # List of scan metadata dicts
        self.scan_colors = {}  # scan_index -> QColor

        # Create UI
        self._create_ui()

    def _create_ui(self):
        """Create scan selector UI."""
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(10)

        # Header
        header_layout = QHBoxLayout()
        header_label = QLabel("Completed Scans")
        header_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        header_layout.addWidget(header_label)

        self.scan_count_label = QLabel("0 scans")
        self.scan_count_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        header_layout.addWidget(self.scan_count_label)
        header_layout.addStretch()

        layout.addLayout(header_layout)

        # Scan list
        self.scan_list = QListWidget()
        self.scan_list.setMaximumHeight(200)
        self.scan_list.setAlternatingRowColors(True)
        self.scan_list.itemChanged.connect(self._on_selection_changed)
        layout.addWidget(self.scan_list)

        # Control buttons
        button_layout = QHBoxLayout()

        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self._select_all)
        button_layout.addWidget(self.select_all_btn)

        self.deselect_all_btn = QPushButton("Deselect All")
        self.deselect_all_btn.clicked.connect(self._deselect_all)
        button_layout.addWidget(self.deselect_all_btn)

        self.change_color_btn = QPushButton("Change Color")
        self.change_color_btn.clicked.connect(self._change_selected_color)
        button_layout.addWidget(self.change_color_btn)

        layout.addLayout(button_layout)

        # Scan info
        info_group = QGroupBox("Selected Scan Info")
        info_layout = QVBoxLayout()

        self.info_label = QLabel("Select a scan to view details")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet(f"color: {TEXT_SECONDARY};")
        info_layout.addWidget(self.info_label)

        info_group.setLayout(info_layout)
        layout.addWidget(info_group)

        layout.addStretch()

        self.setLayout(layout)

        # Connect list selection for info display
        self.scan_list.currentItemChanged.connect(self._on_current_item_changed)

    def add_scan(self, scan_data: Tuple[np.ndarray, np.ndarray], metadata: Dict[str, Any] = None):
        """
        Add a new scan to the selector.

        Args:
            scan_data: Tuple of (x_data, y_data)
            metadata: Optional metadata dict (timestamp, parameters, etc.)
        """
        scan_index = len(self.scans)

        # Store scan
        scan_info = {
            'index': scan_index,
            'data': scan_data,
            'metadata': metadata or {},
            'timestamp': datetime.now()
        }
        self.scans.append(scan_info)

        # Assign color
        color_index = scan_index % len(self.DEFAULT_COLORS)
        self.scan_colors[scan_index] = QColor(self.DEFAULT_COLORS[color_index])

        # Add to list
        self._add_scan_item(scan_info)

        # Update count
        self.scan_count_label.setText(f"{len(self.scans)} scans")

    def _add_scan_item(self, scan_info: Dict[str, Any]):
        """Add scan item to list widget."""
        scan_index = scan_info['index']
        timestamp = scan_info['timestamp'].strftime("%H:%M:%S")

        # Create list item
        item = QListWidgetItem()
        item.setText(f"Scan #{scan_index + 1} ({timestamp})")
        item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
        item.setCheckState(Qt.Checked)  # Auto-select new scans

        # Set color indicator
        color = self.scan_colors[scan_index]
        item.setForeground(color)

        # Store scan index in item data
        item.setData(Qt.UserRole, scan_index)

        self.scan_list.addItem(item)

    def _on_selection_changed(self, item: QListWidgetItem):
        """Handle scan selection change."""
        selected_indices = self.get_selected_scan_indices()
        self.selection_changed.emit(selected_indices)

    def _on_current_item_changed(self, current: QListWidgetItem, previous: QListWidgetItem):
        """Handle current item change (for info display)."""
        if current is None:
            self.info_label.setText("Select a scan to view details")
            return

        scan_index = current.data(Qt.UserRole)
        scan_info = self.scans[scan_index]

        # Format info
        timestamp = scan_info['timestamp'].strftime("%Y-%m-%d %H:%M:%S")
        x_data, y_data = scan_info['data']
        metadata = scan_info['metadata']

        info_text = f"Scan #{scan_index + 1}\n"
        info_text += f"Time: {timestamp}\n"
        info_text += f"Points: {len(x_data)}\n"
        info_text += f"X Range: {x_data[0]:.3e} to {x_data[-1]:.3e}\n"
        info_text += f"Y Range: {y_data.min():.3f} to {y_data.max():.3f}\n"

        if metadata:
            info_text += "\nMetadata:\n"
            for key, value in metadata.items():
                info_text += f"  {key}: {value}\n"

        self.info_label.setText(info_text)

    def _select_all(self):
        """Select all scans."""
        for i in range(self.scan_list.count()):
            item = self.scan_list.item(i)
            item.setCheckState(Qt.Checked)

    def _deselect_all(self):
        """Deselect all scans."""
        for i in range(self.scan_list.count()):
            item = self.scan_list.item(i)
            item.setCheckState(Qt.Unchecked)

    def _change_selected_color(self):
        """Change color of currently selected scan."""
        current_item = self.scan_list.currentItem()
        if current_item is None:
            return

        scan_index = current_item.data(Qt.UserRole)
        current_color = self.scan_colors[scan_index]

        # Open color dialog
        color = QColorDialog.getColor(current_color, self, "Select Scan Color")
        if color.isValid():
            self.scan_colors[scan_index] = color
            current_item.setForeground(color)
            self.color_changed.emit(scan_index, color)

    def get_selected_scan_indices(self) -> List[int]:
        """
        Get indices of selected scans.

        Returns:
            List of selected scan indices
        """
        selected = []
        for i in range(self.scan_list.count()):
            item = self.scan_list.item(i)
            if item.checkState() == Qt.Checked:
                scan_index = item.data(Qt.UserRole)
                selected.append(scan_index)
        return selected

    def get_selected_scans(self) -> List[Dict[str, Any]]:
        """
        Get selected scan data.

        Returns:
            List of scan info dicts for selected scans
        """
        selected_indices = self.get_selected_scan_indices()
        return [self.scans[i] for i in selected_indices]

    def get_scan_color(self, scan_index: int) -> QColor:
        """
        Get color for a scan.

        Args:
            scan_index: Scan index

        Returns:
            QColor for the scan
        """
        return self.scan_colors.get(scan_index, QColor("#ffffff"))

    def clear_scans(self):
        """Clear all scans."""
        self.scans = []
        self.scan_colors = {}
        self.scan_list.clear()
        self.scan_count_label.setText("0 scans")
        self.info_label.setText("Select a scan to view details")

    def remove_selected_scans(self):
        """Remove selected scans."""
        selected_indices = self.get_selected_scan_indices()

        # Remove in reverse order to maintain indices
        for index in sorted(selected_indices, reverse=True):
            # Remove from list widget
            for i in range(self.scan_list.count()):
                item = self.scan_list.item(i)
                if item.data(Qt.UserRole) == index:
                    self.scan_list.takeItem(i)
                    break

            # Remove from data
            del self.scans[index]
            del self.scan_colors[index]

        # Re-index remaining scans
        for i, scan in enumerate(self.scans):
            scan['index'] = i

        # Update count
        self.scan_count_label.setText(f"{len(self.scans)} scans")

        # Emit change
        self.selection_changed.emit(self.get_selected_scan_indices())


# Example usage
if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Apply dark theme
    from gui.main import set_dark_palette, load_stylesheet
    set_dark_palette(app)
    load_stylesheet(app)

    selector = ScanSelectorWidget()

    # Add test scans
    for i in range(5):
        x_data = np.linspace(2.85e9, 2.89e9, 51)
        y_data = 2.5 - 0.5 / (1 + ((x_data - 2.87e9) / 5e6)**2) + 0.05 * np.random.randn(len(x_data))

        metadata = {
            'experiment_type': 'ESR',
            'averages': 1000,
            'power': 8.0
        }

        selector.add_scan((x_data, y_data), metadata)

    def on_selection_changed(selected):
        print(f"Selected scans: {selected}")

    selector.selection_changed.connect(on_selection_changed)

    selector.show()

    sys.exit(app.exec())
