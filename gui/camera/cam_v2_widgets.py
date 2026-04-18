"""
Camera Viewer v2 - Widgets Module
==================================

Dockable widgets for camera viewer:
- ImageDisplayDock: Main image with crosshair and line profiles
- HistogramDock: Vertical histogram with Hi/Lo controls
- IntensityTraceDock: Rolling intensity plot
- ImageFFTDock: 2D FFT with filtering
- ControlDock: Tabbed camera controls

Author: Claude/BP
Date: January 2025
"""

import numpy as np
import pyqtgraph as pg
from typing import Optional, Tuple, List, Dict, Callable
from enum import Enum

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
    QPushButton, QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QTabWidget, QSlider, QDockWidget, QSplitter,
    QFrame, QSizePolicy, QLineEdit, QFormLayout, QProgressBar,
    QRadioButton, QButtonGroup
)
from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtGui import QColor, QPen, QFont

from cam_v2_core import (
    CameraType, TriggerMode, StorageMode, IntensityMode,
    CameraConfig, AcquisitionConfig, DisplayConfig,
    auto_levels, compute_histogram
)


# =============================================================================
# CUSTOM PYQTGRAPH ITEMS
# =============================================================================

class CrosshairROI(pg.GraphicsObject):
    """
    Crosshair cursor with draggable line for arbitrary angle profiles.
    
    Features:
    - Vertical and horizontal lines following mouse
    - Draggable angle line for arbitrary profiles
    - Position/value readout
    """
    
    sigPositionChanged = Signal(float, float)  # x, y position
    sigLineChanged = Signal(object, object)  # start_point, end_point
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Crosshair lines (follow mouse)
        self.vline = pg.InfiniteLine(angle=90, movable=False,
                                      pen=pg.mkPen('#FFAA00', width=1, 
                                                   style=Qt.PenStyle.DashLine))
        self.hline = pg.InfiniteLine(angle=0, movable=False,
                                      pen=pg.mkPen('#FFAA00', width=1,
                                                   style=Qt.PenStyle.DashLine))
        
        # Arbitrary angle line (draggable)
        self.profile_line = pg.LineSegmentROI(
            positions=[[0, 0], [100, 100]],
            pen=pg.mkPen('#00FF00', width=2)
        )
        self.profile_line.sigRegionChanged.connect(self._on_line_changed)
        
        self._visible = True
        self._show_profile_line = False
        
    def set_position(self, x: float, y: float):
        """Set crosshair position."""
        self.vline.setPos(x)
        self.hline.setPos(y)
        self.sigPositionChanged.emit(x, y)
        
    def show_profile_line(self, show: bool):
        """Show/hide arbitrary angle profile line."""
        self._show_profile_line = show
        self.profile_line.setVisible(show)
        
    def get_profile_line_coords(self) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """Get profile line start and end coordinates."""
        handles = self.profile_line.getHandles()
        if len(handles) >= 2:
            p1 = handles[0].pos()
            p2 = handles[1].pos()
            return (p1.x(), p1.y()), (p2.x(), p2.y())
        return (0, 0), (100, 100)
    
    def _on_line_changed(self):
        """Emit signal when profile line changes."""
        p1, p2 = self.get_profile_line_coords()
        self.sigLineChanged.emit(p1, p2)
        
    def paint(self, p, *args):
        pass
    
    def boundingRect(self):
        return pg.QtCore.QRectF()


# =============================================================================
# IMAGE DISPLAY DOCK
# =============================================================================

class ImageDisplayDock(QDockWidget):
    """
    Main image display with:
    - PyQtGraph ImageView (without built-in histogram)
    - Mouse crosshair cursor
    - Horizontal and vertical line profiles
    - Arbitrary angle line profile (draggable)
    - ROI selection
    - Colormap control
    - Click-to-draw mode for line and ROI
    
    Signals:
        roiChanged: ROI selection changed (x, y, w, h)
        cursorMoved: Cursor position changed (x, y, value)
        profileLineChanged: Arbitrary line profile changed (p1, p2)
    """
    
    roiChanged = Signal(int, int, int, int)
    roiApplyToSensor = Signal(int, int, int, int)  # Apply ROI to camera hardware
    cursorMoved = Signal(float, float, float)
    profileLineChanged = Signal(object, object)
    
    def __init__(self, parent=None):
        super().__init__("Image Display", parent)
        self.setObjectName("ImageDisplayDock")
        
        # Enable all dock features
        self.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable |
            QDockWidget.DockWidgetFeature.DockWidgetFloatable |
            QDockWidget.DockWidgetFeature.DockWidgetClosable
        )
        
        # Main widget
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)
        
        # Create splitter for image + profiles
        self.splitter = QSplitter(Qt.Orientation.Vertical)
        
        # === TOP: Image + Right Profile ===
        top_splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Image view (custom, no histogram)
        self.image_widget = pg.PlotWidget()
        self.image_widget.setAspectLocked(True)
        self.image_widget.setBackground('#1a1a1a')
        self.image_widget.hideAxis('left')
        self.image_widget.hideAxis('bottom')
        
        # Image item
        self.image_item = pg.ImageItem()
        self.image_widget.addItem(self.image_item)
        
        # Crosshair
        self.vline = pg.InfiniteLine(angle=90, movable=False,
                                      pen=pg.mkPen('#FFAA00', width=1,
                                                   style=Qt.PenStyle.DashLine))
        self.hline = pg.InfiniteLine(angle=0, movable=False,
                                      pen=pg.mkPen('#FFAA00', width=1,
                                                   style=Qt.PenStyle.DashLine))
        self.vline.setVisible(False)
        self.hline.setVisible(False)
        self.image_widget.addItem(self.vline, ignoreBounds=True)
        self.image_widget.addItem(self.hline, ignoreBounds=True)
        
        # Arbitrary angle line profile ROI
        self.profile_roi = pg.LineSegmentROI(
            positions=[[50, 50], [200, 200]],
            pen=pg.mkPen('#00FF00', width=2)
        )
        self.profile_roi.setVisible(False)
        self.profile_roi.sigRegionChanged.connect(self._on_profile_roi_changed)
        self.image_widget.addItem(self.profile_roi)
        
        # ROI for intensity calculation
        self.intensity_roi = pg.RectROI([0, 0], [100, 100], 
                                         pen=pg.mkPen('#00BFFF', width=2))
        self.intensity_roi.setVisible(False)
        self.intensity_roi.sigRegionChanged.connect(self._on_roi_changed)
        self.image_widget.addItem(self.intensity_roi)
        
        # Cursor label
        self.cursor_label = pg.TextItem(color='#FFAA00', anchor=(0, 1))
        self.cursor_label.setFont(QFont('Consolas', 9))
        self.cursor_label.setVisible(False)
        self.image_widget.addItem(self.cursor_label, ignoreBounds=True)
        
        # Enable mouse tracking
        self.image_widget.scene().sigMouseMoved.connect(self._on_mouse_moved)
        self.image_widget.scene().sigMouseClicked.connect(self._on_mouse_clicked)
        
        top_splitter.addWidget(self.image_widget)
        
        # Vertical profile (right side) - shows column profile at cursor
        self.v_profile_widget = pg.PlotWidget()
        self.v_profile_widget.setBackground('#1a1a1a')
        self.v_profile_widget.setMaximumWidth(150)
        self.v_profile_widget.setLabel('left', 'Y')
        self.v_profile_widget.showGrid(y=True, alpha=0.3)
        self.v_profile_curve = self.v_profile_widget.plot(pen='#00BFFF')
        # Rotate for vertical orientation
        self.v_profile_widget.getPlotItem().invertY(True)
        top_splitter.addWidget(self.v_profile_widget)
        
        top_splitter.setSizes([800, 150])
        self.splitter.addWidget(top_splitter)
        
        # === BOTTOM: Horizontal profile ===
        self.h_profile_widget = pg.PlotWidget()
        self.h_profile_widget.setBackground('#1a1a1a')
        self.h_profile_widget.setMaximumHeight(120)
        self.h_profile_widget.setLabel('bottom', 'X')
        self.h_profile_widget.showGrid(x=True, alpha=0.3)
        self.h_profile_curve = self.h_profile_widget.plot(pen='#FF6B6B')
        
        # Arbitrary line profile
        self.arb_profile_curve = self.h_profile_widget.plot(pen='#00FF00')
        self.arb_profile_curve.setVisible(False)
        
        self.splitter.addWidget(self.h_profile_widget)
        
        self.splitter.setSizes([600, 120])
        layout.addWidget(self.splitter)
        
        # === STATUS BAR ===
        status_layout = QHBoxLayout()
        
        self.pos_label = QLabel("X: -- Y: --")
        self.pos_label.setStyleSheet("color: #888; font-family: Consolas;")
        status_layout.addWidget(self.pos_label)
        
        self.val_label = QLabel("Value: --")
        self.val_label.setStyleSheet("color: #888; font-family: Consolas;")
        status_layout.addWidget(self.val_label)
        
        status_layout.addStretch()
        
        # Crosshair toggle
        self.crosshair_cb = QCheckBox("Crosshair")
        self.crosshair_cb.setChecked(False)
        self.crosshair_cb.toggled.connect(self._toggle_crosshair)
        status_layout.addWidget(self.crosshair_cb)
        
        # Draw Line button
        self.draw_line_btn = QPushButton("Draw Line")
        self.draw_line_btn.setCheckable(True)
        self.draw_line_btn.setMaximumWidth(80)
        self.draw_line_btn.toggled.connect(self._on_draw_line_toggled)
        status_layout.addWidget(self.draw_line_btn)
        
        # Clear Line button
        self.clear_line_btn = QPushButton("×")
        self.clear_line_btn.setMaximumWidth(25)
        self.clear_line_btn.setToolTip("Clear line profile")
        self.clear_line_btn.clicked.connect(self._clear_line_profile)
        status_layout.addWidget(self.clear_line_btn)
        
        # Show ROI button (toggle draggable rect)
        self.show_roi_btn = QPushButton("Show ROI")
        self.show_roi_btn.setCheckable(True)
        self.show_roi_btn.setMaximumWidth(80)
        self.show_roi_btn.toggled.connect(self._on_show_roi_toggled)
        status_layout.addWidget(self.show_roi_btn)
        
        # Apply ROI to sensor (hardware crop for speed)
        self.apply_roi_btn = QPushButton("→ Sensor")
        self.apply_roi_btn.setMaximumWidth(70)
        self.apply_roi_btn.setToolTip("Apply ROI to camera sensor (stops acquisition)")
        self.apply_roi_btn.clicked.connect(self._apply_roi_to_sensor)
        self.apply_roi_btn.setEnabled(False)
        status_layout.addWidget(self.apply_roi_btn)
        
        # Clear ROI button
        self.clear_roi_btn = QPushButton("×")
        self.clear_roi_btn.setMaximumWidth(25)
        self.clear_roi_btn.setToolTip("Clear ROI")
        self.clear_roi_btn.clicked.connect(self._clear_roi)
        status_layout.addWidget(self.clear_roi_btn)
        
        layout.addLayout(status_layout)
        
        self.setWidget(widget)
        
        # State
        self._current_frame: Optional[np.ndarray] = None
        self._levels = (0, 65535)
        self._cursor_pos = (0, 0)
        self._colormap = 'gray'
        
        # Drawing state (only for line - ROI is always draggable)
        self._draw_mode = None  # None or 'line'
        self._drawing = False
        self._draw_start = None
        
        # Profile update debounce (lazy updates - don't recompute every mouse move)
        self._profile_timer = QTimer(self)
        self._profile_timer.setSingleShot(True)
        self._profile_timer.setInterval(50)  # 50ms debounce = max 20 Hz
        self._profile_timer.timeout.connect(self._update_profiles)
        
    def set_image(self, frame: np.ndarray, autoRange: bool = False):
        """Update displayed image."""
        self._current_frame = frame
        self.image_item.setImage(frame, autoLevels=False)
        self.image_item.setLevels(self._levels)
        
        if autoRange:
            self.image_widget.autoRange()
            
        # Update profiles at current cursor position (only if not in line/ROI mode)
        if not self.profile_roi.isVisible() and not self.intensity_roi.isVisible():
            self._profile_timer.start()  # Debounced
        elif self.profile_roi.isVisible():
            self._on_profile_roi_changed()
        
    def set_levels(self, lo: float, hi: float):
        """Set display levels."""
        self._levels = (lo, hi)
        self.image_item.setLevels(self._levels)
        
    def set_colormap(self, name: str):
        """Set colormap."""
        self._colormap = name
        try:
            cmap = pg.colormap.get(name)
            self.image_item.setLookupTable(cmap.getLookupTable())
        except:
            pass  # Use default
            
    def _on_mouse_moved(self, pos):
        """Handle mouse movement."""
        if self._current_frame is None:
            return
            
        # Map to image coordinates
        mouse_point = self.image_widget.plotItem.vb.mapSceneToView(pos)
        x, y = int(mouse_point.x()), int(mouse_point.y())
        
        h, w = self._current_frame.shape[:2]
        if 0 <= x < w and 0 <= y < h:
            self._cursor_pos = (x, y)
            
            # Update crosshair position (even if not visible)
            self.vline.setPos(x)
            self.hline.setPos(y)
            
            # Get pixel value
            if self._current_frame.ndim == 2:
                value = self._current_frame[y, x]
            else:
                value = self._current_frame[y, x, 0]  # First channel
                
            # Update labels (always)
            self.pos_label.setText(f"X: {x} Y: {y}")
            self.val_label.setText(f"Value: {value}")
            
            # Update cursor label
            self.cursor_label.setText(f"({x}, {y}): {value}")
            self.cursor_label.setPos(x + 10, y + 10)
            
            # Only update line profiles if line profile and ROI are NOT visible
            # This prevents messy overlapping profiles
            # Use debounce timer to avoid recomputing every pixel of mouse move
            if not self.profile_roi.isVisible() and not self.intensity_roi.isVisible():
                self._profile_timer.start()  # Restart debounce
            
            # Emit signal
            self.cursorMoved.emit(float(x), float(y), float(value))
            
    def _on_mouse_clicked(self, event):
        """Handle mouse click for drawing line only (ROI is always draggable)."""
        if self._draw_mode != 'line':
            return
            
        # Get click position in image coordinates
        pos = self.image_widget.plotItem.vb.mapSceneToView(event.scenePos())
        x, y = pos.x(), pos.y()
        
        if not self._drawing:
            # First click - start drawing line
            self._draw_start = (x, y)
            self._drawing = True
            self.profile_roi.setVisible(True)
        else:
            # Second click - finish line
            end = (x, y)
            handles = self.profile_roi.getHandles()
            if len(handles) >= 2:
                handles[0].setPos(self._draw_start[0], self._draw_start[1])
                handles[1].setPos(end[0], end[1])
            self.profile_roi.setVisible(True)
            self._on_profile_roi_changed()
            
            # Reset drawing mode
            self._draw_mode = None
            self._drawing = False
            self._draw_start = None
            self.draw_line_btn.setChecked(False)
            self.image_widget.setCursor(Qt.CursorShape.ArrowCursor)
            
    def _on_draw_line_toggled(self, checked: bool):
        """Handle draw line button toggle."""
        if checked:
            self._draw_mode = 'line'
            self._drawing = False
            self.image_widget.setCursor(Qt.CursorShape.CrossCursor)
        else:
            if self._draw_mode == 'line':
                self._draw_mode = None
                self._drawing = False
                self.image_widget.setCursor(Qt.CursorShape.ArrowCursor)
                
    def _on_show_roi_toggled(self, checked: bool):
        """Toggle draggable ROI visibility."""
        self.intensity_roi.setVisible(checked)
        self.apply_roi_btn.setEnabled(checked)
        if checked:
            # Initialize ROI to center quarter of current frame
            if self._current_frame is not None:
                h, w = self._current_frame.shape[:2]
                self.intensity_roi.setPos([w//4, h//4])
                self.intensity_roi.setSize([w//2, h//2])
            self._on_roi_changed()
                
    def _clear_line_profile(self):
        """Clear line profile and its plot data."""
        self.profile_roi.setVisible(False)
        self.arb_profile_curve.setData([], [])  # Actually clear the data
        self.arb_profile_curve.setVisible(False)
        self.draw_line_btn.setChecked(False)
        self._draw_mode = None
        self._drawing = False
        
    def _clear_roi(self):
        """Clear ROI."""
        self.intensity_roi.setVisible(False)
        self.show_roi_btn.setChecked(False)
        self.apply_roi_btn.setEnabled(False)
        self._draw_mode = None
        self._drawing = False
        # Emit signal to clear ROI-based intensity calculation
        self.roiChanged.emit(0, 0, 0, 0)
        
    def _apply_roi_to_sensor(self):
        """Emit signal to apply ROI to camera hardware sensor."""
        if not self.intensity_roi.isVisible():
            return
        pos = self.intensity_roi.pos()
        size = self.intensity_roi.size()
        x, y = int(pos.x()), int(pos.y())
        w, h = int(size.x()), int(size.y())
        if w > 0 and h > 0:
            self.roiApplyToSensor.emit(x, y, w, h)
            
    def _update_profiles(self):
        """Update line profiles at cursor position."""
        if self._current_frame is None:
            return
            
        x, y = self._cursor_pos
        h, w = self._current_frame.shape[:2]
        
        # Horizontal profile (row at y)
        if 0 <= y < h:
            if self._current_frame.ndim == 2:
                h_profile = self._current_frame[y, :]
            else:
                h_profile = self._current_frame[y, :, 0]
            self.h_profile_curve.setData(np.arange(len(h_profile)), h_profile)
            
        # Vertical profile (column at x)
        if 0 <= x < w:
            if self._current_frame.ndim == 2:
                v_profile = self._current_frame[:, x]
            else:
                v_profile = self._current_frame[:, x, 0]
            # Plot with Y as position (vertical)
            self.v_profile_curve.setData(v_profile, np.arange(len(v_profile)))
            
    def _on_profile_roi_changed(self):
        """Handle arbitrary line profile ROI change."""
        if self._current_frame is None:
            return
            
        # Get line coordinates
        handles = self.profile_roi.getHandles()
        if len(handles) < 2:
            return
            
        p1 = handles[0].pos()
        p2 = handles[1].pos()
        
        # Extract profile along line
        profile = self._extract_line_profile(
            self._current_frame, 
            (p1.x(), p1.y()), 
            (p2.x(), p2.y())
        )
        
        if profile is not None:
            self.arb_profile_curve.setData(profile)
            self.arb_profile_curve.setVisible(True)
            
        self.profileLineChanged.emit((p1.x(), p1.y()), (p2.x(), p2.y()))
        
    def _extract_line_profile(self, image: np.ndarray, 
                               p1: Tuple[float, float],
                               p2: Tuple[float, float],
                               num_points: int = None) -> Optional[np.ndarray]:
        """Extract profile along arbitrary line."""
        x1, y1 = p1
        x2, y2 = p2
        
        length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if length < 1:
            return None
            
        if num_points is None:
            num_points = int(length)
            
        # Sample points along line
        x_coords = np.linspace(x1, x2, num_points)
        y_coords = np.linspace(y1, y2, num_points)
        
        # Clip to image bounds
        h, w = image.shape[:2]
        valid = (x_coords >= 0) & (x_coords < w) & (y_coords >= 0) & (y_coords < h)
        
        if not np.any(valid):
            return None
            
        x_coords = x_coords[valid].astype(int)
        y_coords = y_coords[valid].astype(int)
        
        if image.ndim == 2:
            return image[y_coords, x_coords]
        else:
            return image[y_coords, x_coords, 0]
            
    def _on_roi_changed(self):
        """Handle ROI change."""
        pos = self.intensity_roi.pos()
        size = self.intensity_roi.size()
        x, y = int(pos.x()), int(pos.y())
        w, h = int(size.x()), int(size.y())
        self.roiChanged.emit(x, y, w, h)
        
    def _toggle_crosshair(self, enabled: bool):
        """Toggle crosshair visibility."""
        self.vline.setVisible(enabled)
        self.hline.setVisible(enabled)
        self.cursor_label.setVisible(enabled)
            
    def get_roi(self) -> Optional[Tuple[int, int, int, int]]:
        """Get current ROI if visible."""
        if not self.intensity_roi.isVisible():
            return None
        pos = self.intensity_roi.pos()
        size = self.intensity_roi.size()
        return (int(pos.x()), int(pos.y()), int(size.x()), int(size.y()))


# =============================================================================
# HISTOGRAM DOCK (Vertical orientation with RGB support)
# =============================================================================

class HistogramDock(QDockWidget):
    """
    Histogram with level controls (vertical orientation).
    
    Features:
    - Rotated histogram aligned with image Y-axis
    - Hi/Lo level sliders with numeric input
    - Auto Hi/Lo button (percentile-based)
    - Statistics display
    - Colormap selection
    - RGB histogram for color cameras
    
    Signals:
        levelsChanged: Manual level change (lo, hi)
        colormapChanged: Colormap selection changed (name)
    """
    
    levelsChanged = Signal(float, float)
    colormapChanged = Signal(str)
    
    def __init__(self, parent=None):
        super().__init__("Histogram", parent)
        self.setObjectName("HistogramDock")
        
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        
        # === HISTOGRAM PLOT (vertical) ===
        self.hist_widget = pg.PlotWidget()
        self.hist_widget.setBackground('#1a1a1a')
        self.hist_widget.setMaximumWidth(200)
        self.hist_widget.showGrid(y=True, alpha=0.3)
        
        # Rotate: X is counts, Y is intensity
        self.hist_widget.setLabel('left', 'Intensity')
        self.hist_widget.setLabel('bottom', 'Count')
        
        # Histogram curves - grayscale and RGB
        self.hist_item = pg.PlotCurveItem(pen='#888888', fillLevel=0, 
                                           brush=(100, 100, 100, 150))
        self.hist_widget.addItem(self.hist_item)
        
        # RGB histogram curves (hidden by default)
        self.hist_r = pg.PlotCurveItem(pen=pg.mkPen('#FF4444', width=1.5), fillLevel=0,
                                        brush=(255, 68, 68, 80))
        self.hist_g = pg.PlotCurveItem(pen=pg.mkPen('#44FF44', width=1.5), fillLevel=0,
                                        brush=(68, 255, 68, 80))
        self.hist_b = pg.PlotCurveItem(pen=pg.mkPen('#4444FF', width=1.5), fillLevel=0,
                                        brush=(68, 68, 255, 80))
        self.hist_r.setVisible(False)
        self.hist_g.setVisible(False)
        self.hist_b.setVisible(False)
        self.hist_widget.addItem(self.hist_r)
        self.hist_widget.addItem(self.hist_g)
        self.hist_widget.addItem(self.hist_b)
        
        # Level lines
        self.lo_line = pg.InfiniteLine(angle=0, movable=True,
                                        pen=pg.mkPen('#00BFFF', width=2))
        self.hi_line = pg.InfiniteLine(angle=0, movable=True,
                                        pen=pg.mkPen('#FF6B6B', width=2))
        self.lo_line.sigPositionChanged.connect(self._on_lo_line_moved)
        self.hi_line.sigPositionChanged.connect(self._on_hi_line_moved)
        self.hist_widget.addItem(self.lo_line)
        self.hist_widget.addItem(self.hi_line)
        
        layout.addWidget(self.hist_widget, stretch=1)
        
        # === RGB MODE TOGGLE ===
        rgb_layout = QHBoxLayout()
        self.rgb_cb = QCheckBox("RGB Mode")
        self.rgb_cb.setChecked(False)
        self.rgb_cb.toggled.connect(self._on_rgb_toggled)
        rgb_layout.addWidget(self.rgb_cb)
        rgb_layout.addStretch()
        layout.addLayout(rgb_layout)
        
        # === LEVEL CONTROLS ===
        level_group = QGroupBox("Levels")
        level_layout = QGridLayout(level_group)
        level_layout.setContentsMargins(4, 4, 4, 4)
        
        # Lo level
        level_layout.addWidget(QLabel("Lo:"), 0, 0)
        self.lo_spin = QDoubleSpinBox()
        self.lo_spin.setRange(0, 1000000)
        self.lo_spin.setValue(0)
        self.lo_spin.valueChanged.connect(self._on_lo_spin_changed)
        level_layout.addWidget(self.lo_spin, 0, 1)
        
        # Hi level
        level_layout.addWidget(QLabel("Hi:"), 1, 0)
        self.hi_spin = QDoubleSpinBox()
        self.hi_spin.setRange(0, 1000000)
        self.hi_spin.setValue(65535)
        self.hi_spin.valueChanged.connect(self._on_hi_spin_changed)
        level_layout.addWidget(self.hi_spin, 1, 1)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.auto_btn = QPushButton("Auto")
        self.auto_btn.clicked.connect(self._on_auto_clicked)
        btn_layout.addWidget(self.auto_btn)
        
        self.reset_btn = QPushButton("Reset")
        self.reset_btn.clicked.connect(self._on_reset_clicked)
        btn_layout.addWidget(self.reset_btn)
        
        level_layout.addLayout(btn_layout, 2, 0, 1, 2)
        
        # Auto levels checkbox
        self.auto_cb = QCheckBox("Auto Levels")
        self.auto_cb.setChecked(True)
        level_layout.addWidget(self.auto_cb, 3, 0, 1, 2)
        
        layout.addWidget(level_group)
        
        # === STATISTICS ===
        stats_group = QGroupBox("Statistics")
        stats_layout = QFormLayout(stats_group)
        stats_layout.setContentsMargins(4, 4, 4, 4)
        
        self.mean_label = QLabel("--")
        stats_layout.addRow("Mean:", self.mean_label)
        
        self.std_label = QLabel("--")
        stats_layout.addRow("Std:", self.std_label)
        
        self.min_label = QLabel("--")
        stats_layout.addRow("Min:", self.min_label)
        
        self.max_label = QLabel("--")
        stats_layout.addRow("Max:", self.max_label)
        
        layout.addWidget(stats_group)
        
        # === COLORMAP ===
        cmap_group = QGroupBox("Display")
        cmap_layout = QVBoxLayout(cmap_group)
        cmap_layout.setContentsMargins(4, 4, 4, 4)
        
        cmap_row = QHBoxLayout()
        cmap_row.addWidget(QLabel("Colormap:"))
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(['gray', 'viridis', 'inferno', 'plasma', 
                                   'magma', 'hot', 'cool', 'jet'])
        self.cmap_combo.currentTextChanged.connect(self._on_cmap_changed)
        cmap_row.addWidget(self.cmap_combo)
        cmap_layout.addLayout(cmap_row)
        
        # Scale options
        scale_row = QHBoxLayout()
        self.linear_rb = QRadioButton("Linear")
        self.linear_rb.setChecked(True)
        self.log_rb = QRadioButton("Log")
        self.gamma_rb = QRadioButton("Gamma")
        scale_row.addWidget(self.linear_rb)
        scale_row.addWidget(self.log_rb)
        scale_row.addWidget(self.gamma_rb)
        cmap_layout.addLayout(scale_row)
        
        layout.addWidget(cmap_group)
        
        layout.addStretch()
        
        self.setWidget(widget)
        
        # State
        self._updating = False
        self._bit_depth = 16
        self._max_value = 65535
        self._is_rgb = False
        
    def update_histogram(self, counts: np.ndarray, edges: np.ndarray,
                         lo: float = None, hi: float = None):
        """Update histogram display (grayscale mode)."""
        if len(counts) == 0:
            return
            
        # For vertical orientation: Y is bin centers, X is counts
        bin_centers = (edges[:-1] + edges[1:]) / 2
        
        # Plot as step function (vertical bars)
        self.hist_item.setData(x=counts, y=bin_centers)
        
        # Update level lines
        self._updating = True
        if lo is not None:
            self.lo_line.setValue(lo)
            self.lo_spin.setValue(lo)
        if hi is not None:
            self.hi_line.setValue(hi)
            self.hi_spin.setValue(hi)
        self._updating = False
        
    def update_histogram_rgb(self, frame: np.ndarray):
        """
        Update RGB histogram from color frame.
        
        Args:
            frame: 3D array (H, W, 3) with RGB channels
        """
        if frame.ndim != 3 or frame.shape[2] < 3:
            return
            
        bins = 256
        
        # Determine range based on dtype
        if frame.dtype == np.uint8:
            range_ = (0, 255)
        elif frame.dtype == np.uint16:
            range_ = (0, 65535)
        else:
            range_ = (float(frame.min()), float(frame.max()))
        
        # Compute histograms for each channel
        counts_r, edges = np.histogram(frame[:, :, 0].ravel(), bins=bins, range=range_)
        counts_g, _ = np.histogram(frame[:, :, 1].ravel(), bins=bins, range=range_)
        counts_b, _ = np.histogram(frame[:, :, 2].ravel(), bins=bins, range=range_)
        
        bin_centers = (edges[:-1] + edges[1:]) / 2
        
        # Update RGB curves
        self.hist_r.setData(x=counts_r, y=bin_centers)
        self.hist_g.setData(x=counts_g, y=bin_centers)
        self.hist_b.setData(x=counts_b, y=bin_centers)
        
        # Hide grayscale histogram in RGB mode
        if self._is_rgb:
            self.hist_item.setVisible(False)
            self.hist_r.setVisible(True)
            self.hist_g.setVisible(True)
            self.hist_b.setVisible(True)
        
    def update_statistics(self, frame: np.ndarray):
        """Update statistics from frame."""
        if frame.ndim == 3:
            # Color frame - show luminance stats
            luminance = 0.299 * frame[:,:,0] + 0.587 * frame[:,:,1] + 0.114 * frame[:,:,2]
            self.mean_label.setText(f"{np.mean(luminance):.1f}")
            self.std_label.setText(f"{np.std(luminance):.1f}")
            self.min_label.setText(f"{np.min(luminance):.0f}")
            self.max_label.setText(f"{np.max(luminance):.0f}")
        else:
            self.mean_label.setText(f"{np.mean(frame):.1f}")
            self.std_label.setText(f"{np.std(frame):.1f}")
            self.min_label.setText(f"{np.min(frame)}")
            self.max_label.setText(f"{np.max(frame)}")
        
    def set_bit_depth(self, bits: int):
        """Set bit depth for max value."""
        self._bit_depth = bits
        self._max_value = 2**bits - 1
        self.hi_spin.setMaximum(self._max_value)
        self.lo_spin.setMaximum(self._max_value)
        
    def set_rgb_mode(self, enabled: bool):
        """Enable/disable RGB histogram mode."""
        self._is_rgb = enabled
        self.rgb_cb.setChecked(enabled)
        self._on_rgb_toggled(enabled)
        
    def get_levels(self) -> Tuple[float, float]:
        """Get current levels."""
        return (self.lo_spin.value(), self.hi_spin.value())
    
    def is_auto_levels(self) -> bool:
        """Check if auto levels is enabled."""
        return self.auto_cb.isChecked()
    
    def _on_rgb_toggled(self, enabled: bool):
        """Handle RGB mode toggle."""
        self._is_rgb = enabled
        self.hist_item.setVisible(not enabled)
        self.hist_r.setVisible(enabled)
        self.hist_g.setVisible(enabled)
        self.hist_b.setVisible(enabled)
        
        # Disable colormap selection in RGB mode
        self.cmap_combo.setEnabled(not enabled)
        
    def _on_lo_line_moved(self):
        if self._updating:
            return
        lo = self.lo_line.value()
        self._updating = True
        self.lo_spin.setValue(lo)
        self._updating = False
        self.auto_cb.setChecked(False)
        self.levelsChanged.emit(lo, self.hi_spin.value())
        
    def _on_hi_line_moved(self):
        if self._updating:
            return
        hi = self.hi_line.value()
        self._updating = True
        self.hi_spin.setValue(hi)
        self._updating = False
        self.auto_cb.setChecked(False)
        self.levelsChanged.emit(self.lo_spin.value(), hi)
        
    def _on_lo_spin_changed(self, value):
        if self._updating:
            return
        self._updating = True
        self.lo_line.setValue(value)
        self._updating = False
        self.auto_cb.setChecked(False)
        self.levelsChanged.emit(value, self.hi_spin.value())
        
    def _on_hi_spin_changed(self, value):
        if self._updating:
            return
        self._updating = True
        self.hi_line.setValue(value)
        self._updating = False
        self.auto_cb.setChecked(False)
        self.levelsChanged.emit(self.lo_spin.value(), value)
        
    def _on_auto_clicked(self):
        """Request auto levels calculation."""
        self.auto_cb.setChecked(True)
        
    def _on_reset_clicked(self):
        """Reset to full range."""
        self._updating = True
        self.lo_spin.setValue(0)
        self.hi_spin.setValue(self._max_value)
        self.lo_line.setValue(0)
        self.hi_line.setValue(self._max_value)
        self._updating = False
        self.levelsChanged.emit(0, self._max_value)
        
    def _on_cmap_changed(self, name):
        self.colormapChanged.emit(name)


# =============================================================================
# INTENSITY TRACE DOCK
# =============================================================================

class IntensityTraceDock(QDockWidget):
    """
    Rolling intensity time trace.
    
    Features:
    - Scrolling plot of intensity vs time
    - Adjustable time window
    - Clear button
    - Statistics
    """
    
    def __init__(self, parent=None):
        super().__init__("Intensity Trace", parent)
        self.setObjectName("IntensityTraceDock")
        
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(2, 2, 2, 2)
        
        # Plot widget
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('#1a1a1a')
        self.plot_widget.setLabel('left', 'Intensity')
        self.plot_widget.setLabel('bottom', 'Time', 's')
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        
        self.curve = self.plot_widget.plot(pen='#00BFFF')
        
        layout.addWidget(self.plot_widget, stretch=1)
        
        # Controls
        ctrl_layout = QHBoxLayout()
        
        ctrl_layout.addWidget(QLabel("Window:"))
        self.window_spin = QDoubleSpinBox()
        self.window_spin.setRange(1, 3600)
        self.window_spin.setValue(10)
        self.window_spin.setSuffix(" s")
        ctrl_layout.addWidget(self.window_spin)
        
        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self._on_clear)
        ctrl_layout.addWidget(self.clear_btn)
        
        ctrl_layout.addStretch()
        
        # Stats
        self.fps_label = QLabel("FPS: --")
        ctrl_layout.addWidget(self.fps_label)
        
        layout.addLayout(ctrl_layout)
        
        self.setWidget(widget)
        
        # State
        self._start_time = None
        
    def update_trace(self, times: np.ndarray, intensities: np.ndarray):
        """Update intensity trace."""
        if len(times) == 0:
            return
            
        # Convert to relative time if needed
        if self._start_time is None:
            self._start_time = times[0]
            
        rel_times = times - self._start_time
        
        self.curve.setData(rel_times, intensities)
        
        # Auto-scroll
        window = self.window_spin.value()
        if len(rel_times) > 0:
            t_max = rel_times[-1]
            t_min = max(0, t_max - window)
            self.plot_widget.setXRange(t_min, t_max, padding=0)
            
    def update_fps(self, fps: float):
        """Update FPS display."""
        self.fps_label.setText(f"FPS: {fps:.1f}")
        
    def _on_clear(self):
        """Clear trace."""
        self.curve.setData([], [])
        self._start_time = None


# =============================================================================
# IMAGE FFT DOCK
# =============================================================================

class ImageFFTDock(QDockWidget):
    """
    2D FFT display with filtering.
    
    Features:
    - 2D FFT magnitude display
    - Log/linear scale toggle
    - Filter mask drawing
    - Apply filter to get filtered image
    
    Signals:
        filteredImageReady: Filtered image available (np.ndarray)
        filterMaskChanged: Filter mask changed (np.ndarray)
    """
    
    filteredImageReady = Signal(object)
    filterMaskChanged = Signal(object)
    
    def __init__(self, parent=None):
        super().__init__("Image FFT", parent)
        self.setObjectName("ImageFFTDock")
        
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(2, 2, 2, 2)
        
        # FFT display
        self.fft_widget = pg.PlotWidget()
        self.fft_widget.setBackground('#1a1a1a')
        self.fft_widget.setAspectLocked(True)
        
        self.fft_image = pg.ImageItem()
        self.fft_widget.addItem(self.fft_image)
        
        layout.addWidget(self.fft_widget, stretch=1)
        
        # Controls
        ctrl_layout = QHBoxLayout()
        
        self.enable_cb = QCheckBox("Enable FFT")
        self.enable_cb.setChecked(False)
        ctrl_layout.addWidget(self.enable_cb)
        
        self.log_cb = QCheckBox("Log Scale")
        self.log_cb.setChecked(True)
        ctrl_layout.addWidget(self.log_cb)
        
        ctrl_layout.addWidget(QLabel("Window:"))
        self.window_combo = QComboBox()
        self.window_combo.addItems(['None', 'Hann', 'Hamming', 'Blackman'])
        self.window_combo.setCurrentText('Hann')
        ctrl_layout.addWidget(self.window_combo)
        
        ctrl_layout.addStretch()
        
        layout.addLayout(ctrl_layout)
        
        # Filter controls
        filter_layout = QHBoxLayout()
        
        filter_layout.addWidget(QLabel("Filter:"))
        self.filter_combo = QComboBox()
        self.filter_combo.addItems(['None', 'Lowpass', 'Highpass', 'Bandpass'])
        filter_layout.addWidget(self.filter_combo)
        
        filter_layout.addWidget(QLabel("Radius:"))
        self.radius_spin = QSpinBox()
        self.radius_spin.setRange(1, 1000)
        self.radius_spin.setValue(50)
        filter_layout.addWidget(self.radius_spin)
        
        self.apply_btn = QPushButton("Apply Filter")
        self.apply_btn.clicked.connect(self._on_apply_filter)
        filter_layout.addWidget(self.apply_btn)
        
        filter_layout.addStretch()
        
        layout.addLayout(filter_layout)
        
        self.setWidget(widget)
        
        # State
        self._fft_processor = None
        self._current_fft = None
        
    def set_fft_processor(self, processor):
        """Set FFT processor reference."""
        self._fft_processor = processor
        
    def update_fft(self, magnitude: np.ndarray):
        """Update FFT display."""
        self.fft_image.setImage(magnitude)
        
    def is_enabled(self) -> bool:
        """Check if FFT is enabled."""
        return self.enable_cb.isChecked()
        
    def _on_apply_filter(self):
        """Apply current filter settings."""
        if self._fft_processor is None:
            return
            
        filter_type = self.filter_combo.currentText().lower()
        radius = self.radius_spin.value()
        
        if filter_type == 'none':
            return
            
        # Get FFT shape from processor
        if self._fft_processor._last_magnitude is None:
            return
            
        shape = self._fft_processor._last_magnitude.shape
        
        # Create mask
        if filter_type == 'lowpass':
            mask = self._fft_processor.create_circular_mask(shape, radius=radius)
        elif filter_type == 'highpass':
            mask = self._fft_processor.create_circular_mask(shape, radius=radius, invert=True)
        elif filter_type == 'bandpass':
            inner = radius // 2
            mask = self._fft_processor.create_bandpass_mask(shape, inner, radius)
        else:
            return
            
        # Apply filter
        filtered = self._fft_processor.apply_filter(mask)
        if filtered is not None:
            self.filteredImageReady.emit(filtered)
            self.filterMaskChanged.emit(mask)


# =============================================================================
# CONTROL DOCK (Tabbed Interface)
# =============================================================================

class ControlDock(QDockWidget):
    """
    Main control panel with tabs.
    
    Tabs:
    - Capture: Live view, exposure, gain, binning, trigger
    - Sequence: Time lapse, high-speed, autosave, RAM management
    - Devices: Camera selection, info, EMCCD controls
    
    Signals:
        startClicked, stopClicked: Acquisition control
        exposureChanged, gainChanged, etc.: Parameter changes
        cameraSelected: Camera selection changed
        emccdGainChanged: EMCCD gain changed (for Andor)
        coolingChanged: Cooling target changed
    """
    
    startClicked = Signal()
    stopClicked = Signal()
    recordClicked = Signal(str, str)  # storage_mode ('none','ram','hdf5','tiff'), save_path
    exposureChanged = Signal(float)
    gainChanged = Signal(float)
    binningChanged = Signal(int)
    triggerModeChanged = Signal(object)
    roiChanged = Signal(int, int, int, int)
    cameraSelected = Signal(str, int)  # type, index
    emccdGainChanged = Signal(int)
    coolingChanged = Signal(bool, float)  # enabled, target_temp
    displayRateChanged = Signal(float)
    
    def __init__(self, parent=None):
        super().__init__("Controls", parent)
        self.setObjectName("ControlDock")
        
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(4, 4, 4, 4)
        
        # Tab widget
        self.tabs = QTabWidget()
        
        # === CAPTURE TAB ===
        capture_tab = QWidget()
        capture_layout = QVBoxLayout(capture_tab)
        
        # Live button
        self.live_btn = QPushButton("LIVE")
        self.live_btn.setCheckable(True)
        self.live_btn.setStyleSheet("""
            QPushButton { 
                font-size: 16px; font-weight: bold; 
                padding: 10px; background: #2d5a27; 
            }
            QPushButton:checked { background: #8b0000; }
        """)
        self.live_btn.toggled.connect(self._on_live_toggled)
        capture_layout.addWidget(self.live_btn)
        
        # Exposure
        exp_group = QGroupBox("Exposure")
        exp_layout = QFormLayout(exp_group)
        
        self.exposure_spin = QDoubleSpinBox()
        self.exposure_spin.setRange(0.001, 10000)
        self.exposure_spin.setValue(10)
        self.exposure_spin.setSuffix(" ms")
        self.exposure_spin.valueChanged.connect(lambda v: self.exposureChanged.emit(v))
        exp_layout.addRow("Exposure:", self.exposure_spin)
        
        self.gain_spin = QDoubleSpinBox()
        self.gain_spin.setRange(0, 100)
        self.gain_spin.setValue(1)
        self.gain_spin.valueChanged.connect(lambda v: self.gainChanged.emit(v))
        exp_layout.addRow("Gain:", self.gain_spin)
        
        capture_layout.addWidget(exp_group)
        
        # Binning
        bin_group = QGroupBox("Binning")
        bin_layout = QHBoxLayout(bin_group)
        
        self.bin_group = QButtonGroup()
        for i, b in enumerate([1, 2, 4]):
            rb = QRadioButton(f"{b}x{b}")
            if b == 1:
                rb.setChecked(True)
            self.bin_group.addButton(rb, b)
            bin_layout.addWidget(rb)
        self.bin_group.idClicked.connect(lambda b: self.binningChanged.emit(b))
        
        capture_layout.addWidget(bin_group)
        
        # Trigger
        trig_group = QGroupBox("Trigger")
        trig_layout = QFormLayout(trig_group)
        
        self.trigger_combo = QComboBox()
        self.trigger_combo.addItems(['Internal', 'External', 'Software'])
        self.trigger_combo.currentTextChanged.connect(self._on_trigger_changed)
        trig_layout.addRow("Mode:", self.trigger_combo)
        
        capture_layout.addWidget(trig_group)
        
        # Display rate control
        disp_group = QGroupBox("Display")
        disp_layout = QFormLayout(disp_group)
        
        self.disp_rate_spin = QDoubleSpinBox()
        self.disp_rate_spin.setRange(1, 120)
        self.disp_rate_spin.setValue(30)
        self.disp_rate_spin.setSuffix(" Hz")
        self.disp_rate_spin.valueChanged.connect(lambda v: self.displayRateChanged.emit(v))
        disp_layout.addRow("Rate:", self.disp_rate_spin)
        
        capture_layout.addWidget(disp_group)
        
        capture_layout.addStretch()
        self.tabs.addTab(capture_tab, "Capture")
        
        # === SEQUENCE TAB ===
        sequence_tab = QWidget()
        seq_layout = QVBoxLayout(sequence_tab)
        
        # Sequence type
        type_group = QGroupBox("Sequence Type")
        type_layout = QVBoxLayout(type_group)
        
        self.timelapse_rb = QRadioButton("Time Lapse")
        self.timelapse_rb.setChecked(True)
        type_layout.addWidget(self.timelapse_rb)
        
        self.highspeed_rb = QRadioButton("High Speed Streaming")
        type_layout.addWidget(self.highspeed_rb)
        
        seq_layout.addWidget(type_group)
        
        # Progress
        prog_group = QGroupBox("Progress")
        prog_layout = QFormLayout(prog_group)
        
        self.frame_count_label = QLabel("0")
        prog_layout.addRow("Frames:", self.frame_count_label)
        
        self.fps_label = QLabel("0.0")
        prog_layout.addRow("FPS:", self.fps_label)
        
        self.elapsed_label = QLabel("00:00:00")
        prog_layout.addRow("Elapsed:", self.elapsed_label)
        
        seq_layout.addWidget(prog_group)
        
        # Storage
        storage_group = QGroupBox("Storage")
        storage_layout = QVBoxLayout(storage_group)
        
        self.storage_none_rb = QRadioButton("None (display only)")
        self.storage_none_rb.setChecked(True)
        storage_layout.addWidget(self.storage_none_rb)
        
        self.storage_ram_rb = QRadioButton("To RAM")
        storage_layout.addWidget(self.storage_ram_rb)
        
        self.storage_hdf5_rb = QRadioButton("To Disk (HDF5)")
        storage_layout.addWidget(self.storage_hdf5_rb)
        
        self.storage_tiff_rb = QRadioButton("To Disk (TIFF)")
        storage_layout.addWidget(self.storage_tiff_rb)
        
        # RAM ceiling
        ram_layout = QHBoxLayout()
        ram_layout.addWidget(QLabel("RAM Ceiling:"))
        self.ram_spin = QDoubleSpinBox()
        self.ram_spin.setRange(100, 100000)
        self.ram_spin.setValue(2000)
        self.ram_spin.setSuffix(" MB")
        ram_layout.addWidget(self.ram_spin)
        storage_layout.addLayout(ram_layout)
        
        # Save path
        path_layout = QHBoxLayout()
        path_layout.addWidget(QLabel("Path:"))
        self.save_path_edit = QLineEdit()
        self.save_path_edit.setPlaceholderText("Select save location...")
        path_layout.addWidget(self.save_path_edit)
        self.browse_btn = QPushButton("...")
        self.browse_btn.setMaximumWidth(30)
        self.browse_btn.clicked.connect(self._on_browse_save_path)
        path_layout.addWidget(self.browse_btn)
        storage_layout.addLayout(path_layout)
        
        # Record button
        self.record_btn = QPushButton("RECORD")
        self.record_btn.setCheckable(True)
        self.record_btn.setStyleSheet("""
            QPushButton { 
                font-size: 14px; font-weight: bold;
                padding: 8px; background: #8b0000;
            }
            QPushButton:checked { background: #cc0000; }
            QPushButton:disabled { background: #3c3c3c; color: #666; }
        """)
        self.record_btn.toggled.connect(self._on_record_toggled)
        storage_layout.addWidget(self.record_btn)
        
        seq_layout.addWidget(storage_group)
        
        seq_layout.addStretch()
        self.tabs.addTab(sequence_tab, "Sequence")
        
        # === DEVICES TAB ===
        devices_tab = QWidget()
        dev_layout = QVBoxLayout(devices_tab)
        
        # Camera selection
        cam_group = QGroupBox("Camera")
        cam_layout = QFormLayout(cam_group)
        
        self.camera_type_combo = QComboBox()
        self.camera_type_combo.addItems(['Simulated', 'Hamamatsu', 'Thorlabs', 'UC480', 'Andor'])
        cam_layout.addRow("Type:", self.camera_type_combo)
        
        self.camera_index_spin = QSpinBox()
        self.camera_index_spin.setRange(0, 10)
        cam_layout.addRow("Index:", self.camera_index_spin)
        
        self.connect_btn = QPushButton("Connect")
        self.connect_btn.clicked.connect(self._on_connect_clicked)
        cam_layout.addRow(self.connect_btn)
        
        dev_layout.addWidget(cam_group)
        
        # Camera info
        info_group = QGroupBox("Camera Info")
        info_layout = QFormLayout(info_group)
        
        self.model_label = QLabel("--")
        info_layout.addRow("Model:", self.model_label)
        
        self.serial_label = QLabel("--")
        info_layout.addRow("Serial:", self.serial_label)
        
        self.sensor_label = QLabel("--")
        info_layout.addRow("Sensor:", self.sensor_label)
        
        self.temp_label = QLabel("--")
        info_layout.addRow("Temp:", self.temp_label)
        
        dev_layout.addWidget(info_group)
        
        # === EMCCD / Cooling Controls (for Andor) ===
        self.emccd_group = QGroupBox("EMCCD / Cooling")
        emccd_layout = QFormLayout(self.emccd_group)
        
        self.emccd_gain_spin = QSpinBox()
        self.emccd_gain_spin.setRange(1, 300)
        self.emccd_gain_spin.setValue(1)
        self.emccd_gain_spin.valueChanged.connect(lambda v: self.emccdGainChanged.emit(v))
        emccd_layout.addRow("EM Gain:", self.emccd_gain_spin)
        
        self.cooling_cb = QCheckBox("Enable Cooling")
        self.cooling_cb.setChecked(True)
        emccd_layout.addRow(self.cooling_cb)
        
        self.cooling_temp_spin = QDoubleSpinBox()
        self.cooling_temp_spin.setRange(-100, 20)
        self.cooling_temp_spin.setValue(-70)
        self.cooling_temp_spin.setSuffix(" °C")
        self.cooling_temp_spin.valueChanged.connect(self._on_cooling_changed)
        emccd_layout.addRow("Target:", self.cooling_temp_spin)
        
        self.current_temp_label = QLabel("-- °C")
        self.current_temp_label.setStyleSheet("font-weight: bold;")
        emccd_layout.addRow("Current:", self.current_temp_label)
        
        self.emccd_group.setVisible(False)  # Hidden by default, shown for Andor
        dev_layout.addWidget(self.emccd_group)
        
        dev_layout.addStretch()
        self.tabs.addTab(devices_tab, "Devices")
        
        layout.addWidget(self.tabs)
        
        self.setWidget(widget)
        
        # Temperature update timer
        self._temp_timer = QTimer(self)
        self._temp_timer.timeout.connect(self._update_temperature)
        self._camera_ref = None
        
    def _on_live_toggled(self, checked):
        if checked:
            self.live_btn.setText("STOP")
            self.startClicked.emit()
        else:
            self.live_btn.setText("LIVE")
            self.stopClicked.emit()
            
    def _on_trigger_changed(self, text):
        mode_map = {
            'Internal': TriggerMode.INTERNAL,
            'External': TriggerMode.EXTERNAL,
            'Software': TriggerMode.SOFTWARE
        }
        self.triggerModeChanged.emit(mode_map.get(text, TriggerMode.INTERNAL))
        
    def _on_connect_clicked(self):
        cam_type = self.camera_type_combo.currentText().lower()
        cam_idx = self.camera_index_spin.value()
        self.cameraSelected.emit(cam_type, cam_idx)
        
    def _on_cooling_changed(self):
        enabled = self.cooling_cb.isChecked()
        temp = self.cooling_temp_spin.value()
        self.coolingChanged.emit(enabled, temp)
    
    def _on_browse_save_path(self):
        """Browse for save path."""
        from PySide6.QtWidgets import QFileDialog
        if self.storage_tiff_rb.isChecked():
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Location", "recording.tif",
                "TIFF Files (*.tif *.tiff)")
        elif self.storage_hdf5_rb.isChecked():
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Location", "recording.h5",
                "HDF5 Files (*.h5 *.hdf5)")
        else:
            path = QFileDialog.getExistingDirectory(
                self, "Save Location")
        if path:
            self.save_path_edit.setText(path)
    
    def _on_record_toggled(self, checked: bool):
        """Handle record button toggle."""
        if checked:
            # Determine storage mode
            if self.storage_hdf5_rb.isChecked():
                mode = 'hdf5'
            elif self.storage_tiff_rb.isChecked():
                mode = 'tiff'
            elif self.storage_ram_rb.isChecked():
                mode = 'ram'
            else:
                mode = 'none'
            
            path = self.save_path_edit.text()
            if mode in ('hdf5', 'tiff') and not path:
                self._on_browse_save_path()
                path = self.save_path_edit.text()
                if not path:
                    self.record_btn.setChecked(False)
                    return
            
            self.record_btn.setText("STOP REC")
            self.recordClicked.emit(mode, path)
        else:
            self.record_btn.setText("RECORD")
            self.recordClicked.emit('none', '')
        
    def update_camera_info(self, info: Dict):
        """Update camera info display."""
        self.model_label.setText(str(info.get('model', '--')))
        self.serial_label.setText(str(info.get('serial', '--')))
        sensor = info.get('sensor_size', ('--', '--'))
        self.sensor_label.setText(f"{sensor[0]} x {sensor[1]}")
        temp = info.get('temperature')
        self.temp_label.setText(f"{temp:.1f} °C" if temp else "--")
        
        # Show EMCCD controls for Andor cameras
        is_andor = 'andor' in str(info.get('model', '')).lower() or \
                   'ixon' in str(info.get('model', '')).lower()
        self.emccd_group.setVisible(is_andor)
        
        # Update EMCCD gain if available
        emccd_gain = info.get('emccd_gain')
        if emccd_gain is not None:
            self.emccd_gain_spin.setValue(emccd_gain)
            
    def set_camera_reference(self, camera):
        """Set camera reference for temperature updates."""
        self._camera_ref = camera
        if camera is not None and hasattr(camera, 'get_temperature'):
            self._temp_timer.start(2000)  # Update every 2 seconds
        else:
            self._temp_timer.stop()
            
    def _update_temperature(self):
        """Update temperature display."""
        if self._camera_ref is not None:
            temp = self._camera_ref.get_temperature()
            if temp is not None:
                self.current_temp_label.setText(f"{temp:.1f} °C")
                self.temp_label.setText(f"{temp:.1f} °C")
        
    def update_progress(self, frame_count: int, fps: float, elapsed: float):
        """Update progress display."""
        self.frame_count_label.setText(str(frame_count))
        self.fps_label.setText(f"{fps:.1f}")
        
        hours = int(elapsed // 3600)
        mins = int((elapsed % 3600) // 60)
        secs = int(elapsed % 60)
        self.elapsed_label.setText(f"{hours:02d}:{mins:02d}:{secs:02d}")
        
    def set_live_state(self, running: bool):
        """Set live button state without triggering signal."""
        self.live_btn.blockSignals(True)
        self.live_btn.setChecked(running)
        self.live_btn.setText("STOP" if running else "LIVE")
        self.live_btn.blockSignals(False)
