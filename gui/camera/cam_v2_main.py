"""
Camera Viewer v2 - Main Window
===============================

Main application window with dockable interface.

Usage:
    from cam_v2_main import CameraViewerV2
    app = QApplication(sys.argv)
    viewer = CameraViewerV2()
    viewer.show()
    app.exec()

Author: Claude/BP
Date: January 2025
"""

import sys
import numpy as np
import logging
import json
from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QMenuBar, QMenu, QStatusBar, QLabel, QMessageBox, QFileDialog, QDockWidget
)
from PySide6.QtCore import Qt, QTimer, Slot, QSettings, QByteArray
from PySide6.QtGui import QAction, QKeySequence, QShortcut, QPalette, QColor

from cam_v2_core import (
    CameraType, TriggerMode, StorageMode, IntensityMode,
    CameraConfig, AcquisitionConfig, DisplayConfig, StorageConfig,
    LatestFrameBuffer, IntensityRingBuffer, FrameRingBuffer,
    auto_levels, compute_histogram
)
from cam_v2_backends import create_camera, list_available_cameras
from cam_v2_acquisition import (
    AcquisitionThread, DisplayManager, DisplayManagerConfig,
    ImageFFTProcessor, create_acquisition_system, HDF5Writer, TIFFWriter
)
from cam_v2_widgets import (
    ImageDisplayDock, HistogramDock, IntensityTraceDock,
    ImageFFTDock, ControlDock
)


class CameraViewerV2(QMainWindow):
    """
    Camera Viewer v2 Main Window.
    
    Features:
    - Multi-camera backend support
    - Decoupled acquisition and display rates
    - Live image FFT with filtering
    - Intensity time trace
    - Histogram with level controls
    - Dockable widget interface
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.setWindowTitle("Camera Viewer v2")
        self.setMinimumSize(1200, 800)
        
        # Camera and acquisition
        self.camera = None
        self.acq_thread: Optional[AcquisitionThread] = None
        self.display_manager: Optional[DisplayManager] = None
        self.frame_buffer: Optional[LatestFrameBuffer] = None
        self.intensity_buffer: Optional[IntensityRingBuffer] = None
        self.storage_buffer: Optional[FrameRingBuffer] = None
        
        # FFT processor
        self.fft_processor = ImageFFTProcessor()
        self._fft_enabled = False
        
        # Configs
        self.camera_config = CameraConfig()
        self.acquisition_config = AcquisitionConfig()
        self.display_config = DisplayConfig()
        self.storage_config = StorageConfig()
        
        # Setup UI
        self._setup_ui()
        self._setup_menus()
        self._setup_status_bar()
        self._connect_signals()
        
        # Initialize with simulated camera
        self._connect_camera(CameraType.SIMULATED, 0)
        
        # Apply dark theme
        self._apply_dark_theme()
        
        # Setup keyboard shortcuts
        self._setup_shortcuts()
        
        # Load settings (after UI setup)
        self._load_settings()
        
    def _setup_shortcuts(self):
        """Setup keyboard shortcuts."""
        # Space - Start/Stop
        QShortcut(QKeySequence(Qt.Key.Key_Space), self, self._toggle_acquisition)
        
        # S - Save frame
        QShortcut(QKeySequence("Ctrl+S"), self, self._on_save_frame)
        
        # A - Auto levels
        QShortcut(QKeySequence("A"), self, self._trigger_auto_levels)
        
        # F - Toggle FFT
        QShortcut(QKeySequence("F"), self, self._toggle_fft)
        
        # C - Toggle crosshair
        QShortcut(QKeySequence("C"), self, lambda: self.image_dock.crosshair_cb.toggle())
        
        # R - Reset view
        QShortcut(QKeySequence("R"), self, lambda: self.image_dock.image_widget.autoRange())
        
        # Escape - Stop acquisition
        QShortcut(QKeySequence(Qt.Key.Key_Escape), self, self._stop_acquisition)
        
    def _toggle_acquisition(self):
        """Toggle acquisition on/off."""
        if self.acq_thread and self.acq_thread.isRunning():
            self._stop_acquisition()
        else:
            self._start_acquisition()
            self.control_dock.set_live_state(True)
            
    def _trigger_auto_levels(self):
        """Trigger auto levels calculation."""
        self.histogram_dock.auto_cb.setChecked(True)
        
    def _toggle_fft(self):
        """Toggle FFT display."""
        self.fft_dock.enable_cb.toggle()
        
    def _setup_ui(self):
        """Setup dockable widget interface."""
        # Central widget (empty - all content in docks)
        central = QWidget()
        central.setMaximumWidth(0)
        self.setCentralWidget(central)
        
        # Enable dock nesting
        self.setDockNestingEnabled(True)
        
        # === CREATE DOCKS ===
        # Layout: Controls (left) | Image (center) | Histogram (right)
        #         Intensity Trace / FFT (bottom, tabbed)
        
        # Controls (LEFT - narrow)
        self.control_dock = ControlDock(self)
        self.control_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable |
            QDockWidget.DockWidgetFeature.DockWidgetFloatable |
            QDockWidget.DockWidgetFeature.DockWidgetClosable
        )
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.control_dock)
        
        # Image Display (CENTER - wide)
        self.image_dock = ImageDisplayDock(self)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.image_dock)
        
        # Histogram (RIGHT)
        self.histogram_dock = HistogramDock(self)
        self.histogram_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable |
            QDockWidget.DockWidgetFeature.DockWidgetFloatable |
            QDockWidget.DockWidgetFeature.DockWidgetClosable
        )
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.histogram_dock)
        
        # Intensity Trace (BOTTOM)
        self.intensity_dock = IntensityTraceDock(self)
        self.intensity_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable |
            QDockWidget.DockWidgetFeature.DockWidgetFloatable |
            QDockWidget.DockWidgetFeature.DockWidgetClosable
        )
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.intensity_dock)
        
        # FFT (BOTTOM, tabbed with intensity)
        self.fft_dock = ImageFFTDock(self)
        self.fft_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable |
            QDockWidget.DockWidgetFeature.DockWidgetFloatable |
            QDockWidget.DockWidgetFeature.DockWidgetClosable
        )
        self.tabifyDockWidget(self.intensity_dock, self.fft_dock)
        self.intensity_dock.raise_()  # Show intensity by default
        
        # Set FFT processor reference
        self.fft_dock.set_fft_processor(self.fft_processor)
        
        # Set initial dock sizes: Controls narrow, Image wide
        # This needs to be called after show() in practice, but we try here
        self.resizeDocks(
            [self.control_dock, self.image_dock], 
            [250, 800], 
            Qt.Orientation.Horizontal
        )
        
    def _setup_menus(self):
        """Setup menu bar."""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        
        save_action = QAction("&Save Frame...", self)
        save_action.setShortcut(QKeySequence.StandardKey.Save)
        save_action.triggered.connect(self._on_save_frame)
        file_menu.addAction(save_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut(QKeySequence.StandardKey.Quit)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # View menu
        view_menu = menubar.addMenu("&View")
        
        view_menu.addAction(self.image_dock.toggleViewAction())
        view_menu.addAction(self.histogram_dock.toggleViewAction())
        view_menu.addAction(self.control_dock.toggleViewAction())
        view_menu.addAction(self.intensity_dock.toggleViewAction())
        view_menu.addAction(self.fft_dock.toggleViewAction())
        
        view_menu.addSeparator()
        
        reset_layout_action = QAction("Reset Layout", self)
        reset_layout_action.triggered.connect(self._reset_layout)
        view_menu.addAction(reset_layout_action)
        
        # Camera menu
        camera_menu = menubar.addMenu("&Camera")
        
        # Add camera type submenu
        for cam_type in CameraType:
            action = QAction(cam_type.value.capitalize(), self)
            action.triggered.connect(
                lambda checked, ct=cam_type: self._connect_camera(ct, 0)
            )
            camera_menu.addAction(action)
        
        camera_menu.addSeparator()
        
        disconnect_action = QAction("&Disconnect Camera", self)
        disconnect_action.triggered.connect(self._disconnect_camera)
        camera_menu.addAction(disconnect_action)
        
        # Help menu
        help_menu = menubar.addMenu("&Help")
        
        about_action = QAction("&About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)
        
    def _setup_status_bar(self):
        """Setup status bar."""
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        
        # Acquisition rate
        self.acq_rate_label = QLabel("Acq: -- Hz")
        self.status_bar.addPermanentWidget(self.acq_rate_label)
        
        # Display rate
        self.disp_rate_label = QLabel("Disp: -- Hz")
        self.status_bar.addPermanentWidget(self.disp_rate_label)
        
        # Frame count
        self.frame_count_label = QLabel("Frames: 0")
        self.status_bar.addPermanentWidget(self.frame_count_label)
        
        # Camera status
        self.camera_label = QLabel("Camera: --")
        self.status_bar.addWidget(self.camera_label)
        
    def _connect_signals(self):
        """Connect widget signals."""
        # Control dock
        self.control_dock.startClicked.connect(self._start_acquisition)
        self.control_dock.stopClicked.connect(self._stop_acquisition)
        self.control_dock.exposureChanged.connect(self._on_exposure_changed)
        self.control_dock.gainChanged.connect(self._on_gain_changed)
        self.control_dock.binningChanged.connect(self._on_binning_changed)
        self.control_dock.triggerModeChanged.connect(self._on_trigger_changed)
        self.control_dock.cameraSelected.connect(self._on_camera_selected)
        self.control_dock.emccdGainChanged.connect(self._on_emccd_gain_changed)
        self.control_dock.coolingChanged.connect(self._on_cooling_changed)
        self.control_dock.displayRateChanged.connect(self._on_display_rate_changed)
        self.control_dock.recordClicked.connect(self._on_record_clicked)
        
        # Histogram dock
        self.histogram_dock.levelsChanged.connect(self._on_levels_changed)
        self.histogram_dock.colormapChanged.connect(self._on_colormap_changed)
        
        # Image dock
        self.image_dock.roiChanged.connect(self._on_roi_changed)
        self.image_dock.roiApplyToSensor.connect(self._on_roi_apply_to_sensor)
        
        # FFT dock
        self.fft_dock.filteredImageReady.connect(self._on_filtered_image)
        self.fft_dock.enable_cb.toggled.connect(self._on_fft_toggled)
        
    def _connect_camera(self, cam_type: CameraType, index: int):
        """Connect to camera. Guards against re-initializing the same camera."""
        # Guard: if same camera type and already connected, skip
        if (self.camera is not None and self.camera.is_open and
                self.camera_config.camera_type == cam_type and
                self.camera_config.device_index == index):
            self.status_bar.showMessage(
                f"{cam_type.value.capitalize()} already connected", 3000)
            return
        
        # Stop any running acquisition
        self._stop_acquisition()
        
        # Close existing camera
        if self.camera is not None:
            self.camera.close()
            self.camera = None
        
        # Update config
        self.camera_config.camera_type = cam_type
        self.camera_config.device_index = index
        
        # Create camera
        try:
            self.camera = create_camera(cam_type, device_index=index)
            
            if not self.camera.open():
                raise RuntimeError("Failed to open camera")
            
            # Update UI
            info = self.camera.camera_info
            self.control_dock.update_camera_info(info)
            self.camera_label.setText(f"Camera: {info.get('model', 'Unknown')}")
            
            # Set bit depth for histogram
            bit_depth = info.get('bit_depth', 16)
            self.histogram_dock.set_bit_depth(bit_depth)
            
            # Create acquisition system
            self._setup_acquisition_system()
            
            self.status_bar.showMessage(f"Connected to {cam_type.value} camera", 3000)
            logging.info(f"Connected to {cam_type.value} camera: {info}")
            
            # Set camera reference for temperature updates (Andor)
            self.control_dock.set_camera_reference(self.camera)
            
        except Exception as e:
            logging.exception(f"Failed to connect to camera: {e}")
            QMessageBox.critical(self, "Camera Error", 
                                 f"Failed to connect to camera:\n{e}")
            self.camera_label.setText("Camera: Not connected")
            
    def _disconnect_camera(self):
        """Disconnect current camera and clean up."""
        self._stop_acquisition()
        
        if self.camera is not None:
            self.camera.close()
            self.camera = None
            
        self.acq_thread = None
        self.display_manager = None
        self.frame_buffer = None
        self.intensity_buffer = None
        
        self.control_dock.set_camera_reference(None)
        self.camera_label.setText("Camera: Disconnected")
        self.status_bar.showMessage("Camera disconnected", 3000)
        logging.info("Camera disconnected")
            
    def _setup_acquisition_system(self):
        """Setup acquisition thread and display manager."""
        if self.camera is None:
            return
            
        # Create buffers
        self.frame_buffer = LatestFrameBuffer()
        self.intensity_buffer = IntensityRingBuffer(
            capacity=self.acquisition_config.intensity_buffer_size
        )
        
        # Create acquisition thread
        self.acq_thread = AcquisitionThread(
            camera=self.camera,
            frame_buffer=self.frame_buffer,
            intensity_buffer=self.intensity_buffer,
            config=self.acquisition_config
        )
        
        # Connect acquisition signals
        self.acq_thread.frameAcquired.connect(self._on_frame_acquired)
        self.acq_thread.errorOccurred.connect(self._on_acq_error)
        self.acq_thread.statusChanged.connect(self._on_acq_status)
        self.acq_thread.acquisitionStats.connect(self._on_acq_stats)
        
        # Create display manager
        disp_config = DisplayManagerConfig(
            display_rate_hz=self.acquisition_config.display_rate_hz,
            fft_enabled=self.fft_dock.is_enabled()
        )
        
        self.display_manager = DisplayManager(
            frame_buffer=self.frame_buffer,
            intensity_buffer=self.intensity_buffer,
            config=disp_config
        )
        
        # Connect display signals
        self.display_manager.imageReady.connect(self._on_image_ready)
        self.display_manager.histogramReady.connect(self._on_histogram_ready)
        self.display_manager.intensityTraceReady.connect(self._on_intensity_ready)
        self.display_manager.fftReady.connect(self._on_fft_ready)
        self.display_manager.statsReady.connect(self._on_display_stats)
        self.display_manager.levelsComputed.connect(self._on_auto_levels)
        
    # === ACQUISITION CONTROL ===
    
    @Slot()
    def _start_acquisition(self):
        """Start acquisition."""
        if self.camera is None or self.acq_thread is None:
            return
            
        # Apply current settings to camera
        self.camera.set_exposure(self.camera_config.exposure_ms)
        if hasattr(self.camera, 'set_gain'):
            self.camera.set_gain(self.camera_config.gain)
        self.camera.set_binning(self.camera_config.binning)
        self.camera.set_trigger_mode(self.camera_config.trigger_mode)
        
        # Start acquisition thread
        self.acq_thread.start()
        
        # Start display manager
        self.display_manager.start()
        
        self.status_bar.showMessage("Acquisition started", 2000)
        
    @Slot()
    def _stop_acquisition(self):
        """Stop acquisition."""
        if self.acq_thread is not None and self.acq_thread.isRunning():
            self.acq_thread.request_stop()
            self.acq_thread.wait(2000)  # Wait up to 2 seconds
            
            # Close disk writer if active
            if self.acq_thread.disk_writer is not None:
                try:
                    self.acq_thread.disk_writer.close()
                except Exception:
                    pass
                self.acq_thread.disk_writer = None
            
        if self.display_manager is not None:
            self.display_manager.stop()
            
        if self.camera is not None:
            self.camera.stop_acquisition()
            
        self.control_dock.set_live_state(False)
        self.control_dock.record_btn.setChecked(False)
        self.status_bar.showMessage("Acquisition stopped", 2000)
        
    # === SIGNAL HANDLERS ===
    
    @Slot(int, float, float)
    def _on_frame_acquired(self, frame_num: int, timestamp: float, intensity: float):
        """Handle frame acquired signal (called for every frame)."""
        pass  # Stats handled by _on_acq_stats
        
    @Slot(object, float)
    def _on_image_ready(self, frame: np.ndarray, timestamp: float):
        """Handle image ready from display manager."""
        self.image_dock.set_image(frame)
        
        # Statistics use the (possibly downsampled) display frame
        self.histogram_dock.update_statistics(frame)
            
    @Slot(object, object, float, float)
    def _on_histogram_ready(self, counts: np.ndarray, edges: np.ndarray,
                            lo: float, hi: float):
        """Handle histogram data ready."""
        self.histogram_dock.update_histogram(counts, edges, lo, hi)
        
        # When fast_8bit_display is on, the DisplayManager already scales
        # the image to 0-255 using these levels, so we set image levels
        # to the full 8-bit range. Otherwise pass through raw levels.
        disp_cfg = self.display_manager.config if self.display_manager else None
        if disp_cfg and disp_cfg.fast_8bit_display:
            self.image_dock.set_levels(0, 255)
        else:
            self.image_dock.set_levels(lo, hi)
        
    @Slot(object, object)
    def _on_intensity_ready(self, times: np.ndarray, intensities: np.ndarray):
        """Handle intensity trace ready."""
        self.intensity_dock.update_trace(times, intensities)
        
    @Slot(object)
    def _on_fft_ready(self, magnitude: np.ndarray):
        """Handle FFT magnitude ready."""
        self.fft_dock.update_fft(magnitude)
        
    @Slot(float, float)
    def _on_auto_levels(self, lo: float, hi: float):
        """Handle auto levels computed."""
        if self.histogram_dock.is_auto_levels():
            disp_cfg = self.display_manager.config if self.display_manager else None
            if disp_cfg and disp_cfg.fast_8bit_display:
                # DisplayManager already applied levels during 8-bit conversion
                self.image_dock.set_levels(0, 255)
            else:
                self.image_dock.set_levels(lo, hi)
            
    @Slot(dict)
    def _on_acq_stats(self, stats: dict):
        """Handle acquisition statistics."""
        fps = stats.get('fps', 0)
        count = stats.get('frame_count', 0)
        
        self.acq_rate_label.setText(f"Acq: {fps:.1f} Hz")
        self.frame_count_label.setText(f"Frames: {count}")
        
        self.control_dock.update_progress(count, fps, stats.get('elapsed_time', 0))
        self.intensity_dock.update_fps(fps)
        
    @Slot(dict)
    def _on_display_stats(self, stats: dict):
        """Handle display statistics."""
        rate = stats.get('display_rate', 0)
        self.disp_rate_label.setText(f"Disp: {rate:.1f} Hz")
        
    @Slot(str)
    def _on_acq_error(self, error: str):
        """Handle acquisition error."""
        logging.error(f"Acquisition error: {error}")
        self.status_bar.showMessage(f"Error: {error}", 5000)
        
    @Slot(str)
    def _on_acq_status(self, status: str):
        """Handle acquisition status change."""
        self.status_bar.showMessage(status, 2000)
        
    # === PARAMETER CHANGES ===
    
    @Slot(float)
    def _on_exposure_changed(self, value: float):
        """Handle exposure change."""
        self.camera_config.exposure_ms = value
        if self.camera is not None:
            self.camera.set_exposure(value)
            
    @Slot(float)
    def _on_gain_changed(self, value: float):
        """Handle gain change."""
        self.camera_config.gain = value
        if self.camera is not None:
            self.camera.set_gain(value)
            
    @Slot(int)
    def _on_binning_changed(self, value: int):
        """Handle binning change."""
        self.camera_config.binning = value
        if self.camera is not None:
            self.camera.set_binning(value)
            
    @Slot(object)
    def _on_trigger_changed(self, mode: TriggerMode):
        """Handle trigger mode change."""
        self.camera_config.trigger_mode = mode
        if self.camera is not None:
            self.camera.set_trigger_mode(mode)
            
    @Slot(str, int)
    def _on_camera_selected(self, cam_type: str, index: int):
        """Handle camera selection."""
        type_map = {
            'simulated': CameraType.SIMULATED,
            'hamamatsu': CameraType.HAMAMATSU,
            'thorlabs': CameraType.THORLABS,
            'uc480': CameraType.UC480,
            'andor': CameraType.ANDOR
        }
        ct = type_map.get(cam_type.lower(), CameraType.SIMULATED)
        self._connect_camera(ct, index)
        
    @Slot(float, float)
    def _on_levels_changed(self, lo: float, hi: float):
        """Handle manual level change."""
        self.image_dock.set_levels(lo, hi)
        if self.display_manager is not None:
            self.display_manager.set_manual_levels(lo, hi)
            
    @Slot(str)
    def _on_colormap_changed(self, name: str):
        """Handle colormap change."""
        self.image_dock.set_colormap(name)
        
    @Slot(int, int, int, int)
    def _on_roi_changed(self, x: int, y: int, w: int, h: int):
        """Handle ROI change from image dock (intensity calculation only)."""
        if w == 0 and h == 0:
            # ROI cleared - revert to full-frame mean
            if self.acq_thread is not None:
                self.acq_thread.set_intensity_mode(IntensityMode.MEAN, None)
        else:
            if self.acq_thread is not None:
                self.acq_thread.set_intensity_mode(IntensityMode.ROI_MEAN, (x, y, w, h))
    
    @Slot(int, int, int, int)
    def _on_roi_apply_to_sensor(self, x: int, y: int, w: int, h: int):
        """Apply ROI to camera hardware sensor (stops/restarts acquisition)."""
        if self.camera is None:
            return
            
        was_running = (self.acq_thread is not None and self.acq_thread.isRunning())
        
        # Stop acquisition (required to change sensor ROI)
        if was_running:
            self._stop_acquisition()
        
        # Apply ROI to camera hardware
        try:
            self.camera.set_roi(x, y, w, h)
            actual_roi = self.camera.get_roi()
            self.status_bar.showMessage(
                f"Sensor ROI set: {actual_roi[2]}x{actual_roi[3]} at ({actual_roi[0]},{actual_roi[1]})",
                5000)
            logging.info(f"Sensor ROI applied: {actual_roi}")
            
            # Rebuild acquisition system with new frame shape
            self._setup_acquisition_system()
            
        except Exception as e:
            logging.warning(f"Failed to set sensor ROI: {e}")
            self.status_bar.showMessage(f"ROI error: {e}", 5000)
        
        # Restart if was running
        if was_running:
            self._start_acquisition()
            self.control_dock.set_live_state(True)
            
    @Slot(int)
    def _on_emccd_gain_changed(self, gain: int):
        """Handle EMCCD gain change."""
        if self.camera is not None and hasattr(self.camera, 'set_emccd_gain'):
            self.camera.set_emccd_gain(gain)
            
    @Slot(bool, float)
    def _on_cooling_changed(self, enabled: bool, target_temp: float):
        """Handle cooling change."""
        if self.camera is not None and hasattr(self.camera, 'set_cooling'):
            self.camera.set_cooling(enabled, target_temp)
            
    @Slot(float)
    def _on_display_rate_changed(self, rate: float):
        """Handle display rate change."""
        self.acquisition_config.display_rate_hz = rate
        if self.display_manager is not None:
            self.display_manager.set_display_rate(rate)
    
    @Slot(str, str)
    def _on_record_clicked(self, mode: str, path: str):
        """Handle record start/stop from Sequence tab."""
        if mode == 'none' or not mode:
            # Stop recording
            if self.acq_thread is not None and self.acq_thread.disk_writer is not None:
                try:
                    self.acq_thread.disk_writer.close()
                except Exception:
                    pass
                self.acq_thread.disk_writer = None
            self.status_bar.showMessage("Recording stopped", 3000)
            return
        
        if mode == 'ram':
            # Setup RAM storage buffer
            if self.acq_thread is not None:
                ram_mb = self.control_dock.ram_spin.value()
                # Estimate frames: assume 2048x2048 uint16 = 8 MB/frame
                frame_bytes = 2048 * 2048 * 2  # conservative estimate
                max_frames = int(ram_mb * 1024 * 1024 / frame_bytes)
                self.storage_buffer = FrameRingBuffer(max_frames=max_frames)
                self.acq_thread.set_storage_buffer(self.storage_buffer)
            self.status_bar.showMessage(f"Recording to RAM (max {max_frames} frames)", 3000)
            return
        
        # Disk modes (HDF5 or TIFF)
        if not path:
            self.status_bar.showMessage("No save path specified", 3000)
            return
        
        try:
            if mode == 'hdf5':
                writer = HDF5Writer(filepath=path)
            elif mode == 'tiff':
                from cam_v2_acquisition import TIFFWriter
                writer = TIFFWriter(filepath=path)
            else:
                return
            
            writer.open()
            
            if self.acq_thread is not None:
                self.acq_thread.disk_writer = writer
            
            self.status_bar.showMessage(f"Recording to {path}", 5000)
            logging.info(f"Started recording: mode={mode}, path={path}")
            
        except Exception as e:
            logging.exception(f"Failed to start recording: {e}")
            QMessageBox.critical(self, "Record Error", f"Failed to start recording:\n{e}")
            
    @Slot(bool)
    def _on_fft_toggled(self, enabled: bool):
        """Handle FFT enable/disable toggle."""
        if self.display_manager is not None:
            self.display_manager.set_fft_enabled(enabled)
        # Also set on the standalone processor used in _on_image_ready
        self._fft_enabled = enabled
    
    @Slot(object)
    def _on_filtered_image(self, filtered: np.ndarray):
        """Handle filtered image from FFT dock."""
        self.image_dock.set_image(filtered)
        
    # === SETTINGS PERSISTENCE ===
    
    def _load_settings(self):
        """Load settings from QSettings."""
        settings = QSettings("BPLab", "CameraViewerV2")
        
        # Window geometry
        geometry = settings.value("geometry")
        if geometry:
            self.restoreGeometry(geometry)
            
        # Window state (dock positions)
        state = settings.value("windowState")
        if state:
            self.restoreState(state)
            
        # Last camera type
        last_cam_type = settings.value("lastCameraType", "simulated")
        idx = self.control_dock.camera_type_combo.findText(
            last_cam_type.capitalize())
        if idx >= 0:
            self.control_dock.camera_type_combo.setCurrentIndex(idx)
            
        # Exposure
        exposure = settings.value("exposure", 10.0, type=float)
        self.control_dock.exposure_spin.setValue(exposure)
        
        # Display rate
        disp_rate = settings.value("displayRate", 30.0, type=float)
        self.control_dock.disp_rate_spin.setValue(disp_rate)
        
        # Colormap
        colormap = settings.value("colormap", "gray")
        idx = self.histogram_dock.cmap_combo.findText(colormap)
        if idx >= 0:
            self.histogram_dock.cmap_combo.setCurrentIndex(idx)
            
        logging.info("Settings loaded")
        
    def _save_settings(self):
        """Save settings to QSettings."""
        settings = QSettings("BPLab", "CameraViewerV2")
        
        # Window geometry and state
        settings.setValue("geometry", self.saveGeometry())
        settings.setValue("windowState", self.saveState())
        
        # Camera type
        settings.setValue("lastCameraType", 
                         self.control_dock.camera_type_combo.currentText().lower())
        
        # Parameters
        settings.setValue("exposure", self.control_dock.exposure_spin.value())
        settings.setValue("displayRate", self.control_dock.disp_rate_spin.value())
        settings.setValue("colormap", self.histogram_dock.cmap_combo.currentText())
        
        logging.info("Settings saved")
        
    # === ACTIONS ===
    
    def _on_save_frame(self):
        """Save current frame."""
        frame, _, _, _ = self.frame_buffer.peek() if self.frame_buffer else (None, 0, 0, 0)
        
        if frame is None:
            QMessageBox.warning(self, "No Frame", "No frame available to save")
            return
            
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Frame", "frame.tiff",
            "TIFF Files (*.tiff *.tif);;NumPy Files (*.npy);;All Files (*)"
        )
        
        if not path:
            return
            
        try:
            if path.endswith('.npy'):
                np.save(path, frame)
            else:
                import tifffile
                tifffile.imwrite(path, frame)
            self.status_bar.showMessage(f"Saved to {path}", 3000)
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save:\n{e}")
            
    def _reset_layout(self):
        """Reset dock layout to default."""
        for dock in [self.image_dock, self.histogram_dock, self.control_dock,
                     self.intensity_dock, self.fft_dock]:
            self.removeDockWidget(dock)
            
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.control_dock)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.image_dock)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.histogram_dock)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.intensity_dock)
        self.tabifyDockWidget(self.intensity_dock, self.fft_dock)
        self.intensity_dock.raise_()
        
        for dock in [self.image_dock, self.histogram_dock, self.control_dock,
                     self.intensity_dock, self.fft_dock]:
            dock.show()
            
    def _show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self, "About Camera Viewer v2",
            "<h2>Camera Viewer v2</h2>"
            "<p>Multi-camera streaming viewer with:</p>"
            "<ul>"
            "<li>Multiple camera backend support</li>"
            "<li>Decoupled acquisition and display rates</li>"
            "<li>Live image FFT with filtering</li>"
            "<li>Intensity time trace</li>"
            "<li>Histogram with level controls</li>"
            "</ul>"
            "<p>Author: Claude/BP</p>"
            "<p>Date: January 2025</p>"
        )
        
    def _apply_dark_theme(self):
        """Apply dark theme stylesheet."""
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background-color: #1e1e1e;
                color: #cccccc;
            }
            QDockWidget {
                background-color: #252526;
            }
            QDockWidget::title {
                background-color: #333333;
                padding: 4px;
            }
            QGroupBox {
                border: 1px solid #3c3c3c;
                border-radius: 4px;
                margin-top: 8px;
                padding-top: 8px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 8px;
            }
            QPushButton {
                background-color: #0e639c;
                border: none;
                padding: 6px 12px;
                border-radius: 2px;
            }
            QPushButton:hover {
                background-color: #1177bb;
            }
            QPushButton:pressed {
                background-color: #094771;
            }
            QPushButton:disabled {
                background-color: #3c3c3c;
                color: #666666;
            }
            QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit {
                background-color: #3c3c3c;
                border: 1px solid #555555;
                border-radius: 2px;
                padding: 4px;
            }
            QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover {
                border-color: #007acc;
            }
            QTabWidget::pane {
                border: 1px solid #3c3c3c;
            }
            QTabBar::tab {
                background-color: #2d2d2d;
                padding: 8px 16px;
                border: 1px solid #3c3c3c;
            }
            QTabBar::tab:selected {
                background-color: #1e1e1e;
                border-bottom: 2px solid #007acc;
            }
            QSlider::groove:horizontal {
                height: 4px;
                background: #3c3c3c;
            }
            QSlider::handle:horizontal {
                background: #007acc;
                width: 16px;
                margin: -6px 0;
                border-radius: 8px;
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
            QMenuBar {
                background-color: #333333;
            }
            QMenuBar::item:selected {
                background-color: #094771;
            }
            QMenu {
                background-color: #252526;
                border: 1px solid #3c3c3c;
            }
            QMenu::item:selected {
                background-color: #094771;
            }
        """)
        
    def closeEvent(self, event):
        """Handle window close."""
        # Save settings
        self._save_settings()
        
        # Stop acquisition
        self._stop_acquisition()
        
        # Close camera
        if self.camera is not None:
            self.camera.close()
        event.accept()

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

def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    app = QApplication(sys.argv)
    set_dark_theme(app)
    app.setApplicationName("Camera Viewer v2")
    app.setOrganizationName("BP Lab")
    
    viewer = CameraViewerV2()
    viewer.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
