"""
Camera Viewer v2 - Multi-camera streaming viewer
==================================================

A PySide6 + pyqtgraph camera viewer supporting multiple backends
with decoupled acquisition and display rates.

Usage:
    from camera_viewer_v2 import CameraViewerV2
    from PySide6.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)
    viewer = CameraViewerV2()
    viewer.show()
    app.exec()

Backends:
    - SimulatedBackend: Test patterns (no hardware needed)
    - HamamatsuBackend: ORCA series via dcam.py
    - ThorlabsBackend: CS165CU etc. via pylablib TLCamera SDK
    - UC480Backend: DCC1545M, DCC1645C via pylablib uc480
    - AndorBackend: iXon Ultra 888 via pylablib SDK2

Author: Claude/BP
Date: January 2025
"""

# Core types and configurations
from cam_v2_core import (
    CameraType, TriggerMode, TriggerPolarity,
    StorageMode, IntensityMode, FFTWindow, FFTDisplayMode,
    CameraConfig, AcquisitionConfig, StorageConfig,
    DisplayConfig, FFTConfig,
    CameraInterface,
    LatestFrameBuffer, IntensityRingBuffer, FrameRingBuffer,
    ObservableCameraConfig,
    compute_intensity, compute_histogram, auto_levels, apply_transform
)

# Camera backends
from cam_v2_backends import (
    SimulatedBackend, HamamatsuBackend, ThorlabsBackend,
    UC480Backend, AndorBackend,
    create_camera, list_available_cameras
)

# Acquisition and display
from cam_v2_acquisition import (
    AcquisitionThread, DisplayManager, DisplayManagerConfig,
    ImageFFTProcessor,
    DiskWriter, HDF5Writer, TIFFWriter,
    create_acquisition_system, create_disk_writer
)

# Widgets
from cam_v2_widgets import (
    ImageDisplayDock, HistogramDock, IntensityTraceDock,
    ImageFFTDock, ControlDock
)

# Main application
from cam_v2_main import CameraViewerV2

__all__ = [
    # Core
    'CameraType', 'TriggerMode', 'TriggerPolarity',
    'StorageMode', 'IntensityMode', 'FFTWindow', 'FFTDisplayMode',
    'CameraConfig', 'AcquisitionConfig', 'StorageConfig',
    'DisplayConfig', 'FFTConfig',
    'CameraInterface',
    'LatestFrameBuffer', 'IntensityRingBuffer', 'FrameRingBuffer',
    'ObservableCameraConfig',
    'compute_intensity', 'compute_histogram', 'auto_levels', 'apply_transform',
    # Backends
    'SimulatedBackend', 'HamamatsuBackend', 'ThorlabsBackend',
    'UC480Backend', 'AndorBackend',
    'create_camera', 'list_available_cameras',
    # Acquisition
    'AcquisitionThread', 'DisplayManager', 'DisplayManagerConfig',
    'ImageFFTProcessor',
    'DiskWriter', 'HDF5Writer', 'TIFFWriter',
    'create_acquisition_system', 'create_disk_writer',
    # Widgets
    'ImageDisplayDock', 'HistogramDock', 'IntensityTraceDock',
    'ImageFFTDock', 'ControlDock',
    # Main
    'CameraViewerV2',
]

__version__ = '2.0.0'
