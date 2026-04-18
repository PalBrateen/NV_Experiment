"""
Camera Viewer v2 - Core Module
==============================

Configuration classes, enums, buffers, and abstract camera interface.

Follows Scope_v2 architecture patterns for consistency.

Author: Claude/BP
Date: January 2025
"""

import numpy as np
import logging
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any
from enum import Enum, auto
from collections import deque

from PySide6.QtCore import QObject, Signal, QMutex, QMutexLocker


# =============================================================================
# ENUMS
# =============================================================================

class CameraType(Enum):
    """Supported camera types"""
    HAMAMATSU = "hamamatsu"
    THORLABS = "thorlabs"        # TLCamera SDK (CS165CU, etc.)
    UC480 = "uc480"              # UC480/IDS cameras (DCC1545M, DCC1645C)
    ANDOR = "andor"
    SIMULATED = "simulated"


class TriggerMode(Enum):
    """Camera trigger modes"""
    INTERNAL = "internal"      # Free-running
    EXTERNAL = "external"      # Hardware trigger
    SOFTWARE = "software"      # Software trigger


class TriggerPolarity(Enum):
    """Trigger edge polarity"""
    POSITIVE = "positive"
    NEGATIVE = "negative"


class StorageMode(Enum):
    """Frame storage modes"""
    NONE = "none"              # Display only
    RAM = "ram"                # Store in RAM buffer
    DISK_HDF5 = "hdf5"         # Stream to HDF5
    DISK_TIFF = "tiff"         # Stream to TIFF stack
    DISK_NPY = "npy"           # Stream to NumPy files


class IntensityMode(Enum):
    """Intensity calculation modes"""
    SUM = "sum"
    MEAN = "mean"
    MAX = "max"
    ROI_SUM = "roi_sum"
    ROI_MEAN = "roi_mean"


class FFTWindow(Enum):
    """2D FFT window functions"""
    NONE = "none"
    HANN = "hann"
    HAMMING = "hamming"
    BLACKMAN = "blackman"


class FFTDisplayMode(Enum):
    """FFT display modes"""
    MAGNITUDE = "magnitude"
    LOG_MAGNITUDE = "log_magnitude"
    PHASE = "phase"
    POWER = "power"


# =============================================================================
# CONFIGURATION CLASSES
# =============================================================================

@dataclass
class CameraConfig:
    """Camera hardware configuration"""
    camera_type: CameraType = CameraType.SIMULATED
    device_index: int = 0
    
    # Basic settings
    exposure_ms: float = 10.0
    gain: float = 1.0
    
    # ROI (None = full frame)
    roi: Optional[Tuple[int, int, int, int]] = None  # (x, y, w, h)
    binning: int = 1
    
    # Trigger
    trigger_mode: TriggerMode = TriggerMode.INTERNAL
    trigger_polarity: TriggerPolarity = TriggerPolarity.POSITIVE
    trigger_delay_us: float = 0.0
    
    # Cooling (for cooled cameras like Hamamatsu, Andor)
    cooling_enabled: bool = True
    target_temperature_c: float = -10.0
    
    # Andor-specific
    emccd_gain: int = 1
    fan_mode: str = "on"  # "on", "low", "off"


@dataclass
class AcquisitionConfig:
    """Acquisition and display configuration"""
    # Display rate (decoupled from acquisition)
    display_rate_hz: float = 30.0
    
    # Intensity calculation
    intensity_mode: IntensityMode = IntensityMode.MEAN
    intensity_roi: Optional[Tuple[int, int, int, int]] = None
    
    # Buffer sizes
    intensity_buffer_size: int = 1000  # Rolling intensity trace points
    
    # RAM ceiling (for RAM storage mode)
    ram_ceiling_mb: float = 2000.0
    max_ram_frames: int = 1000


@dataclass
class StorageConfig:
    """Frame storage configuration"""
    mode: StorageMode = StorageMode.NONE
    save_path: str = ""
    
    # File naming
    prefix: str = "frame"
    start_number: int = 0
    use_leading_zeros: bool = True
    
    # HDF5 options
    compression: bool = True
    chunk_size: int = 100
    
    # What to save
    save_timestamps: bool = True
    save_intensities: bool = True


@dataclass
class DisplayConfig:
    """Image display configuration"""
    # Colormap
    colormap: str = "gray"  # gray, viridis, inferno, etc.
    
    # Levels (None = auto)
    level_min: Optional[float] = None
    level_max: Optional[float] = None
    auto_levels: bool = True
    auto_levels_percentile: Tuple[float, float] = (1.0, 99.0)
    
    # Transform
    flip_horizontal: bool = False
    flip_vertical: bool = False
    rotation_degrees: int = 0  # 0, 90, 180, 270
    
    # Scaling
    gamma: float = 1.0
    log_scale: bool = False


@dataclass
class FFTConfig:
    """Image FFT configuration"""
    enabled: bool = False
    live_update: bool = True  # Update FFT on live stream
    
    # Display
    window: FFTWindow = FFTWindow.HANN
    display_mode: FFTDisplayMode = FFTDisplayMode.LOG_MAGNITUDE
    
    # Filtering
    filter_enabled: bool = False
    filter_type: str = "none"  # lowpass, highpass, bandpass, notch
    filter_params: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# CAMERA INTERFACE (ABSTRACT BASE CLASS)
# =============================================================================

class CameraInterface(ABC):
    """
    Abstract base class for camera backends.
    
    Implementations:
    - HamamatsuBackend (using dcam.py)
    - ThorlabsBackend (using pylablib)
    - AndorBackend (using pylablib SDK2)
    - SimulatedBackend (for testing)
    """
    
    def __init__(self):
        self._is_open = False
        self._is_acquiring = False
        self._frame_count = 0
        
    # === Connection ===
    
    @abstractmethod
    def open(self) -> bool:
        """
        Open connection to camera.
        
        Returns:
            True if successful, False otherwise
        """
        pass
    
    @abstractmethod
    def close(self):
        """Close connection to camera."""
        pass
    
    @property
    def is_open(self) -> bool:
        """Check if camera is connected."""
        return self._is_open
    
    # === Acquisition Control ===
    
    @abstractmethod
    def start_acquisition(self):
        """Start continuous acquisition."""
        pass
    
    @abstractmethod
    def stop_acquisition(self):
        """Stop acquisition."""
        pass
    
    @property
    def is_acquiring(self) -> bool:
        """Check if currently acquiring."""
        return self._is_acquiring
    
    @abstractmethod
    def get_frame(self, timeout_ms: int = 1000) -> Optional[np.ndarray]:
        """
        Get next frame from camera.
        
        Args:
            timeout_ms: Timeout in milliseconds
            
        Returns:
            Frame as numpy array, or None if timeout/error
        """
        pass
    
    def snap(self) -> Optional[np.ndarray]:
        """Acquire single frame (convenience method)."""
        was_acquiring = self._is_acquiring
        if not was_acquiring:
            self.start_acquisition()
        
        frame = self.get_frame()
        
        if not was_acquiring:
            self.stop_acquisition()
        
        return frame
    
    # === Settings ===
    
    @abstractmethod
    def set_exposure(self, exposure_ms: float):
        """Set exposure time in milliseconds."""
        pass
    
    @abstractmethod
    def get_exposure(self) -> float:
        """Get current exposure time in milliseconds."""
        pass
    
    @abstractmethod
    def set_roi(self, x: int, y: int, width: int, height: int):
        """Set region of interest."""
        pass
    
    @abstractmethod
    def get_roi(self) -> Tuple[int, int, int, int]:
        """Get current ROI (x, y, width, height)."""
        pass
    
    @abstractmethod
    def clear_roi(self):
        """Reset ROI to full frame."""
        pass
    
    @abstractmethod
    def set_binning(self, binning: int):
        """Set pixel binning (1, 2, 4, etc.)."""
        pass
    
    @abstractmethod
    def get_binning(self) -> int:
        """Get current binning."""
        pass
    
    @abstractmethod
    def set_trigger_mode(self, mode: TriggerMode):
        """Set trigger mode."""
        pass
    
    @abstractmethod
    def get_trigger_mode(self) -> TriggerMode:
        """Get current trigger mode."""
        pass
    
    # === Properties ===
    
    @property
    @abstractmethod
    def frame_shape(self) -> Tuple[int, int]:
        """Get frame shape (height, width)."""
        pass
    
    @property
    @abstractmethod
    def pixel_dtype(self) -> np.dtype:
        """Get pixel data type (e.g., uint16)."""
        pass
    
    @property
    @abstractmethod
    def is_color(self) -> bool:
        """Check if camera is color (RGB)."""
        pass
    
    @property
    @abstractmethod
    def camera_info(self) -> Dict[str, Any]:
        """
        Get camera information.
        
        Returns dict with keys like:
        - model, serial, firmware, driver_version
        - sensor_size, pixel_size, bit_depth
        - min_exposure, max_exposure
        - temperature (if cooled)
        """
        pass
    
    @property
    def frame_count(self) -> int:
        """Get total frames acquired since start."""
        return self._frame_count
    
    # === Optional: Cooling (for cooled cameras) ===
    
    def set_cooling(self, enabled: bool, target_temp_c: float = -10.0):
        """Set cooling (override in cooled camera backends)."""
        pass
    
    def get_temperature(self) -> Optional[float]:
        """Get current sensor temperature (override in cooled cameras)."""
        return None
    
    # === Optional: Gain (for cameras with gain control) ===
    
    def set_gain(self, gain: float):
        """Set gain (override in cameras with gain control)."""
        pass
    
    def get_gain(self) -> float:
        """Get current gain."""
        return 1.0
    
    # === Optional: EMCCD Gain (for EMCCDs) ===
    
    def set_emccd_gain(self, gain: int):
        """Set EMCCD gain (override in EMCCD backends)."""
        pass
    
    def get_emccd_gain(self) -> int:
        """Get current EMCCD gain."""
        return 1


# =============================================================================
# BUFFERS
# =============================================================================

class LatestFrameBuffer:
    """
    Thread-safe buffer for display - keeps only the latest frame.
    
    For decoupling acquisition (fast) from display (slower).
    Acquisition writes at camera rate, display reads at display rate.
    """
    
    def __init__(self):
        self._frame: Optional[np.ndarray] = None
        self._timestamp: float = 0.0
        self._intensity: float = 0.0
        self._frame_number: int = 0
        self._mutex = QMutex()
        self._new_frame_available = False
        
    def write(self, frame: np.ndarray, timestamp: float, 
              intensity: float, frame_number: int):
        """Write new frame (overwrites previous)."""
        with QMutexLocker(self._mutex):
            self._frame = frame.copy()
            self._timestamp = timestamp
            self._intensity = intensity
            self._frame_number = frame_number
            self._new_frame_available = True
            
    def get_latest(self) -> Tuple[Optional[np.ndarray], float, float, int]:
        """
        Get latest frame if available.
        
        Returns:
            (frame, timestamp, intensity, frame_number) or (None, 0, 0, 0)
        """
        with QMutexLocker(self._mutex):
            if self._new_frame_available and self._frame is not None:
                self._new_frame_available = False
                return (self._frame.copy(), self._timestamp, 
                        self._intensity, self._frame_number)
            return None, 0.0, 0.0, 0
    
    def peek(self) -> Tuple[Optional[np.ndarray], float, float, int]:
        """Get latest frame without clearing flag."""
        with QMutexLocker(self._mutex):
            if self._frame is not None:
                return (self._frame.copy(), self._timestamp,
                        self._intensity, self._frame_number)
            return None, 0.0, 0.0, 0
    
    def clear(self):
        """Clear buffer."""
        with QMutexLocker(self._mutex):
            self._frame = None
            self._new_frame_available = False
            
    @property
    def has_new_frame(self) -> bool:
        with QMutexLocker(self._mutex):
            return self._new_frame_available


class IntensityRingBuffer:
    """
    Thread-safe ring buffer for intensity time trace.
    
    Stores (timestamp, intensity) pairs for rolling plot.
    """
    
    def __init__(self, capacity: int = 1000):
        self.capacity = capacity
        self._timestamps = np.zeros(capacity, dtype=np.float64)
        self._intensities = np.zeros(capacity, dtype=np.float64)
        self._write_idx = 0
        self._count = 0
        self._mutex = QMutex()
        
    def write(self, timestamp: float, intensity: float):
        """Add new data point."""
        with QMutexLocker(self._mutex):
            self._timestamps[self._write_idx] = timestamp
            self._intensities[self._write_idx] = intensity
            self._write_idx = (self._write_idx + 1) % self.capacity
            self._count = min(self._count + 1, self.capacity)
            
    def write_batch(self, timestamps: np.ndarray, intensities: np.ndarray):
        """Add multiple data points."""
        with QMutexLocker(self._mutex):
            n = len(timestamps)
            if n >= self.capacity:
                # Just keep the last capacity points
                self._timestamps[:] = timestamps[-self.capacity:]
                self._intensities[:] = intensities[-self.capacity:]
                self._write_idx = 0
                self._count = self.capacity
                return
            
            end_idx = self._write_idx + n
            if end_idx <= self.capacity:
                self._timestamps[self._write_idx:end_idx] = timestamps
                self._intensities[self._write_idx:end_idx] = intensities
            else:
                first = self.capacity - self._write_idx
                self._timestamps[self._write_idx:] = timestamps[:first]
                self._intensities[self._write_idx:] = intensities[:first]
                self._timestamps[:n-first] = timestamps[first:]
                self._intensities[:n-first] = intensities[first:]
            
            self._write_idx = end_idx % self.capacity
            self._count = min(self._count + n, self.capacity)
    
    def get_latest(self, n: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get latest n points (or all available).
        
        Returns:
            (timestamps, intensities) arrays in chronological order
        """
        with QMutexLocker(self._mutex):
            if n is None:
                n = self._count
            n = min(n, self._count)
            
            if n == 0:
                return np.array([]), np.array([])
            
            # Read in chronological order
            start = (self._write_idx - n) % self.capacity
            
            if start + n <= self.capacity:
                times = self._timestamps[start:start+n].copy()
                values = self._intensities[start:start+n].copy()
            else:
                first = self.capacity - start
                times = np.concatenate([
                    self._timestamps[start:],
                    self._timestamps[:n-first]
                ])
                values = np.concatenate([
                    self._intensities[start:],
                    self._intensities[:n-first]
                ])
            
            return times, values
    
    def clear(self):
        """Clear buffer."""
        with QMutexLocker(self._mutex):
            self._write_idx = 0
            self._count = 0
            
    @property
    def count(self) -> int:
        with QMutexLocker(self._mutex):
            return self._count


class FrameRingBuffer:
    """
    Thread-safe ring buffer for frame history (RAM storage mode).
    
    Stores frames with timestamps and intensities.
    """
    
    def __init__(self, max_frames: int = 100, 
                 frame_shape: Tuple[int, int] = (512, 512),
                 dtype: np.dtype = np.uint16):
        self.max_frames = max_frames
        self._frames = np.zeros((max_frames,) + frame_shape, dtype=dtype)
        self._timestamps = np.zeros(max_frames, dtype=np.float64)
        self._intensities = np.zeros(max_frames, dtype=np.float64)
        self._write_idx = 0
        self._count = 0
        self._mutex = QMutex()
        self._initialized = False
        
    def initialize(self, frame_shape: Tuple[int, int], dtype: np.dtype):
        """Initialize buffer with correct frame shape."""
        with QMutexLocker(self._mutex):
            self._frames = np.zeros((self.max_frames,) + frame_shape, dtype=dtype)
            self._initialized = True
            
    def write(self, frame: np.ndarray, timestamp: float, intensity: float) -> bool:
        """
        Add frame to buffer.
        
        Returns:
            True if stored, False if buffer full (non-circular mode)
        """
        with QMutexLocker(self._mutex):
            if not self._initialized:
                self._frames = np.zeros((self.max_frames,) + frame.shape, 
                                        dtype=frame.dtype)
                self._initialized = True
            
            self._frames[self._write_idx] = frame
            self._timestamps[self._write_idx] = timestamp
            self._intensities[self._write_idx] = intensity
            
            self._write_idx = (self._write_idx + 1) % self.max_frames
            self._count = min(self._count + 1, self.max_frames)
            return True
    
    def get_frame(self, index: int) -> Optional[Tuple[np.ndarray, float, float]]:
        """Get frame by index (0 = oldest available)."""
        with QMutexLocker(self._mutex):
            if index < 0 or index >= self._count:
                return None
            
            if self._count < self.max_frames:
                actual_idx = index
            else:
                actual_idx = (self._write_idx + index) % self.max_frames
            
            return (self._frames[actual_idx].copy(),
                    self._timestamps[actual_idx],
                    self._intensities[actual_idx])
    
    def get_all(self) -> Dict[str, np.ndarray]:
        """Get all stored frames in chronological order."""
        with QMutexLocker(self._mutex):
            if self._count == 0:
                return {'frames': np.array([]), 
                        'timestamps': np.array([]),
                        'intensities': np.array([])}
            
            if self._count < self.max_frames:
                return {
                    'frames': self._frames[:self._count].copy(),
                    'timestamps': self._timestamps[:self._count].copy(),
                    'intensities': self._intensities[:self._count].copy()
                }
            else:
                # Reorder to chronological
                start = self._write_idx
                return {
                    'frames': np.concatenate([
                        self._frames[start:],
                        self._frames[:start]
                    ]),
                    'timestamps': np.concatenate([
                        self._timestamps[start:],
                        self._timestamps[:start]
                    ]),
                    'intensities': np.concatenate([
                        self._intensities[start:],
                        self._intensities[:start]
                    ])
                }
    
    def clear(self):
        """Clear buffer."""
        with QMutexLocker(self._mutex):
            self._write_idx = 0
            self._count = 0
            
    @property
    def count(self) -> int:
        with QMutexLocker(self._mutex):
            return self._count
    
    @property
    def memory_usage_mb(self) -> float:
        """Estimate memory usage in MB."""
        with QMutexLocker(self._mutex):
            if not self._initialized:
                return 0.0
            return self._frames.nbytes / (1024 * 1024)


# =============================================================================
# OBSERVABLE CONFIG (for IPython/GUI sync)
# =============================================================================

class ObservableCameraConfig(QObject):
    """
    Camera configuration with signals for property changes.
    
    Allows IPython console to modify camera settings and
    have the GUI automatically update.
    """
    
    exposureChanged = Signal(float)
    gainChanged = Signal(float)
    roiChanged = Signal(object)  # Tuple or None
    triggerModeChanged = Signal(object)
    displayRateChanged = Signal(float)
    acquisitionStateChanged = Signal(bool)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._camera_config = CameraConfig()
        self._acquisition_config = AcquisitionConfig()
        self._display_config = DisplayConfig()
        self._is_running = False
        
    # === Camera Config Properties ===
    
    @property
    def exposure_ms(self) -> float:
        return self._camera_config.exposure_ms
    
    @exposure_ms.setter
    def exposure_ms(self, value: float):
        if value != self._camera_config.exposure_ms:
            self._camera_config.exposure_ms = value
            self.exposureChanged.emit(value)
            
    @property
    def gain(self) -> float:
        return self._camera_config.gain
    
    @gain.setter
    def gain(self, value: float):
        if value != self._camera_config.gain:
            self._camera_config.gain = value
            self.gainChanged.emit(value)
            
    @property
    def roi(self) -> Optional[Tuple[int, int, int, int]]:
        return self._camera_config.roi
    
    @roi.setter
    def roi(self, value: Optional[Tuple[int, int, int, int]]):
        if value != self._camera_config.roi:
            self._camera_config.roi = value
            self.roiChanged.emit(value)
            
    @property
    def trigger_mode(self) -> TriggerMode:
        return self._camera_config.trigger_mode
    
    @trigger_mode.setter
    def trigger_mode(self, value: TriggerMode):
        if value != self._camera_config.trigger_mode:
            self._camera_config.trigger_mode = value
            self.triggerModeChanged.emit(value)
            
    # === Acquisition Config Properties ===
    
    @property
    def display_rate_hz(self) -> float:
        return self._acquisition_config.display_rate_hz
    
    @display_rate_hz.setter
    def display_rate_hz(self, value: float):
        if value != self._acquisition_config.display_rate_hz:
            self._acquisition_config.display_rate_hz = value
            self.displayRateChanged.emit(value)
            
    # === State ===
    
    @property
    def is_running(self) -> bool:
        return self._is_running
    
    @is_running.setter
    def is_running(self, value: bool):
        if value != self._is_running:
            self._is_running = value
            self.acquisitionStateChanged.emit(value)
            
    # === Config Access ===
    
    def get_camera_config(self) -> CameraConfig:
        return self._camera_config
    
    def get_acquisition_config(self) -> AcquisitionConfig:
        return self._acquisition_config
    
    def get_display_config(self) -> DisplayConfig:
        return self._display_config


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def compute_intensity(frame: np.ndarray, 
                      mode: IntensityMode = IntensityMode.SUM,
                      roi: Optional[Tuple[int, int, int, int]] = None) -> float:
    """
    Compute intensity value from frame.
    
    Args:
        frame: 2D image array
        mode: Intensity calculation mode
        roi: Optional ROI (x, y, w, h) for ROI modes
        
    Returns:
        Intensity value
    """
    if roi is not None and mode in [IntensityMode.ROI_SUM, IntensityMode.ROI_MEAN]:
        x, y, w, h = roi
        frame = frame[y:y+h, x:x+w]
        mode = IntensityMode.SUM if mode == IntensityMode.ROI_SUM else IntensityMode.MEAN
    
    if mode == IntensityMode.SUM:
        return float(np.sum(frame))
    elif mode == IntensityMode.MEAN:
        return float(np.mean(frame))
    elif mode == IntensityMode.MAX:
        return float(np.max(frame))
    else:
        return float(np.sum(frame))


def compute_histogram(frame: np.ndarray, 
                      bins: int = 256,
                      range_: Optional[Tuple[float, float]] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute histogram of frame.
    
    Args:
        frame: 2D image array
        bins: Number of histogram bins
        range_: (min, max) range, or None for auto
        
    Returns:
        (counts, bin_edges) arrays
    """
    if range_ is None:
        if frame.dtype == np.uint8:
            range_ = (0, 255)
        elif frame.dtype == np.uint16:
            range_ = (0, 65535)
        else:
            range_ = (float(frame.min()), float(frame.max()))
    
    counts, edges = np.histogram(frame.ravel(), bins=bins, range=range_)
    return counts, edges


def auto_levels(frame: np.ndarray, 
                percentile: Tuple[float, float] = (1.0, 99.0)) -> Tuple[float, float]:
    """
    Compute auto-levels based on percentiles.
    
    Args:
        frame: 2D image array
        percentile: (low, high) percentiles
        
    Returns:
        (level_min, level_max)
    """
    lo = np.percentile(frame, percentile[0])
    hi = np.percentile(frame, percentile[1])
    return float(lo), float(hi)


def apply_transform(frame: np.ndarray, 
                    config: DisplayConfig) -> np.ndarray:
    """
    Apply display transforms to frame.
    
    Args:
        frame: Input frame
        config: Display configuration
        
    Returns:
        Transformed frame
    """
    result = frame
    
    # Flip
    if config.flip_horizontal:
        result = np.fliplr(result)
    if config.flip_vertical:
        result = np.flipud(result)
    
    # Rotate
    if config.rotation_degrees != 0:
        k = config.rotation_degrees // 90
        result = np.rot90(result, k)
    
    return result
