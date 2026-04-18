"""
Camera Viewer v2 - Acquisition Module
======================================

Acquisition thread and display manager with decoupled rates.

Architecture:
- AcquisitionThread: Runs at camera rate, writes all frames to disk/RAM
- DisplayManager: QTimer at user-defined rate (e.g., 30 Hz), reads latest frame

Author: Claude/BP
Date: January 2025
"""

import numpy as np
import time
import logging
from typing import Optional, Dict, Tuple, Callable
from dataclasses import dataclass, field

from PySide6.QtCore import QThread, Signal, QObject, QTimer, QMutex, QMutexLocker

from cam_v2_core import (
    CameraInterface, CameraConfig, AcquisitionConfig, StorageConfig, DisplayConfig,
    LatestFrameBuffer, IntensityRingBuffer, FrameRingBuffer,
    IntensityMode, StorageMode,
    compute_intensity, compute_histogram, auto_levels
)


# =============================================================================
# ACQUISITION THREAD
# =============================================================================

class AcquisitionThread(QThread):
    """
    Background acquisition thread.
    
    Runs at camera frame rate (as fast as camera delivers frames).
    - Writes ALL frames to disk if streaming enabled
    - Stores ALL intensities in ring buffer for complete trace
    - Updates latest frame buffer for display (display reads at its own rate)
    
    Signals:
        frameAcquired: Emitted for each frame (frame_number, timestamp, intensity)
        errorOccurred: Emitted on error (error_message)
        statusChanged: Emitted on status change (status_message)
        acquisitionStats: Emitted periodically with stats (dict)
    """
    
    frameAcquired = Signal(int, float, float)  # frame_number, timestamp, intensity
    errorOccurred = Signal(str)
    statusChanged = Signal(str)
    acquisitionStats = Signal(dict)  # fps, frame_count, dropped, etc.
    
    def __init__(self, 
                 camera: CameraInterface,
                 frame_buffer: LatestFrameBuffer,
                 intensity_buffer: IntensityRingBuffer,
                 config: AcquisitionConfig,
                 parent=None):
        super().__init__(parent)
        
        self.camera: CameraInterface = camera
        self.frame_buffer = frame_buffer
        self.intensity_buffer = intensity_buffer
        self.config = config
        
        # Storage (optional)
        self.storage_buffer: Optional[FrameRingBuffer] = None
        self.disk_writer = None  # Set externally if disk streaming
        
        # Control
        self._stop_requested = False
        self._pause_requested = False
        self._mutex = QMutex()
        
        # Stats
        self._frame_count = 0
        self._dropped_count = 0
        self._start_time = 0.0
        self._last_stats_time = 0.0
        self._stats_interval = 1.0  # Report stats every 1 second
        
        # Intensity calculation
        self._intensity_mode = config.intensity_mode
        self._intensity_roi = config.intensity_roi
        
    def request_stop(self):
        """Request thread to stop."""
        with QMutexLocker(self._mutex):
            self._stop_requested = True
            
    def request_pause(self, pause: bool):
        """Request pause/resume."""
        with QMutexLocker(self._mutex):
            self._pause_requested = pause
            
    def set_intensity_mode(self, mode: IntensityMode, roi: Optional[Tuple[int, int, int, int]] = None):
        """Update intensity calculation mode."""
        with QMutexLocker(self._mutex):
            self._intensity_mode = mode
            self._intensity_roi = roi
            
    def set_storage_buffer(self, buffer: Optional[FrameRingBuffer]):
        """Set RAM storage buffer."""
        self.storage_buffer = buffer
        
    def run(self):
        """Main acquisition loop."""
        self._stop_requested = False
        self._pause_requested = False
        self._frame_count = 0
        self._dropped_count = 0
        self._start_time = time.time()
        self._last_stats_time = self._start_time
        
        self.statusChanged.emit("Starting acquisition...")
        # print("ACQ THREAD: Starting run loop")

        try:
            # Start camera acquisition if not already running
            if not self.camera.is_acquiring:
                self.camera.start_acquisition()
                # print("start??")
            
            # print("ACQ THREAD: Camera acquisition started")
            # print(f"camera = {self.camera}")

            self.statusChanged.emit("Acquiring")
            # i=0
            while not self._stop_requested:
                # Check for pause
                if self._pause_requested:
                    self.statusChanged.emit("Paused")
                    time.sleep(0.1)
                    continue
                # i+=1
                # if i==5: break
                # Get frame from camera
                frame = self.camera.get_frame(timeout_ms=1000)
                # print(f"ACQ THREAD: Got frame: {type(frame)}, shape={frame.shape if frame is not None else None}")

                if frame is None:
                    # Timeout or no frame - could be waiting for trigger
                    self._dropped_count += 1
                    continue
                
                timestamp = time.time()
                self._frame_count += 1
                
                # Compute intensity
                with QMutexLocker(self._mutex):
                    mode = self._intensity_mode
                    roi = self._intensity_roi
                intensity = compute_intensity(frame, mode, roi)
                
                # === ALWAYS: Update intensity buffer (for complete trace) ===
                self.intensity_buffer.write(timestamp, intensity)
                
                # === ALWAYS: Update latest frame buffer (for display) ===
                # Display thread reads this at its own rate
                self.frame_buffer.write(frame, timestamp, intensity, self._frame_count)
                
                # === OPTIONAL: Store in RAM buffer ===
                if self.storage_buffer is not None:
                    self.storage_buffer.write(frame, timestamp, intensity)
                
                # === OPTIONAL: Write to disk ===
                if self.disk_writer is not None:
                    try:
                        self.disk_writer.write_frame(frame, timestamp, intensity, 
                                                     self._frame_count)
                    except Exception as e:
                        logging.warning(f"Disk write error: {e}")
                
                # Emit frame acquired signal (lightweight, just numbers)
                self.frameAcquired.emit(self._frame_count, timestamp, intensity)
                
                # Periodic stats reporting
                now = time.time()
                if now - self._last_stats_time >= self._stats_interval:
                    self._report_stats(now)
                    self._last_stats_time = now
                
                    
        except Exception as e:
            self.errorOccurred.emit(f"Acquisition error: {e}")
            logging.exception(f"Acquisition thread error: {e}")
            print(f"ACQ THREAD: Exception: {e}")
            import traceback
            traceback.print_exc()
            
        finally:
            self.statusChanged.emit("Stopped")
            self._report_stats(time.time())  # Final stats
    
    def _report_stats(self, now: float):
        """Emit acquisition statistics."""
        elapsed = now - self._start_time
        fps = self._frame_count / elapsed if elapsed > 0 else 0.0
        
        stats = {
            'frame_count': self._frame_count,
            'dropped_count': self._dropped_count,
            'elapsed_time': elapsed,
            'fps': fps,
            'intensity_buffer_count': self.intensity_buffer.count,
        }
        
        if self.storage_buffer is not None:
            stats['storage_count'] = self.storage_buffer.count
            stats['storage_memory_mb'] = self.storage_buffer.memory_usage_mb
            
        self.acquisitionStats.emit(stats)


# =============================================================================
# DISPLAY MANAGER
# =============================================================================

@dataclass
class DisplayManagerConfig:
    """Display manager configuration."""
    display_rate_hz: float = 30.0
    intensity_trace_points: int = 500  # Points to show in trace
    auto_levels: bool = True
    auto_levels_percentile: Tuple[float, float] = (1.0, 99.0)
    fft_enabled: bool = False
    histogram_bins: int = 256
    
    # --- Performance optimizations ---
    # Skip histogram: compute every Nth displayed frame (1 = every frame)
    histogram_skip: int = 3
    # Downsample frames larger than this for display (pixels per side)
    downsample_threshold: int = 1024
    # Convert 16-bit to 8-bit before display (reduces setImage cost)
    fast_8bit_display: bool = True


class DisplayManager(QObject):
    """
    Display manager - decoupled from acquisition rate.
    
    Uses QTimer to update displays at a fixed rate (e.g., 30 Hz).
    Reads latest frame from buffer (skips intermediate frames).
    
    Signals:
        imageReady: New image for display (frame, timestamp)
        histogramReady: Histogram data (counts, edges, levels)
        intensityTraceReady: Intensity time series (times, values)
        fftReady: FFT data (magnitude_2d, phase_2d) - if enabled
        statsReady: Display statistics (dict)
        levelsComputed: Auto-computed levels (lo, hi)
    """
    
    imageReady = Signal(object, float)  # frame, timestamp
    histogramReady = Signal(object, object, float, float)  # counts, edges, lo, hi
    intensityTraceReady = Signal(object, object)  # times, values
    fftReady = Signal(object)  # magnitude_2d (complex fft stored internally)
    statsReady = Signal(dict)
    levelsComputed = Signal(float, float)  # lo, hi
    
    def __init__(self,
                 frame_buffer: LatestFrameBuffer,
                 intensity_buffer: IntensityRingBuffer,
                 config: DisplayManagerConfig = None,
                 parent=None):
        super().__init__(parent)
        
        self.frame_buffer = frame_buffer
        self.intensity_buffer = intensity_buffer
        self.config = config or DisplayManagerConfig()
        
        # Display timer
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._update_display)
        
        # Stats
        self._running = False
        self._display_count = 0
        self._last_stats_time = 0.0
        self._actual_display_rate = 0.0
        
        # Levels (manual or auto)
        self._level_lo = 0.0
        self._level_hi = 65535.0
        self._manual_levels = False
        
        # FFT processor (lazy init)
        self._fft_processor = None
        self._last_fft_frame = None
        
        # Rate limiter for stats
        self._stats_timer = QTimer(self)
        self._stats_timer.timeout.connect(self._report_stats)
        
    def start(self):
        """Start display updates."""
        if self._running:
            return
            
        self._running = True
        self._display_count = 0
        self._last_stats_time = time.time()
        
        # Start display timer
        interval_ms = int(1000.0 / self.config.display_rate_hz)
        self._timer.start(interval_ms)
        
        # Start stats timer (every 1 second)
        self._stats_timer.start(1000)
        
    def stop(self):
        """Stop display updates."""
        self._running = False
        self._timer.stop()
        self._stats_timer.stop()
        
    def set_display_rate(self, rate_hz: float):
        """Update display rate."""
        self.config.display_rate_hz = max(1.0, min(120.0, rate_hz))
        if self._running:
            interval_ms = int(1000.0 / self.config.display_rate_hz)
            self._timer.setInterval(interval_ms)
            
    def set_manual_levels(self, lo: float, hi: float):
        """Set manual display levels."""
        self._level_lo = lo
        self._level_hi = hi
        self._manual_levels = True
        
    def set_auto_levels(self, enabled: bool, percentile: Tuple[float, float] = None):
        """Enable/disable auto levels."""
        self.config.auto_levels = enabled
        self._manual_levels = not enabled
        if percentile:
            self.config.auto_levels_percentile = percentile
            
    def set_fft_enabled(self, enabled: bool):
        """Enable/disable FFT computation."""
        self.config.fft_enabled = enabled
        if enabled and self._fft_processor is None:
            self._fft_processor = ImageFFTProcessor()
            
    def _update_display(self):
        """Called by timer - update all displays with performance optimizations."""
        if not self._running:
            return
        
        # Get latest frame (returns None if no new frame)
        frame, timestamp, intensity, frame_num = self.frame_buffer.get_latest()
        
        if frame is None:
            return  # No new frame available
            
        self._display_count += 1
        
        # --- Optimization: downsample large frames for display ---
        display_frame = frame
        downsampled_frame = frame  # Keep pre-8bit version for FFT
        threshold = self.config.downsample_threshold
        h, w = frame.shape[:2]
        if h > threshold or w > threshold:
            # Pick stride that brings both dims below threshold
            stride = max(h // threshold, w // threshold, 1) + 1
            if frame.ndim == 2:
                display_frame = frame[::stride, ::stride]
            else:
                display_frame = frame[::stride, ::stride, :]
            downsampled_frame = display_frame
        
        # --- Optimization: 8-bit display for 16-bit cameras ---
        if self.config.fast_8bit_display and display_frame.dtype == np.uint16:
            # Scale to 8-bit using current levels for best contrast
            lo, hi = self._level_lo, self._level_hi
            if hi > lo:
                scaled = np.clip(display_frame, lo, hi)
                scaled = ((scaled - lo) * 255.0 / (hi - lo)).astype(np.uint8)
                display_frame = scaled
        
        # === IMAGE ===
        self.imageReady.emit(display_frame, timestamp)
        
        # === HISTOGRAM + LEVELS (skip N-1 of N frames) ===
        do_histogram = (self._display_count % self.config.histogram_skip == 0)
        
        if do_histogram:
            # Use full-res frame for accurate histogram/levels
            if self.config.auto_levels and not self._manual_levels:
                lo, hi = auto_levels(frame, self.config.auto_levels_percentile)
                self._level_lo = lo
                self._level_hi = hi
                self.levelsComputed.emit(lo, hi)
                
            counts, edges = compute_histogram(frame, self.config.histogram_bins)
            self.histogramReady.emit(counts, edges, self._level_lo, self._level_hi)
        
        # === INTENSITY TRACE ===
        times, values = self.intensity_buffer.get_latest(
            self.config.intensity_trace_points
        )
        if len(times) > 0:
            self.intensityTraceReady.emit(times, values)
        
        # === FFT (if enabled) ===
        if self.config.fft_enabled and self._fft_processor is not None:
            # Use downsampled but NOT 8-bit converted frame for FFT accuracy
            magnitude = self._fft_processor.compute_magnitude(downsampled_frame)
            self.fftReady.emit(magnitude)
            
    def _report_stats(self):
        """Report display statistics."""
        now = time.time()
        elapsed = now - self._last_stats_time
        
        if elapsed > 0:
            self._actual_display_rate = self._display_count / elapsed
            
        stats = {
            'display_rate': self._actual_display_rate,
            'display_count': self._display_count,
            'target_rate': self.config.display_rate_hz,
        }
        self.statsReady.emit(stats)
        
        # Reset counter
        self._display_count = 0
        self._last_stats_time = now
        
    def get_actual_display_rate(self) -> float:
        """Get measured display rate."""
        return self._actual_display_rate


# =============================================================================
# IMAGE FFT PROCESSOR
# =============================================================================

class ImageFFTProcessor:
    """
    2D FFT processor for camera images.
    
    Features:
    - 2D FFT with window functions
    - Log magnitude display
    - Filtering (lowpass, highpass, bandpass, notch)
    - Apply filter to get filtered image
    """
    
    def __init__(self):
        self._window_cache: Dict[Tuple, np.ndarray] = {}
        self._last_fft: Optional[np.ndarray] = None
        self._last_magnitude: Optional[np.ndarray] = None
        self._current_mask: Optional[np.ndarray] = None
        
        # Settings
        self.window_type = "hann"
        self.log_scale = True
        
    def compute_fft(self, image: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute 2D FFT.
        
        Args:
            image: 2D grayscale image
            
        Returns:
            (fft_complex, magnitude) - complex FFT and magnitude for display
        """
        import scipy.fft as fft
        
        # Convert to float
        img = image.astype(np.float64)
        
        # Apply window
        window = self._get_2d_window(img.shape, self.window_type)
        windowed = img * window
        
        # Compute FFT
        fft_result = fft.fft2(windowed)
        fft_shifted = fft.fftshift(fft_result)
        
        # Store for filtering
        self._last_fft = fft_shifted
        
        # Compute magnitude
        magnitude = np.abs(fft_shifted)
        if self.log_scale:
            magnitude = np.log1p(magnitude)
            
        self._last_magnitude = magnitude
        
        return fft_shifted, magnitude
    
    def compute_magnitude(self, image: np.ndarray) -> np.ndarray:
        """Compute FFT and return magnitude only (for display)."""
        _, magnitude = self.compute_fft(image)
        return magnitude
    
    def create_circular_mask(self, shape: Tuple[int, int],
                              center: Tuple[int, int] = None,
                              radius: float = None,
                              invert: bool = False) -> np.ndarray:
        """
        Create circular mask for FFT filtering.
        
        Args:
            shape: Image shape (height, width)
            center: Center of mask (default: image center)
            radius: Radius in pixels
            invert: If True, mask outside circle (highpass)
            
        Returns:
            Binary mask (float)
        """
        h, w = shape
        if center is None:
            center = (h // 2, w // 2)
        if radius is None:
            radius = min(h, w) // 4
            
        y, x = np.ogrid[:h, :w]
        dist = np.sqrt((x - center[1])**2 + (y - center[0])**2)
        
        mask = dist <= radius
        if invert:
            mask = ~mask
            
        return mask.astype(np.float64)
    
    def create_bandpass_mask(self, shape: Tuple[int, int],
                              inner_radius: float,
                              outer_radius: float) -> np.ndarray:
        """Create bandpass filter mask."""
        h, w = shape
        cy, cx = h // 2, w // 2
        y, x = np.ogrid[:h, :w]
        dist = np.sqrt((x - cx)**2 + (y - cy)**2)
        
        mask = (dist >= inner_radius) & (dist <= outer_radius)
        return mask.astype(np.float64)
    
    def create_notch_mask(self, shape: Tuple[int, int],
                           freq_points: list,  # List of (y, x) positions
                           radius: float = 5) -> np.ndarray:
        """
        Create notch filter mask to remove specific frequencies.
        
        Args:
            shape: Image shape
            freq_points: List of (y, x) frequency positions to notch out
            radius: Radius of each notch
            
        Returns:
            Mask with notches (1 everywhere except notch regions)
        """
        h, w = shape
        mask = np.ones((h, w), dtype=np.float64)
        cy, cx = h // 2, w // 2
        
        y, x = np.ogrid[:h, :w]
        
        for fy, fx in freq_points:
            # Notch at (fy, fx) relative to center
            dist1 = np.sqrt((x - cx - fx)**2 + (y - cy - fy)**2)
            # Also notch conjugate
            dist2 = np.sqrt((x - cx + fx)**2 + (y - cy + fy)**2)
            
            mask *= (dist1 > radius).astype(np.float64)
            mask *= (dist2 > radius).astype(np.float64)
            
        return mask
    
    def apply_filter(self, mask: np.ndarray = None) -> Optional[np.ndarray]:
        """
        Apply filter mask to last FFT and inverse transform.
        
        Args:
            mask: Filter mask (or use stored mask if None)
            
        Returns:
            Filtered image
        """
        import scipy.fft as fft
        
        if self._last_fft is None:
            return None
            
        if mask is None:
            mask = self._current_mask
        if mask is None:
            return None
            
        # Apply mask
        filtered_fft = self._last_fft * mask
        
        # Inverse FFT
        result = fft.ifft2(fft.ifftshift(filtered_fft))
        
        return np.abs(result)
    
    def set_mask(self, mask: np.ndarray):
        """Store mask for repeated filtering."""
        self._current_mask = mask
        
    def _get_2d_window(self, shape: Tuple[int, int], 
                        window_type: str) -> np.ndarray:
        """Get cached 2D separable window."""
        key = (shape, window_type)
        
        if key not in self._window_cache:
            h, w = shape
            
            if window_type == "hann":
                win_h = np.hanning(h)
                win_w = np.hanning(w)
            elif window_type == "hamming":
                win_h = np.hamming(h)
                win_w = np.hamming(w)
            elif window_type == "blackman":
                win_h = np.blackman(h)
                win_w = np.blackman(w)
            else:  # none/rectangular
                win_h = np.ones(h)
                win_w = np.ones(w)
                
            self._window_cache[key] = np.outer(win_h, win_w)
            
        return self._window_cache[key]
    
    def get_last_magnitude(self) -> Optional[np.ndarray]:
        """Get last computed magnitude."""
        return self._last_magnitude
    
    def clear_cache(self):
        """Clear window cache."""
        self._window_cache.clear()


# =============================================================================
# DISK WRITER (Abstract base + HDF5 implementation)
# =============================================================================

class DiskWriter:
    """Abstract base class for disk streaming."""
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self._is_open = False
        
    def open(self):
        """Open file for writing."""
        raise NotImplementedError
        
    def close(self):
        """Close file."""
        raise NotImplementedError
        
    def write_frame(self, frame: np.ndarray, timestamp: float, 
                    intensity: float, frame_number: int):
        """Write single frame."""
        raise NotImplementedError
        
    @property
    def is_open(self) -> bool:
        return self._is_open
        
    def __enter__(self):
        self.open()
        return self
        
    def __exit__(self, *args):
        self.close()


class HDF5Writer(DiskWriter):
    """
    HDF5 disk writer for frame streaming.
    
    Creates datasets:
    - /frames: (N, H, W) uint16 frames
    - /timestamps: (N,) float64 timestamps
    - /intensities: (N,) float64 intensities
    """
    
    def __init__(self, filepath: str, 
                 frame_shape: Tuple[int, int] = None,
                 dtype: np.dtype = np.uint16,
                 compression: bool = True,
                 chunk_frames: int = 100):
        super().__init__(filepath)
        self.frame_shape = frame_shape
        self.dtype = dtype
        self.compression = compression
        self.chunk_frames = chunk_frames
        
        self._file = None
        self._frames_dset = None
        self._timestamps_dset = None
        self._intensities_dset = None
        self._write_index = 0
        
    def open(self):
        """Open HDF5 file and create datasets."""
        import h5py
        
        self._file = h5py.File(self.filepath, 'w')
        
        # Create datasets with unlimited first dimension
        if self.frame_shape is not None:
            maxshape = (None,) + self.frame_shape
            chunks = (self.chunk_frames,) + self.frame_shape
            
            compression = 'gzip' if self.compression else None
            
            self._frames_dset = self._file.create_dataset(
                'frames',
                shape=(0,) + self.frame_shape,
                maxshape=maxshape,
                dtype=self.dtype,
                chunks=chunks,
                compression=compression
            )
        
        self._timestamps_dset = self._file.create_dataset(
            'timestamps',
            shape=(0,),
            maxshape=(None,),
            dtype=np.float64,
            chunks=(self.chunk_frames,)
        )
        
        self._intensities_dset = self._file.create_dataset(
            'intensities',
            shape=(0,),
            maxshape=(None,),
            dtype=np.float64,
            chunks=(self.chunk_frames,)
        )
        
        self._write_index = 0
        self._is_open = True
        
    def close(self):
        """Close HDF5 file."""
        if self._file is not None:
            self._file.close()
            self._file = None
        self._is_open = False
        
    def write_frame(self, frame: np.ndarray, timestamp: float,
                    intensity: float, frame_number: int):
        """Write frame to HDF5."""
        if not self._is_open:
            return
            
        # Initialize frames dataset if needed
        if self._frames_dset is None and frame is not None:
            self.frame_shape = frame.shape
            maxshape = (None,) + self.frame_shape
            chunks = (self.chunk_frames,) + self.frame_shape
            compression = 'gzip' if self.compression else None
            
            self._frames_dset = self._file.create_dataset(
                'frames',
                shape=(0,) + self.frame_shape,
                maxshape=maxshape,
                dtype=frame.dtype,
                chunks=chunks,
                compression=compression
            )
        
        # Resize and write
        n = self._write_index + 1
        
        if self._frames_dset is not None:
            self._frames_dset.resize((n,) + self.frame_shape)
            self._frames_dset[self._write_index] = frame
            
        self._timestamps_dset.resize((n,))
        self._timestamps_dset[self._write_index] = timestamp
        
        self._intensities_dset.resize((n,))
        self._intensities_dset[self._write_index] = intensity
        
        self._write_index = n


class TIFFWriter(DiskWriter):
    """
    TIFF stack writer for frame streaming.
    
    Uses tifffile for BigTIFF support (>4GB files).
    """
    
    def __init__(self, filepath: str, 
                 bigtiff: bool = True,
                 compression: str = None):
        super().__init__(filepath)
        self.bigtiff = bigtiff
        self.compression = compression
        
        self._tiff = None
        self._timestamps = []
        self._intensities = []
        
    def open(self):
        """Open TIFF file."""
        try:
            import tifffile
            self._tiff = tifffile.TiffWriter(
                self.filepath, 
                bigtiff=self.bigtiff
            )
            self._timestamps = []
            self._intensities = []
            self._is_open = True
        except ImportError:
            logging.error("tifffile not available")
            raise
        
    def close(self):
        """Close TIFF and save metadata."""
        if self._tiff is not None:
            self._tiff.close()
            self._tiff = None
            
            # Save timestamps/intensities to companion file
            metadata_path = self.filepath.replace('.tif', '_metadata.npz')
            np.savez(metadata_path,
                     timestamps=np.array(self._timestamps),
                     intensities=np.array(self._intensities))
                     
        self._is_open = False
        
    def write_frame(self, frame: np.ndarray, timestamp: float,
                    intensity: float, frame_number: int):
        """Write frame to TIFF stack."""
        if not self._is_open or self._tiff is None:
            return
            
        self._tiff.write(frame, compression=self.compression)
        self._timestamps.append(timestamp)
        self._intensities.append(intensity)


# =============================================================================
# CONVENIENCE FACTORY FUNCTIONS
# =============================================================================

def create_acquisition_system(camera: CameraInterface,
                               acquisition_config: AcquisitionConfig = None,
                               display_config: DisplayManagerConfig = None):
    """
    Create complete acquisition system with buffers.
    
    Args:
        camera: Camera backend instance
        acquisition_config: Acquisition configuration
        display_config: Display manager configuration
        
    Returns:
        (acq_thread, display_manager, frame_buffer, intensity_buffer)
    """
    acq_config = acquisition_config or AcquisitionConfig()
    disp_config = display_config or DisplayManagerConfig()
    
    # Create buffers
    frame_buffer = LatestFrameBuffer()
    intensity_buffer = IntensityRingBuffer(
        capacity=acq_config.intensity_buffer_size
    )
    
    # Create acquisition thread
    acq_thread = AcquisitionThread(
        camera=camera,
        frame_buffer=frame_buffer,
        intensity_buffer=intensity_buffer,
        config=acq_config
    )
    
    # Create display manager
    display_manager = DisplayManager(
        frame_buffer=frame_buffer,
        intensity_buffer=intensity_buffer,
        config=disp_config
    )
    
    return acq_thread, display_manager, frame_buffer, intensity_buffer


def create_disk_writer(storage_config: StorageConfig,
                        frame_shape: Tuple[int, int] = None) -> Optional[DiskWriter]:
    """
    Create disk writer based on storage config.
    
    Args:
        storage_config: Storage configuration
        frame_shape: Expected frame shape
        
    Returns:
        DiskWriter instance or None
    """
    if storage_config.mode == StorageMode.NONE:
        return None
        
    if storage_config.mode == StorageMode.DISK_HDF5:
        return HDF5Writer(
            filepath=storage_config.save_path,
            frame_shape=frame_shape,
            compression=storage_config.compression,
            chunk_frames=storage_config.chunk_size
        )
        
    if storage_config.mode == StorageMode.DISK_TIFF:
        return TIFFWriter(
            filepath=storage_config.save_path
        )
        
    return None
