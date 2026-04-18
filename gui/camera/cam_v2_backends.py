"""
Camera Viewer v2 - Backend Module
=================================

Camera interface implementations:
- SimulatedBackend: For testing without hardware
- HamamatsuBackend: Using dcam.py
- ThorlabsBackend: Using pylablib
- AndorBackend: Using pylablib SDK2 (for iXon Ultra 888)

Author: Claude/BP
Date: January 2025
"""

import numpy as np
import time
import logging
from typing import Optional, Dict, Tuple, Any, List
from dataclasses import dataclass

from cam_v2_core import (
    CameraInterface, CameraType, TriggerMode, TriggerPolarity
)


# =============================================================================
# SIMULATED CAMERA BACKEND
# =============================================================================

class SimulatedBackend(CameraInterface):
    """
    Simulated camera for testing without hardware.
    
    Generates test patterns with configurable noise and features.
    """
    
    def __init__(self, width: int = 512, height: int = 512, 
                 dtype: np.dtype = np.uint16, device_index: int = 0):
        super().__init__()
        self._width = width
        self._height = height
        self._dtype = dtype
        self._device_index = device_index  # Ignored but accepted for API consistency
        self._exposure_ms = 10.0
        self._gain = 1.0
        self._roi = None
        self._binning = 1
        self._trigger_mode = TriggerMode.INTERNAL
        
        # Test pattern parameters
        self._spot_x = width // 2
        self._spot_y = height // 2
        self._spot_sigma = 30
        self._spot_intensity = 50000 if dtype == np.uint16 else 200
        self._noise_level = 100 if dtype == np.uint16 else 10
        self._background = 500 if dtype == np.uint16 else 20
        
        # Animation
        self._time_offset = 0.0
        self._move_spot = True
        
    def open(self) -> bool:
        self._is_open = True
        logging.info("SimulatedBackend: Opened")
        return True
    
    def close(self):
        self.stop_acquisition()
        self._is_open = False
        logging.info("SimulatedBackend: Closed")
        
    def start_acquisition(self):
        if not self._is_open:
            return
        self._is_acquiring = True
        self._frame_count = 0
        self._time_offset = time.time()
        logging.info("SimulatedBackend: Started acquisition")
        
    def stop_acquisition(self):
        self._is_acquiring = False
        logging.info("SimulatedBackend: Stopped acquisition")
        
    def get_frame(self, timeout_ms: int = 1000) -> Optional[np.ndarray]:
        if not self._is_acquiring:
            return None
        
        # Simulate exposure time
        time.sleep(self._exposure_ms / 1000.0)
        
        # Generate frame
        frame = self._generate_frame()
        self._frame_count += 1
        
        return frame
    
    def _generate_frame(self) -> np.ndarray:
        """Generate test pattern with moving Gaussian spot."""
        h, w = self._height, self._width
        
        # Apply binning
        if self._binning > 1:
            h = h // self._binning
            w = w // self._binning
        
        # Apply ROI
        if self._roi is not None:
            x0, y0, roi_w, roi_h = self._roi
            # Clip to valid range
            roi_w = min(roi_w, w - x0)
            roi_h = min(roi_h, h - y0)
        else:
            x0, y0 = 0, 0
            roi_w, roi_h = w, h
        
        # Create coordinate grids
        y_coords, x_coords = np.ogrid[:roi_h, :roi_w]
        x_coords = x_coords + x0
        y_coords = y_coords + y0
        
        # Moving spot position
        t = time.time() - self._time_offset
        if self._move_spot:
            cx = w // 2 + int(100 * np.sin(t * 0.5))
            cy = h // 2 + int(100 * np.cos(t * 0.7))
        else:
            cx, cy = self._spot_x, self._spot_y
        
        # Gaussian spot
        sigma = self._spot_sigma / self._binning
        r2 = (x_coords - cx)**2 + (y_coords - cy)**2
        spot = self._spot_intensity * np.exp(-r2 / (2 * sigma**2))
        
        # Background + spot + noise
        frame = self._background + spot
        frame += np.random.normal(0, self._noise_level, (roi_h, roi_w))
        
        # Clip and convert
        if self._dtype == np.uint16:
            frame = np.clip(frame, 0, 65535)
        else:
            frame = np.clip(frame, 0, 255)
        
        return frame.astype(self._dtype)
    
    def set_exposure(self, exposure_ms: float):
        self._exposure_ms = max(0.1, exposure_ms)
        
    def get_exposure(self) -> float:
        return self._exposure_ms
    
    def set_roi(self, x: int, y: int, width: int, height: int):
        self._roi = (x, y, width, height)
        
    def get_roi(self) -> Tuple[int, int, int, int]:
        if self._roi is not None:
            return self._roi
        return (0, 0, self._width, self._height)
    
    def clear_roi(self):
        self._roi = None
        
    def set_binning(self, binning: int):
        self._binning = max(1, binning)
        
    def get_binning(self) -> int:
        return self._binning
    
    def set_trigger_mode(self, mode: TriggerMode):
        self._trigger_mode = mode
        
    def get_trigger_mode(self) -> TriggerMode:
        return self._trigger_mode
    
    def set_gain(self, gain: float):
        self._gain = gain
        # Adjust noise/spot based on gain
        base_noise = 100 if self._dtype == np.uint16 else 10
        self._noise_level = base_noise * gain
        
    def get_gain(self) -> float:
        return self._gain
    
    @property
    def frame_shape(self) -> Tuple[int, int]:
        h, w = self._height, self._width
        if self._binning > 1:
            h = h // self._binning
            w = w // self._binning
        if self._roi is not None:
            return (self._roi[3], self._roi[2])  # (height, width)
        return (h, w)
    
    @property
    def pixel_dtype(self) -> np.dtype:
        return self._dtype
    
    @property
    def is_color(self) -> bool:
        return False
    
    @property
    def camera_info(self) -> Dict[str, Any]:
        return {
            'model': 'Simulated Camera',
            'serial': 'SIM-001',
            'firmware': '1.0.0',
            'driver_version': '1.0.0',
            'sensor_size': (self._width, self._height),
            'pixel_size': 6.5,  # µm
            'bit_depth': 16 if self._dtype == np.uint16 else 8,
            'min_exposure': 0.1,
            'max_exposure': 10000.0,
            'temperature': None
        }
    
    # === Simulation Controls ===
    
    def set_spot_position(self, x: int, y: int):
        """Set spot center position."""
        self._spot_x = x
        self._spot_y = y
        
    def set_spot_size(self, sigma: float):
        """Set spot sigma in pixels."""
        self._spot_sigma = sigma
        
    def set_noise_level(self, level: float):
        """Set noise standard deviation."""
        self._noise_level = level
        
    def enable_spot_motion(self, enable: bool):
        """Enable/disable spot motion."""
        self._move_spot = enable


# =============================================================================
# HAMAMATSU BACKEND (using dcam.py)
# =============================================================================

# Try to import dcam
DCAM_AVAILABLE = False
dcam = None
dcamapi4 = None

try:
    import sys
    import os
    # Add paths where dcam.py might be
    sys.path.append(os.path.abspath(r'D:\Brateen\NV_Experiment'))
    
    import dcam as _dcam
    import dcamapi4 as _dcamapi4
    dcam = _dcam
    dcamapi4 = _dcamapi4
    DCAM_AVAILABLE = True
except (ImportError, OSError) as e:
    logging.info(f"DCAM libraries not available: {e}. HamamatsuBackend disabled.")


class HamamatsuBackend(CameraInterface):
    """
    Hamamatsu camera backend using dcam.py/dcamapi4.py.
    
    Supports ORCA series cameras like C13440-20CU.
    """
    
    def __init__(self, device_index: int = 0):
        super().__init__()
        self._device_index = device_index
        self._dcam = None
        self._exposure_ms = 10.0
        self._roi = None
        self._binning = 1
        self._trigger_mode = TriggerMode.INTERNAL
        self._api_initialized = False
        
    def open(self) -> bool:
        if not DCAM_AVAILABLE:
            logging.error("DCAM not available")
            return False
        
        try:
            # Initialize API if not already done
            if not self._api_initialized:
                if not dcam.Dcamapi.init():
                    logging.error("Failed to initialize DCAM API")
                    return False
                self._api_initialized = True
            
            # Get number of cameras
            n_cameras = dcam.Dcamapi.get_devicecount()
            if n_cameras == 0:
                logging.error("No Hamamatsu cameras found")
                return False
            
            if self._device_index >= n_cameras:
                logging.error(f"Camera index {self._device_index} out of range (found {n_cameras})")
                return False
            
            # Open camera
            self._dcam = dcam.Dcam(self._device_index)
            if not self._dcam.dev_open():
                logging.error("Failed to open Hamamatsu camera")
                return False
            
            self._is_open = True
            logging.info(f"HamamatsuBackend: Opened camera {self._device_index}")
            
            # Set default exposure
            self.set_exposure(self._exposure_ms)
            
            return True
            
        except Exception as e:
            logging.exception(f"Error opening Hamamatsu camera: {e}")
            return False
    
    def close(self):
        if self._dcam is not None:
            try:
                self.stop_acquisition()
                self._dcam.dev_close()
            except:
                pass
            self._dcam = None
        self._is_open = False
        logging.info("HamamatsuBackend: Closed")
        
    def start_acquisition(self):
        if not self._is_open or self._dcam is None:
            return
        
        try:
            # Allocate buffers
            self._dcam.buf_alloc(10)  # 10 frame buffer
            
            # Start capture
            self._dcam.cap_start()
            
            self._is_acquiring = True
            self._frame_count = 0
            logging.info("HamamatsuBackend: Started acquisition")
            
        except Exception as e:
            logging.exception(f"Error starting acquisition: {e}")
    
    def stop_acquisition(self):
        if self._dcam is None:
            return
        
        try:
            self._dcam.cap_stop()
            self._dcam.buf_release()
        except:
            pass
        
        self._is_acquiring = False
        logging.info("HamamatsuBackend: Stopped acquisition")
        
    def get_frame(self, timeout_ms: int = 1000) -> Optional[np.ndarray]:
        if not self._is_acquiring or self._dcam is None:
            return None
        
        try:
            # Wait for frame
            if self._dcam.wait_capevent_frameready(timeout_ms):
                frame = self._dcam.buf_getlastframedata()
                if frame is not None:
                    self._frame_count += 1
                    return frame.copy()
        except Exception as e:
            logging.debug(f"Error getting frame: {e}")
        
        return None
    
    def set_exposure(self, exposure_ms: float):
        self._exposure_ms = exposure_ms
        if self._dcam is not None:
            try:
                # DCAM uses seconds
                self._dcam.prop_setvalue(
                    dcamapi4.DCAM_IDPROP.EXPOSURETIME,
                    exposure_ms / 1000.0
                )
            except Exception as e:
                logging.warning(f"Error setting exposure: {e}")
                
    def get_exposure(self) -> float:
        if self._dcam is not None:
            try:
                exp_s = self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.EXPOSURETIME)
                return exp_s * 1000.0
            except:
                pass
        return self._exposure_ms
    
    def set_roi(self, x: int, y: int, width: int, height: int):
        self._roi = (x, y, width, height)
        if self._dcam is not None:
            try:
                # Hamamatsu uses subarray properties
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SUBARRAYHPOS, x)
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SUBARRAYVPOS, y)
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SUBARRAYHSIZE, width)
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SUBARRAYVSIZE, height)
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SUBARRAYMODE, 2)  # ON
            except Exception as e:
                logging.warning(f"Error setting ROI: {e}")
                
    def get_roi(self) -> Tuple[int, int, int, int]:
        if self._dcam is not None:
            try:
                x = int(self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.SUBARRAYHPOS))
                y = int(self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.SUBARRAYVPOS))
                w = int(self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.SUBARRAYHSIZE))
                h = int(self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.SUBARRAYVSIZE))
                return (x, y, w, h)
            except:
                pass
        if self._roi is not None:
            return self._roi
        return (0, 0, 2048, 2048)  # Default for ORCA
    
    def clear_roi(self):
        self._roi = None
        if self._dcam is not None:
            try:
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SUBARRAYMODE, 1)  # OFF
            except:
                pass
                
    def set_binning(self, binning: int):
        self._binning = binning
        if self._dcam is not None:
            try:
                self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.BINNING, binning)
            except Exception as e:
                logging.warning(f"Error setting binning: {e}")
                
    def get_binning(self) -> int:
        if self._dcam is not None:
            try:
                return int(self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.BINNING))
            except:
                pass
        return self._binning
    
    def set_trigger_mode(self, mode: TriggerMode):
        self._trigger_mode = mode
        if self._dcam is not None:
            try:
                if mode == TriggerMode.INTERNAL:
                    self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.TRIGGERSOURCE, 1)  # Internal
                elif mode == TriggerMode.EXTERNAL:
                    self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.TRIGGERSOURCE, 2)  # External
                elif mode == TriggerMode.SOFTWARE:
                    self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.TRIGGERSOURCE, 3)  # Software
            except Exception as e:
                logging.warning(f"Error setting trigger mode: {e}")
                
    def get_trigger_mode(self) -> TriggerMode:
        return self._trigger_mode
    
    def set_cooling(self, enabled: bool, target_temp_c: float = -10.0):
        if self._dcam is not None:
            try:
                if enabled:
                    self._dcam.prop_setvalue(dcamapi4.DCAM_IDPROP.SENSORTEMPERATURETARGET, 
                                             target_temp_c)
            except Exception as e:
                logging.warning(f"Error setting cooling: {e}")
                
    def get_temperature(self) -> Optional[float]:
        if self._dcam is not None:
            try:
                return self._dcam.prop_getvalue(dcamapi4.DCAM_IDPROP.SENSORTEMPERATURE)
            except:
                pass
        return None
    
    @property
    def frame_shape(self) -> Tuple[int, int]:
        roi = self.get_roi()
        return (roi[3], roi[2])  # (height, width)
    
    @property
    def pixel_dtype(self) -> np.dtype:
        return np.uint16
    
    @property
    def is_color(self) -> bool:
        return False
    
    @property
    def camera_info(self) -> Dict[str, Any]:
        info = {
            'model': 'Hamamatsu ORCA',
            'serial': 'Unknown',
            'firmware': 'Unknown',
            'driver_version': 'Unknown',
            'sensor_size': (2048, 2048),
            'pixel_size': 6.5,
            'bit_depth': 16,
            'min_exposure': 0.01,
            'max_exposure': 10000.0,
            'temperature': self.get_temperature()
        }
        
        if self._dcam is not None:
            try:
                info['model'] = self._dcam.dev_getstring(dcamapi4.DCAM_IDSTR.MODEL)
                info['serial'] = self._dcam.dev_getstring(dcamapi4.DCAM_IDSTR.CAMERAID)
            except:
                pass
        
        return info


# =============================================================================
# THORLABS BACKEND (using pylablib)
# =============================================================================

try:
    from pylablib.devices import Thorlabs
    THORLABS_AVAILABLE = True
except ImportError:
    THORLABS_AVAILABLE = False
    logging.info("pylablib Thorlabs not available. ThorlabsBackend disabled.")


class ThorlabsBackend(CameraInterface):
    """
    Thorlabs camera backend using pylablib TLCamera SDK.
    
    Supports scientific cameras like CS165CU.
    For DCx/UC480 cameras (DCC1545M, DCC1645C), use UC480Backend instead.
    """
    
    def __init__(self, serial: Optional[str] = None, device_index: int = 0):
        super().__init__()
        self._serial = serial
        self._device_index = device_index
        self._cam = None
        self._exposure_ms = 10.0
        self._roi = None
        self._binning = 1
        self._trigger_mode = TriggerMode.INTERNAL
        
    def open(self) -> bool:
        if not THORLABS_AVAILABLE:
            logging.error("pylablib Thorlabs not available")
            return False
        
        try:
            # List available cameras
            cameras = Thorlabs.list_cameras_tlcam()
            if not cameras:
                logging.error("No Thorlabs TLCamera cameras found")
                return False
            
            # Open camera by serial, index, or first available
            if self._serial:
                self._cam = Thorlabs.ThorlabsTLCamera(serial=self._serial)
            elif self._device_index < len(cameras):
                self._cam = Thorlabs.ThorlabsTLCamera(serial=cameras[self._device_index])
            else:
                self._cam = Thorlabs.ThorlabsTLCamera()
            
            self._is_open = True
            logging.info(f"ThorlabsBackend: Opened camera {self._device_index}")
            
            # Set default exposure
            self.set_exposure(self._exposure_ms)
            
            return True
            
        except Exception as e:
            logging.exception(f"Error opening Thorlabs camera: {e}")
            return False
    
    def close(self):
        if self._cam is not None:
            try:
                self.stop_acquisition()
                self._cam.close()
            except:
                pass
            self._cam = None
        self._is_open = False
        logging.info("ThorlabsBackend: Closed")
        
    def start_acquisition(self):
        if not self._is_open or self._cam is None:
            return
        
        try:
            self._cam.start_acquisition()
            self._is_acquiring = True
            self._frame_count = 0
            logging.info("ThorlabsBackend: Started acquisition")
        except Exception as e:
            logging.exception(f"Error starting acquisition: {e}")
    
    def stop_acquisition(self):
        if self._cam is None:
            return
        
        try:
            self._cam.stop_acquisition()
        except:
            pass
        
        self._is_acquiring = False
        logging.info("ThorlabsBackend: Stopped acquisition")
        
    def get_frame(self, timeout_ms: int = 1000) -> Optional[np.ndarray]:
        if not self._is_acquiring or self._cam is None:
            return None
        
        try:
            # Wait for frame
            self._cam.wait_for_frame(timeout=timeout_ms/1000.0)
            frame = self._cam.read_newest_image()
            if frame is not None:
                self._frame_count += 1
                return frame
        except Exception as e:
            logging.debug(f"Error getting frame: {e}")
        
        return None
    
    def set_exposure(self, exposure_ms: float):
        self._exposure_ms = exposure_ms
        if self._cam is not None:
            try:
                self._cam.set_exposure(exposure_ms / 1000.0)  # pylablib uses seconds
            except Exception as e:
                logging.warning(f"Error setting exposure: {e}")
                
    def get_exposure(self) -> float:
        if self._cam is not None:
            try:
                return self._cam.get_exposure() * 1000.0
            except:
                pass
        return self._exposure_ms
    
    def set_roi(self, x: int, y: int, width: int, height: int):
        self._roi = (x, y, width, height)
        if self._cam is not None:
            try:
                self._cam.set_roi(x, x+width, y, y+height)
            except Exception as e:
                logging.warning(f"Error setting ROI: {e}")
                
    def get_roi(self) -> Tuple[int, int, int, int]:
        if self._cam is not None:
            try:
                roi = self._cam.get_roi()
                return (roi[0], roi[2], roi[1]-roi[0], roi[3]-roi[2])
            except:
                pass
        if self._roi is not None:
            return self._roi
        return (0, 0, 1024, 1024)
    
    def clear_roi(self):
        self._roi = None
        if self._cam is not None:
            try:
                self._cam.set_roi()  # Reset to full
            except:
                pass
                
    def set_binning(self, binning: int):
        self._binning = binning
        # Note: Thorlabs cameras may not support binning through this interface
                
    def get_binning(self) -> int:
        return self._binning
    
    def set_trigger_mode(self, mode: TriggerMode):
        self._trigger_mode = mode
        # Note: Trigger configuration varies by Thorlabs model
                
    def get_trigger_mode(self) -> TriggerMode:
        return self._trigger_mode
    
    def set_gain(self, gain: float):
        if self._cam is not None:
            try:
                self._cam.set_gain(gain)
            except:
                pass
    
    def get_gain(self) -> float:
        if self._cam is not None:
            try:
                return self._cam.get_gain()
            except:
                pass
        return 1.0
    
    @property
    def frame_shape(self) -> Tuple[int, int]:
        if self._cam is not None:
            try:
                size = self._cam.get_detector_size()
                return (size[1], size[0])  # (height, width)
            except:
                pass
        roi = self.get_roi()
        return (roi[3], roi[2])
    
    @property
    def pixel_dtype(self) -> np.dtype:
        # Most Thorlabs cameras are 8 or 12 bit
        return np.uint16
    
    @property
    def is_color(self) -> bool:
        if self._cam is not None:
            try:
                # Check if color camera
                size = self._cam.get_detector_size()
                return len(size) > 2 and size[2] > 1
            except:
                pass
        return False
    
    @property
    def camera_info(self) -> Dict[str, Any]:
        info = {
            'model': 'Thorlabs Camera',
            'serial': self._serial or 'Unknown',
            'firmware': 'Unknown',
            'driver_version': 'Unknown',
            'sensor_size': self.frame_shape[::-1],
            'pixel_size': 5.2,
            'bit_depth': 12,
            'min_exposure': 0.01,
            'max_exposure': 10000.0,
            'temperature': None
        }
        
        if self._cam is not None:
            try:
                info['model'] = self._cam.get_device_info().get('model', 'Unknown')
                info['serial'] = self._cam.get_device_info().get('serial_number', 'Unknown')
            except:
                pass
        
        return info


# =============================================================================
# UC480 BACKEND (for DCC1545M, DCC1645C, and similar IDS cameras)
# =============================================================================

try:
    from pylablib.devices import uc480
    UC480_AVAILABLE = True
except ImportError:
    UC480_AVAILABLE = False
    logging.info("pylablib uc480 not available. UC480Backend disabled.")


class UC480Backend(CameraInterface):
    """
    UC480/IDS camera backend using pylablib.
    
    Supports Thorlabs DCx series and IDS uEye cameras:
    - DCC1545M (mono)
    - DCC1645C (color)
    - Other UC480-compatible cameras
    
    Note: For Thorlabs scientific cameras (CS165CU, etc.), use ThorlabsBackend.
    """
    
    def __init__(self, device_index: int = 0):
        super().__init__()
        self._device_index = device_index
        self._cam = None
        self._exposure_ms = 10.0
        self._roi = None
        self._binning = 1
        self._trigger_mode = TriggerMode.INTERNAL
        self._is_color = False
        
    def open(self) -> bool:
        if not UC480_AVAILABLE:
            logging.error("pylablib uc480 not available")
            return False
        
        try:
            # List available cameras
            cameras = uc480.list_cameras()
            if not cameras:
                logging.error("No UC480 cameras found")
                return False
            
            if self._device_index >= len(cameras):
                logging.error(f"Camera index {self._device_index} out of range (found {len(cameras)})")
                return False
            
            # Get camera ID
            cam_id = cameras[self._device_index].cam_id
            # if isinstance(cam_id, dict):
            #     cam_id = cam_id.get('cam_id', cam_id.get('id', 0))
            # print(f"Camera id = {cam_id}")
            # Open camera
            self._cam = uc480.UC480Camera(cam_id=cam_id)
            
            # Check if color camera
            try:
                sensor_info = self._cam.get_sensor_info()
                color_mode = sensor_info.get('color_mode', 'mono')
                self._is_color = 'color' in str(color_mode).lower() or 'bayer' in str(color_mode).lower()
            except:
                self._is_color = False
            
            self._is_open = True
            logging.info(f"UC480Backend: Opened camera {self._device_index} (color={self._is_color})")
            
            # Set default exposure
            self.set_exposure(self._exposure_ms)
            
            return True
            
        except Exception as e:
            logging.exception(f"Error opening UC480 camera: {e}")
            return False
    
    def close(self):
        if self._cam is not None:
            try:
                self.stop_acquisition()
                self._cam.close()
            except:
                pass
            self._cam = None
        self._is_open = False
        logging.info("UC480Backend: Closed")
        
    def start_acquisition(self):
        if not self._is_open or self._cam is None:
            return
        
        try:
            self._cam.start_acquisition()
            self._is_acquiring = True
            self._frame_count = 0
            logging.info("UC480Backend: Started acquisition")
        except Exception as e:
            logging.exception(f"Error starting acquisition: {e}")
    
    def stop_acquisition(self):
        if self._cam is None:
            return
        
        try:
            self._cam.stop_acquisition()
        except:
            pass
        
        self._is_acquiring = False
        logging.info("UC480Backend: Stopped acquisition")
        
    def get_frame(self, timeout_ms: int = 1000) -> Optional[np.ndarray]:
        if not self._is_acquiring or self._cam is None:
            return None
        
        try:
            # Wait for frame
            self._cam.wait_for_frame(timeout=timeout_ms/1000.0)
            frame = self._cam.read_newest_image()
            
            if frame is not None:
                self._frame_count += 1
                
                # Transpose: UC480 returns (W,H) or (W,H,C), we need (H,W) or (H,W,C)
                if frame.ndim == 2:
                    frame = frame.T
                elif frame.ndim == 3:
                    frame = np.transpose(frame, (1, 0, 2))
                
                return frame
                
        except Exception as e:
            logging.debug(f"Error getting frame: {e}")
        
        return None
    
    def set_exposure(self, exposure_ms: float):
        self._exposure_ms = exposure_ms
        if self._cam is not None:
            try:
                self._cam.set_exposure(exposure_ms / 1000.0)  # pylablib uses seconds
            except Exception as e:
                logging.warning(f"Error setting exposure: {e}")
                
    def get_exposure(self) -> float:
        if self._cam is not None:
            try:
                return self._cam.get_exposure() * 1000.0
            except:
                pass
        return self._exposure_ms
    
    def set_roi(self, x: int, y: int, width: int, height: int):
        self._roi = (x, y, width, height)
        if self._cam is not None:
            try:
                self._cam.set_roi(x, x+width, y, y+height)
            except Exception as e:
                logging.warning(f"Error setting ROI: {e}")
    
    def get_roi(self) -> Tuple[int, int, int, int]:
        if self._cam is not None:
            try:
                roi = self._cam.get_roi()
                return (roi[0], roi[2], roi[1]-roi[0], roi[3]-roi[2])
            except:
                pass
        return self._roi or (0, 0, 1280, 1024)
    
    def clear_roi(self):
        self._roi = None
        if self._cam is not None:
            try:
                self._cam.set_roi()  # Reset to full
            except:
                pass
    
    def set_binning(self, binning: int):
        self._binning = binning
        # UC480 binning is often limited or not available
        if self._cam is not None:
            try:
                # Try subsampling instead of binning
                self._cam.set_subsampling(binning, binning)
            except Exception as e:
                logging.warning(f"Binning not supported on UC480: {e}")
                
    def get_binning(self) -> int:
        return self._binning
    
    def set_trigger_mode(self, mode: TriggerMode):
        self._trigger_mode = mode
        if self._cam is not None:
            try:
                if mode == TriggerMode.INTERNAL:
                    self._cam.set_trigger_mode("int")
                elif mode == TriggerMode.EXTERNAL:
                    self._cam.set_trigger_mode("ext")
                elif mode == TriggerMode.SOFTWARE:
                    self._cam.set_trigger_mode("soft")
            except Exception as e:
                logging.warning(f"Error setting trigger mode: {e}")
                
    def get_trigger_mode(self) -> TriggerMode:
        return self._trigger_mode
    
    def set_gain(self, gain: float):
        if self._cam is not None:
            try:
                # UC480 uses master gain (0-100 typically)
                self._cam.set_gain(int(gain))
            except Exception as e:
                logging.warning(f"Error setting gain: {e}")
                
    def get_gain(self) -> float:
        if self._cam is not None:
            try:
                return float(self._cam.get_gain())
            except:
                pass
        return 1.0
    
    def get_temperature(self) -> Optional[float]:
        # UC480 cameras typically don't have temperature readout
        return None
    
    @property
    def frame_shape(self) -> Tuple[int, int]:
        roi = self.get_roi()
        h, w = roi[3], roi[2]
        if self._is_color:
            return (h, w, 3)
        return (h, w)
    
    @property
    def pixel_dtype(self) -> np.dtype:
        # UC480 cameras are typically 8 or 10/12 bit
        return np.uint8 if not self._is_color else np.uint8
    
    @property
    def is_color(self) -> bool:
        return self._is_color
    
    @property
    def camera_info(self) -> Dict[str, Any]:
        info = {
            'model': 'UC480 Camera',
            'serial': 'Unknown',
            'firmware': 'Unknown',
            'driver_version': 'Unknown',
            'sensor_size': (1280, 1024),
            'pixel_size': 5.2,
            'bit_depth': 8,
            'min_exposure': 0.01,
            'max_exposure': 10000.0,
            'temperature': None,
            'is_color': self._is_color
        }
        
        if self._cam is not None:
            try:
                sensor_info = self._cam.get_sensor_info()
                info['model'] = sensor_info.get('sensor_name', 'UC480')
                info['sensor_size'] = (
                    sensor_info.get('max_width', 1280),
                    sensor_info.get('max_height', 1024)
                )
                info['pixel_size'] = sensor_info.get('pixel_size', 5.2)
            except:
                pass
            
            try:
                cam_info = self._cam.get_camera_info()
                info['serial'] = cam_info.get('serial_number', 'Unknown')
            except:
                pass
        
        return info


# =============================================================================
# ANDOR BACKEND (using pylablib SDK2 for iXon)
# =============================================================================

try:
    from pylablib.devices import Andor
    ANDOR_AVAILABLE = True
except ImportError:
    ANDOR_AVAILABLE = False
    logging.info("pylablib Andor not available. AndorBackend disabled.")


class AndorBackend(CameraInterface):
    """
    Andor EMCCD camera backend using pylablib SDK2.
    
    Supports iXon series including Ultra 888.
    """
    
    def __init__(self, device_index: int = 0, 
                 temperature: float = -70.0,
                 fan_mode: str = "on"):
        super().__init__()
        self._device_index = device_index
        self._target_temp = temperature
        self._fan_mode = fan_mode
        self._cam = None
        self._exposure_ms = 10.0
        self._roi = None
        self._binning = 1
        self._trigger_mode = TriggerMode.INTERNAL
        self._emccd_gain = 1
        
    def open(self) -> bool:
        if not ANDOR_AVAILABLE:
            logging.error("pylablib Andor not available")
            return False
        
        try:
            # Check number of cameras
            n_cameras = Andor.get_cameras_number_SDK2()
            if n_cameras == 0:
                logging.error("No Andor cameras found")
                return False
            
            if self._device_index >= n_cameras:
                logging.error(f"Camera index {self._device_index} out of range")
                return False
            
            # Open camera with initial settings
            self._cam = Andor.AndorSDK2Camera(
                idx=self._device_index,
                temperature=self._target_temp,
                fan_mode=self._fan_mode
            )
            
            # Set shutter to open for acquisition
            self._cam.set_shutter("open")
            
            self._is_open = True
            logging.info(f"AndorBackend: Opened camera {self._device_index}")
            
            # Set default exposure
            self.set_exposure(self._exposure_ms)
            
            return True
            
        except Exception as e:
            logging.exception(f"Error opening Andor camera: {e}")
            return False
    
    def close(self):
        if self._cam is not None:
            try:
                self.stop_acquisition()
                # Important: must close to prevent DLL resource leak
                self._cam.close()
            except:
                pass
            self._cam = None
        self._is_open = False
        logging.info("AndorBackend: Closed")
        
    def start_acquisition(self):
        if not self._is_open or self._cam is None:
            return
        
        try:
            # Setup continuous acquisition
            self._cam.set_acquisition_mode("cont")
            self._cam.setup_acquisition()
            self._cam.start_acquisition()
            
            self._is_acquiring = True
            self._frame_count = 0
            logging.info("AndorBackend: Started acquisition")
        except Exception as e:
            logging.exception(f"Error starting acquisition: {e}")
    
    def stop_acquisition(self):
        if self._cam is None:
            return
        
        try:
            self._cam.stop_acquisition()
        except:
            pass
        
        self._is_acquiring = False
        logging.info("AndorBackend: Stopped acquisition")
        
    def get_frame(self, timeout_ms: int = 1000) -> Optional[np.ndarray]:
        if not self._is_acquiring or self._cam is None:
            return None
        
        try:
            # Wait for frame
            self._cam.wait_for_frame(timeout=timeout_ms/1000.0)
            frame = self._cam.read_newest_image()
            if frame is not None:
                self._frame_count += 1
                return frame
        except Exception as e:
            logging.debug(f"Error getting frame: {e}")
        
        return None
    
    def set_exposure(self, exposure_ms: float):
        self._exposure_ms = exposure_ms
        if self._cam is not None:
            try:
                self._cam.set_exposure(exposure_ms / 1000.0)  # seconds
            except Exception as e:
                logging.warning(f"Error setting exposure: {e}")
                
    def get_exposure(self) -> float:
        if self._cam is not None:
            try:
                return self._cam.get_exposure() * 1000.0
            except:
                pass
        return self._exposure_ms
    
    def set_roi(self, x: int, y: int, width: int, height: int):
        self._roi = (x, y, width, height)
        if self._cam is not None:
            try:
                self._cam.set_roi(x, x+width, y, y+height, 
                                  hbin=self._binning, vbin=self._binning)
            except Exception as e:
                logging.warning(f"Error setting ROI: {e}")
                
    def get_roi(self) -> Tuple[int, int, int, int]:
        if self._cam is not None:
            try:
                roi = self._cam.get_roi()
                return (roi[0], roi[2], roi[1]-roi[0], roi[3]-roi[2])
            except:
                pass
        if self._roi is not None:
            return self._roi
        return (0, 0, 1024, 1024)  # iXon Ultra 888 is 1024x1024
    
    def clear_roi(self):
        self._roi = None
        if self._cam is not None:
            try:
                self._cam.set_roi()  # Reset to full
            except:
                pass
                
    def set_binning(self, binning: int):
        self._binning = binning
        if self._cam is not None and self._roi is not None:
            # Re-apply ROI with new binning
            self.set_roi(*self._roi)
                
    def get_binning(self) -> int:
        return self._binning
    
    def set_trigger_mode(self, mode: TriggerMode):
        self._trigger_mode = mode
        if self._cam is not None:
            try:
                if mode == TriggerMode.INTERNAL:
                    self._cam.set_trigger_mode("int")
                elif mode == TriggerMode.EXTERNAL:
                    self._cam.set_trigger_mode("ext")
                elif mode == TriggerMode.SOFTWARE:
                    self._cam.set_trigger_mode("software")
            except Exception as e:
                logging.warning(f"Error setting trigger mode: {e}")
                
    def get_trigger_mode(self) -> TriggerMode:
        return self._trigger_mode
    
    def set_cooling(self, enabled: bool, target_temp_c: float = -70.0):
        if self._cam is not None:
            try:
                if enabled:
                    self._cam.set_temperature(target_temp_c)
                self._target_temp = target_temp_c
            except Exception as e:
                logging.warning(f"Error setting cooling: {e}")
                
    def get_temperature(self) -> Optional[float]:
        if self._cam is not None:
            try:
                return self._cam.get_temperature()
            except:
                pass
        return None
    
    def set_emccd_gain(self, gain: int):
        self._emccd_gain = gain
        if self._cam is not None:
            try:
                self._cam.set_EMCCD_gain(gain)
            except Exception as e:
                logging.warning(f"Error setting EMCCD gain: {e}")
    
    def get_emccd_gain(self) -> int:
        if self._cam is not None:
            try:
                return self._cam.get_EMCCD_gain()
            except:
                pass
        return self._emccd_gain
    
    @property
    def frame_shape(self) -> Tuple[int, int]:
        if self._cam is not None:
            try:
                size = self._cam.get_detector_size()
                return (size[1], size[0])  # (height, width)
            except:
                pass
        roi = self.get_roi()
        h = roi[3] // self._binning
        w = roi[2] // self._binning
        return (h, w)
    
    @property
    def pixel_dtype(self) -> np.dtype:
        return np.uint16
    
    @property
    def is_color(self) -> bool:
        return False
    
    @property
    def camera_info(self) -> Dict[str, Any]:
        info = {
            'model': 'Andor iXon',
            'serial': 'Unknown',
            'firmware': 'Unknown',
            'driver_version': 'Unknown',
            'sensor_size': (1024, 1024),
            'pixel_size': 13.0,  # µm for Ultra 888
            'bit_depth': 16,
            'min_exposure': 0.00001,  # Very short possible with EMCCD
            'max_exposure': 10000.0,
            'temperature': self.get_temperature(),
            'emccd_gain': self.get_emccd_gain()
        }
        
        if self._cam is not None:
            try:
                full_info = self._cam.get_full_info()
                info.update({
                    'model': full_info.get('camera_model', 'Andor iXon'),
                    'serial': full_info.get('serial_number', 'Unknown'),
                })
            except:
                pass
        
        return info


# =============================================================================
# CAMERA FACTORY
# =============================================================================

def create_camera(camera_type: CameraType, **kwargs) -> CameraInterface:
    """
    Factory function to create camera backend.
    
    Args:
        camera_type: Type of camera
        **kwargs: Camera-specific arguments
        
    Returns:
        Camera backend instance
    """
    if camera_type == CameraType.SIMULATED:
        return SimulatedBackend(**kwargs)
    
    elif camera_type == CameraType.HAMAMATSU:
        if not DCAM_AVAILABLE:
            logging.warning("DCAM not available, falling back to simulated")
            return SimulatedBackend()
        return HamamatsuBackend(**kwargs)
    
    elif camera_type == CameraType.THORLABS:
        if not THORLABS_AVAILABLE:
            logging.warning("Thorlabs not available, falling back to simulated")
            return SimulatedBackend()
        return ThorlabsBackend(**kwargs)
    
    elif camera_type == CameraType.UC480:
        if not UC480_AVAILABLE:
            logging.warning("UC480 not available, falling back to simulated")
            return SimulatedBackend()
        return UC480Backend(**kwargs)
    
    elif camera_type == CameraType.ANDOR:
        if not ANDOR_AVAILABLE:
            logging.warning("Andor not available, falling back to simulated")
            return SimulatedBackend()
        return AndorBackend(**kwargs)
    
    else:
        logging.warning(f"Unknown camera type {camera_type}, using simulated")
        return SimulatedBackend()


def list_available_cameras() -> Dict[CameraType, List[str]]:
    """
    List all available cameras by type.
    
    Returns:
        Dict mapping camera type to list of identifiers
    """
    available = {}
    
    # Simulated always available
    available[CameraType.SIMULATED] = ["Simulated Camera"]
    
    # Hamamatsu
    if DCAM_AVAILABLE:
        try:
            dcam.Dcamapi.init()
            n = dcam.Dcamapi.get_devicecount()
            available[CameraType.HAMAMATSU] = [f"Camera {i}" for i in range(n)]
        except:
            pass
    
    # Thorlabs TLCamera
    if THORLABS_AVAILABLE:
        try:
            cameras = Thorlabs.list_cameras_tlcam()
            available[CameraType.THORLABS] = cameras if cameras else []
        except:
            pass
    
    # UC480 / IDS cameras
    if UC480_AVAILABLE:
        try:
            cameras = uc480.list_cameras()
            available[CameraType.UC480] = [f"Camera {i}" for i in range(len(cameras))]
        except:
            pass
    
    # Andor
    if ANDOR_AVAILABLE:
        try:
            n = Andor.get_cameras_number_SDK2()
            available[CameraType.ANDOR] = [f"Camera {i}" for i in range(n)]
        except:
            pass
    
    return available
