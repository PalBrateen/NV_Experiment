"""
Scope Viewer v2 - Core Module
=============================

Configuration classes, device enumeration, and utility functions.
"""

# # Optional nidaqmx
# try:
import nidaqmx, nidaqmx.system, numpy as np, logging
from nidaqmx.constants import AcquisitionType, Edge, TerminalConfiguration, Level
HAS_NIDAQMX = True
# except ImportError:
#     HAS_NIDAQMX = False

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any
from enum import Enum

from PySide6.QtCore import QObject, Signal, QMutex, QMutexLocker


# =============================================================================
# ENUMS FOR TRIGGER/CLOCK CONFIGURATION
# =============================================================================

class TriggerMode(Enum):
    """Trigger mode options"""
    NONE = "none"
    START = "start"
    START_PAUSE = "start_pause"


class TriggerEdge(Enum):
    """Trigger edge options"""
    RISING = "rising"
    FALLING = "falling"


class PauseWhen(Enum):
    """Pause trigger level options"""
    HIGH = "high"
    LOW = "low"


class ClockSource(Enum):
    """Sample clock source options"""
    INTERNAL = "internal"
    EXTERNAL = "external"

# =============================================================================
# DEVICE ENUMERATION
# =============================================================================

@dataclass
class DAQDeviceInfo:
    """Cached information about a DAQ device"""
    name: str
    product_type: str
    ai_channels: List[str]
    ai_voltage_ranges: List[Tuple[float, float]]
    ai_max_single_rate: float
    ai_max_multi_rate: float
    ai_min_rate: float
    ai_simultaneous_sampling: bool
    terminal_configs: List[str]
    pfi_terminals: List[str] = field(default_factory=list)  # Available PFI lines

def check_for_ai_devices(devices) -> list:  #DeviceCollection
    devs = []
    for dev in devices:
        if dev.product_category.name in ['M_SERIES_DAQ', 'X_SERIES_DAQ', 'MIODAQ', 'E_SERIES_DAQ', 'S_SERIES_DAQ']:
            devs.append(dev)
    
    return devs

def enumerate_devices() -> List[str]:
    """Get list of available DAQ devices"""
    if not HAS_NIDAQMX:
        return ["Simulated"]
        
    try:
        system = nidaqmx.system.System.local()
        return [dev.name for dev in check_for_ai_devices(system.devices)]
    except Exception as e:
        logging.exception(f"Failed to enumerate devices: {e}")
        return []


def get_device_info(device_name: str) -> Optional[DAQDeviceInfo]:
    """Get detailed device information"""
    if not HAS_NIDAQMX:
        # Return simulated device info
        return DAQDeviceInfo(
            name="Simulated",
            product_type="Simulated Device",
            ai_channels=[f"ai{i}" for i in range(8)],
            ai_voltage_ranges=[(-10, 10), (-5, 5), (-1, 1)],
            ai_max_single_rate=2e6,
            ai_max_multi_rate=1e6,
            ai_min_rate=1.0,
            ai_simultaneous_sampling=True,
            terminal_configs=["RSE", "NRSE", "DIFF"],
            pfi_terminals=[f"PFI{i}" for i in range(16)]
        )
        
    try:
        system = nidaqmx.system.System.local()
        for dev in check_for_ai_devices(system.devices):
            if dev.product_category.name in ['M_SERIES_DAQ', 'X_SERIES_DAQ', 'MIODAQ', 'E_SERIES_DAQ', 'S_SERIES_DAQ']:
                if dev.name == device_name:
                    # Get voltage ranges
                    ranges = []
                    try:
                        for i in range(0, len(dev.ai_voltage_rngs), 2):
                            ranges.append((dev.ai_voltage_rngs[i], dev.ai_voltage_rngs[i+1]))
                    except:
                        ranges = [(-10., 10.), (-5., 5.), (-1., 1.)]
                        
                    # Get terminal configurations
                    term_configs = []
                    try:
                        for tc in dev.ai_term_cfgs:
                            term_configs.append(tc.name)
                    except:
                        term_configs = ["RSE", "NRSE", "DIFF"]
                    
                    # Get PFI terminals
                    pfi_terminals = []
                    try:
                        # Get terminals from device
                        terminals = dev.terminals
                        for term in terminals:
                            # Extract just the terminal name (e.g., "/Dev1/PFI0" -> "PFI0")
                            term_name = term.split('/')[-1]
                            if term_name.startswith('PFI'):
                                pfi_terminals.append(term_name)
                        pfi_terminals = sorted(pfi_terminals, key=lambda x: int(x[3:]) if x[3:].isdigit() else 0)
                    except:
                        # Default PFI lines for most X-series devices
                        pfi_terminals = [f"PFI{i}" for i in range(16)]
                        
                    return DAQDeviceInfo(
                        name=dev.name,
                        product_type=dev.product_type,
                        ai_channels=[ch.name for ch in dev.ai_physical_chans],
                        ai_voltage_ranges=ranges,
                        ai_max_single_rate=dev.ai_max_single_chan_rate,
                        ai_max_multi_rate=dev.ai_max_multi_chan_rate,
                        ai_min_rate=dev.ai_min_rate,
                        ai_simultaneous_sampling=dev.ai_simultaneous_sampling_supported,
                        terminal_configs=term_configs,
                        pfi_terminals=pfi_terminals
                    )
    except Exception as e:
        logging.exception(f"Failed to get device info: {e}")
    return None


# =============================================================================
# CHANNEL & ACQUISITION CONFIG
# =============================================================================

@dataclass
class ChannelConfig:
    """Configuration for a single analog input channel"""
    physical_channel: str = "ai0"
    enabled: bool = True
    min_voltage: float = -10.0
    max_voltage: float = 10.0
    terminal_config: str = "RSE"
    color: str = "#00BFFF"
    name: str = "Channel 1"


@dataclass
class TriggerConfig:
    """Trigger configuration for acquisition"""
    mode: TriggerMode = TriggerMode.NONE
    
    # Start trigger settings
    start_source: str = "PFI0"
    start_edge: TriggerEdge = TriggerEdge.RISING
    retriggerable: bool = False
    
    # Pause trigger settings
    pause_source: str = "PFI1"
    pause_when: PauseWhen = PauseWhen.HIGH


@dataclass
class ClockConfig:
    """Sample clock configuration"""
    source: ClockSource = ClockSource.INTERNAL
    
    # Internal clock settings
    rate: float = 100000.0
    
    # External clock settings
    external_source: str = "PFI0"
    external_edge: TriggerEdge = TriggerEdge.RISING


@dataclass 
class AcquisitionConfig:
    """Master acquisition configuration"""
    device: str = "Dev1"
    sample_rate: float = 100000.0
    display_length: int = 10000
    channels: List[ChannelConfig] = field(default_factory=list)
    buffer_size: int = 524288
    chunk_size: int = 4096
    
    # Trigger and clock configuration
    trigger: TriggerConfig = field(default_factory=TriggerConfig)
    clock: ClockConfig = field(default_factory=ClockConfig)


# =============================================================================
# OBSERVABLE CONFIG (for IPython sync)
# =============================================================================

class ObservableConfig(QObject):
    """Configuration with signals for property changes"""
    
    sampleRateChanged = Signal(float)
    deviceChanged = Signal(str)
    channelChanged = Signal(int, object)
    acquisitionStateChanged = Signal(bool)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._config = AcquisitionConfig()
        self._is_running = False
        
    @property
    def sample_rate(self) -> float:
        return self._config.sample_rate
        
    @sample_rate.setter
    def sample_rate(self, value: float):
        if value != self._config.sample_rate:
            self._config.sample_rate = value
            self._update_chunk_size()
            self.sampleRateChanged.emit(value)
            
    @property
    def device(self) -> str:
        return self._config.device
        
    @device.setter
    def device(self, value: str):
        if value != self._config.device:
            self._config.device = value
            self.deviceChanged.emit(value)
            
    @property
    def channels(self) -> List[ChannelConfig]:
        return self._config.channels
        
    @property
    def is_running(self) -> bool:
        return self._is_running
        
    @is_running.setter
    def is_running(self, value: bool):
        if value != self._is_running:
            self._is_running = value
            self.acquisitionStateChanged.emit(value)
            
    def _update_chunk_size(self):
        """Auto-calculate chunk size for ~30 Hz updates"""
        chunk = int(self._config.sample_rate / 30)
        self._config.chunk_size = max(256, min(16384, chunk))
        
    def get_config(self) -> AcquisitionConfig:
        return self._config
        
    def add_channel(self, config: ChannelConfig):
        self._config.channels.append(config)
        self.channelChanged.emit(len(self._config.channels) - 1, config)
        
    def update_channel(self, index: int, config: ChannelConfig):
        if 0 <= index < len(self._config.channels):
            self._config.channels[index] = config
            self.channelChanged.emit(index, config)


# =============================================================================
# CIRCULAR BUFFER
# =============================================================================

class CircularBuffer:
    """Thread-safe circular buffer"""
    
    def __init__(self, capacity: int, dtype=np.float32):
        self.capacity = capacity
        self.buffer = np.zeros(capacity, dtype=dtype)
        self.write_idx = 0
        self.count = 0
        self._mutex = QMutex()
        
    def write(self, data: np.ndarray):
        with QMutexLocker(self._mutex):
            n = len(data)
            if n >= self.capacity:
                self.buffer[:] = data[-self.capacity:]
                self.write_idx = 0
                self.count = self.capacity
                return
                
            end = self.write_idx + n
            if end <= self.capacity:
                self.buffer[self.write_idx:end] = data
            else:
                first = self.capacity - self.write_idx
                self.buffer[self.write_idx:] = data[:first]
                self.buffer[:n - first] = data[first:]
                
            self.write_idx = end % self.capacity
            self.count = min(self.count + n, self.capacity)
            
    def read_latest(self, n: int) -> np.ndarray:
        with QMutexLocker(self._mutex):
            n = min(n, self.count)
            if n == 0:
                return np.array([], dtype=self.buffer.dtype)
                
            start = (self.write_idx - n) % self.capacity
            if start + n <= self.capacity:
                return self.buffer[start:start + n].copy()
            else:
                first = self.capacity - start
                return np.concatenate([self.buffer[start:], self.buffer[:n - first]])
                
    def clear(self):
        with QMutexLocker(self._mutex):
            self.write_idx = 0
            self.count = 0
            
    @property
    def available(self) -> int:
        with QMutexLocker(self._mutex):
            return self.count
