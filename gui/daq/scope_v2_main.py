"""
Scope Viewer v2 - Main Application (Refactored)
===============================================

Updates:
- Multi-channel FFT support
- Averaging for scope and FFT
- Spectral density / Power display modes
- Cursors with status bar display
- Frequency unit selection

Author: BP Lab
Date: 2025
"""
# TODO: channel-wise offset inputs needed

import sys, numpy as np, gc, pyqtgraph as pg, ctypes, psutil, json, h5py, os
from typing import Optional, Dict, List
from pathlib import Path
from datetime import datetime

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QSpinBox, QDoubleSpinBox,
    QGroupBox, QComboBox, QTabWidget, QDockWidget, QStatusBar,
    QScrollArea, QToolButton, QSlider, QCheckBox, QFileDialog
)
from PySide6.QtCore import QTimer, Qt, Signal
from PySide6.QtGui import QColor, QPalette, QPainter, QFontMetrics, QIcon

from scope_v2_core import (
    ObservableConfig, ChannelConfig, AcquisitionConfig,
    CircularBuffer, enumerate_devices, get_device_info, DAQDeviceInfo,
    TriggerConfig, ClockConfig, TriggerMode, TriggerEdge, PauseWhen, ClockSource
)
from scope_v2_acquisition import AcquisitionThread, DisplayManager, DisplayConfig
from scope_v2_widgets import (
    ToggleSwitch, SegmentedButton, ScopeWidget, TrendsWidget, 
    FFTWidget, ChannelConfigWidget, FFTDisplayMode, FrequencyUnit
)

try:
    from qtconsole.rich_ipython_widget import RichIPythonWidget
    from qtconsole.inprocess import QtInProcessKernelManager
    HAS_IPYTHON = True
except ImportError:
    HAS_IPYTHON = False

expt_dir = os.path.abspath(r'D:\Brateen\NV_Experiment')     # absolute path to the directory
sys.path.append(expt_dir)   # Add to sys.path

DARK_STYLE = """
QMainWindow,QWidget{background-color:#1e1e1e;color:#d4d4d4;font-family:'Segoe UI',Arial,sans-serif;}
QGroupBox { border: 1px solid #3d3d3d; border-radius: 4px; margin-top: 4px; padding-top: 4px; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; color: #a0a0a0; }
QTabWidget::pane { border: 1px solid #3c3c3c; background-color: #1e1e1e; }
QTabBar::tab { background-color: #2d2d2d; color: #a0a0a0; padding: 4px 10px;
    border: 1px solid #3c3c3c; padding:4px 10px;
    margin-right:2px;border-top-left-radius:3px; border-top-right-radius:3px}
QTabBar::tab:selected { background-color: #1e1e1e; color: #fff; border-bottom-color:#1e1e1e}
QPushButton { background-color: #3d3d3d; border: 1px solid #4d4d4d; border-radius: 4px; padding: 6px 12px; color: #e0e0e0; }
QPushButton:hover {
    background-color: #1aa0d9; color: #ffffff;
}
QPushButton:checked { background-color: #0078d4; }
QLineEdit, QSpinBox, QComboBox, QDoubleSpinBox { background-color: #484848; border: 1px solid #3f3f46; border-radius: 3px;
    padding: 4px; color: #ffffff; selection-background-color: #264f78;
    height: 15px}
QLineEdit:focus, QSpinBox:focus, QComboBox:focus {
    border: 1px solid #009de0;
}
/*QComboBox{padding-right: 20px;}*/
QComboBox::drop-down {
    border: none;
    width: 20px;
}
QComboBox::down-arrow {
    width: 0px; height: 0px;
    border-left: 4px solid #484848;
    border-right: 4px solid #484848;
    border-top: 5px solid #ffffff;
    margin-right: 5px;
}
QComboBox QAbstractItemView {
    background-color: #2d2d30; color: #ffffff;
    selection-background-color: #094771;
    border: 1px solid #3f3f46;
}
/* --- SCROLLBAR STYLING --- */
QComboBox QAbstractItemView QScrollBar:vertical {
    border: none; background: #737475;
    width: 10px; margin: 0px 0px 0px 0px;
}
/* Handle */
QComboBox QAbstractItemView QScrollBar::handle:vertical {
    background: #888; min-height: 20px; border-radius: 5px;
}
/* Handle hover */
QComboBox QAbstractItemView QScrollBar::handle:vertical:hover {
    background: #555;
}
/* Remove arrows */
QComboBox QAbstractItemView QScrollBar::add-line:vertical, 
QComboBox QAbstractItemView QScrollBar::sub-line:vertical {
    height: 0px;
}
/* The Main Box */
QSpinBox, QDoubleSpinBox {
    background-color: #2d2d2d; color: #ffffff;
    border: 1px solid #555555; border-radius: 4px;
    padding-right: 5px; /* Leave space for buttons */
    selection-background-color: #444444;
}
/* The Buttons Container */
QSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button  {
    background-color: #3d3d3d;
    border-left: 1px solid #555555;
    width: 20px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover, QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {
    background-color: #4d4d4d;
}
/* The Arrows (Triangle Hack) */
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {
    width: 0px; height: 0px; border-left: 4px solid #3d3d3d;
    border-right: 4px solid #3d3d3d; border-bottom: 5px solid #ffffff;
}
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {
    width: 0px; height: 0px; border-left: 4px solid #3d3d3d;
    border-right: 4px solid #3d3d3d; border-top: 5px solid #ffffff;
}
QStatusBar{background:#252526}
QDockWidget::title { background-color: #2d2d2d; padding: 6px; }
QScrollArea { border: none; }
QSlider::groove:horizontal { height: 6px; background: #3d3d3d; border-radius: 3px; }
QSlider::handle:horizontal { background: #0078d4; width: 16px; margin: -5px 0; border-radius: 8px; }
"""


class StatusBarLongLabel(QLabel):
    """Status bar label with elided text and tooltip for long messages"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(100)

    def paintEvent(self, event):
        """Custom paint event to draw elided text"""
        painter = QPainter(self)
        metrics = QFontMetrics(self.font())
        elided_text = metrics.elidedText(self.text(), Qt.TextElideMode.ElideRight, self.width())
        painter.drawText(self.rect(), self.alignment(), elided_text)
        painter.end()

    def set_msg(self, message: str, error: bool = False):
        """Update text and set the full message as a tooltip"""
        self.setText(message)
        if error:
            self.setStyleSheet("color: #ff6b6b;")  # Light red for errors
        else:
            self.setStyleSheet("color: #4ae34a;")  # Light green for messages
        self.setToolTip(message)


class ControlTab(QWidget):
    """Main control tab with acquisition and display settings"""
    deviceChanged = Signal(object)  # Emits DAQDeviceInfo when device changes
    sampleRateChanged = Signal(float)  # Emits when user changes sample rate
    chunkSizeChanged = Signal(int)  # Emits when user changes chunk size
    channelRemoved = Signal(int)  # Emits channel index when removed
    inputTypeChanged = Signal(str)  # Emits 'analog' or 'counter'
    
    def __init__(self, obs_config: ObservableConfig, parent=None):
        super().__init__(parent)
        self.obs_config = obs_config
        self.device_info: Optional[DAQDeviceInfo] = None
        self.channel_widgets: List[ChannelConfigWidget] = []
        self._setup_ui()
        self._connect_signals()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        dev_layout = QHBoxLayout()
        dev_layout.addWidget(QLabel("Device"))
        self.device_combo = QComboBox()
        self.device_combo.currentTextChanged.connect(self._on_device_changed)
        dev_layout.addWidget(self.device_combo, stretch=1)
        
        self.refresh_btn = QPushButton()
        self.refresh_btn.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; padding: 2px 2px;}
            QPushButton:hover {background-color: #009de0; color: #ffffff;}
        """)
        self.refresh_btn.setFixedSize(25, 28)
        self.refresh_btn.setText("↻")
        self.refresh_btn.clicked.connect(self._refresh_devices)
        dev_layout.addWidget(self.refresh_btn)
        layout.addLayout(dev_layout)
        
        acq_group = QGroupBox("Acquisition")
        acq_layout = QGridLayout(acq_group)
        
        # Input type selector (Analog/Counter)
        acq_layout.addWidget(QLabel("Input Type"), 0, 0)
        self.input_type_combo = QComboBox()
        self.input_type_combo.addItems(["Analog (Voltage)", "Counter (Edges)"])
        self.input_type_combo.currentIndexChanged.connect(self._on_input_type_changed)
        acq_layout.addWidget(self.input_type_combo, 0, 1)
        
        acq_layout.addWidget(QLabel("Sample Rate"), 1, 0)
        self.rate_combo = QComboBox()
        self.rate_combo.setEditable(True)
        self.rate_combo.addItems(["2M", "1M", "500k", "200k", "100k", "50k", "20k", "10k", "1k"])
        self.rate_combo.setCurrentText("100k")
        self.rate_combo.currentTextChanged.connect(self._on_rate_changed)
        acq_layout.addWidget(self.rate_combo, 1, 1)
        acq_layout.addWidget(QLabel("Chunk Size"), 2, 0)
        self.chunk_size_combo = QComboBox()
        self.chunk_size_combo.addItems(["1024", "2048", "4096", "8192", "16384", "32768", "65536"])
        self.chunk_size_combo.setCurrentText("8192")
        self.chunk_size_combo.currentTextChanged.connect(self._on_chunk_size_changed)
        acq_layout.addWidget(self.chunk_size_combo, 2, 1)
        acq_layout.addWidget(QLabel("Buffer Length"), 3, 0)
        self.buffer_combo = QComboBox()
        self.buffer_combo.addItems(["0.5 s", "1 s", "2 s", "5 s", "10 s", "30 s"])
        self.buffer_combo.setCurrentText("2 s")
        self.buffer_combo.currentTextChanged.connect(self._on_buffer_changed)
        acq_layout.addWidget(self.buffer_combo, 3, 1)
        self.buffer_info = QLabel("Buffer: -- samples (-- MB)")
        self.buffer_info.setStyleSheet("color: #888; font-size: 10px;")
        acq_layout.addWidget(self.buffer_info, 4, 0, 1, 2)
        layout.addWidget(acq_group)
        
        chan_group = QGroupBox("Channels")
        self.chan_layout = QVBoxLayout(chan_group)
        chan_btn_layout = QHBoxLayout()
        self.add_chan_btn = QPushButton("+ Add")
        self.add_chan_btn.clicked.connect(self._add_channel)
        chan_btn_layout.addWidget(self.add_chan_btn)
        self.remove_chan_btn = QPushButton("- Remove")
        self.remove_chan_btn.clicked.connect(self._remove_channel)
        chan_btn_layout.addWidget(self.remove_chan_btn)
        # chan_btn_layout.addStretch()
        self.chan_layout.addLayout(chan_btn_layout)

        self.chan_scroll = QScrollArea()
        self.chan_scroll.setWidgetResizable(True)
        self.chan_scroll.setMinimumHeight(100)
        self.chan_scroll.setMaximumHeight(250)
        self.chan_container = QWidget()
        self.chan_container_layout = QVBoxLayout(self.chan_container)
        self.chan_container_layout.setContentsMargins(0, 0, 0, 0)
        self.chan_scroll.setWidget(self.chan_container)
        self.chan_layout.addWidget(self.chan_scroll)
        layout.addWidget(chan_group)
        
        layout.addStretch()
        display_group = QGroupBox("Display")
        display_layout = QGridLayout(display_group)
        display_layout.addWidget(QLabel("Scope"), 0, 0)
        self.scope_toggle = ToggleSwitch()
        self.scope_toggle.setChecked(True)
        display_layout.addWidget(self.scope_toggle, 0, 1)
        display_layout.addWidget(QLabel("Trends"), 0, 2)
        self.trends_toggle = ToggleSwitch()
        display_layout.addWidget(self.trends_toggle, 0, 3)
        display_layout.addWidget(QLabel("History:"), 1, 0)
        self.trends_length = QComboBox()
        self.trends_length.addItems(["10 s", "30 s", "60 s"])
        self.trends_length.setCurrentText("60 s")
        display_layout.addWidget(self.trends_length, 1, 1)
        self.clear_trends_btn = QPushButton("Clear")
        self.clear_trends_btn.setStyleSheet("background-color: #5a3d2d;")
        display_layout.addWidget(self.clear_trends_btn, 1, 2, 1, 2)
        layout.addWidget(display_group)
        
        
        self._refresh_devices()
        self._add_channel()
        self._update_buffer_info()
        
    def _connect_signals(self):
        self.obs_config.sampleRateChanged.connect(self._on_external_rate)
        self.obs_config.acquisitionStateChanged.connect(self._on_external_state)
        
    def _on_external_rate(self, rate: float):
        s = self._rate_to_str(rate)
        if self.rate_combo.currentText() != s:
            self.rate_combo.blockSignals(True)
            self.rate_combo.setCurrentText(s)
            self.rate_combo.blockSignals(False)
        self._update_buffer_info()
            
    def _on_external_state(self, running: bool):
        """Update UI state based on running status - controls are now disabled/enabled"""
        # Disable device and rate changes while running
        self.device_combo.setEnabled(not running)
        self.rate_combo.setEnabled(not running)
        self.add_chan_btn.setEnabled(not running)
        self.remove_chan_btn.setEnabled(not running)
        
    def _rate_to_str(self, rate: float) -> str:
        if rate >= 1e6: return f"{rate/1e6:.0f}M"
        elif rate >= 1e3: return f"{rate/1e3:.0f}k"
        return f"{rate:.0f}"
            
    def _refresh_devices(self):
        self.device_combo.clear()
        devices = enumerate_devices()
        self.device_combo.addItems(devices if devices else ["No devices"])
        if devices: self._on_device_changed(devices[0])
            
    def _on_device_changed(self, name: str):
        self.obs_config.device = name
        self.device_info = get_device_info(name)
        for w in self.channel_widgets: w.update_device_info(self.device_info)
        # Emit signal with device info for other components (like trigger tab)
        if self.device_info:
            self.deviceChanged.emit(self.device_info)
                
    def _on_rate_changed(self, text: str):
        rate = self._parse_rate(text)
        if rate:
            self.obs_config.sample_rate = rate
            self._update_buffer_info()
            self.sampleRateChanged.emit(rate)  # Notify other components
    
    def set_sample_rate(self, rate: float):
        """Set sample rate from external source (e.g., TriggerClockTab)"""
        self.rate_combo.blockSignals(True)
        self.rate_combo.setCurrentText(self._rate_to_str(rate))
        self.rate_combo.blockSignals(False)
        self._update_buffer_info()
    
    def _on_chunk_size_changed(self, text: str):
        """Handle chunk size change from UI"""
        try:
            size = int(text)
            self.chunkSizeChanged.emit(size)
        except:
            pass
    
    def set_chunk_size(self, size: int):
        """Set chunk size from external source (sync with FFT tab)"""
        self.chunk_size_combo.blockSignals(True)
        self.chunk_size_combo.setCurrentText(str(size))
        self.chunk_size_combo.blockSignals(False)
    
    def get_chunk_size(self) -> int:
        """Get current chunk size"""
        try:
            return int(self.chunk_size_combo.currentText())
        except:
            return 8192
    
    def _on_input_type_changed(self, index: int):
        """Handle input type change (Analog/Counter)"""
        input_type = 'analog' if index == 0 else 'counter'
        self.inputTypeChanged.emit(input_type)
    
    def get_input_type(self) -> str:
        """Get current input type: 'analog' or 'counter'"""
        return 'analog' if self.input_type_combo.currentIndex() == 0 else 'counter'
            
    def _on_buffer_changed(self, text: str):
        self._update_buffer_info()
            
    def _parse_rate(self, text: str) -> Optional[float]:
        text = text.strip().upper()
        mult = {'K': 1e3, 'M': 1e6, 'G': 1e9}
        try:
            for s, m in mult.items():
                if text.endswith(s): return float(text[:-1]) * m
            return float(text)
        except ValueError: return None
            
    def _update_buffer_info(self):
        rate = self.obs_config.sample_rate
        duration = self.get_buffer_duration()
        samples = int(rate * duration)
        mb = samples * 4 / 1024 / 1024
        self.buffer_info.setText(f"Buffer: {samples:,} samples ({mb:.1f} MB/ch)")
        
    def get_buffer_duration(self) -> float:
        text = self.buffer_combo.currentText()
        try: return float(text.replace('s', '').strip())
        except: return 2.0
            
    def _add_channel(self):
        idx = len(self.channel_widgets)
        w = ChannelConfigWidget(idx, self.device_info)
        w.configChanged.connect(self._on_channel_changed)
        self.channel_widgets.append(w)
        self.chan_container_layout.addWidget(w)
        self.obs_config.add_channel(w.get_config())
        
    def _remove_channel(self):
        if len(self.channel_widgets) > 1:
            w = self.channel_widgets.pop()
            removed_idx = len(self.channel_widgets)  # Index of removed channel
            w.deleteLater()
            self.obs_config._config.channels.pop()
            self.channelRemoved.emit(removed_idx)  # Notify to remove curves
            
    def _on_channel_changed(self, idx: int, config: ChannelConfig):
        self.obs_config.update_channel(idx, config)
        
    def get_all_configs(self) -> List[ChannelConfig]:
        return [w.get_config() for w in self.channel_widgets]
        
    def set_running(self, running: bool):
        """Update UI state based on running status"""
        self._on_external_state(running)

class ScopeTab(QWidget):
    """Scope display settings with averaging and cursors"""
    chunkSizeChanged = Signal(int)  # Synced with FFT Size
    pointsChanged = Signal(int)
    histogramChanged = Signal(bool)  # Histogram toggle
    averagingChanged = Signal(bool, int)
    cursorsChanged = Signal(bool)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Acquisition group - chunk size control
        acq_group = QGroupBox("Acquisition")
        acq_layout = QGridLayout(acq_group)
        
        acq_layout.addWidget(QLabel("Chunk Size"), 0, 0)
        self.chunk_size_combo = QComboBox()
        self.chunk_size_combo.addItems(["1024", "2048", "4096", "8192", "16384", "32768", "65536"])
        self.chunk_size_combo.setCurrentText("8192")
        self.chunk_size_combo.currentTextChanged.connect(self._on_chunk_size)
        acq_layout.addWidget(self.chunk_size_combo, 0, 1)
        
        self.time_info_label = QLabel("Duration: -- ms")
        self.time_info_label.setStyleSheet("color: #4ec9b0;")
        acq_layout.addWidget(self.time_info_label, 1, 0, 1, 2)
        
        layout.addWidget(acq_group)
        
        # Display settings
        disp_group = QGroupBox("Display")
        disp_layout = QGridLayout(disp_group)
        disp_layout.addWidget(QLabel("Max Points"), 0, 0)
        self.points_combo = QComboBox()
        self.points_combo.addItems(["1000", "2000", "5000", "10000", "20000", "50000"])
        self.points_combo.setCurrentText("5000")
        self.points_combo.currentTextChanged.connect(self._on_points)
        disp_layout.addWidget(self.points_combo, 0, 1)
        self.decimation_label = QLabel("Decimation: 1x")
        self.decimation_label.setStyleSheet("color: #888; font-size: 10px;")
        disp_layout.addWidget(self.decimation_label, 1, 0, 1, 2)
        
        # Histogram toggle
        disp_layout.addWidget(QLabel("Histogram"), 2, 0)
        self.histogram_toggle = ToggleSwitch()
        self.histogram_toggle.toggled.connect(lambda h: self.histogramChanged.emit(h))
        disp_layout.addWidget(self.histogram_toggle, 2, 1)
        
        layout.addWidget(disp_group)
        
        avg_group = QGroupBox("Averaging")
        avg_layout = QGridLayout(avg_group)
        avg_layout.addWidget(QLabel("Enable"), 0, 0)
        self.avg_toggle = ToggleSwitch()
        self.avg_toggle.toggled.connect(self._on_averaging_changed)
        avg_layout.addWidget(self.avg_toggle, 0, 1)
        avg_layout.addWidget(QLabel("Count"), 1, 0)
        self.avg_count = QComboBox()
        self.avg_count.addItems(["2", "4", "8", "16", "32", "64", "128"])
        self.avg_count.setCurrentText("8")
        self.avg_count.currentTextChanged.connect(self._on_averaging_changed)
        avg_layout.addWidget(self.avg_count, 1, 1)
        self.clear_avg_btn = QPushButton("Clear Avg")
        self.clear_avg_btn.clicked.connect(lambda: self.averagingChanged.emit(self.avg_toggle.isChecked(), 0))
        avg_layout.addWidget(self.clear_avg_btn, 2, 0, 1, 2)
        layout.addWidget(avg_group)
        
        cursor_group = QGroupBox("Cursors")
        cursor_layout = QGridLayout(cursor_group)
        cursor_layout.addWidget(QLabel("Show"), 0, 0)
        self.cursor_toggle = ToggleSwitch()
        self.cursor_toggle.toggled.connect(lambda c: self.cursorsChanged.emit(c))
        cursor_layout.addWidget(self.cursor_toggle, 0, 1)
        layout.addWidget(cursor_group)
        layout.addStretch()
    
    def _on_chunk_size(self, text: str):
        try: 
            self.chunkSizeChanged.emit(int(text))
        except: 
            pass
        
    def _on_points(self, text: str):
        try: 
            self.pointsChanged.emit(int(text))
        except: 
            pass
            
    def _on_averaging_changed(self):
        enabled = self.avg_toggle.isChecked()
        count = int(self.avg_count.currentText()) if enabled else 1
        self.averagingChanged.emit(enabled, count)
            
    def get_points(self) -> int:
        try: 
            return int(self.points_combo.currentText())
        except: 
            return 5000
            
    def update_chunk_info(self, sample_rate: float, chunk_size: int):
        """Update info labels with current chunk parameters"""
        duration_ms = (chunk_size / sample_rate) * 1000
        
        if duration_ms >= 1000:
            self.time_info_label.setText(f"Duration: {duration_ms/1000:.2f} s")
        else:
            self.time_info_label.setText(f"Duration: {duration_ms:.2f} ms")
        
        # Update decimation info
        points = self.get_points()
        if chunk_size > points:
            dec = chunk_size // points
            self.decimation_label.setText(f"Decimation: {dec}x ({chunk_size} → {points})")
        else:
            self.decimation_label.setText(f"No decimation ({chunk_size} samples)")
    
    def set_chunk_size(self, size: int):
        """Set chunk size from external source (sync with FFT tab)"""
        self.chunk_size_combo.blockSignals(True)
        self.chunk_size_combo.setCurrentText(str(size))
        self.chunk_size_combo.blockSignals(False)
    
    def get_chunk_size(self) -> int:
        try:
            return int(self.chunk_size_combo.currentText())
        except:
            return 8192


class FFTTab(QWidget):
    """FFT settings tab with spectral modes and frequency units"""
    modeChanged = Signal(int)
    xLogChanged = Signal(bool)
    yLogChanged = Signal(bool)
    useDbChanged = Signal(bool)
    displayModeChanged = Signal(object)
    freqUnitChanged = Signal(object)
    fftSizeChanged = Signal(int)
    windowChanged = Signal(str)
    averagingChanged = Signal(bool, int)
    cursorsChanged = Signal(bool)
    correctionChanged = Signal(bool)  # True = Power, False = Amplitude
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        self.fft_group = QGroupBox("FFT Settings")
        fft_layout = QGridLayout(self.fft_group)
        fft_layout.addWidget(QLabel("Window"), 0, 0)
        self.window_combo = QComboBox()
        self.window_combo.addItems(["Hann", "Hamming", "Blackman", "Blackman-Harris", "Flat Top", "Rectangular"])
        self.window_combo.currentTextChanged.connect(lambda w: self.windowChanged.emit(w))
        fft_layout.addWidget(self.window_combo, 0, 1)
        fft_layout.addWidget(QLabel("FFT Size"), 1, 0)
        self.fft_size_combo = QComboBox()
        self.fft_size_combo.addItems(["1024", "2048", "4096", "8192", "16384", "32768", "65536"])
        self.fft_size_combo.setCurrentText("8192")
        self.fft_size_combo.currentTextChanged.connect(self._on_fft_size)
        fft_layout.addWidget(self.fft_size_combo, 1, 1)
        fft_layout.addWidget(QLabel("Correction"), 2, 0)
        self.correction_combo = QComboBox()
        self.correction_combo.addItems(["Amplitude", "Power"])
        self.correction_combo.currentTextChanged.connect(self._on_correction_changed)
        fft_layout.addWidget(self.correction_combo, 2, 1)
        self.resolution_label = QLabel("Resolution: -- Hz")
        self.resolution_label.setStyleSheet("color: #4ec9b0; font-size: 11px;")
        fft_layout.addWidget(self.resolution_label, 3, 0, 1, 2)
        layout.addWidget(self.fft_group)
        
        h_group = QGroupBox("Horizontal")
        h_layout = QGridLayout(h_group)
        h_layout.addWidget(QLabel("Units"), 0, 0)
        self.freq_unit_combo = QComboBox()
        self.freq_unit_combo.addItems(["Hz", "kHz", "MHz"])
        self.freq_unit_combo.currentTextChanged.connect(self._on_freq_unit)
        h_layout.addWidget(self.freq_unit_combo, 0, 1)
        h_layout.addWidget(QLabel("Scale"), 1, 0)
        self.x_scale = SegmentedButton(["Lin", "Log"])
        self.x_scale.valueChanged.connect(lambda i: self.xLogChanged.emit(i == 1))
        h_layout.addWidget(self.x_scale, 1, 1)
        layout.addWidget(h_group)
        
        v_group = QGroupBox("Vertical")
        v_layout = QGridLayout(v_group)
        v_layout.addWidget(QLabel("Spectral Density"), 0, 0)
        self.spectral_density_toggle = ToggleSwitch()
        self.spectral_density_toggle.toggled.connect(self._on_display_mode_changed)
        v_layout.addWidget(self.spectral_density_toggle, 0, 1)
        v_layout.addWidget(QLabel("Power"), 1, 0)
        self.power_toggle = ToggleSwitch()
        self.power_toggle.toggled.connect(self._on_display_mode_changed)
        v_layout.addWidget(self.power_toggle, 1, 1)
        v_layout.addWidget(QLabel("dB"), 2, 0)
        self.db_toggle = ToggleSwitch()
        self.db_toggle.setChecked(True)
        self.db_toggle.toggled.connect(lambda c: self.useDbChanged.emit(c))
        v_layout.addWidget(self.db_toggle, 2, 1)
        v_layout.addWidget(QLabel("Log Scale"), 3, 0)
        self.y_log_toggle = ToggleSwitch()
        self.y_log_toggle.toggled.connect(lambda c: self.yLogChanged.emit(c))
        v_layout.addWidget(self.y_log_toggle, 3, 1)
        layout.addWidget(v_group)
        
        avg_group = QGroupBox("Averaging")
        avg_layout = QGridLayout(avg_group)
        avg_layout.addWidget(QLabel("Enable"), 0, 0)
        self.avg_toggle = ToggleSwitch()
        self.avg_toggle.toggled.connect(self._on_averaging_changed)
        avg_layout.addWidget(self.avg_toggle, 0, 1)
        avg_layout.addWidget(QLabel("Count"), 1, 0)
        self.avg_count = QComboBox()
        self.avg_count.addItems(["2", "4", "8", "16", "32", "64", "128"])
        self.avg_count.setCurrentText("8")
        self.avg_count.currentTextChanged.connect(self._on_averaging_changed)
        avg_layout.addWidget(self.avg_count, 1, 1)
        layout.addWidget(avg_group)
        
        opt_group = QGroupBox("Options")
        opt_layout = QGridLayout(opt_group)
        opt_layout.addWidget(QLabel("Peaks"), 0, 0)
        self.peaks_toggle = ToggleSwitch()
        self.peaks_toggle.setChecked(True)
        opt_layout.addWidget(self.peaks_toggle, 0, 1)
        opt_layout.addWidget(QLabel("Cursors"), 1, 0)
        self.cursor_toggle = ToggleSwitch()
        self.cursor_toggle.toggled.connect(lambda c: self.cursorsChanged.emit(c))
        opt_layout.addWidget(self.cursor_toggle, 1, 1)
        layout.addWidget(opt_group)
        
        self.fft_group.setEnabled(False)
        h_group.setEnabled(False)
        v_group.setEnabled(False)
        avg_group.setEnabled(False)
        opt_group.setEnabled(False)
        self._fft_controls = [self.fft_group, h_group, v_group, avg_group, opt_group]
        layout.addStretch()
    
    def set_fft_mode_enabled(self, enabled: bool):
        """Enable/disable FFT controls based on display mode"""
        for ctrl in self._fft_controls: 
            ctrl.setEnabled(enabled)
        
    def _on_fft_size(self, text: str):
        try: self.fftSizeChanged.emit(int(text))
        except: pass
    
    def _on_correction_changed(self, text: str):
        """Emit correction mode change (Amplitude=False, Power=True)"""
        self.correctionChanged.emit(text == "Power")
            
    def _on_freq_unit(self, text: str):
        unit_map = {"Hz": FrequencyUnit.HZ, "kHz": FrequencyUnit.KHZ, "MHz": FrequencyUnit.MHZ}
        if text in unit_map: self.freqUnitChanged.emit(unit_map[text])
            
    def _on_display_mode_changed(self):
        spectral = self.spectral_density_toggle.isChecked()
        power = self.power_toggle.isChecked()
        if spectral and power: mode = FFTDisplayMode.PSD
        elif spectral: mode = FFTDisplayMode.SPECTRAL_DENSITY
        elif power: mode = FFTDisplayMode.POWER
        else: mode = FFTDisplayMode.NORMAL
        self.displayModeChanged.emit(mode)
        
    def _on_averaging_changed(self):
        enabled = self.avg_toggle.isChecked()
        count = int(self.avg_count.currentText()) if enabled else 1
        self.averagingChanged.emit(enabled, count)
        
    def get_fft_size(self) -> int: return int(self.fft_size_combo.currentText())
    def get_window(self) -> str: return self.window_combo.currentText()
    def is_power_correction(self) -> bool: return self.correction_combo.currentText() == "Power"
        
    def update_resolution(self, sample_rate: float):
        fft_size = self.get_fft_size()
        resolution = sample_rate / fft_size
        self.resolution_label.setText(f"Resolution: {resolution:.2f} Hz")
    
    def set_fft_size(self, size: int):
        """Set FFT size from external source (sync with Scope tab)"""
        self.fft_size_combo.blockSignals(True)
        self.fft_size_combo.setCurrentText(str(size))
        self.fft_size_combo.blockSignals(False)


class TriggerClockTab(QWidget):
    """Trigger and sample clock configuration tab"""
    
    # Signals to notify config changes
    triggerConfigChanged = Signal(object)  # TriggerConfig
    clockConfigChanged = Signal(object)    # ClockConfig
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._pfi_terminals = [f"PFI{i}" for i in range(16)]  # Default
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        
        # === TRIGGER SECTION ===
        # trigger_group = QGroupBox("Trigger")
        trigger_layout = QGridLayout()
        
        # Trigger Mode
        trigger_layout.addWidget(QLabel("Mode"), 0, 0)
        self.trigger_mode = QComboBox()
        self.trigger_mode.addItems(["None", "Start", "Start + Pause"])
        self.trigger_mode.currentIndexChanged.connect(self._on_trigger_mode_changed)
        trigger_layout.addWidget(self.trigger_mode, 0, 1)
        
        # Start Trigger Settings
        self.start_trigger_group = QGroupBox("Start Trigger")
        start_layout = QGridLayout(self.start_trigger_group)
        
        start_layout.addWidget(QLabel("Source"), 0, 0)
        self.start_source = QComboBox()
        self.start_source.addItems(self._pfi_terminals)
        self.start_source.currentTextChanged.connect(self._emit_trigger_config)
        start_layout.addWidget(self.start_source, 0, 1)
        
        start_layout.addWidget(QLabel("Edge"), 1, 0)
        self.start_edge = QComboBox()
        self.start_edge.addItems(["Rising", "Falling"])
        self.start_edge.currentTextChanged.connect(self._emit_trigger_config)
        start_layout.addWidget(self.start_edge, 1, 1)
        
        start_layout.addWidget(QLabel("Retriggerable"), 2, 0)
        self.retriggerable = ToggleSwitch()
        self.retriggerable.toggled.connect(self._emit_trigger_config)
        start_layout.addWidget(self.retriggerable, 2, 1)
        
        trigger_layout.addWidget(self.start_trigger_group, 1, 0, 1, 2)
        
        # Pause Trigger Settings
        self.pause_trigger_group = QGroupBox("Pause Trigger")
        pause_layout = QGridLayout(self.pause_trigger_group)
        
        pause_layout.addWidget(QLabel("Source"), 0, 0)
        self.pause_source = QComboBox()
        self.pause_source.addItems(self._pfi_terminals)
        self.pause_source.currentTextChanged.connect(self._emit_trigger_config)
        pause_layout.addWidget(self.pause_source, 0, 1)
        
        pause_layout.addWidget(QLabel("Pause When"), 1, 0)
        self.pause_when = QComboBox()
        self.pause_when.addItems(["High", "Low"])
        self.pause_when.currentTextChanged.connect(self._emit_trigger_config)
        pause_layout.addWidget(self.pause_when, 1, 1)
        
        trigger_layout.addWidget(self.pause_trigger_group, 2, 0, 1, 2)
        
        layout.addLayout(trigger_layout)
        
        # === SAMPLE CLOCK SECTION ===
        clock_group = QGroupBox("Sample Clock")
        clock_layout = QGridLayout(clock_group)
        
        # Clock Source
        clock_layout.addWidget(QLabel("Source"), 0, 0)
        self.clock_source = QComboBox()
        self.clock_source.addItems(["Internal", "External"])
        self.clock_source.currentIndexChanged.connect(self._on_clock_source_changed)
        clock_layout.addWidget(self.clock_source, 0, 1)
        
        # Internal Clock Settings
        self.internal_clock_group = QGroupBox("Internal Clock")
        internal_layout = QGridLayout(self.internal_clock_group)
        
        internal_layout.addWidget(QLabel("Rate"), 0, 0)
        self.clock_rate = QComboBox()
        self.clock_rate.setEditable(True)
        self.clock_rate.addItems(["10000", "50000", "100000", "250000", "500000", "1000000", "2000000"])
        self.clock_rate.setCurrentText("100000")
        self.clock_rate.currentTextChanged.connect(self._emit_clock_config)
        internal_layout.addWidget(self.clock_rate, 0, 1)
        internal_layout.addWidget(QLabel("Sa/s"), 0, 2)
        
        clock_layout.addWidget(self.internal_clock_group, 1, 0, 1, 2)
        
        # External Clock Settings
        self.external_clock_group = QGroupBox("External Clock")
        external_layout = QGridLayout(self.external_clock_group)
        
        external_layout.addWidget(QLabel("Source"), 0, 0)
        self.ext_clock_source = QComboBox()
        self.ext_clock_source.addItems(self._pfi_terminals)
        self.ext_clock_source.currentTextChanged.connect(self._emit_clock_config)
        external_layout.addWidget(self.ext_clock_source, 0, 1)
        
        external_layout.addWidget(QLabel("Edge"), 1, 0)
        self.ext_clock_edge = QComboBox()
        self.ext_clock_edge.addItems(["Rising", "Falling"])
        self.ext_clock_edge.currentTextChanged.connect(self._emit_clock_config)
        external_layout.addWidget(self.ext_clock_edge, 1, 1)
        
        external_layout.addWidget(QLabel("Rate"), 2, 0)
        self.ext_clock_rate = QComboBox()
        self.ext_clock_rate.setEditable(True)
        self.ext_clock_rate.addItems(["1000", "10000", "50000", "100000", "500000", "1000000", "2000000"])
        self.ext_clock_rate.setCurrentText("100000")
        self.ext_clock_rate.currentTextChanged.connect(self._emit_clock_config)
        external_layout.addWidget(self.ext_clock_rate, 2, 1)
        external_layout.addWidget(QLabel("Sa/s"), 2, 2)
        
        # Note about external clock rate
        ext_note = QLabel("Note: Enter the expected external clock frequency. Used for FFT resolution and time axis scaling.")
        ext_note.setStyleSheet("color: #888; font-size: 10px;")
        ext_note.setWordWrap(True)
        external_layout.addWidget(ext_note, 3, 0, 1, 3)
        
        clock_layout.addWidget(self.external_clock_group, 2, 0, 1, 2)
        
        layout.addWidget(clock_group)
        
        # Info label
        self.info_label = QLabel("")
        self.info_label.setStyleSheet("color: #4ec9b0; font-size: 11px;")
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)
        
        layout.addStretch()
        
        # Initialize visibility
        self._on_trigger_mode_changed(0)
        self._on_clock_source_changed(0)
        
    def update_pfi_terminals(self, terminals: List[str]):
        """Update available PFI terminals from device info"""
        self._pfi_terminals = terminals if terminals else [f"PFI{i}" for i in range(16)]
        
        # Update all PFI combo boxes
        for combo in [self.start_source, self.pause_source, self.ext_clock_source]:
            current = combo.currentText()
            combo.clear()
            combo.addItems(self._pfi_terminals)
            # Try to restore previous selection
            idx = combo.findText(current)
            if idx >= 0:
                combo.setCurrentIndex(idx)
    
    def _on_trigger_mode_changed(self, idx: int):
        """Handle trigger mode change"""
        mode = ["None", "Start", "Start + Pause"][idx]
        
        self.start_trigger_group.setEnabled(idx >= 1)
        self.pause_trigger_group.setEnabled(idx >= 2)
        
        # Update info
        if idx == 0:
            self.info_label.setText("Free-running acquisition, no trigger required.")
        elif idx == 1:
            self.info_label.setText("Acquisition starts on trigger edge, runs continuously.")
        else:
            self.info_label.setText("Acquisition starts on trigger, pauses when pause condition met.")
        
        self._emit_trigger_config()
    
    def _on_clock_source_changed(self, idx: int):
        """Handle clock source change"""
        self.internal_clock_group.setEnabled(idx == 0)
        self.external_clock_group.setEnabled(idx == 1)
        self._emit_clock_config()
    
    def _emit_trigger_config(self):
        """Emit current trigger configuration"""
        mode_map = {0: TriggerMode.NONE, 1: TriggerMode.START, 2: TriggerMode.START_PAUSE}
        edge_map = {"Rising": TriggerEdge.RISING, "Falling": TriggerEdge.FALLING}
        when_map = {"High": PauseWhen.HIGH, "Low": PauseWhen.LOW}
        
        config = TriggerConfig(
            mode=mode_map.get(self.trigger_mode.currentIndex(), TriggerMode.NONE),
            start_source=self.start_source.currentText(),
            start_edge=edge_map.get(self.start_edge.currentText(), TriggerEdge.RISING),
            retriggerable=self.retriggerable.isChecked(),
            pause_source=self.pause_source.currentText(),
            pause_when=when_map.get(self.pause_when.currentText(), PauseWhen.HIGH)
        )
        self.triggerConfigChanged.emit(config)
    
    def _emit_clock_config(self):
        """Emit current clock configuration"""
        source_map = {0: ClockSource.INTERNAL, 1: ClockSource.EXTERNAL}
        edge_map = {"Rising": TriggerEdge.RISING, "Falling": TriggerEdge.FALLING}
        
        try:
            internal_rate = float(self.clock_rate.currentText())
        except:
            internal_rate = 100000.0
        
        try:
            external_rate = float(self.ext_clock_rate.currentText())
        except:
            external_rate = 100000.0
        
        # Use the appropriate rate based on clock source
        is_external = self.clock_source.currentIndex() == 1
        
        config = ClockConfig(
            source=source_map.get(self.clock_source.currentIndex(), ClockSource.INTERNAL),
            rate=external_rate if is_external else internal_rate,
            external_source=self.ext_clock_source.currentText(),
            external_edge=edge_map.get(self.ext_clock_edge.currentText(), TriggerEdge.RISING)
        )
        self.clockConfigChanged.emit(config)
    
    def get_trigger_config(self) -> TriggerConfig:
        """Get current trigger configuration"""
        mode_map = {0: TriggerMode.NONE, 1: TriggerMode.START, 2: TriggerMode.START_PAUSE}
        edge_map = {"Rising": TriggerEdge.RISING, "Falling": TriggerEdge.FALLING}
        when_map = {"High": PauseWhen.HIGH, "Low": PauseWhen.LOW}
        
        return TriggerConfig(
            mode=mode_map.get(self.trigger_mode.currentIndex(), TriggerMode.NONE),
            start_source=self.start_source.currentText(),
            start_edge=edge_map.get(self.start_edge.currentText(), TriggerEdge.RISING),
            retriggerable=self.retriggerable.isChecked(),
            pause_source=self.pause_source.currentText(),
            pause_when=when_map.get(self.pause_when.currentText(), PauseWhen.HIGH)
        )
    
    def get_clock_config(self) -> ClockConfig:
        """Get current clock configuration"""
        source_map = {0: ClockSource.INTERNAL, 1: ClockSource.EXTERNAL}
        edge_map = {"Rising": TriggerEdge.RISING, "Falling": TriggerEdge.FALLING}
        
        try:
            internal_rate = float(self.clock_rate.currentText())
        except:
            internal_rate = 100000.0
        
        try:
            external_rate = float(self.ext_clock_rate.currentText())
        except:
            external_rate = 100000.0
        
        # Use the appropriate rate based on clock source
        is_external = self.clock_source.currentIndex() == 1
        
        return ClockConfig(
            source=source_map.get(self.clock_source.currentIndex(), ClockSource.INTERNAL),
            rate=external_rate if is_external else internal_rate,
            external_source=self.ext_clock_source.currentText(),
            external_edge=edge_map.get(self.ext_clock_edge.currentText(), TriggerEdge.RISING)
        )
    
    def set_sample_rate(self, rate: float):
        """Set sample rate from external source (e.g., ControlTab)"""
        rate_str = str(int(rate)) if rate == int(rate) else str(rate)
        
        # Update the appropriate field based on current clock source
        if self.clock_source.currentIndex() == 0:  # Internal
            self.clock_rate.blockSignals(True)
            self.clock_rate.setCurrentText(rate_str)
            self.clock_rate.blockSignals(False)
        else:  # External
            self.ext_clock_rate.blockSignals(True)
            self.ext_clock_rate.setCurrentText(rate_str)
            self.ext_clock_rate.blockSignals(False)


class IPythonWidget(QWidget):
    """IPython console widget"""
    def __init__(self, namespace: dict, parent=None):
        super().__init__(parent)
        self.namespace = namespace
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        if not HAS_IPYTHON:
            layout.addWidget(QLabel("IPython not available\npip install qtconsole"))
            return
        self.kernel_mgr = QtInProcessKernelManager()
        self.kernel_mgr.start_kernel()
        self.kernel_client = self.kernel_mgr.client()
        self.kernel_client.start_channels()
        self.console = RichIPythonWidget()
        self.console.kernel_manager = self.kernel_mgr
        self.console.kernel_client = self.kernel_client
        self.console.setStyleSheet("QPlainTextEdit, QTextEdit { background-color: #002b36; color: #839496; font-family: 'Consolas', monospace; font-size: 10pt; }")
        self.console.syntax_style = 'monokai'
        layout.addWidget(self.console)
        self.kernel_mgr.kernel.shell.push(self.namespace)
        self.console.execute('print("Scope Viewer v2 Console\\nObjects: scope, config, display, np")')
        
    def update_namespace(self, key: str, value):
        if HAS_IPYTHON: self.kernel_mgr.kernel.shell.push({key: value})


class ScopeViewerV2(QMainWindow):
    """Main application"""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("NI-DAQ Scope Viewer v2")
        self.setGeometry(100, 100, 1600, 900)
        self.setStyleSheet(DARK_STYLE)
        self.obs_config = ObservableConfig()
        self.buffers: Dict[int, CircularBuffer] = {}
        self.acq_thread: Optional[AcquisitionThread] = None
        self.display_mgr: Optional[DisplayManager] = None
        self.fft_mode = False
        self._hdf_file = None
        self._input_type = 'analog'  # 'analog' or 'counter'
        
        # Single acquisition mode state
        self._single_mode = False
        self._single_avg_target = 1
        self._single_avg_count = 0
        
        self._setup_menus()
        self._setup_ui()
        self._connect_signals()
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.cursor_label = QLabel("")
        self.cursor_label.setStyleSheet("color: #FFAA00;")
        self.status_bar.addWidget(self.cursor_label, stretch=1)
        # Add error/message label with elided text and tooltip
        self.error_label = StatusBarLongLabel()
        self.error_label.setMinimumWidth(200)
        self.status_bar.addWidget(self.error_label, stretch=2)
        self.mem_label = QLabel("Mem: -- MB")
        self.status_bar.addPermanentWidget(self.mem_label)
        self.mem_timer = QTimer()
        self.mem_timer.timeout.connect(self._update_mem)
        self.mem_timer.start(2000)
        # Try to load last state
        self._try_load_last_state()
        
    def _setup_menus(self):
        """Setup menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        
        # State management
        save_state_action = file_menu.addAction("Save State...")
        save_state_action.setShortcut("Ctrl+S")
        save_state_action.triggered.connect(lambda: self.save_state())
        
        load_state_action = file_menu.addAction("Load State...")
        load_state_action.setShortcut("Ctrl+O")
        load_state_action.triggered.connect(lambda: self.load_state())
        
        file_menu.addSeparator()
        
        # Data saving
        save_scope_action = file_menu.addAction("Save Scope Data...")
        save_scope_action.triggered.connect(lambda: self.save_current_scope_data())
        
        save_fft_action = file_menu.addAction("Save FFT Data...")
        save_fft_action.triggered.connect(lambda: self.save_current_fft_data())
        
        file_menu.addSeparator()
        
        # HDF streaming
        self.start_stream_action = file_menu.addAction("Start HDF5 Stream...")
        self.start_stream_action.triggered.connect(lambda: self.start_hdf_stream())
        
        self.stop_stream_action = file_menu.addAction("Stop HDF5 Stream")
        self.stop_stream_action.triggered.connect(lambda: self.stop_hdf_stream())
        self.stop_stream_action.setEnabled(False)
        
        file_menu.addSeparator()
        
        exit_action = file_menu.addAction("Exit")
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        
        # View menu
        view_menu = menubar.addMenu("&View")
        
        # Reset View (at top for easy access)
        reset_view_action = view_menu.addAction("Reset View")
        reset_view_action.setShortcut("Ctrl+R")
        reset_view_action.triggered.connect(self._reset_view)
        
        view_menu.addSeparator()
        
        # Dock visibility toggles (will be populated after docks are created)
        self._view_menu = view_menu
        
    def _populate_view_menu(self):
        """Populate View menu with dock toggle actions (called after docks are created)"""
        # Scope dock
        self.scope_action = self._view_menu.addAction("Scope")
        self.scope_action.setCheckable(True)
        self.scope_action.setChecked(True)
        self.scope_action.triggered.connect(lambda checked: self.scope_dock.setVisible(checked))
        self.scope_dock.visibilityChanged.connect(self.scope_action.setChecked)
        
        # FFT dock
        self.fft_action = self._view_menu.addAction("FFT")
        self.fft_action.setCheckable(True)
        self.fft_action.setChecked(False)
        self.fft_action.triggered.connect(lambda checked: self.fft_dock.setVisible(checked))
        self.fft_dock.visibilityChanged.connect(self.fft_action.setChecked)
        
        # Trends dock
        self.trends_action = self._view_menu.addAction("Trends")
        self.trends_action.setCheckable(True)
        self.trends_action.setChecked(False)
        self.trends_action.triggered.connect(lambda checked: self.trends_dock.setVisible(checked))
        self.trends_dock.visibilityChanged.connect(self.trends_action.setChecked)
        
        # Control dock
        self.control_action = self._view_menu.addAction("Control Panel")
        self.control_action.setCheckable(True)
        self.control_action.setChecked(True)
        self.control_action.triggered.connect(lambda checked: self.control_dock.setVisible(checked))
        self.control_dock.visibilityChanged.connect(self.control_action.setChecked)
        
        # Console dock
        self.console_action = self._view_menu.addAction("Console")
        self.console_action.setCheckable(True)
        self.console_action.setChecked(True)
        self.console_action.triggered.connect(lambda checked: self.console_dock.setVisible(checked))
        self.console_dock.visibilityChanged.connect(self.console_action.setChecked)
    
    def _reset_view(self):
        """Reset all docks to default layout"""
        # Show main docks
        self.scope_dock.show()
        self.control_dock.show()
        self.console_dock.show()
        
        # Hide optional docks
        self.fft_dock.hide()
        self.trends_dock.hide()
        
        # Reset dock positions
        self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, self.scope_dock)
        self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, self.fft_dock)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.control_dock)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.trends_dock)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.console_dock)
        
        # Tabify bottom docks
        self.tabifyDockWidget(self.trends_dock, self.console_dock)
        
        # Set control dock width
        self.control_dock.setMaximumWidth(400)
        
        # Update toggle states in Control tab
        self.control_tab.scope_toggle.setChecked(True)
        self.control_tab.trends_toggle.setChecked(False)
        
        self.status_bar.showMessage("View reset to default", 2000)
        print("[MAIN] View reset to default layout")
        
    def _try_load_last_state(self):
        """Try to load the last saved state on startup"""
        # Always reset view first to ensure consistent starting state
        self._reset_view()
        
        try:
            last_state = self.get_state_directory() / "last_state.json"
            if last_state.exists():
                self.load_state(str(last_state))
        except Exception as e:
            print(f"Could not load last state: {e}")
        
    def _setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        QVBoxLayout(central)
        
        self.scope_dock = QDockWidget("Scope", self)
        self.scope_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.scope_widget = ScopeWidget()
        self.scope_dock.setWidget(self.scope_widget)
        self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, self.scope_dock)
        
        self.fft_dock = QDockWidget("FFT", self)
        self.fft_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.fft_widget = FFTWidget()
        self.fft_dock.setWidget(self.fft_widget)
        self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, self.fft_dock)
        self.fft_dock.hide()
        
        self.trends_dock = QDockWidget("Trends", self)
        self.trends_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.trends_widget = TrendsWidget(max_seconds=60.0)
        self.trends_dock.setWidget(self.trends_widget)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.trends_dock)
        self.trends_dock.hide()
        
        self.control_dock = QDockWidget("Control", self)
        self.control_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        
        # Main container for control dock
        control_container = QWidget()
        control_layout = QVBoxLayout(control_container)
        control_layout.setContentsMargins(4, 4, 4, 4)
        control_layout.setSpacing(6)
        
        # Common controls at top (Run/Stop, Single, Mode selector)
        common_group = QGroupBox()
        common_layout = QVBoxLayout(common_group)
        common_layout.setContentsMargins(4, 4, 4, 4)
        
        # Run/Stop and Single buttons in a row
        btn_layout = QHBoxLayout()
        
        self.run_btn = QPushButton("▷ Run")
        self.run_btn.setCheckable(True)
        self.run_btn.setStyleSheet("background-color: #2d5a2d; font-size: 14px; padding: 8px;")
        self.run_btn.clicked.connect(self._on_run_clicked)
        btn_layout.addWidget(self.run_btn)
        
        self.single_btn = QPushButton("◇ Single")
        self.single_btn.setStyleSheet("background-color: #2d4a5a; font-size: 14px; padding: 8px;")
        self.single_btn.clicked.connect(self._on_single_clicked)
        btn_layout.addWidget(self.single_btn)
        
        common_layout.addLayout(btn_layout)
        
        # Display Mode selector (Time/FFT)
        mode_layout = QHBoxLayout()
        mode_layout.addWidget(QLabel("Display:"))
        self.mode_selector = SegmentedButton(["Time", "FFT"])
        self.mode_selector.valueChanged.connect(self._on_mode_changed)
        mode_layout.addWidget(self.mode_selector, stretch=1)
        common_layout.addLayout(mode_layout)
        
        control_layout.addWidget(common_group)
        
        # Tabs below common controls
        tabs = QTabWidget()
        self.control_tab = ControlTab(self.obs_config)
        tabs.addTab(self.control_tab, "Control")
        self.scope_tab = ScopeTab()
        tabs.addTab(self.scope_tab, "Scope")
        self.fft_tab = FFTTab()
        tabs.addTab(self.fft_tab, "FFT")
        self.trigger_tab = TriggerClockTab()
        tabs.addTab(self.trigger_tab, "Trigger")
        control_layout.addWidget(tabs)
        
        self.control_dock.setWidget(control_container)
        self.control_dock.setMaximumWidth(400)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.control_dock)
        
        self.console_dock = QDockWidget("Console", self)
        self.console_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetFloatable | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.console_widget = IPythonWidget({'scope': self, 'config': self.obs_config, 'np': np})
        self.console_dock.setWidget(self.console_widget)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.console_dock)
        self.tabifyDockWidget(self.trends_dock, self.console_dock)
        
        # Populate View menu now that docks exist
        self._populate_view_menu()
        
    def _connect_signals(self):
        # Control tab signals
        self.control_tab.scope_toggle.toggled.connect(lambda c: self.scope_dock.setVisible(c))
        self.control_tab.trends_toggle.toggled.connect(self._on_trends_toggled)
        self.control_tab.clear_trends_btn.clicked.connect(self.trends_widget.clear_and_restart)
        self.control_tab.trends_length.currentTextChanged.connect(self._on_trends_length)
        self.control_tab.deviceChanged.connect(self._on_device_changed)
        self.control_tab.sampleRateChanged.connect(self._on_control_rate_changed)
        self.control_tab.chunkSizeChanged.connect(self._on_chunk_size_changed)
        self.control_tab.channelRemoved.connect(self._on_channel_removed)
        self.control_tab.inputTypeChanged.connect(self._on_input_type_changed)
        # Scope tab signals
        self.scope_tab.chunkSizeChanged.connect(self._on_chunk_size_changed)
        self.scope_tab.pointsChanged.connect(self._on_scope_points)
        self.scope_tab.histogramChanged.connect(self._on_histogram_changed)
        self.scope_tab.averagingChanged.connect(self._on_scope_averaging)
        self.scope_tab.cursorsChanged.connect(self._on_scope_cursors)
        # FFT tab signals
        self.fft_tab.xLogChanged.connect(self.fft_widget.set_x_log)
        self.fft_tab.yLogChanged.connect(self.fft_widget.set_y_log)
        self.fft_tab.useDbChanged.connect(self.fft_widget.set_use_db)
        self.fft_tab.displayModeChanged.connect(self.fft_widget.set_display_mode)
        self.fft_tab.freqUnitChanged.connect(self.fft_widget.set_freq_unit)
        self.fft_tab.correctionChanged.connect(self._on_correction_changed)
        self.fft_tab.fftSizeChanged.connect(self._on_chunk_size_changed)  # FFT size = chunk size
        self.fft_tab.windowChanged.connect(self._on_fft_window)
        self.fft_tab.peaks_toggle.toggled.connect(self.fft_widget.toggle_peaks)
        self.fft_tab.averagingChanged.connect(self._on_fft_averaging)
        self.fft_tab.cursorsChanged.connect(self._on_fft_cursors)
        self.trigger_tab.triggerConfigChanged.connect(self._on_trigger_config_changed)
        self.trigger_tab.clockConfigChanged.connect(self._on_clock_config_changed)
        self.obs_config.sampleRateChanged.connect(self._on_sample_rate_changed)
        self.scope_widget.cursorMoved.connect(self._on_scope_cursor_moved)
        self.fft_widget.cursorMoved.connect(self._on_fft_cursor_moved)
    
    def _on_run_clicked(self):
        """Handle Run/Stop button click"""
        if self.run_btn.isChecked():
            self._single_mode = False  # Continuous mode
            self.run_btn.setText("■ Stop")
            self.run_btn.setStyleSheet("background-color: #5a2d2d; font-size: 14px; padding: 8px;")
            self.start()
        else:
            self.run_btn.setText("▷ Run")
            self.run_btn.setStyleSheet("background-color: #2d5a2d; font-size: 14px; padding: 8px;")
            self.stop()
    
    def _on_single_clicked(self):
        """Handle Single acquisition button click"""
        if self.obs_config.is_running:
            # If already running, stop first
            self.stop()
            return
        
        # Get averaging count based on current mode
        if self.fft_mode:
            avg_enabled = self.fft_tab.avg_toggle.isChecked()
            avg_count = int(self.fft_tab.avg_count.currentText()) if avg_enabled else 1
        else:
            avg_enabled = self.scope_tab.avg_toggle.isChecked()
            avg_count = int(self.scope_tab.avg_count.currentText()) if avg_enabled else 1
        
        # Setup single acquisition mode
        self._single_mode = True
        self._single_avg_target = avg_count
        self._single_avg_count = 0
        
        # Update button appearance
        self.single_btn.setText("■ Stop")
        self.single_btn.setStyleSheet("background-color: #5a4a2d; font-size: 14px; padding: 8px;")
        
        # Clear averaging buffers for fresh start
        self.scope_widget.clear_averaging()
        self.fft_widget.clear_averaging()
        
        self.error_label.set_msg(f"Single: Acquiring {avg_count} chunk(s)...", error=False)
        self.start()
    
    def _set_running_ui(self, running: bool):
        """Update Run/Single button UI state"""
        self.run_btn.setChecked(running and not self._single_mode)
        
        if running:
            if self._single_mode:
                self.single_btn.setText("■ Stop")
                self.single_btn.setStyleSheet("background-color: #5a4a2d; font-size: 14px; padding: 8px;")
                self.run_btn.setEnabled(False)
            else:
                self.run_btn.setText("■ Stop")
                self.run_btn.setStyleSheet("background-color: #5a2d2d; font-size: 14px; padding: 8px;")
                self.single_btn.setEnabled(False)
        else:
            self.run_btn.setText("▷ Run")
            self.run_btn.setStyleSheet("background-color: #2d5a2d; font-size: 14px; padding: 8px;")
            self.single_btn.setText("◇ Single")
            self.single_btn.setStyleSheet("background-color: #2d4a5a; font-size: 14px; padding: 8px;")
            self.run_btn.setEnabled(True)
            self.single_btn.setEnabled(True)
            self._single_mode = False
    
    def _on_mode_changed(self, idx: int):
        """Handle display mode change (Time/FFT)"""
        self.fft_mode = idx == 1
        if self.fft_mode:
            self.scope_dock.hide()
            self.fft_dock.show()
        else:
            self.fft_dock.hide()
            self.scope_dock.show()
        if self.display_mgr:
            self.display_mgr.set_fft_enabled(self.fft_mode)
        # Enable/disable FFT controls
        self.fft_tab.set_fft_mode_enabled(self.fft_mode)
        self.error_label.set_msg(f"Display mode: {'FFT' if self.fft_mode else 'Time'}", error=False)
    
    def _on_correction_changed(self, power_correction: bool):
        """Handle FFT correction mode change (Amplitude/Power)"""
        if self.display_mgr:
            self.display_mgr.set_power_correction(power_correction)
        self.fft_widget.set_power_correction(power_correction)
        mode_str = "Power" if power_correction else "Amplitude"
        self.error_label.set_msg(f"FFT correction: {mode_str}", error=False)
        
    def start(self):
        if self.acq_thread and self.acq_thread.isRunning(): 
            return
            
        configs = self.control_tab.get_all_configs()
        self.obs_config._config.channels = configs
        
        # Get trigger and clock config
        self.obs_config._config.trigger = self.trigger_tab.get_trigger_config()
        self.obs_config._config.clock = self.trigger_tab.get_clock_config()
        
        # Use clock rate as sample rate
        clock_cfg = self.obs_config._config.clock
        if clock_cfg.source == ClockSource.INTERNAL:
            self.obs_config._config.sample_rate = clock_cfg.rate
        else:
            # For external clock, use the expected rate setting
            self.obs_config._config.sample_rate = clock_cfg.rate
        
        # Get chunk size = fft_size
        chunk_size = self.fft_tab.get_fft_size()
        sample_rate = self.obs_config.sample_rate
        
        # Calculate natural update rate
        natural_rate = sample_rate / chunk_size
        actual_rate = min(natural_rate, 30.0)
        
        self.error_label.set_msg(
            f"Starting: {sample_rate/1e3:.1f} kSa/s | Chunk: {chunk_size} | Update: {actual_rate:.1f} Hz",
            error=False
        )
        
        # Create circular buffers for trends history
        buffer_duration = self.control_tab.get_buffer_duration()
        buffer_samples = int(sample_rate * buffer_duration)
        n_channels = sum(1 for c in configs if c.enabled)
        total_mb = buffer_samples * 4 * n_channels / 1024 / 1024
        
        self.buffers.clear()
        for i, ch in enumerate(configs):
            if ch.enabled: 
                self.buffers[i] = CircularBuffer(buffer_samples)
        
        # Create display manager
        self.display_mgr = DisplayManager(self.obs_config, self.buffers)
        self.display_mgr.scopeDataReady.connect(self._on_scope_data)
        self.display_mgr.fftDataReady.connect(self._on_fft_data)
        self.display_mgr.statsReady.connect(self._on_stats)
        self.display_mgr.displayRateUpdated.connect(self._on_display_rate_updated)
        self.display_mgr.hangDetected.connect(self._on_hang_detected)
        self.display_mgr.chunkProcessed.connect(self._on_chunk_processed)
        self.display_mgr.set_scope_points(self.scope_tab.get_points())
        self.display_mgr.set_fft_size(chunk_size)
        self.display_mgr.set_fft_window(self.fft_tab.get_window())
        self.display_mgr.set_fft_enabled(self.fft_mode)
        
        self.fft_widget.set_sample_rate(sample_rate)
        self.fft_widget.set_fft_size(chunk_size)
        
        # Set power correction mode
        power_correction = self.fft_tab.is_power_correction()
        self.display_mgr.set_power_correction(power_correction)
        self.fft_widget.set_power_correction(power_correction)
        
        if HAS_IPYTHON: 
            self.console_widget.update_namespace('display', self.display_mgr)
        
        if self.control_tab.trends_toggle.isChecked(): 
            self.trends_widget.start()
        
        # Create acquisition thread with chunk_size = fft_size
        self.acq_thread = AcquisitionThread(self.obs_config, self.buffers, chunk_size=chunk_size)
        self.acq_thread.set_input_type(self._input_type)  # Set input type (analog/counter)
        self.acq_thread.errorOccurred.connect(self._on_error)
        self.acq_thread.statusChanged.connect(self._on_status)
        # Use QueuedConnection for cross-thread signal
        self.acq_thread.chunkReady.connect(self.display_mgr.on_chunk_ready, Qt.ConnectionType.QueuedConnection)
        self.acq_thread.start()
        
        self.display_mgr.start()
        self.obs_config.is_running = True
        self._set_running_ui(True)
        self.control_tab.set_running(True)
        self._update_info_labels()
        
    def stop(self):
        if self.display_mgr: 
            self.display_mgr.stop()
            self.display_mgr = None
        if self.acq_thread: 
            self.acq_thread.request_stop()
            # Wait up to 5 seconds for thread to finish
            if not self.acq_thread.wait(5000):
                self.error_label.set_msg("Warning: Acquisition thread did not stop cleanly", error=True)
            self.acq_thread = None
        self.trends_widget.stop()
        self.obs_config.is_running = False
        self._set_running_ui(False)
        self.control_tab.set_running(False)
        self.error_label.set_msg("Stopped", error=False)
        gc.collect()
    
    def _on_display_rate_updated(self, rate: float):
        """Update status bar with actual display rate"""
        sample_rate = self.obs_config.sample_rate
        chunk_size = self.fft_tab.get_fft_size()
        self.error_label.set_msg(
            f"Running: {sample_rate/1e3:.1f} kSa/s | Chunk: {chunk_size} | Display: {rate:.1f} Hz",
            error=False
        )
    
    def _on_chunk_processed(self, chunk_count: int):
        """Handle chunk processed - check for single mode completion"""
        if not self._single_mode:
            return
        
        self._single_avg_count = chunk_count
        
        # Update status
        self.error_label.set_msg(
            f"Single: {chunk_count}/{self._single_avg_target} chunks",
            error=False
        )
        
        # Check if we've reached the target
        if chunk_count >= self._single_avg_target:
            # Use QTimer.singleShot to stop after current event processing
            QTimer.singleShot(0, self._complete_single_acquisition)
    
    def _complete_single_acquisition(self):
        """Complete single acquisition - stop and show result"""
        self.stop()
        self.error_label.set_msg(
            f"Single complete: {self._single_avg_target} chunk(s) averaged",
            error=False
        )
        
    def _on_sample_rate_changed(self, rate: float):
        self._update_info_labels()
        self.fft_widget.set_sample_rate(rate)
        if self.display_mgr:
            self.status_bar.showMessage(f"Restarting with {rate/1e3:.1f} kSa/s", 3000)
            self.stop()
            self.start()
    
    def _on_device_changed(self, device_info: DAQDeviceInfo):
        """Handle device change - update trigger tab with PFI terminals"""
        if device_info and device_info.pfi_terminals:
            self.trigger_tab.update_pfi_terminals(device_info.pfi_terminals)
    
    def _on_trigger_config_changed(self, config: TriggerConfig):
        """Handle trigger config change"""
        self.obs_config._config.trigger = config
        mode_str = config.mode.value
        self.status_bar.showMessage(f"Trigger: {mode_str}", 2000)
    
    def _on_clock_config_changed(self, config: ClockConfig):
        """Handle clock config change from TriggerClockTab - sync to ControlTab"""
        self.obs_config._config.clock = config
        # Sync rate to ControlTab
        self.control_tab.set_sample_rate(config.rate)
        self.obs_config.sample_rate = config.rate
        
        if config.source == ClockSource.INTERNAL:
            self.status_bar.showMessage(f"Clock: Internal {config.rate/1e3:.1f} kHz", 2000)
        else:
            self.status_bar.showMessage(f"Clock: External {config.external_source} @ {config.rate/1e3:.1f} kHz", 2000)
    
    def _on_control_rate_changed(self, rate: float):
        """Handle sample rate change from ControlTab - sync to TriggerClockTab"""
        self.trigger_tab.set_sample_rate(rate)
        self.obs_config.sample_rate = rate
        self._update_info_labels()
    
    def _on_channel_removed(self, ch_idx: int):
        """Handle channel removal - remove curves from plots"""
        self.scope_widget.remove_channel(ch_idx)
        self.fft_widget.remove_channel(ch_idx)
        self.status_bar.showMessage(f"Channel {ch_idx + 1} removed", 2000)
    
    def _on_input_type_changed(self, input_type: str):
        """Handle input type change (analog/counter)"""
        self._input_type = input_type
        
        # Update Y-axis labels based on input type
        if input_type == 'counter':
            self.scope_widget.setLabel('left', 'Count Rate', 'Hz')
            self.trends_widget.setLabel('left', 'Count Rate', 'Hz')
            self.error_label.set_msg("Input: Counter (Edges)", error=False)
        else:
            self.scope_widget.setLabel('left', 'Voltage', 'V')
            self.trends_widget.setLabel('left', 'Value', 'µV')
            self.error_label.set_msg("Input: Analog (Voltage)", error=False)
            
    def _update_info_labels(self):
        """Update info labels with current settings"""
        rate = self.obs_config.sample_rate
        chunk_size = self.fft_tab.get_fft_size()
        self.scope_tab.update_chunk_info(rate, chunk_size)
        self.fft_tab.update_resolution(rate)
        
    def _on_scope_data(self, ch_idx: int, time_data, volt_data):
        self.scope_widget.update_data(ch_idx, time_data, volt_data)
        
    def _on_fft_data(self, ch_idx: int, freq, mag):
        self.fft_widget.update_spectrum(ch_idx, freq, mag)
        
    def _on_stats(self, stats: dict):
        """Update status bar with multi-channel statistics"""
        if stats['mean'] and not self.scope_widget._cursors_enabled and not self.fft_widget._cursors_enabled:
            parts = []
            for i, (mean, std) in enumerate(zip(stats['mean'], stats['std'])):
                parts.append(f"Ch{i+1}: {mean*1e3:.2f}±{std*1e3:.2f} mV")
            self.cursor_label.setText(" | ".join(parts))
        if self.control_tab.trends_toggle.isChecked(): 
            self.trends_widget.add_stats(stats)
        # Stream to HDF if enabled
        if hasattr(self, '_hdf_file') and self._hdf_file is not None:
            self._stream_stats_to_hdf(stats)
            
    def _on_error(self, msg: str): 
        self.error_label.set_msg(f"Error: {msg}", error=True)
        
    def _on_status(self, status: str):
        self.error_label.set_msg(status, error=False)
    
    def _on_hang_detected(self):
        """Handle acquisition hang - auto-restart"""
        self.error_label.set_msg("Hang detected - restarting acquisition...", error=True)
        # Stop and restart
        self.stop()
        QTimer.singleShot(500, self.start)  # Restart after 500ms
        
    def _on_scope_points(self, points: int):
        if self.display_mgr: 
            self.display_mgr.set_scope_points(points)
        self._update_info_labels()
    def _on_scope_averaging(self, enabled: bool, count: int):
        if count == 0: self.scope_widget.clear_averaging()
        else: self.scope_widget.set_averaging(enabled, count)
    def _on_scope_cursors(self, enabled: bool):
        self.scope_widget.set_cursors_enabled(enabled)
        if not enabled: self.cursor_label.setText("")
    def _on_histogram_changed(self, enabled: bool):
        """Handle histogram toggle"""
        self.scope_widget.set_histogram_enabled(enabled)
        self.error_label.set_msg(f"Histogram: {'On' if enabled else 'Off'}", error=False)
    def _on_scope_cursor_moved(self, cursor_data: dict):
        """Handle mouse cursor movement on scope"""
        if not cursor_data: 
            return
        x = cursor_data.get('x', 0)
        y = cursor_data.get('y', 0)
        channel_values = cursor_data.get('channel_values', {})
        
        # Format time
        if abs(x) >= 1:
            t_str = f"T: {x:.4f} s"
        elif abs(x) >= 1e-3:
            t_str = f"T: {x*1e3:.3f} ms"
        else:
            t_str = f"T: {x*1e6:.2f} µs"
        
        # Format voltage at cursor
        v_str = f"Cursor: {y*1e3:.2f} mV"
        
        # Format channel values
        ch_parts = []
        for ch_idx, val in sorted(channel_values.items()):
            ch_parts.append(f"Ch{ch_idx+1}: {val*1e3:.2f} mV")
        
        text = f"{t_str} | {v_str}"
        if ch_parts:
            text += " | " + " | ".join(ch_parts)
        self.cursor_label.setText(text)
        
    def _on_chunk_size_changed(self, size: int):
        """Handle chunk/FFT size change - syncs all tabs and restarts if needed"""
        # Sync all three locations
        self.control_tab.set_chunk_size(size)
        self.scope_tab.set_chunk_size(size)
        self.fft_tab.set_fft_size(size)
        
        # Update widget
        self.fft_widget.set_fft_size(size)
        self._update_info_labels()
        
        # Restart acquisition if running (chunk size changed)
        if self.obs_config.is_running:
            self.error_label.set_msg(f"Restarting with chunk size {size}...", error=False)
            self.stop()
            self.start()
        else:
            self.error_label.set_msg(f"Chunk size: {size}", error=False)
            
    def _on_fft_window(self, window: str):
        if self.display_mgr: 
            self.display_mgr.set_fft_window(window)
        self.error_label.set_msg(f"FFT window: {window}", error=False)
        
    def _on_fft_averaging(self, enabled: bool, count: int):
        self.fft_widget.set_averaging(enabled, count)
        if enabled:
            self.error_label.set_msg(f"FFT averaging: {count}x", error=False)
        else:
            self.error_label.set_msg("FFT averaging: Off", error=False)
            
    def _on_fft_cursors(self, enabled: bool):
        self.fft_widget.set_cursors_enabled(enabled)
        if not enabled: 
            self.cursor_label.setText("")
            
    def _on_fft_cursor_moved(self, cursor_data: dict):
        """Handle mouse cursor movement on FFT"""
        if not cursor_data: 
            return
        freq = cursor_data.get('freq', 0)
        mag_display = cursor_data.get('mag_display', 0)
        channel_mags = cursor_data.get('channel_mags', {})
        
        # Format frequency
        if freq >= 1e6:
            f_str = f"F: {freq/1e6:.3f} MHz"
        elif freq >= 1e3:
            f_str = f"F: {freq/1e3:.3f} kHz"
        else:
            f_str = f"F: {freq:.1f} Hz"
        
        # Format magnitude for each channel
        ch_parts = []
        for ch_idx, mag in sorted(channel_mags.items()):
            mag_db = 20 * np.log10(mag + 1e-12)
            ch_parts.append(f"Ch{ch_idx+1}: {mag_db:.1f} dB")
        
        text = f_str
        if ch_parts:
            text += " | " + " | ".join(ch_parts)
        self.cursor_label.setText(text)
    def _on_trends_toggled(self, checked: bool):
        self.trends_dock.setVisible(checked)
        if checked and self.obs_config.is_running: self.trends_widget.start()
        else: self.trends_widget.stop()
    def _on_trends_length(self, text: str):
        try: self.trends_widget.set_max_history(float(text.replace('s', '').strip()))
        except: pass
        
    def _update_mem(self):
        try:
            mb = psutil.Process().memory_info().rss / 1024 / 1024
            self.mem_label.setText(f"Mem: {mb:.0f} MB")
            color = "#ff6666" if mb > 800 else "#ffcc66" if mb > 500 else "#66ff66"
            self.mem_label.setStyleSheet(f"color: {color};")
        except: pass
    
    # =========================================================================
    # STATE MANAGEMENT
    # =========================================================================
    
    def get_state_directory(self) -> Path:
        """Get the directory for state files"""
        state_dir = Path.home() / "NVExperiment" / "ScopeStates"
        state_dir.mkdir(parents=True, exist_ok=True)
        return state_dir
    
    def save_state(self, filepath: Optional[str] = None):
        """Save current GUI state to JSON file"""
        if filepath is None:
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Save State", str(self.get_state_directory()),
                "JSON Files (*.json)"
            )
            if not filepath:
                return
        
        filepath = Path(filepath)
        if not filepath.suffix:
            filepath = filepath.with_suffix('.json')
        
        state = {
            'version': '2.0',
            'timestamp': datetime.now().isoformat(),
            'acquisition': {
                'device': self.obs_config.device,
                'sample_rate': self.obs_config.sample_rate,
                'buffer_duration': self.control_tab.get_buffer_duration(),
            },
            'channels': [
                {
                    'physical_channel': ch.physical_channel,
                    'enabled': ch.enabled,
                    'min_voltage': ch.min_voltage,
                    'max_voltage': ch.max_voltage,
                    'terminal_config': ch.terminal_config,
                }
                for ch in self.control_tab.get_all_configs()
            ],
            'scope': {
                'points': self.scope_tab.get_points(),
                'averaging_enabled': self.scope_tab.avg_toggle.isChecked(),
                'averaging_count': int(self.scope_tab.avg_count.currentText()),
                'cursors_enabled': self.scope_tab.cursor_toggle.isChecked(),
            },
            'fft': {
                'mode': self.mode_selector.value(),  # Now on main window
                'window': self.fft_tab.get_window(),
                'fft_size': self.fft_tab.get_fft_size(),
                'correction': self.fft_tab.correction_combo.currentText(),
                'x_log': self.fft_tab.x_scale.value() == 1,
                'y_log': self.fft_tab.y_log_toggle.isChecked(),
                'use_db': self.fft_tab.db_toggle.isChecked(),
                'spectral_density': self.fft_tab.spectral_density_toggle.isChecked(),
                'power': self.fft_tab.power_toggle.isChecked(),
                'freq_unit': self.fft_tab.freq_unit_combo.currentText(),
                'averaging_enabled': self.fft_tab.avg_toggle.isChecked(),
                'averaging_count': int(self.fft_tab.avg_count.currentText()),
                'show_peaks': self.fft_tab.peaks_toggle.isChecked(),
                'cursors_enabled': self.fft_tab.cursor_toggle.isChecked(),
            },
            'display': {
                'scope_visible': self.control_tab.scope_toggle.isChecked(),
                'trends_visible': self.control_tab.trends_toggle.isChecked(),
                'trends_length': self.control_tab.trends_length.currentText(),
            },
            'window': {
                'width': self.width(),
                'height': self.height(),
                'x': self.x(),
                'y': self.y(),
            }
        }
        
        with open(filepath, 'w') as f:
            json.dump(state, f, indent=2)
        
        self.status_bar.showMessage(f"State saved to {filepath.name}", 3000)
        print(f"✓ State saved to {filepath}")
    
    def load_state(self, filepath: Optional[str] = None):
        """Load GUI state from JSON file"""
        if filepath is None:
            filepath, _ = QFileDialog.getOpenFileName(
                self, "Load State", str(self.get_state_directory()),
                "JSON Files (*.json)"
            )
            if not filepath:
                return
        
        filepath = Path(filepath)
        if not filepath.exists():
            self.status_bar.showMessage(f"File not found: {filepath}", 5000)
            return
        
        try:
            with open(filepath, 'r') as f:
                state = json.load(f)
            
            # Acquisition settings
            if 'acquisition' in state:
                acq = state['acquisition']
                # Set device if available
                idx = self.control_tab.device_combo.findText(acq.get('device', ''))
                if idx >= 0:
                    self.control_tab.device_combo.setCurrentIndex(idx)
                self.control_tab.rate_combo.setCurrentText(
                    self.control_tab._rate_to_str(acq.get('sample_rate', 100000))
                )
                # Buffer duration
                buf_dur = acq.get('buffer_duration', 2.0)
                self.control_tab.buffer_combo.setCurrentText(f"{buf_dur} s")
            
            # Channels - recreate to match saved count
            if 'channels' in state:
                # Remove extra channels
                while len(self.control_tab.channel_widgets) > len(state['channels']):
                    self.control_tab._remove_channel()
                # Add missing channels
                while len(self.control_tab.channel_widgets) < len(state['channels']):
                    self.control_tab._add_channel()
                # Configure channels
                for i, ch_state in enumerate(state['channels']):
                    if i < len(self.control_tab.channel_widgets):
                        w = self.control_tab.channel_widgets[i]
                        w.enable_cb.setChecked(ch_state.get('enabled', True))
                        idx = w.phys_combo.findText(ch_state.get('physical_channel', f'ai{i}'))
                        if idx >= 0:
                            w.phys_combo.setCurrentIndex(idx)
                        idx = w.term_combo.findText(ch_state.get('terminal_config', 'RSE'))
                        if idx >= 0:
                            w.term_combo.setCurrentIndex(idx)
                        w.min_spin.setValue(ch_state.get('min_voltage', -10))
                        w.max_spin.setValue(ch_state.get('max_voltage', 10))
            
            # Scope settings
            if 'scope' in state:
                sc = state['scope']
                dur = sc.get('duration', 0.05)
                self.scope_tab.points_combo.setCurrentText(str(sc.get('points', 5000)))
                self.scope_tab.avg_toggle.setChecked(sc.get('averaging_enabled', False))
                self.scope_tab.avg_count.setCurrentText(str(sc.get('averaging_count', 8)))
                self.scope_tab.cursor_toggle.setChecked(sc.get('cursors_enabled', False))
            
            # FFT settings
            if 'fft' in state:
                fft = state['fft']
                self.mode_selector.setValue(fft.get('mode', 0))  # Now on main window
                self._on_mode_changed(fft.get('mode', 0))  # Apply the mode change
                self.fft_tab.window_combo.setCurrentText(fft.get('window', 'Hann'))
                self.fft_tab.fft_size_combo.setCurrentText(str(fft.get('fft_size', 8192)))
                self.fft_tab.correction_combo.setCurrentText(fft.get('correction', 'Amplitude'))
                self.fft_tab.x_scale.setValue(1 if fft.get('x_log', False) else 0)
                self.fft_tab.y_log_toggle.setChecked(fft.get('y_log', False))
                self.fft_tab.db_toggle.setChecked(fft.get('use_db', True))
                self.fft_tab.spectral_density_toggle.setChecked(fft.get('spectral_density', False))
                self.fft_tab.power_toggle.setChecked(fft.get('power', False))
                self.fft_tab.freq_unit_combo.setCurrentText(fft.get('freq_unit', 'Hz'))
                self.fft_tab.avg_toggle.setChecked(fft.get('averaging_enabled', False))
                self.fft_tab.avg_count.setCurrentText(str(fft.get('averaging_count', 8)))
                self.fft_tab.peaks_toggle.setChecked(fft.get('show_peaks', True))
                self.fft_tab.cursor_toggle.setChecked(fft.get('cursors_enabled', False))
            
            # Display settings
            if 'display' in state:
                disp = state['display']
                self.control_tab.scope_toggle.setChecked(disp.get('scope_visible', True))
                self.control_tab.trends_toggle.setChecked(disp.get('trends_visible', False))
                self.control_tab.trends_length.setCurrentText(disp.get('trends_length', '60 s'))
            
            # Window geometry
            if 'window' in state:
                win = state['window']
                self.resize(win.get('width', 1600), win.get('height', 900))
                self.move(win.get('x', 100), win.get('y', 100))
            
            self.status_bar.showMessage(f"State loaded from {filepath.name}", 3000)
            print(f"✓ State loaded from {filepath}")
            
        except Exception as e:
            self.status_bar.showMessage(f"Error loading state: {e}", 5000)
            print(f"Error loading state: {e}")
    
    # =========================================================================
    # HDF5 DATA STREAMING
    # =========================================================================
    
    def start_hdf_stream(self, filepath: Optional[str] = None):
        """Start streaming data to HDF5 file"""
        if filepath is None:
            default_name = f"scope_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.h5"
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Save Data Stream", str(Path.home() / "NVExperiment" / "ScopeData" / default_name),
                "HDF5 Files (*.h5 *.hdf5)"
            )
            if not filepath:
                return
        
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        # Create HDF5 file
        self._hdf_file = h5py.File(filepath, 'w')
        self._hdf_filepath = filepath
        
        # Store metadata
        meta = self._hdf_file.create_group('metadata')
        meta.attrs['timestamp'] = datetime.now().isoformat()
        meta.attrs['sample_rate'] = self.obs_config.sample_rate
        meta.attrs['device'] = self.obs_config.device
        
        # Create groups for different data types
        self._hdf_scope_grp = self._hdf_file.create_group('scope')
        self._hdf_fft_grp = self._hdf_file.create_group('fft')
        self._hdf_trends_grp = self._hdf_file.create_group('trends')
        
        # Initialize datasets (extendable)
        self._hdf_trends_idx = 0
        n_channels = sum(1 for c in self.control_tab.get_all_configs() if c.enabled)
        self._hdf_trends_grp.create_dataset(
            'timestamps', shape=(0,), maxshape=(None,), dtype='f8'
        )
        self._hdf_trends_grp.create_dataset(
            'mean', shape=(0, n_channels), maxshape=(None, n_channels), dtype='f4'
        )
        self._hdf_trends_grp.create_dataset(
            'std', shape=(0, n_channels), maxshape=(None, n_channels), dtype='f4'
        )
        
        # Update menu items
        self.start_stream_action.setEnabled(False)
        self.stop_stream_action.setEnabled(True)
        
        self.status_bar.showMessage(f"Started HDF5 stream: {filepath.name}", 3000)
        print(f"✓ HDF5 streaming started: {filepath}")
    
    def stop_hdf_stream(self):
        """Stop HDF5 streaming and close file"""
        if hasattr(self, '_hdf_file') and self._hdf_file is not None:
            self._hdf_file.close()
            self._hdf_file = None
            
            # Update menu items
            self.start_stream_action.setEnabled(True)
            self.stop_stream_action.setEnabled(False)
            
            self.status_bar.showMessage("HDF5 stream stopped", 3000)
            print(f"✓ HDF5 streaming stopped: {self._hdf_filepath}")
    
    def _stream_stats_to_hdf(self, stats: dict):
        """Stream statistics to HDF5 file"""
        if not hasattr(self, '_hdf_file') or self._hdf_file is None:
            return
        
        try:
            # Extend datasets
            ts_ds = self._hdf_trends_grp['timestamps']
            mean_ds = self._hdf_trends_grp['mean']
            std_ds = self._hdf_trends_grp['std']
            
            # Resize
            new_size = ts_ds.shape[0] + 1
            ts_ds.resize(new_size, axis=0)
            mean_ds.resize(new_size, axis=0)
            std_ds.resize(new_size, axis=0)
            
            # Write data
            ts_ds[-1] = stats.get('timestamp', datetime.now().timestamp())
            mean_ds[-1, :len(stats['mean'])] = stats['mean']
            std_ds[-1, :len(stats['std'])] = stats['std']
            
        except Exception as e:
            print(f"Error streaming to HDF5: {e}")
    
    def save_current_scope_data(self, filepath: Optional[str] = None):
        """Save current scope buffer to HDF5 file"""
        if not self.buffers:
            self.status_bar.showMessage("No data to save", 3000)
            return
        
        if filepath is None:
            default_name = f"scope_snapshot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.h5"
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Save Scope Data", str(Path.home() / "NVExperiment" / "ScopeData" / default_name),
                "HDF5 Files (*.h5 *.hdf5)"
            )
            if not filepath:
                return
        
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with h5py.File(filepath, 'w') as f:
            # Metadata
            meta = f.create_group('metadata')
            meta.attrs['timestamp'] = datetime.now().isoformat()
            meta.attrs['sample_rate'] = self.obs_config.sample_rate
            meta.attrs['device'] = self.obs_config.device
            
            # Channel data
            for ch_idx, buffer in self.buffers.items():
                data = buffer.read_latest(buffer.available)
                ch_grp = f.create_group(f'channel_{ch_idx}')
                ch_grp.create_dataset('data', data=data, compression='gzip')
                ch_grp.attrs['samples'] = len(data)
                
                # Time axis
                t = np.arange(len(data)) / self.obs_config.sample_rate
                ch_grp.create_dataset('time', data=t, compression='gzip')
        
        self.status_bar.showMessage(f"Scope data saved: {filepath.name}", 3000)
        print(f"✓ Scope data saved to {filepath}")
    
    def save_current_fft_data(self, filepath: Optional[str] = None):
        """Save current FFT data to HDF5 file"""
        if not self.fft_widget._channel_data:
            self.status_bar.showMessage("No FFT data to save", 3000)
            return
        
        if filepath is None:
            default_name = f"fft_snapshot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.h5"
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Save FFT Data", str(Path.home() / "NVExperiment" / "ScopeData" / default_name),
                "HDF5 Files (*.h5 *.hdf5)"
            )
            if not filepath:
                return
        
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with h5py.File(filepath, 'w') as f:
            # Metadata
            meta = f.create_group('metadata')
            meta.attrs['timestamp'] = datetime.now().isoformat()
            meta.attrs['sample_rate'] = self.obs_config.sample_rate
            meta.attrs['fft_size'] = self.fft_tab.get_fft_size()
            meta.attrs['window'] = self.fft_tab.get_window()
            meta.attrs['display_mode'] = self.fft_widget._display_mode.value
            meta.attrs['freq_unit'] = self.fft_widget._freq_unit.label
            
            # FFT data per channel
            for ch_idx, (freq, mag) in self.fft_widget._channel_data.items():
                ch_grp = f.create_group(f'channel_{ch_idx}')
                ch_grp.create_dataset('frequency', data=freq, compression='gzip')
                ch_grp.create_dataset('magnitude', data=mag, compression='gzip')
        
        self.status_bar.showMessage(f"FFT data saved: {filepath.name}", 3000)
        print(f"✓ FFT data saved to {filepath}")
        
    def closeEvent(self, event): 
        self.stop()
        self.stop_hdf_stream()
        # Auto-save state on close
        try:
            auto_save_path = self.get_state_directory() / "last_state.json"
            self.save_state(str(auto_save_path))
        except:
            pass
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
    try:
        myappid = 'aglab.scope'
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    except: pass
    pg.setConfigOptions(antialias=False, useOpenGL=True)
    app = QApplication(sys.argv)
    set_dark_theme(app)

    window = ScopeViewerV2()

    app_icon = QIcon(expt_dir + r"\gui\daq\gui_icon1.png") # .ico is preferred for Windows
    window.setWindowIcon(app_icon)
    app.setWindowIcon(app_icon) # Sets it for the whole application

    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
