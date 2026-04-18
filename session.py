# session.py — NV Experiment Session Manager
# =============================================================================
# Layout: Tabbed left panel (Instruments | Experiment | Config) + plots right
#
# Features:
#   - Instrument init / per-instrument eject / Init PB Only
#   - Load Config: preview in metadata tree without running
#   - Experiment launch via mainControl_*.run() in QThread
#   - Real-time pyqtgraph with:
#       * Per-run overlay (semi-transparent + running mean)
#       * Per-outer-combo separate color group
#       * Shuffled x-axis: pre-allocated from param values, scatter as data arrives
#       * Click legend to show/hide traces
#   - View Sequence: matplotlib popup (needs PB only)
#   - Play Sequence: program + start PB at last scan point (needs PB only)
#   - RUN: experiment in QThread with real-time pyqtgraph
#   - Save: Saved_Data/YYYY-MM-DD/meas_NNN/ with data + metadata
#   - Default autosave ON, per-run overwrite (crash-safe)
# =============================================================================

import sys, time, json, glob, importlib, traceback, subprocess, os
from typing import Optional, List, Dict, Tuple
import numpy as np
from pathlib import Path
from typing import Optional
from datetime import datetime

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QComboBox, QCheckBox, QStatusBar,
    QLineEdit, QMessageBox, QTabWidget, QTreeWidget, QSizePolicy,
    QTreeWidgetItem, QHeaderView, QSpinBox, QGroupBox, QButtonGroup, QGridLayout
)
from PySide6.QtCore import Qt, QThread, Signal, Slot, QTimer, QPropertyAnimation, Property, QEasingCurve
from PySide6.QtGui import QColor, QPainter, QBrush, QPalette
import pyqtgraph as pg

from SGcontrol import SignalGenerator, SignalGenerator_sim
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogOutputTask
import connectionConfig as concfg

# ── IPython in-process console (optional — degrades gracefully) ──────
try:
    from qtconsole.rich_ipython_widget import RichIPythonWidget
    from qtconsole.inprocess import QtInProcessKernelManager
    _HAS_QTCONSOLE = True
except ImportError:
    _HAS_QTCONSOLE = False

# ── palette ──────────────────────────────────────────────────────────────
_PAL = ['#4ec9b0', '#569cd6', '#dcdcaa', '#ce9178', '#c586c0',
        '#9cdcfe', '#d7ba7d', '#b5cea8', '#f44747', '#6a9955']

STATEFILE_DIRECTORY = r'D:\Brateen\Saved_Data\SavedStates\SessionStates'  # Default working directory for state files
DEFAULT_STATE_FILE = 'last_state.json'
SAVE_DIRECTORY = r'D:\Brateen\Saved_Data'
EXPERIMENT_FILE_DIRECTORY = r'D:\Brateen\NV_Experiment'

def _pen(idx, alpha=255, width=1.5, dash=False):
    c = pg.mkColor(_PAL[idx % len(_PAL)])
    c.setAlpha(alpha)
    style = Qt.PenStyle.DashLine if dash else Qt.PenStyle.SolidLine
    return pg.mkPen(c, width=width, style=style)

# ── clickable legend ────────────────────────────────────────────────────────
class ClickableLegend(pg.LegendItem):
    """LegendItem that toggles curve visibility on click."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def mousePressEvent(self, ev):
        pos = ev.pos()
        for sample, label in self.items:
            if sample.boundingRect().translated(
                    sample.pos()).contains(pos) or \
               label.boundingRect().translated(
                    label.pos()).contains(pos):
                # Find the PlotDataItem
                for item in sample.item.parentItem().items \
                        if hasattr(sample.item, 'parentItem') else []:
                    pass
                # Toggle via the sample's item reference
                pdi = sample.item
                if hasattr(pdi, 'curve'):
                    pdi = pdi  # PlotDataItem
                vis = pdi.isVisible()
                pdi.setVisible(not vis)
                label.setOpacity(0.3 if vis else 1.0)
                ev.accept()
                return
        super().mousePressEvent(ev)

# ── metadata tree helper ────────────────────────────────────────────────
def _populate_tree(parent, data):
    if isinstance(data, dict):
        for key, value in data.items():
            item = QTreeWidgetItem(parent, [str(key), ""])
            if isinstance(value, (dict, list)):
                _populate_tree(item, value)
                item.setExpanded(True)
            else:
                item.setText(1, str(value))
    elif isinstance(data, list):
        for i, value in enumerate(data):
            item = QTreeWidgetItem(parent, [f"[{i}]", ""])
            if isinstance(value, (dict, list)):
                _populate_tree(item, value)
                item.setExpanded(True)
            else:
                item.setText(1, str(value))

class ToggleSwitch(QCheckBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedSize(38, 21)
        self._knob_position = 4  # starting position
        self.animation = QPropertyAnimation(self, b"knob_position")
        self.animation.setDuration(150)
        self.animation.setEasingCurve(QEasingCurve.Type.InOutCubic)
        self.stateChanged.connect(self._start_animation)
    
    def hitButton(self, pos):
        return self.rect().contains(pos)
    
    def _start_animation(self, state):
        self.animation.setStartValue(self._knob_position)
        self.animation.setEndValue(20 if state else 4)
        self.animation.start()

    def get_knob_position(self):
        return self._knob_position

    def set_knob_position(self, pos):
        self._knob_position = pos
        self.update()

    knob_position = Property(float, get_knob_position, set_knob_position)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Track
        track_color = QColor("#2196F3") if self.isChecked() else QColor("#CCCCCC")
        p.setBrush(QBrush(track_color))
        p.setPen(Qt.PenStyle.NoPen)
        p.drawRoundedRect(0, 0, self.width(), self.height(), 10, 10)
        
        # Knob
        p.setBrush(QBrush(QColor("white")))
        p.drawEllipse(int(self._knob_position), 3, 15, 15)

class SegmentedButton(QWidget):
    valueChanged = Signal(int)
    textChanged = Signal(str) # Added for parity with ComboBox

    def __init__(self, options: List[str], parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        
        self.buttons = []
        self._group = QButtonGroup(self)
        self._group.setExclusive(True) # Ensures only one is checked
        
        for i, opt in enumerate(options):
            btn = QPushButton(opt)
            btn.setCheckable(True)
            
            # Styling logic remains the same
            if i == 0:
                btn.setStyleSheet("border-top-right-radius:0;border-bottom-right-radius:0;")
            elif i == len(options) - 1:
                btn.setStyleSheet("border-top-left-radius:0;border-bottom-left-radius:0;")
            else:
                btn.setStyleSheet("border-radius:0;")
            
            layout.addWidget(btn)
            self.buttons.append(btn)
            self._group.addButton(btn, i) # Assign ID as the index

        self._group.idClicked.connect(self._on_id_clicked)
        self.setValue(0)

    def _on_id_clicked(self, idx: int):
        self.valueChanged.emit(idx)
        self.textChanged.emit(self.currentText())

    def value(self) -> int:
        return self._group.checkedId()

    def setValue(self, idx: int):
        btn = self._group.button(idx)
        if btn:
            btn.setChecked(True)
            self.valueChanged.emit(idx)

    def currentText(self) -> str:
        """Mimics QComboBox.currentText()"""
        checked_btn = self._group.checkedButton()
        return checked_btn.text() if checked_btn else ""

    def setCurrentText(self, text: str):
        """Mimics QComboBox.setCurrentText()"""
        for i, btn in enumerate(self.buttons):
            if btn.text() == text:
                self.setValue(i)
                break

# ── ExperimentThread ────────────────────────────────────────────────────────
class ExperimentThread(QThread):
    # inner_idx, param_val, i_run, oc_idx, oc_vals_str, processed, raw
    point_signal = Signal(int, float, int, int, str, object, object)
    finished_signal = Signal(object, object)   # (exp, data)
    error_signal = Signal(str)

    def __init__(self, run_func, instruments, config,
                 save_path=None, folder_number=None):
        super().__init__()
        self.run_func = run_func
        self.instruments = instruments
        self.config = config
        self.save_path = save_path
        self.folder_number = folder_number
        self._stop = False

    def run(self):
        try:
            exp, data = self.run_func(
                instruments=self.instruments,
                config=self.config,
                callback=self._cb,
                stop_check=lambda: self._stop,
                save_path=self.save_path,
                folder_number=self.folder_number,
            )
            self.finished_signal.emit(exp, data)
        except Exception:
            self.error_signal.emit(traceback.format_exc())

    def _cb(self, inner_idx, val, i_run, oc_idx, oc_vals, proc, raw):
        oc_str = "  ".join(f"{k}={v:.4g}" for k, v in oc_vals.items()) \
                 if oc_vals else ""
        self.point_signal.emit(inner_idx, val, i_run, oc_idx, oc_str,
                               proc, raw)

    def request_stop(self):
        self._stop = True


# ── Session Manager ─────────────────────────────────────────────────────────
pg.setConfigOptions(antialias=True)

# QPushButton { background-color: #3d3d3d; border: 1px solid #4d4d4d; border-radius: 4px; padding: 6px 12px; color: #e0e0e0; }
# QPushButton:hover { background-color: #4d4d4d; }
                   
_SS = """
QMainWindow,QWidget{background-color:#1e1e1e;color:#d4d4d4;font-family:'Segoe UI',Arial,sans-serif;}
QGroupBox{border:1px solid #3c3c3c;border-radius:4px;margin-top:8px;
    padding-top:14px;font-weight:bold}
QGroupBox::title{subcontrol-origin:margin;left:10px}
QTabWidget::pane{border:1px solid #3c3c3c;background:#1e1e1e}
QTabBar::tab{background:#2d2d2d; border:1px solid #3c3c3c; color: #a0a0a0;
    padding:4px 10px; margin-right:2px; border-top-left-radius:3px;
    border-top-right-radius:3px}
QTabBar::tab:selected{background:#1e1e1e;border-bottom-color:#1e1e1e}
QPushButton{
    background-color: #484848; color: #ffffff;
    border: 1px solid #3c3c3c; border-radius: 4px;
    padding: 5px 14px
}
QPushButton:hover {
    background-color: #1aa0d9; color: #ffffff;
}
QPushButton:disabled {
    color: #555;
}
QPushButton:checked {
    background-color: #0078d4;
}
QLineEdit, QSpinBox, QComboBox {
    background: #484848;
    color: #ffffff;
    border: 1px solid #3f3f46;
    border-radius: 3px;
    padding: 4px;
    selection-background-color: #264f78;
}
QLineEdit:focus, QSpinBox:focus, QComboBox:focus {
    border: 1px solid #009de0;
}
QLineEdit:disabled, QSpinBox:disabled {
    background-color: #252526;
    color: #6d6d6d;
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
    border: none; background: #737475; max-height: 60px;
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
QSpinBox {
    background-color: #2d2d2d; color: #ffffff;
    border: 1px solid #555555; border-radius: 4px;
    padding-right: 5px; /* Leave space for buttons */
    selection-background-color: #444444;
}
/* The Buttons Container */
QSpinBox::up-button, QSpinBox::down-button {
    background-color: #3d3d3d;
    border-left: 1px solid #555555;
    width: 20px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {
    background-color: #4d4d4d;
}
/* The Arrows (Triangle Hack) */
QSpinBox::up-arrow {
    width: 0px; height: 0px; border-left: 4px solid #3d3d3d;
    border-right: 4px solid #3d3d3d; border-bottom: 5px solid #ffffff;
}
QSpinBox::down-arrow {
    width: 0px; height: 0px; border-left: 4px solid #3d3d3d;
    border-right: 4px solid #3d3d3d; border-top: 5px solid #ffffff;
}
QStatusBar{background:#252526}

QTreeWidget{background:#252526;border:none;color:#d4d4d4}
QTreeWidget::item:alternate{background:#2a2a2a}
QCheckBox {
    color: #ffffff;
    spacing: 8px;
}
QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 1px solid #555555;
    border-radius: 3px;
    background-color: #2d2d30;
}
QCheckBox::indicator:checked {
    background-color: #1397f0;
    border-color: #859199;
}
"""


class SessionManager(QMainWindow):
    _upd = Signal(int, float, int, int, str, object, object)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("NV Session Manager")
        self.setMinimumSize(1100, 720)
        self.setStyleSheet(_SS)

        # Working directory for state files
        self.statefile_directory = Path(STATEFILE_DIRECTORY)
        self.statefile_directory.mkdir(parents=True, exist_ok=True)

        self.sg = None; self.pb = None; self.ao_task = None
        self.instruments_connected = False
        self._pb_ready = False

        self.experiment_thread: Optional[ExperimentThread] = None
        self.is_running = False
        self.last_exp = None  # DiodeExperiment object from last run

        # Per-(outer, run) curve storage
        self._curves: dict[str, pg.PlotDataItem] = {}

        self._loaded_config = None
        self._loaded_config_name = ""
        self._ipy_kernel = None
        self._ipy_widget = None
        self._camera_viewer = None       # In-process CameraViewerV2 window
        self._hci_process = None          # subprocess.Popen for HCImageLive
        self._hci_watcher = None          # QTimer polling HCImageLive exit
        self._last_save_path = None       # Path to last meas_NNN/ directory
        self._last_folder_number = None   # Last folder_number string
        self._deferred_comment = ''       # Comment deferred to next measurement
        self._contrast_vb = None          # Right-axis ViewBox for diode sweep
        self._cam_contrast_vb = None      # Right-axis ViewBox for camera sweep
        self._cam_extra_plots = []        # Camera-specific dynamic plot widgets
        self.last_used_configs = ['esr_config', 'esr_camera_config']      # pos 0 - diode, pos 1 - camera

        self._build_ui()
        self._load_session_state()

        self._upd.connect(self._on_plot, Qt.ConnectionType.QueuedConnection)

    # ── UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        central_widget = QWidget(); self.setCentralWidget(central_widget)
        root = QHBoxLayout(central_widget); root.setContentsMargins(4, 4, 4, 4)

        # Left: tabs
        self._tabs = QTabWidget()
        self._tabs.setMaximumWidth(550); self._tabs.setMinimumWidth(260)
        self._tabs.setFixedWidth(340)
        self._tabs.addTab(self._build_instruments_tab(), "Instruments")
        self._tabs.addTab(self._build_experiment_tab(), "Experiment")
        self._tabs.addTab(self._build_metadata_tab(), "Config")
        self._tabs.addTab(self._build_console_tab(), "Console")
        self._tabs.currentChanged.connect(
            lambda: self._tabs.setFixedWidth(500)if self._tabs.currentIndex() == 3
            else self._tabs.setFixedWidth(340)
            )
        root.addWidget(self._tabs)

        # Right: plots
        plot_wg = QWidget()
        
        plot_wg_sub = QVBoxLayout(plot_wg);  plot_wg_sub.setContentsMargins(0, 0, 0, 0)
        self.sweep_plot = pg.PlotWidget(title="Sweep")
        self.sweep_plot.addLegend(offset=(10, 10))
        self.sweep_plot.setLabel('bottom', 'Parameter')
        self.sweep_plot.setLabel('left', 'Signal')
        self.sweep_plot.showGrid(x=True, y=True, alpha=0.3)
        plot_wg_sub.addWidget(self.sweep_plot, stretch=2)

        self.raw_plot = pg.PlotWidget(title="Raw Trace")
        self.raw_plot.setLabel('bottom', 'Sample')
        self.raw_plot.showGrid(x=True, y=True, alpha=0.3)
        self.raw_curve = self.raw_plot.plot(pen=_pen(3, width=1))
        plot_wg_sub.addWidget(self.raw_plot, stretch=1)
        root.addWidget(plot_wg, stretch=1)

        self.status_bar = QStatusBar(); self.setStatusBar(self.status_bar)
        self.status_label = QLabel("Ready"); self.status_bar.addWidget(self.status_label, 1)

    def _build_instruments_tab(self):
        w = QWidget(); vl = QVBoxLayout(w); vl.setContentsMargins(6,6,6,6)
        lbl = QLabel("Instruments")
        lbl.setStyleSheet("color: #888; font-size: 15px")
        vl.addWidget(lbl)

        mono = "font-family:'Consolas','Courier New',monospace;font-size:12px"
        vl.addSpacing(5)
        self.lbl_sg = QLabel("SG:  ○"); self.lbl_sg.setStyleSheet(mono)
        self.lbl_pb = QLabel("PB:  ○"); self.lbl_pb.setStyleSheet(mono)
        self.lbl_ao = QLabel("AO:  ○"); self.lbl_ao.setStyleSheet(mono)
        self.lbl_cam = QLabel("CAM: ○"); self.lbl_cam.setStyleSheet(mono)
        for index, l in enumerate([self.lbl_sg, self.lbl_pb, self.lbl_ao, self.lbl_cam]):
            vl.addWidget(l)
            vl.addSpacing(5) if index !=3 else None

        # ── Reload modules button (Feature 2) ───────────────────────
        cr = QHBoxLayout()
        cr.addWidget(QLabel("Modules:"))
        # cr.addStretch()
        self.btn_reload_modules = QPushButton("🔄 Reload Modules")
        self.btn_reload_modules.setFixedSize(200, 25)
        self.btn_reload_modules.setToolTip(
            "Re-import PBcontrol, DAQcontrol, SGcontrol,\n"
            "sequencecontrol, connectionConfig, experiment_base,\n"
            "camera_experiment, timeseries_experiment, sweep_utils.\n"
            "Use after editing driver code without restarting.")
        self.btn_reload_modules.setStyleSheet(
            "QPushButton{background:#3d3d1d}QPushButton:hover{background:#4d4d2d}")
        self.btn_reload_modules.clicked.connect(self._reload_modules)
        cr.addWidget(self.btn_reload_modules)
        cr.addStretch()
        vl.addLayout(cr)
        
        vl.addSpacing(10)

        instrs = ('ALL', 'SG', 'PB', 'AO', 'CAM')

        r1 = QHBoxLayout()
        r1.addSpacing(10)
        r2 = QVBoxLayout()
        self.simulate_checkboxes = {}

        for nm in instrs:
            r2.addSpacing(4)
            if nm != 'SG':
                b = QLabel("")
            else:
                b = QCheckBox(f"Simulate {nm}")
            
            r2.addWidget(b)
            r2.addSpacing(4)
            self.simulate_checkboxes[nm] = b
        r2.addStretch()
        r1.addLayout(r2)

        r2 = QVBoxLayout()
        self.init_pushbuttons = {}
        for nm in instrs:
            if nm not in ['AO', ]:
                b = QPushButton(f"Init {nm}"); b.setMaximumWidth(90)
                b.clicked.connect(lambda _, n=nm: self._init(n))
                b.setStyleSheet("QPushButton{background:#2d4a2d}QPushButton:hover{background:#3a5f3a}")
            else:
                r2.addSpacing(6)
                b = QLabel("")
                r2.addSpacing(6)
            
            r2.addWidget(b)
            self.init_pushbuttons[nm] = b
        r2.addStretch()
        r1.addLayout(r2)

        r2 = QVBoxLayout()
        self.eject_pushbuttons = {}
        for nm in instrs:
            b = QPushButton(f"Eject {nm}"); b.setMaximumWidth(90)
            b.clicked.connect(lambda _, n=nm: self._eject(n))
            r2.addWidget(b)
            self.eject_pushbuttons[nm] = b
        r2.addStretch()
        r1.addLayout(r2)
        r1.addSpacing(10)
        vl.addLayout(r1)
        
        self.init_pushbuttons['CAM'].setToolTip(
            "Initialize CameraWorker from loaded config.\n"
            "Reads ROI/exposure from config or camera viewer.")

        # ── Tools ──────────────────

        vl.addWidget(QLabel(""))  # spacer
        lbl = QLabel("Tools")
        lbl.setStyleSheet("color: #888; font-size: 15px")
        vl.addWidget(lbl)
        
        vl.addSpacing(2)

        # ── Camera tools ────────────────────────────────────────────────
        cr = QHBoxLayout()

        cr_sub = QVBoxLayout()
        lbl2 = QLabel("Camera (In-process):")
        lbl2.setStyleSheet("color: #888; font-size: 11px")
        cr_sub.addWidget(lbl2)

        self.btn_launch_cam = QPushButton("📷 Camera Viewer")
        self.btn_launch_cam.setToolTip("Open Camera Viewer v2 (in-process, returns settings)")
        self.btn_launch_cam.clicked.connect(self._launch_camera_viewer)
        cr_sub.addWidget(self.btn_launch_cam)
        cr.addLayout(cr_sub)

        cr_sub = QVBoxLayout()
        lbl2 = QLabel("Camera (Windows App):")
        lbl2.setStyleSheet("color: #888; font-size: 11px")
        cr_sub.addWidget(lbl2)

        self.btn_launch_hci = QPushButton("🎥 HCImageLive")
        self.btn_launch_hci.setToolTip(
            "Launch HCImageLive.exe for high-speed streaming\n"
            "(disconnects camera first, reconnects on close)")
        self.btn_launch_hci.clicked.connect(self._launch_hcimagelive)
        cr_sub.addWidget(self.btn_launch_hci)
        cr.addLayout(cr_sub)

        vl.addLayout(cr)
        vl.addSpacing(5)

        # ── Subprocess — no parameter return ─────────
        lbl = QLabel("Separate process:")
        lbl.setStyleSheet("color: #888; font-size: 11px")
        vl.addWidget(lbl)

        tr = QHBoxLayout()
        tr.addStretch()
        tr_sub = QVBoxLayout()
        self.btn_launch_scope = QPushButton("📊 DAQ Scope")
        self.btn_launch_scope.setToolTip("Launch scope_v2_main.py in a new process")
        self.btn_launch_scope.clicked.connect(lambda: self._launch_tool('scope_v2_main.py'))
        tr_sub.addWidget(self.btn_launch_scope)

        self.btn_launch_pb_gui = QPushButton("🔧 PB GUI")
        self.btn_launch_pb_gui.setToolTip("Launch pyPBLV.py in a new process")
        self.btn_launch_pb_gui.clicked.connect(lambda: self._launch_tool('pyPBLV.py'))
        tr_sub.addWidget(self.btn_launch_pb_gui)
        tr.addLayout(tr_sub)

        tr_sub = QVBoxLayout()
        self.btn_launch_sg_gui = QPushButton("🔧 SG GUI")
        self.btn_launch_sg_gui.setToolTip("Launch SG384gui.py in a new process")
        self.btn_launch_sg_gui.clicked.connect(lambda: self._launch_tool('SG384gui.py'))
        tr_sub.addWidget(self.btn_launch_sg_gui)

        self.btn_launch_verdi_gui = QPushButton("🔦 Verdi GUI")
        self.btn_launch_verdi_gui.setToolTip("Launch verdi_gui.py in a new process")
        self.btn_launch_verdi_gui.clicked.connect(lambda: self._launch_tool('verdi_gui.py'))
        tr_sub.addWidget(self.btn_launch_verdi_gui)
        tr.addLayout(tr_sub)
        tr.addStretch()
        vl.addLayout(tr)

        # ── Session state save/load ─────────────────────────────────
        vl.addStretch()
        self._build_state_controls(vl)

        return w
    
    def _refresh_config_file_list(self, script_file:str=''):
        self.config_files = [file.stem
                             for file in list(Path.glob(Path(EXPERIMENT_FILE_DIRECTORY),
                                                        '*_config.py'))
                            ]
        self.config_files.remove('experiment_config')

        if 'camera' in script_file:
            self.config_files = [file if 'camera' in file else '' for file in self.config_files]
            self.config_files = list(filter(lambda x: x != '', self.config_files))
        elif 'diode' in script_file:
            self.config_files = [file if 'camera' not in file else '' for file in self.config_files]
            self.config_files = list(filter(lambda x: x != '', self.config_files))
    
    def _update_config_file_display(self):
        self.cb_cfg.clear()  # Clear old items
        self.cb_cfg.addItems(self.config_files)

        script_file = self.cb_scr.currentText()
        if 'diode' in script_file:
            self.cb_cfg.setCurrentText(self.last_used_configs[0])
        elif 'camera' in script_file:
            self.cb_cfg.setCurrentText(self.last_used_configs[1])

    def _update_config_file_list(self):
        script_file = self.cb_scr.currentText()

        if 'diode' in script_file:
            self.last_used_configs[1] = self.cb_cfg.currentText()
        elif 'camera' in script_file:
            self.last_used_configs[0] = self.cb_cfg.currentText()
        self._refresh_config_file_list(script_file=script_file)
        self._update_config_file_display()
    
    def _refresh_script_file_list(self):
        self.script_files = [file.stem
                             for file in list(Path.glob(Path(EXPERIMENT_FILE_DIRECTORY),
                                                        'mainControl_*.py'))
                            ]
    
    def _update_script_list_display(self):
        self.cb_scr.clear()  # Clear old items
        self.cb_scr.addItems(self.script_files)

    def _build_experiment_tab(self):
        w = QWidget()
        
        vl = QVBoxLayout(w); vl.setContentsMargins(6,6,6,6)

        vl.addSpacing(4)
        config_grp = QGroupBox("Experiment Config")
        config_grp.setStyleSheet(
            "QGroupBox {"
            "border: 1px solid #3c3c3c; border-radius: 3px;"
            "margin-top: 6px; padding-top: 14px }"
            "QGroupBox::title {"
            "subcontrol-origin: margin; subcontrol-position: top left; padding: 0 4px}")
        
        r1 = QVBoxLayout(config_grp)
        
        r1_sub = QHBoxLayout()
        r1_sub.addWidget(QLabel("Config Module:"))
        self.cb_cfg = QComboBox()
        r1_sub.addWidget(self.cb_cfg)

        self.open_config = QPushButton("📝")
        self.open_config.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; padding: 2px 2px;}
            QPushButton:hover {background-color: #009de0; color: #ffffff;}
        """)
        self.open_config.setFixedSize(25, 28)
        self.open_config.setToolTip("Edit Config File")
        self.open_config.clicked.connect(lambda: self._open_in_editor(self.cb_cfg.currentText()))
        r1_sub.addWidget(self.open_config)

        # Refresh button inline with combo
        self.refresh_config_btn = QPushButton('↻')
        self.refresh_config_btn.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; padding: 2px 2px;}
            QPushButton:hover {background-color: #009de0; color: #ffffff;}
        """)
        self.refresh_config_btn.setFixedSize(25, 28)
        self.refresh_config_btn.setToolTip("Refresh config list")
        self.refresh_config_btn.clicked.connect(self._refresh_config_file_list)
        r1_sub.addWidget(self.refresh_config_btn)
        r1.addLayout(r1_sub)
        r1.addSpacing(5)

        r1_sub = QHBoxLayout()
        self.btn_load_cfg = QPushButton("📋 Load")
        self.btn_load_cfg.clicked.connect(self._load_config)
        r1_sub.addWidget(self.btn_load_cfg)

        self.btn_view_seq = QPushButton("🔍 View")
        self.btn_view_seq.setToolTip("View Sequence (needs PB)")
        # self.btn_view_seq.setStyleSheet(
        #     "QPushButton{background:#2d5f2d}QPushButton:hover{background:#3a7a3a}")
        self.btn_view_seq.clicked.connect(self._view_sequence)
        # self.btn_view_seq.setEnabled(False)
        r1_sub.addWidget(self.btn_view_seq)

        self.btn_play_seq = QPushButton("▶ Play")
        self.btn_play_seq.setToolTip("Play Sequence on PB (needs PB)")
        # self.btn_play_seq.setStyleSheet(
        #     "QPushButton{background:#2d4a5f}QPushButton:hover{background:#3a6a7a}")
        self.btn_play_seq.clicked.connect(self._play_sequence)
        self.btn_play_seq.setEnabled(False)
        r1_sub.addWidget(self.btn_play_seq)

        self.btn_stop_seq = QPushButton("🛑 Stop")
        self.btn_stop_seq.setToolTip("Stop Sequence on PB (needs PB)")
        # self.btn_stop_seq.setStyleSheet(
        #     "QPushButton{background:#2d4a5f}QPushButton:hover{background:#3a6a7a}")
        self.btn_stop_seq.clicked.connect(self._stop_sequence)
        self.btn_stop_seq.setEnabled(False)
        r1_sub.addWidget(self.btn_stop_seq)
        r1.addLayout(r1_sub)

        vl.addWidget(config_grp)

        vl.addSpacing(10)
        control_grp = QGroupBox("Experiment Config")
        control_grp.setStyleSheet("""
            QGroupBox {border: 1px solid #3c3c3c; border-radius: 3px;
            margin-top: 6px; padding-top: 14px}
            QGroupBox::title {subcontrol-origin: margin; subcontrol-position: top left; padding: 0 4px}""")
        
        r1 = QVBoxLayout(control_grp)
        r1_sub = QHBoxLayout()
        r1_sub.addWidget(QLabel("Script File:"))
        self.cb_scr = QComboBox()
        self.cb_scr.currentTextChanged.connect(self._update_config_file_list)
        self._refresh_script_file_list()
        self._update_script_list_display(); self._update_config_file_list()
        r1_sub.addWidget(self.cb_scr)

        self.open_script = QPushButton("🖋")
        self.open_script.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; padding: 2px 2px;}
            QPushButton:hover {background-color: #009de0; color: #ffffff;}
        """)
        self.open_script.setFixedSize(25, 28)
        self.open_script.setToolTip("Edit Config File")
        self.open_script.clicked.connect(lambda: self._open_in_editor(self.cb_scr.currentText()))
        r1_sub.addWidget(self.open_script)

        # Refresh button inline with combo
        self.refresh_script_btn = QPushButton('↻')
        self.refresh_script_btn.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; padding: 2px 2px;}
            QPushButton:hover {background-color: #009de0; color: #ffffff;}
        """)
        self.refresh_script_btn.setFixedSize(25, 28)
        self.refresh_script_btn.setToolTip("Refresh config list")
        self.refresh_script_btn.clicked.connect(self._refresh_script_file_list)
        r1_sub.addWidget(self.refresh_script_btn)
        r1.addLayout(r1_sub)
        r1.addSpacing(5)

        r1_sub = QHBoxLayout()
        r1_sub.addWidget(QLabel("Mode:"))
        # self.cb_mode = QComboBox()
        # self.cb_mode.addItems(['Sweep', 'Timeseries'])
        r1_sub.addStretch()
        self.cb_mode = SegmentedButton(['Sweep', 'Timeseries'])
        self.cb_mode.valueChanged.connect(self.toggle_ts_duration)
        r1_sub.addWidget(self.cb_mode)
        r1.addLayout(r1_sub)
        
        # Timeseries duration (0 = manual stop)
        r1_sub = QHBoxLayout()
        self.ts_label = QLabel("TS dur (s):")        
        r1_sub.addWidget(self.ts_label)
        r1_sub.addStretch()

        self.txt_ts_duration = QLineEdit("0")
        self.txt_ts_duration.setMaximumWidth(60)
        self.txt_ts_duration.setToolTip(
            "Timeseries auto-stop duration in seconds.\n"
            "0 = manual stop only.\n"
            "Press Enter to apply mid-acquisition.")
        self.txt_ts_duration.returnPressed.connect(self._apply_ts_duration)
        r1_sub.addWidget(self.txt_ts_duration)
        self.toggle_ts_duration(0)
        
        r1.addLayout(r1_sub)
        # r1.addStretch(50)
        vl.addWidget(control_grp)

        vl.addSpacing(10)
        plot_grp = QGroupBox("Experiment Plot")
        plot_grp.setStyleSheet("""
            QGroupBox {
            border: 1px solid #3c3c3c; border-radius: 3px;"
            margin-top: 6px; padding-top: 14px }
            QGroupBox::title {
            subcontrol-origin: margin; subcontrol-position: top left; padding: 0 4px}""")
        r1 = QVBoxLayout(plot_grp)

        r1_sub = QHBoxLayout()
        r1_sub.addWidget(QLabel("Real-time plot"))
        self.chk_rt = ToggleSwitch()
        self.chk_rt.setChecked(True)
        r1_sub.addWidget(self.chk_rt)
        # r1_sub.addSpacing(10)
        r1_sub.addStretch()
        
        # Contrast operation selector
        r1_sub.addWidget(QLabel("Contrast:"))
        self.cb_contrast = QComboBox()
        self.cb_contrast.addItems(['s/r', 'r/s', '(s-r)/(s+r)',
                                   's-r', 'r-s'])
        self.cb_contrast.setToolTip(
            "Contrast operation for real-time display.\n"
            "s = signal mean, r = reference mean.")
        r1_sub.addWidget(self.cb_contrast)

        r1.addLayout(r1_sub)
        vl.addWidget(plot_grp)
        vl.addSpacing(10)

        # ── Autosave section (HCImageLive-style) ────────────────────        
        save_grp = QGroupBox("Save")
        save_grp.setStyleSheet("""
            QGroupBox{border:1px solid #3c3c3c;border-radius:3px;
            margin-top:6px;padding-top:14px}
            QGroupBox::title{subcontrol-origin:margin;
            subcontrol-position:top left;padding:0 4px}""")
        sg_l = QVBoxLayout(save_grp)

        save_text = QHBoxLayout()
        r_prefix = QVBoxLayout()
        r_prefix.addWidget(QLabel("Prefix"))
        self.txt_prefix = QLineEdit("data")
        self.txt_prefix.setToolTip("Base prefix for saved files")
        r_prefix.addWidget(self.txt_prefix)
        save_text.addLayout(r_prefix)

        r_suffix = QVBoxLayout()
        r_suffix.addWidget(QLabel("Suffix"))
        self.txt_suffix = QLineEdit("")
        self.txt_suffix.setToolTip("Optional suffix appended after number")
        r_suffix.addWidget(self.txt_suffix)
        save_text.addLayout(r_suffix)

        sg_l.addLayout(save_text)

        r_num = QHBoxLayout()
        r_num.addWidget(QLabel("Start #:"))
        self.spin_autosave_num = QSpinBox()
        self.spin_autosave_num.setRange(1, 99999)
        self.spin_autosave_num.setValue(1)
        self.spin_autosave_num.setToolTip("Auto-incrementing measurement number")
        r_num.addWidget(self.spin_autosave_num)
        self._lbl_collision = QLabel("✔")
        self._lbl_collision.setStyleSheet("color:#6a9955;font-weight:bold")
        self._lbl_collision.setMaximumWidth(20)
        r_num.addWidget(self._lbl_collision)
        r_num.addStretch()
        
        self.chk_timestamp = QCheckBox("Timestamp suffix")
        self.chk_timestamp.setToolTip(
            "Append _HHMMSS timestamp after suffix")
        r_num.addWidget(self.chk_timestamp)
        sg_l.addLayout(r_num)
        sg_l.addSpacing(5)

        # Wire collision check to all autosave widget changes
        self.txt_prefix.textChanged.connect(
            lambda: self._check_folder_collision())
        self.spin_autosave_num.valueChanged.connect(
            lambda: self._check_folder_collision())
        self.txt_suffix.textChanged.connect(
            lambda: self._check_folder_collision())        

        # ── Comments section (Phase 6) ──────────────────────────────
        sg_l.addWidget(QLabel("Comments:"))
        cmt_l = QHBoxLayout()
        self.txt_comment = QLineEdit()
        self.txt_comment.setPlaceholderText("Appended to params YAML")
        # self.txt_comment.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.txt_comment.setToolTip(
            "Comment saved into the param file.\n"
            "Cleared after each measurement completes.\n"
            "Use ◀/▶ to annotate previous/next measurement.")
        cmt_l.addWidget(self.txt_comment)

        self.btn_comment_prev = QPushButton("◀")
        self.btn_comment_prev.setFixedSize(28, 28)
        self.btn_comment_prev.setStyleSheet("QPushButton{padding: 2px 2px; font-size: 20px}")
        self.btn_comment_prev.setToolTip(
            "Write comment into the PREVIOUS measurement's param file")
        self.btn_comment_prev.clicked.connect(self._comment_prev)
        cmt_l.addWidget(self.btn_comment_prev)

        self.btn_comment_next = QPushButton("▶")
        self.btn_comment_next.setFixedSize(28, 28)
        self.btn_comment_next.setStyleSheet("QPushButton{padding: 2px 2px; font-size: 20px}")
        self.btn_comment_next.setToolTip(
            "Defer comment — write into the NEXT measurement's param file")
        self.btn_comment_next.clicked.connect(self._comment_next)
        cmt_l.addWidget(self.btn_comment_next)
        sg_l.addLayout(cmt_l)
        sg_l.addSpacing(10)

        chk_save = QHBoxLayout()
        chk_save.addWidget(QLabel("Autosave"))
        self.chk_save = ToggleSwitch(); self.chk_save.setChecked(True)
        chk_save.addWidget(self.chk_save)
        chk_save.addStretch()

        self.save_data_manually = QPushButton("Save Data")
        self.save_data_manually.setFixedSize(100, 45)
        self.save_data_manually.setToolTip("Save Data when Autosave is OFF.\nBut how to save? Autosave saves after each run!")
        # self.save_data_manually.clicked.connect(self._launch_camera_viewer)   # TODO: wire this
        chk_save.addWidget(self.save_data_manually)
        sg_l.addLayout(chk_save)

        vl.addWidget(save_grp)
        vl.addStretch()

        r2 = QHBoxLayout()
        self.btn_run = QPushButton("▶ RUN")
        self.btn_run.setStyleSheet("""
            QPushButton{background:#0e639c;font-weight:bold;
            font-size:14px;padding:8px}
            QPushButton:hover{background:#1177bb;}
            QPushButton:disabled{background:#333;color:#666;}""")
        self.btn_run.clicked.connect(self._run_experiment)
        self.btn_run.setEnabled(False)
        r2.addWidget(self.btn_run)

        self.btn_stop = QPushButton("🛑 STOP")
        self.btn_stop.setStyleSheet("""
            QPushButton{background:#6c1d1d;font-weight:bold;
            font-size:14px;padding:8px}
            QPushButton:hover{background:#8b2222;}""")
        self.btn_stop.clicked.connect(self._stop_experiment)
        self.btn_stop.setEnabled(False)
        r2.addWidget(self.btn_stop)

        self.btn_abort = QPushButton("⚠ ABORT")
        self.btn_abort.setToolTip(
            "Force-abort a stalled experiment.\n"
            "Kills the worker thread, tears down DAQ/PB tasks,\n"
            "and returns session to idle — without closing the UI.")
        self.btn_abort.setStyleSheet("""
            QPushButton{background:#8b0000;font-weight:bold;color:#ffc0c0;
            font-size:14px;padding:8px}
            QPushButton:hover{background:#a00000}""")
        self.btn_abort.clicked.connect(self._abort_experiment)
        self.btn_abort.setEnabled(False)
        r2.addWidget(self.btn_abort)

        vl.addLayout(r2)
        return w

    def toggle_ts_duration(self, index: int):
        """
        Logic: 'Timeseries' is index 1.
        We show the widgets if index == 1, hide otherwise.
        """
        is_timeseries = (index == 1)
        
        # You must set visibility for each widget in the sub-layout
        self.ts_label.setVisible(is_timeseries)
        self.txt_ts_duration.setVisible(is_timeseries)
    
    def on_script_file_changed(self):
        script_file = self.cb_scr.currentText()
        # if 'camera' in script_file:
        self._refresh_config_file_list(script_file=script_file)
        # elif 'diode' in script_file:
        #     self._refresh_config_file_list(filter='diode')
        
    def _open_in_editor(self, file:str='', editor:str='vscode'):
        if editor == 'vscode':
            editor_path = r"%USERPROFILE%\AppData\Local\Programs\Microsoft VS Code\Code.exe"
        elif editor == 'spyder':
            editor_path = r"%USERPROFILE%\AppData\Local\spyder-6\envs\spyder-runtime\Scripts\spyder.exe"
        else:
            editor_path = 'notepad.exe'
        clean_path = os.path.normpath((Path(EXPERIMENT_FILE_DIRECTORY) / file).with_suffix('.py'))
        
        try:
            if editor == 'vscode':
                subprocess.Popen([editor_path, clean_path], shell=True)
            elif editor == 'spyder':
                subprocess.Popen([editor_path, clean_path], creationflags=subprocess.DETACHED_PROCESS,
                                 close_fds=True)
            else:
                subprocess.Popen(['notepad.exe', clean_path])
        except Exception as e:
            print(f"Failed to launch: {e}")
    
    def _build_metadata_tab(self):
        w = QWidget(); vl = QVBoxLayout(w); vl.setContentsMargins(4,4,4,4)
        self._meta_tree = QTreeWidget()
        self._meta_tree.setColumnCount(2)
        self._meta_tree.setHeaderLabels(["Parameter","Value"])
        self._meta_tree.setAlternatingRowColors(True)
        self._meta_tree.setRootIsDecorated(True)
        self._meta_tree.header().setSectionResizeMode(0,QHeaderView.ResizeMode.ResizeToContents)
        self._meta_tree.header().setSectionResizeMode(1,QHeaderView.ResizeMode.Stretch)
        vl.addWidget(self._meta_tree)
        return w

    def _update_metadata_tree(self, config):
        self._meta_tree.clear()
        try:
            _populate_tree(self._meta_tree.invisibleRootItem(),
                           config.to_dict(include_values=False))
        except Exception as e:
            QTreeWidgetItem(self._meta_tree.invisibleRootItem(),
                            ["Error", str(e)])

    # ── IPython CONSOLE TAB ───────────────────────────────────────

    def _build_console_tab(self):
        w = QWidget(); vl = QVBoxLayout(w); vl.setContentsMargins(2, 2, 2, 2)
        self._ipy_widget = None
        self._ipy_kernel = None

        if not _HAS_QTCONSOLE:
            lbl = QLabel("qtconsole not installed.\n"
                         "pip install qtconsole ipykernel")
            lbl.setStyleSheet("color:#888;padding:20px")
            lbl.setAlignment(Qt.AlignCenter)
            vl.addWidget(lbl)
            return w

        try:
            # Kernel (in-process — shares GIL + memory with session)
            km = QtInProcessKernelManager()
            km.start_kernel()
            kc = km.client()
            kc.start_channels()
            self._ipy_kernel = km

            # Console widget
            cw = RichIPythonWidget()
            cw.kernel_manager = km
            cw.kernel_client = kc
            cw.syntax_style = 'solarized-dark'  # Configure pygments style after widget is created # Dark theme
            cw.style_sheet = "" # Clear default style first

            # MANUAL COLOR CONFIGURATION - This should work!
            # Set ANSI colors manually (this is what controls syntax highlighting)
            cw.ansi_codes = True  # Enable ANSI color codes
            cw.kind = 'rich'  # Use rich text formatting

            # Set individual colors for syntax elements
            # These correspond to pygments token types
            cw._ansi_color_names = [
                'black', 'darkred', 'darkgreen', 'brown',
                'darkblue', 'darkviolet', 'steelblue', 'grey',
                'lightgrey', 'red', 'green', 'yellow',
                'blue', 'violet', 'lightblue', 'white'
            ]
            
            cw.style_sheet = ("""
                QPlainTextEdit, QTextEdit {
                    background-color: #002b36;   /* Solarized Base03 (Background) */
                    color: #839496;              /* Solarized Base0 (Foreground) */
                    selection-background-color: #073642; /* Solarized Base02 */
                    selection-color: #93a1a1;    /* Solarized Base1 */
                    font-family: 'Consolas', 'Monaco', monospace;
                    font-size: 11pt;
                    border: none;
                }
                .in-prompt { color: #00ff00; font-weight: normal; }
                .out-prompt { color: #00ff00; font-weight: normal; }
                .error { color: #dc322f; }
            """)
            cw.in_prompt_style = "color: #00ff00; font-weight: bold;"
            cw.out_prompt_style = "color: #dc322f; font-weight: bold;"

            # Manually set color mappings for different text types
            # This is the key part that makes colors work!
            cw._ansi_foreground_colors = {
                0: QColor("#073642"),  # Base02 (Black)
                1: QColor("#dc322f"),  # Red
                2: QColor("#859900"),  # Green
                3: QColor("#b58900"),  # Yellow
                4: QColor("#268bd2"),  # Blue
                5: QColor("#d33682"),  # Magenta
                6: QColor("#2aa198"),  # Cyan
                7: QColor("#eee8d5"),  # Base2 (White)
                8: QColor("#002b36"),  # Base03 (Bright Black)
                9: QColor("#cb4b16"),  # Orange (Bright Red)
                10: QColor("#586e75"), # Base01 (Bright Green)
                11: QColor("#657b83"), # Base00 (Bright Yellow)
                12: QColor("#839496"), # Base0 (Bright Blue)
                13: QColor("#6c71c4"), # Violet (Bright Magenta)
                14: QColor("#93a1a1"), # Base1 (Bright Cyan)
                15: QColor("#fdf6e3"), # Base3 (Bright White)
            }
            cw._ansi_background_colors = {
                0: QColor(0, 0, 0),
                1: QColor(205, 0, 0),
                2: QColor(0, 205, 0),
                3: QColor(205, 205, 0),
                4: QColor(0, 0, 238),
                5: QColor(205, 0, 205),
                6: QColor(0, 205, 205),
                7: QColor(229, 229, 229),
            }

            self._ipy_widget = cw
            vl.addWidget(cw)

            # Push initial namespace (instruments not yet connected)
            self._push_namespace(banner=True)

        except Exception as e:
            lbl = QLabel(f"Console init failed:\n{e}")
            lbl.setStyleSheet("color:#f44;padding:20px")
            lbl.setWordWrap(True)
            vl.addWidget(lbl)

        return w

    def _push_namespace(self, banner: bool = False):
        """Push session state into the IPython kernel namespace.

        Called:
          - Once at console creation (with banner=True)
          - After each sweep experiment completes (_on_done)
          - After each timeseries stops (_ts_stop)
          - After instrument init

        The user always sees the latest state via the pushed references.
        """
        if self._ipy_kernel is None:
            return

        ns = {
            'session': self,
            'np': np,
            'sg': self.sg,
            'pb': self.pb,
            'ao_task': self.ao_task,
            'cfg': self._loaded_config,
            'exp': self.last_exp,
            'data': (self.last_exp.data_array
                     if self.last_exp is not None else None),
        }

        # Timeseries experiment if available
        if hasattr(self, '_ts_exp') and self._ts_exp is not None:
            ns['ts_exp'] = self._ts_exp

        # Camera viewer if available
        if self._camera_viewer is not None:
            ns['cam_viewer'] = self._camera_viewer
            ns['cam_settings'] = self._get_camera_viewer_settings()

        self._ipy_kernel.kernel.shell.push(ns)

        if banner:
            self._ipy_widget.execute(
                'import numpy as np, matplotlib.pyplot as plt, matplotlib\n'
                'matplotlib.use("QtAgg")    # Force the backend to Qt\n'
                # If you want to switch back to inline later:
                # matplotlib.use('module://matplotlib_inline.backend_inline')
                'print("\\033[96m" + "─"*44 + "\\033[0m")\n'
                'print("\\033[96m  NV Session Console\\033[0m")\n'
                'print("\\033[96m" + "─"*44 + "\\033[0m")\n'
                'print("\\033[93mNamespace:\\033[0m")\n'
                'print("  session  — SessionManager (this window)")\n'
                'print("  sg, pb, ao_task — instrument handles")\n'
                'print("  cfg      — loaded ExperimentConfig")\n'
                'print("  exp      — last DiodeExperiment")\n'
                'print("  data     — last exp.data_array")\n'
                'print("  ts_exp   — last TimeseriesExperiment")\n'
                'print("  cam_viewer — Camera Viewer window (if open)")\n'
                'print("  cam_settings — dict from cam_viewer (exposure, ROI…)")\n'
                'print("  np       — numpy")\n'
                'print()\n'
                'print("\\033[90mNamespace refreshed after each run.\\033[0m")\n'
                'print("\\033[90mAvoid heavy computation during acquisition.\\033[0m")\n',
                hidden=True,
            )

    # ── INSTRUMENTS ────────────────────────────────────────

    def _init_instruments(self):
        try:
            self.status_label.setText("Initializing…")
            QApplication.processEvents()

            self._init_sg(simulate=self.simulate_checkboxes['SG'].isChecked())

            if not self._pb_ready:
                self._do_init_pb()
            
            self.ao_task = AnalogOutputTask(
                dev="P6363", channels=[0, 1, 2], coil='small_confocal')
            self.lbl_ao.setText("AO:  ● ready")

            self.instruments_connected = True
            self.btn_run.setEnabled(True)
            self._update_seq_btns()
            self._push_namespace()
            self.status_label.setText("✔ Instruments initialized")
        except Exception as e:
            self.status_label.setText(f"❌ {e}")
            print(traceback.format_exc())

    def _init_sg(self, simulate=False):
        if simulate:
            self.sg = SignalGenerator_sim()
            self.lbl_sg.setText("SG:  ● simulated")
        else:
            self.sg = SignalGenerator()
            self.sg.query_all()
            self.lbl_sg.setText(f"SG:  ● {self.sg.freq/1e9:.4f} GHz")
    
    def _init_pb(self):
        try:
            self.status_label.setText("Initializing PB…"); QApplication.processEvents()
            self._do_init_pb()
            self._update_seq_btns()
            self.status_label.setText("✔ PB ready — View/Play enabled")
        except Exception as e:
            self.status_label.setText(f"❌ PB: {e}"); print(traceback.format_exc())

    def _do_init_pb(self):
        base = {'pb':{'clk_cyc':1e3/concfg.PBclk},'scan':{},'seq':{},'mw':{}}
        self.pb = PulseBlaster(base); self.pb.configure()
        self.lbl_pb.setText("PB:  ● configured")
        self._pb_ready = True

    def _init(self, name: str):
        """Init a single instrument."""
        try:
            if name == 'SG':
                self._init_sg()
            
            elif name == 'PB':
                self._init_pb()
            
            elif name == 'ALL':
                self._init_instruments()
            
            self.status_label.setText(f"✔ {name} Initialized") if name != 'ALL' else None
        except Exception as e:
            self.status_label.setText(f"❌ Init {name} Failed: {e}")
    
    def _update_seq_btns(self):
        # self.btn_view_seq.setEnabled(self._pb_ready)
        self.btn_play_seq.setEnabled(self._pb_ready)

    def _set_running_state(self, running: bool):
        """Toggle UI buttons for running/idle state."""
        self.is_running = running
        self.btn_run.setEnabled(not running)
        self.btn_stop.setEnabled(running)
        self.btn_abort.setEnabled(running)
        self.init_pushbuttons['ALL'].setEnabled(not running)
        if not running:
            self._update_seq_btns()
        else:
            # self.btn_view_seq.setEnabled(False)
            self.btn_play_seq.setEnabled(False)

    def _abort_experiment(self):
        """Force-abort a stalled experiment without closing the UI.

        This is the nuclear option — it terminates the worker thread,
        tears down DAQ/PB tasks, and resets the session to idle.
        Use when STOP doesn't return (e.g. DAQ timeout, PB hang).
        """
        r = QMessageBox.question(
            self, "Abort",
            "Force-abort the running experiment?\n\n"
            "This will kill the acquisition thread and tear down\n"
            "DAQ tasks. Instruments stay connected.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if r == QMessageBox.StandardButton.No:
            return

        self.status_label.setText("⚠ Aborting…")
        QApplication.processEvents()

        # ── Stop timeseries display timers ──
        for attr in ('_ts_display_timer', '_cam_ts_timer'):
            t = getattr(self, attr, None)
            if t is not None:
                t.stop()

        # ── Kill experiment thread ──
        if self.experiment_thread is not None:
            self.experiment_thread.request_stop()
            if not self.experiment_thread.wait(2000):
                self.experiment_thread.terminate()
                self.experiment_thread.wait(1000)
            self.experiment_thread = None

        # ── Tear down timeseries experiments ──
        for attr in ('_ts_exp', '_cam_ts_exp'):
            exp = getattr(self, attr, None)
            if exp is not None:
                try:
                    exp.stop()
                    exp.teardown()
                except Exception as e:
                    print(f"  abort teardown {attr}: {e}")
                setattr(self, attr, None)

        # ── Stop PB (it may be looping) ──
        if self.pb:
            try:
                self.pb.stop_sequence()
            except Exception:
                pass

        self._set_running_state(False)
        if getattr(self, '_is_camera_run', False):
            self._is_camera_run = False
        
        self.status_label.setText("⚠ Aborted — session idle, instruments connected")
        self._push_namespace()

    def _eject(self, name: str):
        """Eject a single instrument.  Object set to None; session stays."""
        try:
            if name == 'SG' and self.sg:
                if hasattr(self.sg, 'uninit'): self.sg.uninit()
                self.sg = None
                self.lbl_sg.setText("SG:  ○ ejected")
            
            elif name == 'PB' and self.pb:
                try: self.pb.stop_sequence()
                except: pass
                self.pb.closePB(); self.pb = None
                self.lbl_pb.setText("PB:  ○ ejected")
                self._pb_ready = False
                self._update_seq_btns()
            
            elif name == 'AO' and self.ao_task:
                self.ao_task.set_outputs_to_constant([0, 0, 0])
                time.sleep(0.2); self.ao_task.__del__()
                self.ao_task = None
                self.lbl_ao.setText("AO:  ○ ejected")
            
            self.status_label.setText(f"✔ {name} ejected")
        except Exception as e:
            self.status_label.setText(f"❌ Eject {name}: {e}")

    def _init_camera(self):
        """Initialize CameraWorker standalone (without running an experiment)."""
        cfg = self._ensure_config()
        if cfg is None:
            self.status_label.setText("❌ Load a config first"); return
        if not self._is_camera_config(cfg):
            self.status_label.setText("❌ Config has no camera block"); return
        try:
            cw = self._get_or_create_camera_worker(cfg)
            if cw is not None:
                self.status_label.setText(f"✔ Camera: {cw.device_title}")
                self._push_namespace()
        except Exception as e:
            self.status_label.setText(f"❌ Camera init: {e}")
            print(traceback.format_exc())

    def _eject_camera(self):
        """Disconnect and release camera resources."""
        try:
            # Close viewer if open
            if self._camera_viewer is not None:
                try:
                    self._camera_viewer._disconnect_camera()
                    self._camera_viewer.close()
                except Exception:
                    pass
                self._camera_viewer = None

            # Uninit camera worker
            if hasattr(self, '_camera_worker') and self._camera_worker is not None:
                try:
                    self._camera_worker.uninit_cam()
                except Exception:
                    pass
                self._camera_worker = None

            self.lbl_cam.setText("Cam: ○ ejected")
            self.status_label.setText("✔ Camera ejected")
        except Exception as e:
            self.status_label.setText(f"❌ Camera eject: {e}")

    def _reload_modules(self):
        """Re-import instrument and experiment modules in-place.

        Use after editing driver source code (PBcontrol, sequencecontrol,
        etc.) without restarting the session manager. Does NOT re-init
        any hardware — just refreshes the Python module objects so that
        subsequent experiments use the updated code.
        """
        modules_to_reload = [
            'connectionConfig',
            'spinapi',
            'sequencecontrol',
            'PBcontrol',
            'SGcontrol',
            'DAQcontrol',
            'sweep_utils',
            'experiment_base',
            'camera_experiment',
            'timeseries_experiment',
            'Camcontrol',
        ]
        reloaded = []
        failed = []
        for name in modules_to_reload:
            if name in sys.modules:
                try:
                    importlib.reload(sys.modules[name])
                    reloaded.append(name)
                except Exception as e:
                    failed.append(f"{name}: {e}")

        msg = f"✔ Reloaded {len(reloaded)} modules"
        if failed:
            msg += f"  ⚠ {len(failed)} failed"
            for f in failed:
                print(f"  ❌ {f}")
        self.status_label.setText(msg)
        print(f"  Reloaded: {', '.join(reloaded)}")

    # ── TOOL LAUNCHERS (subprocess — fire and forget) ──────────────

    def _launch_tool(self, script_name: str):
        """Launch a tool script in a separate process.

        The child process is fully independent — it manages its own
        DAQ tasks, PB connection, etc.  No parameter return needed
        for DAQ Scope or PB GUI (they don't modify session state).

        For camera viewer (future), use in-process launch instead
        so the session can read back exposure/ROI on close.
        """
        script = Path(script_name)
        if not script.exists():
            # Try common locations
            for candidate in [Path('.') / script_name,
                              Path('scope_v2') / script_name,
                              Path('..') / script_name]:
                if candidate.exists():
                    script = candidate
                    break
        try:
            proc = subprocess.Popen(
                [sys.executable, str(script)],
                cwd=str(script.parent) if script.parent != Path('.') else None,
            )
            self.status_label.setText(f"✔ Launched {script_name} (PID {proc.pid})")
        except Exception as e:
            self.status_label.setText(f"❌ Launch {script_name}: {e}")

    # ════════════════════════════════════════════════════════════════
    #  CAMERA VIEWER (in-process — can read back settings)
    # ════════════════════════════════════════════════════════════════

    def _launch_camera_viewer(self):
        """Open Camera Viewer v2 as an in-process window.

        In-process means:
        - Shares the Qt event loop and GIL with session
        - Session can read back exposure/ROI when viewer closes
        - Camera object is owned by viewer, not session
        - Viewer window is independent (not docked)
        """
        # If already open, just raise it
        if self._camera_viewer is not None:
            try:
                if self._camera_viewer.isVisible():
                    self._camera_viewer.raise_()
                    self._camera_viewer.activateWindow()
                    self.status_label.setText("Camera Viewer already open")
                    return
            except RuntimeError:
                # C++ object deleted
                self._camera_viewer = None

        try:
            from cam_v2_main import CameraViewerV2
            self._camera_viewer = CameraViewerV2()
            self._camera_viewer.setWindowTitle("Camera Viewer v2 — Session")
            self._camera_viewer.setAttribute(
                Qt.WidgetAttribute.WA_DeleteOnClose, False)
            self._camera_viewer.show()
            self.lbl_cam.setText("Cam: ● viewer open")
            self.status_label.setText("✔ Camera Viewer opened")

            # Push to console namespace
            if self._ipy_kernel is not None:
                self._ipy_kernel.kernel.shell.push({
                    'cam_viewer': self._camera_viewer})

        except Exception as e:
            self.status_label.setText(f"❌ Camera Viewer: {e}")
            print(traceback.format_exc())

    def _get_camera_viewer_settings(self) -> dict:
        """Read back camera settings from viewer (if open).

        Returns dict with exposure_ms, roi, camera_type etc.
        Useful for setting up experiment after visual inspection.
        """
        if self._camera_viewer is None:
            return {}
        try:
            v = self._camera_viewer
            settings = {
                'exposure_ms': v.camera_config.exposure_ms,
                'gain': v.camera_config.gain,
                'binning': v.camera_config.binning,
                'camera_type': v.camera_config.camera_type.value,
                'roi': v.camera_config.roi,
                'trigger_mode': v.camera_config.trigger_mode.value,
            }
            # Read actual exposure from camera if connected
            if v.camera is not None and v.camera.is_open:
                settings['actual_exposure_ms'] = v.camera.get_exposure()
                settings['actual_roi'] = v.camera.get_roi()
            return settings
        except Exception:
            return {}

    # ── HCIMAGELIVE LAUNCHER (subprocess — high-speed streaming) ─────────

    # Common install paths for HCImageLive
    _HCIMAGELIVE_PATHS = [
        r"C:\Program Files\HCImageLive\HCImageLive.exe"
    ]

    def _find_hcimagelive(self) -> Optional[str]:
        """Find HCImageLive.exe on the system."""
        import shutil
        for p in self._HCIMAGELIVE_PATHS:
            if Path(p).exists():
                return p
        found = shutil.which("HCImageLive")
        return found

    def _launch_hcimagelive(self):
        """Launch HCImageLive.exe for high-speed streaming.

        Workflow:
        1. Disconnect camera from viewer (DCAM allows only one owner)
        2. Launch HCImageLive as subprocess
        3. Poll for exit — reconnect camera viewer on close

        Camera settings (exposure, ROI) set via DCAM properties persist
        on the hardware, so HCImageLive inherits them.
        """
        # Check if already running
        if self._hci_process is not None and self._hci_process.poll() is None:
            self.status_label.setText("HCImageLive already running")
            return

        # Find executable
        exe = self._find_hcimagelive()
        if exe is None:
            self.status_label.setText("❌ HCImageLive.exe not found")
            QMessageBox.warning(
                self, "HCImageLive",
                "Could not find HCImageLive.exe.\n\n"
                "Searched:\n" +
                "\n".join(f"  • {p}" for p in self._HCIMAGELIVE_PATHS) +
                "\n\nInstall from Hamamatsu or add to PATH.")
            return

        # Disconnect camera from viewer (DCAM exclusive access)
        if self._camera_viewer is not None:
            try:
                self._camera_viewer._disconnect_camera()
                self.status_label.setText("Camera released for HCImageLive...")
                QApplication.processEvents()
                time.sleep(0.3)  # Brief pause for DCAM cleanup
            except Exception as e:
                print(f"Camera disconnect warning: {e}")

        # Launch
        try:
            self._hci_process = subprocess.Popen([exe])
            self.lbl_cam.setText("Cam: ● HCImageLive")
            self.status_label.setText(f"✔ HCImageLive launched (PID {self._hci_process.pid})")
            # self.btn_launch_hci.setStyleSheet("QPushButton{background:#2d4a2d}QPushButton:hover{background:#3a5f3a}")
            self.btn_launch_hci.setText("▶ HCI running...")
            self.btn_launch_hci.setEnabled(False)

            # Start polling for exit
            self._hci_watcher = QTimer(self)
            self._hci_watcher.timeout.connect(self._poll_hcimagelive)
            self._hci_watcher.start(1000)  # Check every second

        except Exception as e:
            self.status_label.setText(f"❌ HCImageLive: {e}")
            print(traceback.format_exc())

    def _poll_hcimagelive(self):
        """Poll HCImageLive process — reconnect camera when it exits."""
        if self._hci_process is None or self._hci_process.poll() is not None:
            # Process has exited
            if self._hci_watcher is not None:
                self._hci_watcher.stop()
                self._hci_watcher = None

            self._hci_process = None
            self.btn_launch_hci.setText("🎥 HCImageLive")
            self.btn_launch_hci.setEnabled(True)
            # self.btn_launch_hci.setStyleSheet("QPushButton{background:#2d2d2d}")
            self.lbl_cam.setText("Cam: ○ HCI closed")
            self.status_label.setText("HCImageLive closed — camera available")

            # Offer to reconnect camera in viewer
            if self._camera_viewer is not None and self._camera_viewer.isVisible():
                self.status_label.setText(
                    "HCImageLive closed — reconnect via Camera menu in viewer")

    # ── CONTRAST DISPLAY HELPER ────────────────────────────────────

    def _compute_display_contrast(self, signal, reference):
        """Compute contrast using the currently selected operation.

        Reads self.cb_contrast and applies to signal/reference arrays.
        Used by both diode and camera sweep plot callbacks.
        """
        op = self.cb_contrast.currentText()
        s = np.asarray(signal, dtype=float)
        r = np.asarray(reference, dtype=float)
        ops = {
            's/r':           lambda s, r: np.where(r != 0, s / r, np.nan),
            'r/s':           lambda s, r: np.where(s != 0, r / s, np.nan),
            '(s-r)/(s+r)':   lambda s, r: np.where(
                (s + r) != 0, (s - r) / (s + r), np.nan),
            's-r':           lambda s, r: s - r,
            'r-s':           lambda s, r: r - s,
        }
        fn = ops.get(op, ops['s/r'])
        return fn(s, r)

    # ── AUTOSAVE — HCImageLive-style naming ────────────────────────────────

    def _make_folder_number(self, bump_on_collision: bool = True) -> str:
        """Build the folder_number string from autosave widgets.

        Pattern: {prefix}_{number:03d}[_{suffix}][_{HHMMSS}]
        Also auto-increments spin_autosave_num for next run.

        If bump_on_collision is True and the resulting meas_* folder
        already exists under today's date directory, keeps incrementing
        the number until a free slot is found.
        """
        prefix = self.txt_prefix.text().strip() or "data"
        num = self.spin_autosave_num.value()
        suffix = self.txt_suffix.text().strip()
        ts = datetime.now().strftime('%H%M%S') if self.chk_timestamp.isChecked() else ''

        dd = Path("..") / "Saved_Data" / time.strftime("%Y-%m-%d")

        if bump_on_collision and dd.exists():
            # Scan forward until we find a free number
            for _ in range(1000):
                parts = [prefix, f"{num:03d}"]
                if suffix:
                    parts.append(suffix)
                # Don't include timestamp in collision check (it changes)
                fn_check = '_'.join(parts)
                if not (dd / f"meas_{fn_check}").exists():
                    break
                # Also check with timestamp variants (any HHMMSS)
                pattern = str(dd / f"meas_{fn_check}*")
                if not glob.glob(pattern):
                    break
                num += 1
            # Update spinbox to reflect the resolved number
            self.spin_autosave_num.setValue(num)

        parts = [prefix, f"{num:03d}"]
        if suffix:
            parts.append(suffix)
        if ts:
            parts.append(ts)
        fn = '_'.join(parts)

        # Auto-increment for next run
        self.spin_autosave_num.setValue(num + 1)
        return fn

    def _preview_folder_name(self) -> str:
        """Build a preview of the next folder name WITHOUT incrementing.

        Used for the collision indicator in the UI.
        """
        prefix = self.txt_prefix.text().strip() or "data"
        num = self.spin_autosave_num.value()
        suffix = self.txt_suffix.text().strip()
        parts = [prefix, f"{num:03d}"]
        if suffix:
            parts.append(suffix)
        if self.chk_timestamp.isChecked():
            parts.append("HHMMSS")
        return '_'.join(parts)

    def _check_folder_collision(self):
        """Update the collision indicator label.

        Called when prefix, number, or suffix changes.
        Shows ✔ (free) or ✘ (exists) next to the autosave number.
        """
        prefix = self.txt_prefix.text().strip() or "data"
        num = self.spin_autosave_num.value()
        suffix = self.txt_suffix.text().strip()
        parts = [prefix, f"{num:03d}"]
        if suffix:
            parts.append(suffix)
        fn_check = '_'.join(parts)

        dd = Path("..") / "Saved_Data" / time.strftime("%Y-%m-%d")
        folder = dd / f"meas_{fn_check}"

        if folder.exists():
            self._lbl_collision.setText("✘")
            self._lbl_collision.setStyleSheet("color:#f44;font-weight:bold")
            self._lbl_collision.setToolTip(
                f"Folder exists: meas_{fn_check}\n"
                f"Number will auto-bump at run time.")
        else:
            self._lbl_collision.setText("✔")
            self._lbl_collision.setStyleSheet("color:#6a9955;font-weight:bold")
            self._lbl_collision.setToolTip(f"Free: meas_{fn_check}")

    def _get_save_path_and_fn(self):
        """Return (save_path, folder_number) or (None, None) if autosave off.

        Creates Saved_Data/YYYY-MM-DD/meas_{folder_number}/ directory.
        Auto-bumps the number if a collision is detected.
        """
        if not self.chk_save.isChecked():
            return None, None
        fn = self._make_folder_number(bump_on_collision=True)
        dd = Path("..") / "Saved_Data" / time.strftime("%Y-%m-%d")
        dd.mkdir(parents=True, exist_ok=True)
        sp = dd / f"meas_{fn}"; sp.mkdir(exist_ok=True)
        self._last_save_path = sp
        self._last_folder_number = fn
        return sp, fn

    # ── COMMENTS (Phase 6) ──────────────────────────────────────

    def _get_current_comment(self) -> str:
        """Read and return the comment field text."""
        return self.txt_comment.text().strip()

    def _clear_comment_field(self):
        """Clear comment field after measurement completes.

        Also flushes any deferred 'next' comment into the current params.
        """
        # If there was a deferred comment from a previous ▶ press,
        # append it to the just-completed measurement's param file
        deferred = getattr(self, '_deferred_comment', '')
        if deferred:
            self._append_comment_to_params(
                getattr(self, '_last_save_path', None),
                getattr(self, '_last_folder_number', None),
                deferred)
            self._deferred_comment = ''

        self.txt_comment.clear()

    def _comment_prev(self):
        """Write the current comment into the PREVIOUS measurement's param file."""
        comment = self._get_current_comment()
        if not comment:
            self.status_label.setText("No comment to save"); return

        sp = getattr(self, '_last_save_path', None)
        fn = getattr(self, '_last_folder_number', None)
        if sp is None or fn is None:
            self.status_label.setText("No previous measurement"); return

        self._append_comment_to_params(sp, fn, comment)
        self.txt_comment.clear()
        self.status_label.setText(f"✔ Comment → params_{fn}.yaml")

    def _comment_next(self):
        """Defer comment — it will be written into the NEXT measurement's param file."""
        comment = self._get_current_comment()
        if not comment:
            self.status_label.setText("No comment to defer"); return
        self._deferred_comment = comment
        self.txt_comment.clear()
        self.status_label.setText("Comment deferred → next measurement")

    def _append_comment_to_params(self, save_path, folder_number, comment):
        """Append a comment to an existing params YAML file.

        Appends raw YAML text at the end of the file rather than
        re-serializing. This preserves the original formatting from
        experiment_base.savefile (yaml.dump) exactly.
        """
        if save_path is None or folder_number is None:
            return
        pf = save_path / f"params_{folder_number}.yaml"
        if not pf.exists():
            pf = save_path / f"params_{folder_number}.json"
        if not pf.exists():
            print(f"⚠ Cannot find params file for {folder_number}")
            return
        try:
            ts = datetime.now().strftime('%H:%M:%S')
            line = f"    - '[{ts}] {comment}'\n"

            # Check if _comments key already exists in the file
            text = pf.read_text()
            if '_comments:' in text:
                # Append to existing list
                with open(pf, 'a') as f:
                    f.write(line)
            else:
                # Create the key + first entry
                with open(pf, 'a') as f:
                    f.write("\n_comments:\n")
                    f.write(line)

            print(f"  💬 Comment → {pf.name}")
        except Exception as e:
            print(f"⚠ Comment save error: {e}")

    # ── SESSION STATE SAVE / LOAD (Phase 4) ──────────────────────────────
    def _build_state_controls(self, parent_layout: QVBoxLayout):
        """Add state save/load controls to the instruments tab."""

        parent_layout.addWidget(QLabel(""))  # spacer

        r1 = QHBoxLayout()
        lbl = QLabel("Session State")
        lbl.setStyleSheet("color: #888; font-size: 15px")
        r1.addWidget(lbl)

        btn_save_st = QPushButton('Save')
        btn_save_st.setStyleSheet("""
            QPushButton {background-color: #d3dbde; color: #484848; font-size: 13px}
            QPushButton:hover {background-color: #009de0;color: #ffffff;}
        """)
        btn_save_st.setFixedSize(60, 30)
        btn_save_st.clicked.connect(self._save_session_state)
        r1.addWidget(btn_save_st)

        btn_load_st = QPushButton('Load')
        btn_load_st.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; font-size: 13px}
            QPushButton:hover {background-color: #009de0;}
        """)
        btn_load_st.setFixedSize(60, 30)
        btn_load_st.clicked.connect(self._load_session_state)
        r1.addWidget(btn_load_st)
        parent_layout.addLayout(r1)

        r1 = QHBoxLayout()
        self.state_combo = QComboBox()
        self.state_combo.setEditable(True)
        self.state_combo.setToolTip(
            "Select or type a session state name.\n"
            "Saves/loads: config module, script, mode, autosave settings,\n"
            "contrast op, window geometry, instrument sim flag.")
        r1.addWidget(self.state_combo)

        # Refresh button inline with combo
        refresh_btn = QPushButton('↻')
        refresh_btn.setStyleSheet("""
            QPushButton {background-color: #484848; color: #ffffff; padding: 2px 2px;}
            QPushButton:hover {background-color: #009de0; color: #ffffff;}
        """)
        refresh_btn.setFixedSize(25, 28)
        refresh_btn.setToolTip("Refresh file list")
        refresh_btn.clicked.connect(self._refresh_state_files)
        r1.addWidget(refresh_btn)

        parent_layout.addLayout(r1)
        parent_layout.addSpacing(10)
        self._refresh_state_files()

    def _refresh_state_files(self):
        """Refresh the list of state files (.json) in the state save directory."""
        # Save current text if user typed something
        current_text = self.state_combo.currentText()

        self.state_combo.clear()

        try:
            self.state_combo.addItem(Path(DEFAULT_STATE_FILE).stem, DEFAULT_STATE_FILE)
            if self.statefile_directory.exists():
                json_files = sorted(self.statefile_directory.glob("*.json"))
                # json_files.remove(self.statefile_directory / DEFAULT_STATE_FILE)
                for f in json_files:
                    # if f.name != CHANNEL_CONFIG_FILE:
                    self.state_combo.addItem(f.stem, str(f))
        except Exception as e:
            self.status_label.setText(f'Error refreshing state files: {e}')
            # pass

    def _save_session_state(self):
        """Save current session UI state to a JSON file."""
        try:
            # Get filename from combo (typed or selected)
            filename = Path(self.state_combo.currentText().strip())# if filename is None else Path(filename)

            if not filename:
                filename = Path(f"session_state_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            
            # Remove .json extension if user typed it
            filepath = self.statefile_directory / filename.with_suffix('.json')

            state = {
                'window_maximized': self.isMaximized(),
                'window_width': self.width(),
                'window_height': self.height(),
                'script_file': self.cb_scr.currentText(),
                'config_modules': self.last_used_configs,
                'mode': self.cb_mode.currentText(),
                'contrast_op': self.cb_contrast.currentText(),
                'realtime_plot': self.chk_rt.isChecked(),
                'autosave_enabled': self.chk_save.isChecked(),
                'autosave_prefix': self.txt_prefix.text(),
                'autosave_start_num': self.spin_autosave_num.value(),
                'autosave_suffix': self.txt_suffix.text(),
                'autosave_timestamp': self.chk_timestamp.isChecked(),
                'sim_sg': self.simulate_checkboxes['SG'].isChecked(),
                'ts_duration': self.txt_ts_duration.text(),
            }

            with open(filepath, 'w') as f:
                json.dump(state, f, indent=2)

            self.status_label.setText(f"✔ State saved: {filename}")
            self._refresh_state_files()
            idx = self.state_combo.findText(filename.stem)
            if idx >= 0:
                self.state_combo.setCurrentIndex(idx)
        except Exception as e:
            self.status_label.setText(f"❌ State save: {e}")

    def _load_session_state(self):
        """Load session UI state from a JSON file."""
        try:
            filename = Path(self.state_combo.currentText().strip())# if filename is None else Path(filename)

            if not filename:
                self.status_label.setText("No state file specified")
                return
            
            filepath = self.statefile_directory / filename.with_suffix(".json")
            # else:
            #     filepath = Path(filepath)

            if not filepath.exists():
                self.status_label.setText(f"❌ Not found: {filename}.json")
                return

            with open(filepath, 'r') as f:
                state = json.load(f)

            _was_window_maximized = state.get('window_maximized', False)
            if not _was_window_maximized:
                self.resize(state.get('window_width', 1100),
                            state.get('window_height', 720))
            else:
                self.showMaximized()
            
            script_file = state.get('script_file', '')
            self.cb_scr.setCurrentText(script_file)
            self.last_used_configs = state.get('config_modules', '')
            self.cb_cfg.setCurrentText(self.last_used_configs[0] if 'diode' in script_file else self.last_used_configs[1])
            self.cb_mode.setCurrentText(state.get('mode', 'Sweep'))
            self.cb_contrast.setCurrentText(state.get('contrast_op', 's/r'))
            self.chk_rt.setChecked(state.get('realtime_plot', True))
            self.chk_save.setChecked(state.get('autosave_enabled', True))
            self.txt_prefix.setText(state.get('autosave_prefix', 'data'))
            self.spin_autosave_num.setValue(
                state.get('autosave_start_num', 1))
            self.txt_suffix.setText(state.get('autosave_suffix', ''))
            self.chk_timestamp.setChecked(
                state.get('autosave_timestamp', False))
            self.simulate_checkboxes['SG'].setChecked(state.get('sim_sg', True))
            self.txt_ts_duration.setText(state.get('ts_duration', '0'))

            self.status_label.setText(f"✔ State loaded: {filename}")
        except Exception as e:
            self.status_label.setText(f"❌ State load: {e}")

    # ── LOAD CONFIG ────────────────────────────────────────────────

    def _load_config(self):
        cn = self.cb_cfg.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
        except Exception as e:
            self.status_label.setText(f"❌ Import: {e}"); return
        if not hasattr(cm,'config'):
            self.status_label.setText(f"❌ {cn} has no 'config'"); return
        self._loaded_config = cm.config
        self._loaded_config_name = cn
        self._update_metadata_tree(cm.config)
        self._on_config_load()
        self._push_namespace()
        self.status_label.setText(f"✔ Loaded {cn}")
    
    def _on_config_load(self):
        # self._tabs.setCurrentIndex(2)
        self._tabs.tabBar().setTabTextColor(2, QColor("#00ff00"))
        QTimer.singleShot(5000,
                          lambda: self._tabs.tabBar().setTabTextColor(2, QColor("#ffffff")))
    
    def _ensure_config(self):
        cn = self.cb_cfg.currentText()
        if self._loaded_config is None or self._loaded_config_name != cn:
            self._load_config()
        return self._loaded_config

    # ── RUN ─────────────────────────────────────────────────────────────

    def _run_experiment(self):
        if not self.instruments_connected: return
        mode = self.cb_mode.currentText()
        if mode == 'Timeseries':
            # Check if current config is camera-based
            cfg = self._ensure_config()
            if cfg and self._is_camera_config(cfg):
                self._run_camera_timeseries()
            else:
                self._run_timeseries()
        else:
            self._run_sweep()

    def _run_sweep(self):
        cn = self.cb_cfg.currentText(); sn = self.cb_scr.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
            sm = importlib.import_module(sn); importlib.reload(sm)
        except Exception as e:
            self.status_label.setText(f"❌ Import: {e}"); return

        if not hasattr(sm, 'run'):
            self.status_label.setText(f"❌ {sn} has no run()"); return

        # Get the ExperimentConfig object from config module
        if not hasattr(cm, 'config'):
            self.status_label.setText(f"❌ {cn} has no 'config' attribute"); return
        
        exp_config = cm.config
        self._loaded_config = exp_config
        self._loaded_config_name = cn
        self._update_metadata_tree(exp_config)

        # Detect camera experiment
        is_camera = self._is_camera_config(exp_config)

        # Build instruments dict
        instr = {'sg': self.sg, 'pb': self.pb, 'ao_task': self.ao_task}
        
        if is_camera:
            cw = self._get_or_create_camera_worker(exp_config)
            if cw is None:
                return
            instr['camera_worker'] = cw

        sp, fn = self._get_save_path_and_fn()

        # Inject current comment into config extra dict for param save
        comment = self._get_current_comment()
        if comment:
            exp_config.extra = exp_config.extra or {}
            exp_config.extra['_comment'] = (
                f"[{datetime.now().strftime('%H:%M:%S')}] {comment}")

        # Switch plot layout based on experiment type
        if is_camera:
            self._setup_camera_plots(exp_config)
        else:
            self._setup_diode_plots(exp_config)

        self._exp_config = exp_config
        primary_name = exp_config.scan_names[0]
        self._inner_vals = exp_config.scans[primary_name].values
        self._xu = exp_config.plot.x_units
        self._Nruns = exp_config.runtime.Nruns
        self._is_camera_run = is_camera

        # Thread
        self.experiment_thread = ExperimentThread(sm.run, instr, exp_config, sp, fn)
        self.experiment_thread.point_signal.connect(lambda *a: self._upd.emit(*a))
        self.experiment_thread.finished_signal.connect(self._on_done)
        self.experiment_thread.error_signal.connect(self._on_err)
        self.experiment_thread.start()

        self._set_running_state(True)
        self.status_label.setText(f"Running {cn}…")

    # ── TIMESERIES DURATION ──────────────────────────────────

    def _get_ts_duration(self) -> float:
        """Read the TS duration field. Returns 0.0 for manual stop."""
        try:
            val = float(self.txt_ts_duration.text().strip())
            return max(val, 0.0)
        except ValueError:
            return 0.0

    def _apply_ts_duration(self):
        """Called on Enter key in the TS duration field.

        If a timeseries is running, restarts the auto-stop timer
        with the new duration (counted from NOW, not from start).
        If not running, the value is just stored for the next run.
        """
        dur = self._get_ts_duration()

        if not self.is_running:
            self.status_label.setText(
                f"TS duration set: {dur:.1f}s"
                if dur > 0 else "TS duration: manual stop")
            return

        # ── Diode timeseries: manipulate the experiment's QTimer ──
        exp = getattr(self, '_ts_exp', None)
        if exp is not None and exp.is_running:
            # Cancel existing auto-stop timer
            if exp._auto_stop_timer is not None:
                exp._auto_stop_timer.stop()
                exp._auto_stop_timer = None

            if dur > 0:
                # Remaining = dur - elapsed so far
                elapsed = time.perf_counter() - exp._start_time
                remaining = max(dur - elapsed, 0.1)
                exp._auto_stop_timer = QTimer()
                exp._auto_stop_timer.setSingleShot(True)
                exp._auto_stop_timer.timeout.connect(exp.stop)
                exp._auto_stop_timer.start(int(remaining * 1000))
                self.status_label.setText(
                    f"TS: auto-stop updated → {remaining:.1f}s remaining")
            else:
                self.status_label.setText("TS: switched to manual stop")
            return

        # ── Camera timeseries: session-side timer ──
        cam_exp = getattr(self, '_cam_ts_exp', None)
        if cam_exp is not None and cam_exp.is_running:
            # Cancel existing session timer
            t = getattr(self, '_ts_auto_stop', None)
            if t is not None:
                t.stop()
                self._ts_auto_stop = None

            if dur > 0:
                elapsed = time.perf_counter() - cam_exp._start_time
                remaining = max(dur - elapsed, 0.1)
                self._ts_auto_stop = QTimer()
                self._ts_auto_stop.setSingleShot(True)
                self._ts_auto_stop.timeout.connect(self._cam_ts_stop)
                self._ts_auto_stop.start(int(remaining * 1000))
                self.status_label.setText(f"CamTS: auto-stop updated → {remaining:.1f}s remaining")
            else:
                # Also clear the experiment's internal auto-stop
                cam_exp._auto_stop_time = None
                self.status_label.setText("CamTS: switched to manual stop")
            return

    # ── RUN TIMESERIES ────────────────────────────────────────────

    def _run_timeseries(self):
        cn = self.cb_cfg.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
        except Exception as e:
            self.status_label.setText(f"❌ Import: {e}"); return
        if not hasattr(cm, 'config'):
            self.status_label.setText(f"❌ {cn} has no 'config'"); return

        exp_config = cm.config
        self._loaded_config = exp_config
        self._loaded_config_name = cn
        self._update_metadata_tree(exp_config)

        instr = {'sg': self.sg, 'pb': self.pb, 'ao_task': self.ao_task}

        from timeseries_experiment import TimeseriesExperiment

        try:
            self._ts_exp = TimeseriesExperiment(instr, exp_config)
            self._ts_exp.setup()
        except Exception as e:
            self.status_label.setText(f"❌ TS setup: {e}")
            print(traceback.format_exc()); return

        # Switch plots to timeseries mode
        self.sweep_plot.clear()
        self.sweep_plot.setTitle("Contrast Timeseries")
        self.sweep_plot.setLabel('bottom', 'Time', units='s')
        self.sweep_plot.setLabel('left', 'Contrast')
        self._ts_contrast_curve = self.sweep_plot.plot(pen=_pen(0, 220, 1.5))
        self.raw_curve = self.raw_plot.plot(pen=_pen(3, width=1), clear=True)

        # 30 Hz display refresh
        self._ts_display_timer = QTimer()
        self._ts_display_timer.timeout.connect(self._ts_refresh)
        self._ts_display_timer.start(33)

        # Save path
        sp, fn = self._get_save_path_and_fn()
        self._ts_save_path = sp
        self._ts_fn = fn

        self._ts_exp.start(max_duration=self._get_ts_duration())

        self._set_running_state(True)

        self.status_label.setText(
            f"TS: {self._ts_exp.sequence} @ "
            f"{self._ts_exp.scan_name}={self._ts_exp.fixed_value:.4g}")

    def _ts_refresh(self):
        """30 Hz display update for timeseries mode."""
        if not hasattr(self, '_ts_exp') or self._ts_exp is None:
            return
        exp = self._ts_exp

        # Contrast scrolling plot
        if exp.contrast_ring is not None:
            cdata = exp.contrast_ring.get_ordered()
            if len(cdata) > 0:
                sample_rate = exp.cfg.daq_ai.ai_sample_rate
                contrast_rate = sample_rate / exp.chunk_samples
                t = np.arange(-len(cdata), 0) / contrast_rate
                self._ts_contrast_curve.setData(t, cdata)

        # Raw: show latest chunk
        if exp.raw_ring is not None and exp.raw_ring.total_count > 0:
            rdata = exp.raw_ring.get_ordered()
            n_show = min(len(rdata), exp.chunk_samples * 3)
            if n_show > 0:
                self.raw_curve.setData(rdata[-n_show:])

        # Status
        elapsed = time.perf_counter() - exp._start_time
        n_c = exp.contrast_ring.total_count if exp.contrast_ring else 0
        self.status_label.setText(
            f"TS: {elapsed:.1f}s | {n_c} contrast pts | "
            f"{sum(len(c) for c in exp.full_chunks):,} raw samples")

        # Check if auto-stopped
        if not exp.is_running and self.is_running:
            self._ts_stop()

    def _ts_stop(self):
        """Stop timeseries and save."""
        if hasattr(self, '_ts_display_timer'):
            self._ts_display_timer.stop()
        if hasattr(self, '_ts_exp') and self._ts_exp is not None:
            self._ts_exp.stop()
            if self._ts_save_path:
                fn = getattr(self, '_ts_fn', '001')
                self._ts_exp.save(
                    self._ts_save_path / f"ts_{fn}", fmt='npz')
            self._ts_exp.teardown()

        self._set_running_state(False)
        self._is_camera_run = False
        self.status_label.setText("✔ Timeseries done")
        self._push_namespace()

    # ── CAMERA TIMESERIES ────────────────────────────────

    def _run_camera_timeseries(self):
        """Run PB-triggered continuous camera acquisition at a fixed point."""
        cn = self.cb_cfg.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
        except Exception as e:
            self.status_label.setText(f"❌ Import: {e}"); return
        if not hasattr(cm, 'config'):
            self.status_label.setText(f"❌ {cn} has no 'config'"); return

        exp_config = cm.config
        self._loaded_config = exp_config
        self._loaded_config_name = cn
        self._update_metadata_tree(exp_config)

        # Get camera worker
        cw = self._get_or_create_camera_worker(exp_config)
        if cw is None:
            return

        instr = {
            'sg': self.sg, 'pb': self.pb, 'ao_task': self.ao_task,
            'camera_worker': cw,
        }

        from camera_experiment import CameraTimeseriesExperiment

        try:
            self._cam_ts_exp = CameraTimeseriesExperiment(instr, exp_config)
            self._cam_ts_exp.setup()
        except Exception as e:
            self.status_label.setText(f"❌ CamTS setup: {e}")
            print(traceback.format_exc()); return

        # ── 3-panel layout: Image | Contrast + Intensity ─────────
        self._teardown_dynamic_plots()
        plot_wg = self.sweep_plot.parent()
        layout = plot_wg.layout()
        self.sweep_plot.hide()
        self.raw_plot.hide()

        # Top: Image panel
        self._cam_plot = pg.PlotWidget(title="Camera Frame")
        self._cam_plot.setAspectLocked(True)
        self._cam_plot.invertY(True)
        self._cam_plot.hideAxis('bottom')
        self._cam_plot.hideAxis('left')
        self._cam_image = pg.ImageItem()
        self._cam_plot.addItem(self._cam_image)
        self._cam_cbar = pg.ColorBarItem(
            colorMap=pg.colormap.get('inferno'),
            interactive=False, width=15)
        self._cam_cbar.setImageItem(self._cam_image)
        layout.addWidget(self._cam_plot, stretch=2)

        # Bottom row: Contrast scrolling | Intensity scrolling
        bottom_w = QWidget()
        bottom_l = QHBoxLayout(bottom_w)
        bottom_l.setContentsMargins(0, 0, 0, 0)

        self._contrast_plot = pg.PlotWidget(title="Contrast Timeseries")
        self._contrast_plot.setLabel('bottom', 'Time', units='s')
        self._contrast_plot.setLabel('left', 'S / R')
        self._contrast_plot.showGrid(x=True, y=True, alpha=0.3)
        self._contrast_curve = self._contrast_plot.plot(
            pen=_pen(2, 220, 1.5), name='S/R')
        bottom_l.addWidget(self._contrast_plot)

        self._intensity_plot = pg.PlotWidget(title="Intensity Timeseries")
        self._intensity_plot.addLegend(offset=(10, 10))
        self._intensity_plot.setLabel('bottom', 'Time', units='s')
        self._intensity_plot.setLabel('left', 'Pixel mean')
        self._intensity_plot.showGrid(x=True, y=True, alpha=0.3)
        self._intensity_sig = self._intensity_plot.plot(
            pen=_pen(0, 220, 1.5), name='Sig')
        self._intensity_ref = self._intensity_plot.plot(
            pen=_pen(1, 220, 1.5, dash=True), name='Ref')
        bottom_l.addWidget(self._intensity_plot)

        layout.addWidget(bottom_w, stretch=1)
        self._cam_bottom_widget = bottom_w

        # 30 Hz display refresh
        self._cam_ts_timer = QTimer()
        self._cam_ts_timer.timeout.connect(self._cam_ts_refresh)
        self._cam_ts_timer.start(33)

        # Save path
        sp, fn = self._get_save_path_and_fn()
        self._cam_ts_save_path = sp
        self._cam_ts_fn = fn

        self._cam_ts_exp.max_duration = self._get_ts_duration()
        # Ensure _auto_stop_time exists before polling thread checks it
        # (CameraTimeseriesExperiment.start() sets it, but the polling
        #  thread may race ahead and read it before start() finishes)
        if not hasattr(self._cam_ts_exp, '_auto_stop_time'):
            self._cam_ts_exp._auto_stop_time = None
        self._cam_ts_exp.start()
        self._is_camera_run = True

        # Session-side auto-stop timer (works for both diode and camera TS)
        dur = self._get_ts_duration()
        if dur > 0:
            self._ts_auto_stop = QTimer()
            self._ts_auto_stop.setSingleShot(True)
            self._ts_auto_stop.timeout.connect(self._cam_ts_stop)
            self._ts_auto_stop.start(int(dur * 1000))

        self._set_running_state(True)
        self.status_label.setText(
            f"CamTS: {self._cam_ts_exp.sequence} @ "
            f"{self._cam_ts_exp.scan_name}="
            f"{self._cam_ts_exp.fixed_value:.4g}")

    def _cam_ts_refresh(self):
        """30 Hz display update for camera timeseries — all 3 panels."""
        if not hasattr(self, '_cam_ts_exp') or self._cam_ts_exp is None:
            return
        exp = self._cam_ts_exp

        # ── Image panel: latest signal frame ──
        if exp.last_frame is not None and self._cam_image is not None:
            self._cam_image.setImage(exp.last_frame.T, autoLevels=True)

        # ── Contrast scrolling ──
        if (exp.contrast_ring is not None
                and self._contrast_curve is not None):
            cdata = exp.contrast_ring.get_ordered()
            if len(cdata) > 0:
                t = np.arange(-len(cdata), 0) / exp.contrast_rate
                self._contrast_curve.setData(t, cdata)

        # ── Intensity scrolling (sig + ref) ──
        if (exp.sig_ring is not None
                and self._intensity_sig is not None):
            sdata = exp.sig_ring.get_ordered()
            rdata = exp.ref_ring.get_ordered()
            if len(sdata) > 0:
                t = np.arange(-len(sdata), 0) / exp.contrast_rate
                self._intensity_sig.setData(t, sdata)
            if len(rdata) > 0:
                t = np.arange(-len(rdata), 0) / exp.contrast_rate
                self._intensity_ref.setData(t, rdata)

        # ── Status ──
        elapsed = time.perf_counter() - exp._start_time
        n_c = exp.contrast_ring.total_count if exp.contrast_ring else 0
        self.status_label.setText(
            f"CamTS: {elapsed:.1f}s | {n_c} contrast pts | "
            f"{exp.cycle_count} cycles")

        # Check if auto-stopped
        if not exp.is_running and self.is_running:
            self._cam_ts_stop()

    def _cam_ts_stop(self):
        """Stop camera timeseries and save."""
        if hasattr(self, '_cam_ts_timer'):
            self._cam_ts_timer.stop()
        # Cancel session-side auto-stop timer
        t = getattr(self, '_ts_auto_stop', None)
        if t is not None:
            t.stop()
            self._ts_auto_stop = None
        if hasattr(self, '_cam_ts_exp') and self._cam_ts_exp is not None:
            self._cam_ts_exp.stop()
            if getattr(self, '_cam_ts_save_path', None):
                fn = getattr(self, '_cam_ts_fn', '001')
                self._cam_ts_exp.save(
                    self._cam_ts_save_path / f"cam_ts_{fn}")
            self._cam_ts_exp.teardown()

        self._set_running_state(False)
        self._is_camera_run = False
        self.status_label.setText("✔ Camera timeseries done")
        self._push_namespace()

    # ════════════════════════════════════════════════════════════════
    #  VIEW SEQUENCE (PB only)
    # ════════════════════════════════════════════════════════════════

    def _view_sequence(self):   # TODO: should not require PB - this is just plot
        if not self._pb_ready:
            self.status_label.setText("❌ PB not initialized"); return
        cfg = self._ensure_config()
        if cfg is None: return

        instr = {'sg':self.sg,'pb':self.pb,'ao_task':self.ao_task}
        from experiment_base import DiodeExperiment
        import matplotlib.pyplot as plt
        try:
            exp = DiodeExperiment(instr, cfg)
            print("After experiment init...")
            idx = getattr(cfg.runtime,'seq_plot_indices',[0,-1])
            exp.view_sequences(indices=idx, dpi=100)
            plt.show()
            self.status_label.setText("✔ Sequence plotted")
        except Exception as e:
            self.status_label.setText(f"❌ View: {e}"); print(traceback.format_exc())

    # ════════════════════════════════════════════════════════════════
    #  PLAY SEQUENCE (PB only)
    # ════════════════════════════════════════════════════════════════

    def _play_sequence(self):
        """Program + start PB at last seq_plot_indices value."""
        if not self._pb_ready:
            self.status_label.setText("❌ PB not initialized"); return
        cfg = self._ensure_config()
        if cfg is None: return

        try:
            args_names = cfg.seq.args_names
            seq_args = list(cfg.seq.args_values)
            sn = cfg.scan_names[0] if cfg.scan_names else ''
            plot_wg_sub = cfg.scans[sn].values if sn else np.array([0])

            indices = getattr(cfg.runtime,'seq_plot_indices',[0,-1])
            pi = indices[-1]
            if pi < 0: pi = len(plot_wg_sub) + pi
            val = plot_wg_sub[pi]

            is_freq = cfg.seq.name in (
                'esr_dig_mod_seq','esr_seq','pesr_seq','modesr','drift_seq')
            if is_freq:
                sal = seq_args
            elif sn in args_names:
                sal = list(seq_args); sal[args_names.index(sn)] = val
            else:
                sal = [val] + seq_args

            _, the_list = PulseBlaster.PB_program(
                'diode', cfg.seq.name, sal + [cfg.pb.channels])
            self.pb.run_sequence_for_diode(
                [the_list[i][0] for i in range(len(the_list))])
            self.status_label.setText(f"▶ PB: {cfg.seq.name} @ {sn}={val:.4g}")
        except Exception as e:
            self.status_label.setText(f"❌ Play: {e}"); print(traceback.format_exc())

    def _stop_sequence(self):
        if not self._pb_ready:
            self.status_label.setText("❌ PB not initialized"); return
        cfg = self._ensure_config()
        if cfg is None: return

        try:
            self.pb.stop_sequence()
            self.status_label.setText(f"🛑 PB")
        except Exception as e:
            self.status_label.setText(f"❌ Stop: {e}"); print(traceback.format_exc())

    # ── real-time plot (main thread) ────────────────────────────────────

    @Slot(int, float, int, int, str, object, object)
    def _on_plot(self, inner_idx, val, i_run, oc_idx, oc_str, proc, raw):
        if getattr(self, '_is_camera_run', False):
            self._on_camera_plot(inner_idx, val, i_run, oc_idx, oc_str,
                                 proc, raw)
            return

        self.status_label.setText(
            f"Outer {oc_idx+1}  Run {i_run+1}/{self._Nruns}  "
            f"Pt {inner_idx+1}/{len(self._inner_vals)}"
            f"{'  '+oc_str if oc_str else ''}")

        if not self.chk_rt.isChecked():
            return

        # Curve key: unique per (outer_combo, run, data_type)
        alpha = 100 if self._Nruns > 1 else 220
        color_base = oc_idx * 2  # offset colors per outer combo

        if proc is not None and len(proc) >= 1:
            ms = proc[0]
            n = ms.shape[0]

            # For shuffled data: the processed data is a slice of whatever
            # points have been filled so far.  We need to plot at the
            # correct x positions.  The callback sends inner_idx (actual
            # position in the values array), so we can build scatter-style.
            # However, process_data returns contiguous means for filled pts.
            # Simplest correct approach: plot sig for all filled inner points.

            # Get the full data slice for this (outer, run) from the
            # experiment's perspective.  Since we're in the GUI thread and
            # the experiment thread is blocked on DAQ I/O, accessing the
            # shared data_array is safe here (GIL).

            key_sig = f"sig_o{oc_idx}_r{i_run}"
            if key_sig not in self._curves:
                name = f"S o{oc_idx}" if i_run == 0 else None
                self._curves[key_sig] = self.sweep_plot.plot(
                    pen=_pen(color_base, alpha, 1.5),
                    name=name, symbol='o', symbolSize=4,
                    symbolBrush=pg.mkColor(
                        _PAL[color_base % len(_PAL)]))

            # x-coords for filled points: we use inner_values directly
            # since data_array is indexed by actual inner position
            xs_all = self._inner_vals / self._xu
            # Plot whatever points are nonzero in this run's data
            # This handles shuffled order naturally
            self._curves[key_sig].setData(xs_all[:n], ms[:n, 0])

            # Contrast on right y-axis
            if len(proc) >= 3 and getattr(self, '_contrast_vb', None) is not None:
                contrast_vals = self._compute_display_contrast(
                    proc[0][:n, 0], proc[1][:n, 0])
                key_con = f"con_o{oc_idx}_r{i_run}"
                if key_con not in self._curves:
                    nm = f"C o{oc_idx}" if i_run==0 else None
                    crv = pg.PlotDataItem(
                        pen=_pen(color_base + 1, alpha, 1.5, dash=True),
                        name=nm)
                    self._contrast_vb.addItem(crv)
                    self._curves[key_con] = crv
                self._curves[key_con].setData(xs_all[:n], contrast_vals)
                self._contrast_vb.autoRange()
        # Raw
        if raw is not None:
            self.raw_curve.setData(np.ravel(np.array(raw)))

    @Slot(object, object)
    def _on_done(self, exp, data):
        self.last_exp = exp
        self._set_running_state(False)
        sh = " × ".join(str(s) for s in data.shape)
        self.status_label.setText(f"✔ Done — data {sh}")

        # Restore diode plot layout if camera was used
        if getattr(self, '_is_camera_run', False):
            self._is_camera_run = False
        self._push_namespace()
        self._clear_comment_field()

    @Slot(str)
    def _on_err(self, tb):
        self._set_running_state(False)
        self.status_label.setText("❌ Error (see console)")
        print(tb)
        if getattr(self, '_is_camera_run', False):
            self._is_camera_run = False

    def _stop_experiment(self):
        if self.cb_mode.currentText() == 'Timeseries':
            # Check if it's a camera timeseries
            if hasattr(self, '_cam_ts_exp') and self._cam_ts_exp is not None:
                self._cam_ts_stop()
            else:
                self._ts_stop()
        elif self.experiment_thread:
            self.experiment_thread.request_stop()
            self.status_label.setText("⏹ Stopping…")

    # ════════════════════════════════════════════════════════════════
    #  CAMERA EXPERIMENT — plot layout + callback + helpers
    # ════════════════════════════════════════════════════════════════

    def _is_camera_config(self, config) -> bool:
        """Detect if config is for a camera experiment."""
        return (hasattr(config, 'camera')
                and config._is_active('camera'))

    def _get_or_create_camera_worker(self, config):
        """Get or create a CameraWorker for camera experiments.

        Reads settings from camera viewer if open, otherwise creates
        a fresh CameraWorker with settings from config.
        """
        # If camera viewer is open, close it first (DCAM exclusive access)
        if self._camera_viewer is not None:
            try:
                if self._camera_viewer.isVisible():
                    # Read back settings before closing
                    viewer_settings = self._get_camera_viewer_settings()
                    if viewer_settings:
                        # Apply viewer settings to config if different
                        if 'actual_exposure_ms' in viewer_settings:
                            config.camera.exposure_s = (
                                viewer_settings['actual_exposure_ms'] / 1e3)
                        if 'actual_roi' in viewer_settings:
                            config.camera.roi = list(
                                viewer_settings['actual_roi'])
                    self._camera_viewer._disconnect_camera()
                    self._camera_viewer.close()
                    self.lbl_cam.setText("Cam: ○ closed for experiment")
                    QApplication.processEvents()
                    time.sleep(0.5)  # DCAM cleanup
            except Exception as e:
                print(f"⚠ Viewer close: {e}")
            self._camera_viewer = None

        # Check if we already have a camera worker from a previous run
        if hasattr(self, '_camera_worker') and self._camera_worker is not None:
            try:
                self._camera_worker.query_camera_status()
                self.lbl_cam.setText("Cam: ● reusing worker")
                return self._camera_worker
            except Exception:
                self._camera_worker = None

        # Create new CameraWorker
        try:
            from Camcontrol import CameraWorker
            cw = CameraWorker(
                simulate=self.chk_sim.isChecked(),
                roi=config.camera.roi,
                exposure=config.camera.exposure_s,
            )
            cw.init_cam()
            self._camera_worker = cw
            self.lbl_cam.setText(f"Cam: ● {cw.device_title}")
            self.status_label.setText("✔ Camera initialized")
            return cw
        except Exception as e:
            self.status_label.setText(f"❌ Camera init: {e}")
            print(traceback.format_exc())
            return None

    def _setup_camera_plots(self, config):
        """Switch plot area to 3-panel camera layout.

        Layout:
            Top:    Image panel (latest signal frame)
            Middle: Sweep — left Y = intensity (sig/ref), right Y = contrast
            Bottom: Raw trace (latest raw frame pixel values, 1D)

        This mirrors the diode layout (sweep above, raw below) with
        an image panel added on top.
        """
        # ── Tear down any previous camera or diode panels ────────
        self._teardown_dynamic_plots()

        plot_wg = self.sweep_plot.parent()
        layout = plot_wg.layout()

        # Hide default diode plots (stay parented but hidden)
        self.sweep_plot.hide()
        self.raw_plot.hide()

        # ── Image panel ──────────────────────────────────────────
        self._cam_plot = pg.PlotWidget(title="Camera Frame")
        self._cam_plot.setAspectLocked(True)
        self._cam_plot.invertY(True)
        self._cam_plot.hideAxis('bottom')
        self._cam_plot.hideAxis('left')
        self._cam_image = pg.ImageItem()
        self._cam_plot.addItem(self._cam_image)

        # ROI overlay (yellow rectangle matching camera subarray)
        roi = config.camera.roi
        if roi and len(roi) == 4 and roi != [0, 0, 2048, 2048]:
            roi_rect = pg.RectROI(
                [0, 0], [roi[2], roi[3]],
                pen=pg.mkPen('#dcdcaa', width=2),
                movable=False, removable=False)
            roi_rect.removeHandle(0)  # non-interactive
            self._cam_plot.addItem(roi_rect)
            self._cam_roi_rect = roi_rect
        else:
            self._cam_roi_rect = None

        # Colorbar
        self._cam_cbar = pg.ColorBarItem(
            colorMap=pg.colormap.get('inferno'),
            interactive=False, width=15)
        self._cam_cbar.setImageItem(self._cam_image)

        layout.addWidget(self._cam_plot, stretch=2)

        # ── Sweep plot (intensity left axis, contrast right axis) ──
        self._intensity_plot = pg.PlotWidget(title="Sweep")
        self._intensity_plot.addLegend(offset=(10, 10))
        self._intensity_plot.setLabel('bottom', config.plot.x_label)
        self._intensity_plot.setLabel('left', 'Intensity')
        self._intensity_plot.showGrid(x=True, y=True, alpha=0.3)
        # No initial curves — callback creates per-run curves via _curves dict

        # Right axis for contrast
        self._cam_contrast_vb = pg.ViewBox()
        self._intensity_plot.scene().addItem(self._cam_contrast_vb)
        self._intensity_plot.getAxis('right').linkToView(self._cam_contrast_vb)
        self._cam_contrast_vb.setXLink(self._intensity_plot)
        self._intensity_plot.showAxis('right')
        self._intensity_plot.setLabel('right', 'Contrast')

        def _sync_cam_vb():
            self._cam_contrast_vb.setGeometry(
                self._intensity_plot.getViewBox().sceneBoundingRect())
        self._intensity_plot.getViewBox().sigResized.connect(_sync_cam_vb)

        self._curves.clear()  # fresh per experiment

        layout.addWidget(self._intensity_plot, stretch=2)

        # ── Raw trace (bottom) ───────────────────────────────────
        self._cam_raw_plot = pg.PlotWidget(title="Raw Trace")
        self._cam_raw_plot.setLabel('bottom', 'Pixel')
        self._cam_raw_plot.showGrid(x=True, y=True, alpha=0.3)
        self._cam_raw_curve = self._cam_raw_plot.plot(pen=_pen(3, width=1))
        layout.addWidget(self._cam_raw_plot, stretch=1)

        # Bundle the sweep+raw into a container for teardown
        self._cam_bottom_widget = QWidget()  # placeholder for teardown
        # Store actual widgets for cleanup
        self._cam_extra_plots = [self._intensity_plot, self._cam_raw_plot]

    def _teardown_dynamic_plots(self):
        """Remove any dynamically-added camera/timeseries panels.

        Called before setting up a new plot layout to prevent stacking.
        """
        plot_wg = self.sweep_plot.parent()
        layout = plot_wg.layout()

        # Remove camera extra plots (intensity sweep plot + raw trace)
        for w in getattr(self, '_cam_extra_plots', []):
            if w is not None:
                layout.removeWidget(w)
                w.setParent(None)
                w.deleteLater()
        self._cam_extra_plots = []

        for attr in ('_cam_plot', '_cam_bottom_widget'):
            w = getattr(self, attr, None)
            if w is not None:
                layout.removeWidget(w)
                w.setParent(None)
                w.deleteLater()
                setattr(self, attr, None)

        # Clean up right-axis contrast ViewBox (diode dual-axis sweep plot)
        vb = getattr(self, '_contrast_vb', None)
        if vb is not None:
            self.sweep_plot.scene().removeItem(vb)
            self._contrast_vb = None
            self.sweep_plot.hideAxis('right')

        # Clean up camera contrast ViewBox
        vb2 = getattr(self, '_cam_contrast_vb', None)
        if vb2 is not None:
            self._cam_contrast_vb = None

        # Clear dangling refs to curves/items inside destroyed widgets
        for attr in ('_cam_image', '_cam_roi_rect', '_cam_cbar',
                     '_intensity_plot', '_intensity_sig', '_intensity_ref',
                     '_contrast_plot', '_contrast_curve',
                     '_cam_raw_plot', '_cam_raw_curve'):
            if hasattr(self, attr):
                setattr(self, attr, None)

    def _setup_diode_plots(self, config):
        """Restore standard diode 2-panel layout (sweep + raw trace).

        Sweep plot: left Y = signal mean, right Y = contrast (dual axis).
        Raw plot: latest raw DAQ trace.
        """
        self._teardown_dynamic_plots()

        plot_wg = self.sweep_plot.parent()
        layout = plot_wg.layout()

        # Show diode plots
        self.sweep_plot.show()
        self.raw_plot.show()
        if layout.indexOf(self.sweep_plot) < 0:
            layout.addWidget(self.sweep_plot, stretch=2)
        if layout.indexOf(self.raw_plot) < 0:
            layout.addWidget(self.raw_plot, stretch=1)

        # Reset sweep plot
        self.sweep_plot.clear()
        self.sweep_plot.addLegend(offset=(10, 10))
        self.sweep_plot.setLabel('bottom', config.plot.x_label)
        self.sweep_plot.setLabel('left', 'Signal')

        # Right y-axis for contrast
        self._contrast_vb = pg.ViewBox()
        self.sweep_plot.scene().addItem(self._contrast_vb)
        self.sweep_plot.getAxis('right').linkToView(self._contrast_vb)
        self._contrast_vb.setXLink(self.sweep_plot)
        self.sweep_plot.showAxis('right')
        self.sweep_plot.setLabel('right', 'Contrast')

        # Keep right axis in sync on resize
        def _update_contrast_vb():
            self._contrast_vb.setGeometry(
                self.sweep_plot.getViewBox().sceneBoundingRect())
        self.sweep_plot.getViewBox().sigResized.connect(_update_contrast_vb)

        self.raw_curve = self.raw_plot.plot(pen=_pen(3, width=1), clear=True)
        self._curves.clear()

    def _on_camera_plot(self, inner_idx, val, i_run, oc_idx, oc_str,
                        processed, frame):
        """Handle camera experiment callback — update 3 panels.

        Layout: Image (top) | Sweep dual-axis (middle) | Raw trace (bottom)

        Convention from CameraExperiment:
            inner_idx == -1  → focus frame (processed=None)
            inner_idx >= 0   → acquisition point

        Multi-run overlay: previous runs stay on the plot with reduced
        alpha, matching the diode sweep behavior.
        """
        is_focus = (inner_idx < 0)

        if is_focus:
            self.status_label.setText(
                f"🔍 Focus  Run {i_run+1}/{self._Nruns}  "
                f"{'  '+oc_str if oc_str else ''}")
        else:
            self.status_label.setText(
                f"Outer {oc_idx+1}  Run {i_run+1}/{self._Nruns}  "
                f"Pt {inner_idx+1}/{len(self._inner_vals)}"
                f"{'  '+oc_str if oc_str else ''}")

        if not self.chk_rt.isChecked():
            return

        # Image panel — update for both focus and acquisition
        if frame is not None and getattr(self, '_cam_image', None) is not None:
            # pyqtgraph ImageItem wants (width, height) = frame.T
            self._cam_image.setImage(frame.T, autoLevels=True)
            if getattr(self, '_cam_raw_curve', None) is not None:
                self._cam_raw_curve.setData(frame.ravel()[:2000])

        # Intensity + Contrast panels — only during acquisition
        if not is_focus and processed is not None:
            mean_sig, mean_ref, contrast = processed
            n = len(mean_sig)
            xs = self._inner_vals[:n] / self._xu

            display_contrast = self._compute_display_contrast(
                mean_sig, mean_ref)

            alpha = 80 if self._Nruns > 1 else 220
            iplot = getattr(self, '_intensity_plot', None)
            cvb = getattr(self, '_cam_contrast_vb', None)
            if iplot is None:
                return

            # ── Signal curve (left axis) ──
            k_sig = f"cam_sig_o{oc_idx}_r{i_run}"
            if k_sig not in self._curves:
                self._curves[k_sig] = iplot.plot(
                    pen=_pen(0, alpha, 1.5),
                    name=f"S o{oc_idx}" if i_run == 0 else None)
            self._curves[k_sig].setData(xs, mean_sig)

            # ── Reference curve (left axis, dashed) ──
            k_ref = f"cam_ref_o{oc_idx}_r{i_run}"
            if k_ref not in self._curves:
                self._curves[k_ref] = iplot.plot(
                    pen=_pen(1, alpha, 1.5, dash=True),
                    name=f"R o{oc_idx}" if i_run == 0 else None)
            self._curves[k_ref].setData(xs, mean_ref)

            # ── Contrast curve (right axis) ──
            if cvb is not None:
                k_con = f"cam_con_o{oc_idx}_r{i_run}"
                if k_con not in self._curves:
                    crv = pg.PlotDataItem(
                        pen=_pen(2, alpha, 1.5, dash=True),
                        name=f"C o{oc_idx}" if i_run == 0 else None)
                    cvb.addItem(crv)
                    self._curves[k_con] = crv
                self._curves[k_con].setData(xs, display_contrast)
                cvb.autoRange()

    def closeEvent(self, ev):
        self._save_session_state()

        if self.is_running:
            r = QMessageBox.question(self, "Running", "Stop and quit?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if r == QMessageBox.StandardButton.No: ev.ignore(); return
            # Stop whatever is running
            if hasattr(self, '_ts_exp') and self._ts_exp and self._ts_exp.is_running:
                self._ts_stop()
            if self.experiment_thread:
                self.experiment_thread.request_stop()
                self.experiment_thread.wait(3000)
        try:
            # Close camera viewer
            if self._camera_viewer is not None:
                try:
                    self._camera_viewer.close()
                except Exception:
                    pass
                self._camera_viewer = None

            # Uninit camera worker if owned by session
            if hasattr(self, '_camera_worker') and self._camera_worker is not None:
                try:
                    self._camera_worker.uninit_cam()
                except Exception:
                    pass
                self._camera_worker = None

            # Terminate HCImageLive if running
            if self._hci_process is not None and self._hci_process.poll() is None:
                self._hci_process.terminate()
                self._hci_process = None
            if self._hci_watcher is not None:
                self._hci_watcher.stop()

            if self.ao_task:
                self.ao_task.set_outputs_to_constant([0, 0, 0])
                time.sleep(0.2); self.ao_task.__del__()
            if self.sg and hasattr(self.sg, 'uninit'): self.sg.uninit()
            if self.pb:
                try: self.pb.stop_sequence()
                except: pass
                self.pb.closePB()
            # Shutdown IPython kernel
            if self._ipy_kernel is not None:
                try:
                    self._ipy_kernel.shutdown_kernel()
                except Exception:
                    pass
        except Exception as e:
            print(f"Cleanup: {e}")
        ev.accept()

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
    app = QApplication.instance() or QApplication(sys.argv)
    # Apply dark theme to entire application
    set_dark_theme(app)

    win = SessionManager()
    win.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
