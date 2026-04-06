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
#   - Default autosave ON, per-run overwrite (crash-safe)
#   - View Sequence: matplotlib popup (needs PB only)
#   - Play Sequence: program + start PB at last scan point (needs PB only)
#   - RUN: experiment in QThread with real-time pyqtgraph
#   - Save: Saved_Data/YYYY-MM-DD/meas_NNN/ with data + metadata
# =============================================================================

import sys, time, importlib, traceback, subprocess
from typing import Optional, List, Dict, Tuple
import numpy as np
from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QComboBox, QCheckBox, QStatusBar,
    QLineEdit, QMessageBox, QTabWidget, QTreeWidget,
    QTreeWidgetItem, QHeaderView, QButtonGroup, QGridLayout, QGroupBox
)
from PySide6.QtCore import Qt, QThread, Signal, Slot, QTimer, QPropertyAnimation, Property, QEasingCurve
from PySide6.QtGui import QColor, QPainter, QBrush
import pyqtgraph as pg

from SGcontrol import SignalGenerator, SignalGenerator_sim
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogOutputTask
import connectionConfig as concfg

# ── palette ─────────────────────────────────────────────────────────────────
_PAL = ['#4ec9b0', '#569cd6', '#dcdcaa', '#ce9178', '#c586c0',
        '#9cdcfe', '#d7ba7d', '#b5cea8', '#f44747', '#6a9955']

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
        self.setFixedSize(50, 28)
        self._knob_position = 4  # starting position
        self.animation = QPropertyAnimation(self, b"knob_position")
        self.animation.setDuration(150)
        self.animation.setEasingCurve(QEasingCurve.Type.InOutCubic)
        self.stateChanged.connect(self._start_animation)
    
    def hitButton(self, pos):
        return self.rect().contains(pos)
    
    def _start_animation(self, state):
        self.animation.setStartValue(self._knob_position)
        self.animation.setEndValue(26 if state else 4)
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
        p.drawRoundedRect(0, 0, self.width(), self.height(), 14, 14)
        
        # Knob
        p.setBrush(QBrush(QColor("white")))
        p.drawEllipse(int(self._knob_position), 4, 20, 20)


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
QMainWindow,QWidget{background:#1e1e1e;color:#d4d4d4}
QGroupBox{border:1px solid #3c3c3c;border-radius:4px;margin-top:8px;
  padding-top:14px;font-weight:bold}
QGroupBox::title{subcontrol-origin:margin;left:10px}
QPushButton{background:#2d2d2d;border:1px solid #3c3c3c;border-radius:3px;
  padding:5px 14px}
QPushButton:hover{background:#3c3c3c}
QPushButton:disabled{color:#555}
QPushButton:checked { background-color: #0078d4; }
QComboBox,QLineEdit{background:#2d2d2d;border:1px solid #3c3c3c;
  border-radius:3px;padding:3px 6px}
QStatusBar{background:#252526}
QTabWidget::pane{border:1px solid #3c3c3c;background:#1e1e1e}
QTabBar::tab{background:#2d2d2d;border:1px solid #3c3c3c;
  padding:4px 10px;margin-right:2px;border-top-left-radius:3px;
  border-top-right-radius:3px}
QTabBar::tab:selected{background:#1e1e1e;border-bottom-color:#1e1e1e}
QTreeWidget{background:#252526;border:none;color:#d4d4d4}
QTreeWidget::item:alternate{background:#2a2a2a}
"""


class SessionManager(QMainWindow):
    _upd = Signal(int, float, int, int, str, object, object)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("NV Session Manager")
        self.setMinimumSize(1100, 720)
        self.setStyleSheet(_SS)

        self.sg = None; self.pb = None; self.ao_task = None
        self.camera_worker = None
        self.instruments_connected = False
        self._pb_ready = False

        self.experiment_thread: Optional[ExperimentThread] = None
        self.is_running = False
        self.last_exp = None  # DiodeExperiment object from last run

        # Per-(outer, run) curve storage
        self._curves: dict[str, pg.PlotDataItem] = {}

        self._loaded_config = None
        self._loaded_config_name = ""

        self._build_ui()
        self._upd.connect(self._on_plot, Qt.ConnectionType.QueuedConnection)

    # ── UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        c = QWidget(); self.setCentralWidget(c)
        root = QHBoxLayout(c); root.setContentsMargins(4, 4, 4, 4)

        # Left: tabs
        self._tabs = QTabWidget()
        self._tabs.setMaximumWidth(340); self._tabs.setMinimumWidth(260)
        self._tabs.addTab(self._build_instruments_tab(), "Instruments")
        self._tabs.addTab(self._build_experiment_tab(), "Experiment")
        self._tabs.addTab(self._build_metadata_tab(), "Config")
        root.addWidget(self._tabs)

        # Right: plots
        pw = QWidget(); pv = QVBoxLayout(pw)
        pv.setContentsMargins(0, 0, 0, 0)
        self.sweep_plot = pg.PlotWidget(title="Sweep")
        self.sweep_plot.addLegend(offset=(10, 10))
        self.sweep_plot.setLabel('bottom', 'Parameter')
        self.sweep_plot.setLabel('left', 'Signal')
        self.sweep_plot.showGrid(x=True, y=True, alpha=0.3)
        pv.addWidget(self.sweep_plot, stretch=2)
        self.raw_plot = pg.PlotWidget(title="Raw Trace")
        self.raw_plot.setLabel('bottom', 'Sample')
        self.raw_plot.showGrid(x=True, y=True, alpha=0.3)
        self.raw_curve = self.raw_plot.plot(pen=_pen(3, width=1))
        pv.addWidget(self.raw_plot, stretch=1)
        root.addWidget(pw, stretch=1)

        self.sbar = QStatusBar(); self.setStatusBar(self.sbar)
        self.lbl_st = QLabel("Ready"); self.sbar.addWidget(self.lbl_st, 1)

    def _build_instruments_tab(self):
        w = QWidget(); vl = QVBoxLayout(w); vl.setContentsMargins(6,6,6,6)
        mono = "font-family:'Consolas','Courier New',monospace;font-size:12px"
        self.lbl_sg = QLabel("SG:  ○"); self.lbl_sg.setStyleSheet(mono)
        self.lbl_pb = QLabel("PB:  ○"); self.lbl_pb.setStyleSheet(mono)
        self.lbl_ao = QLabel("AO:  ○"); self.lbl_ao.setStyleSheet(mono)
        self.lbl_cam = QLabel("Cam: ○"); self.lbl_cam.setStyleSheet(mono)
        for l in (self.lbl_sg, self.lbl_pb, self.lbl_ao, self.lbl_cam):
            vl.addWidget(l)
        r1 = QHBoxLayout()
        self.btn_init = QPushButton("Init All")
        self.btn_init.clicked.connect(self._init_instruments)
        r1.addWidget(self.btn_init)
        self.chk_sim = QCheckBox("Sim SG"); self.chk_sim.setChecked(True)
        r1.addWidget(self.chk_sim)
        vl.addLayout(r1)
        r2 = QHBoxLayout()
        for nm in ('SG','PB','AO'):
            b = QPushButton(f"Eject {nm}"); b.setMaximumWidth(80)
            b.clicked.connect(lambda _, n=nm: self._eject(n))
            r2.addWidget(b)
        r2.addStretch()
        vl.addLayout(r2)
        self.btn_init_pb = QPushButton("Init PB Only")
        self.btn_init_pb.setStyleSheet(
            "QPushButton{background:#2d4a2d}QPushButton:hover{background:#3a5f3a}")
        self.btn_init_pb.clicked.connect(self._init_pb_only)
        vl.addWidget(self.btn_init_pb)

        # ── Tool launchers (subprocess — no parameter return) ────────
        vl.addWidget(QLabel(""))  # spacer
        lbl = QLabel("Tools (separate process):")
        lbl.setStyleSheet("color:#888;font-size:11px")
        vl.addWidget(lbl)
        tr = QHBoxLayout()
        self.btn_launch_scope = QPushButton("📊 DAQ Scope")
        self.btn_launch_scope.setToolTip("Launch scope_v2_main.py in a new process")
        self.btn_launch_scope.clicked.connect(
            lambda: self._launch_tool('scope_v2_main.py'))
        tr.addWidget(self.btn_launch_scope)
        self.btn_launch_pb_gui = QPushButton("🔧 PB GUI")
        self.btn_launch_pb_gui.setToolTip("Launch pyPBLV.py in a new process")
        self.btn_launch_pb_gui.clicked.connect(
            lambda: self._launch_tool('pyPBLV.py'))
        tr.addWidget(self.btn_launch_pb_gui)
        vl.addLayout(tr)

        vl.addStretch()
        return w

    def _build_experiment_tab(self):
        w = QWidget(); vl = QVBoxLayout(w); vl.setContentsMargins(6,6,6,6)
        vl.addWidget(QLabel("Config module:"))
        self.cb_cfg = QComboBox(); self.cb_cfg.setEditable(True)
        self.cb_cfg.addItems(['esr_config','rabi_config','echo_config',
                              't1_config','t2_config','ramsey_config'])
        vl.addWidget(self.cb_cfg)
        vl.addWidget(QLabel("Script:"))
        self.cb_scr = QComboBox(); self.cb_scr.setEditable(True)
        self.cb_scr.addItems(['mainControl_diode','mainControl_camera'])
        vl.addWidget(self.cb_scr)

        mr = QHBoxLayout()
        mr.addWidget(QLabel("Mode:"))
        # self.cb_mode = QComboBox()
        # self.cb_mode.addItems(['Sweep', 'Timeseries'])
        self.cb_mode = SegmentedButton(['Sweep', 'Timeseries'])
        mr.addWidget(self.cb_mode)
        mr.addStretch()
        vl.addLayout(mr)

        # group = QGroupBox("Processing")
        # group_layout = QGridLayout(group)
        # self.chk_rt = ToggleSwitch()
        # self.chk_rt.setChecked(True)

        # # self.chk_rt = QCheckBox("Real-time plot"); self.chk_rt.setChecked(True)
        # group_layout.addWidget(self.chk_rt, 0, 2)
        # group_layout.addWidget(QLabel("Real-time plot"), 0, 1)
        # vl.addWidget(group)

        chk_rt = QHBoxLayout()
        chk_rt.addWidget(QLabel("Real-time plot"))
        self.chk_rt = ToggleSwitch()
        self.chk_rt.setChecked(True)
        
        chk_rt.addWidget(self.chk_rt)
        vl.addLayout(chk_rt)

        chk_save = QHBoxLayout()
        chk_save.addWidget(QLabel("Auto-save"));# chk_save.addStretch()
        
        self.txt_f = QLineEdit("001"); self.txt_f.setMaximumWidth(50)
        chk_save.addWidget(self.txt_f); chk_save.addStretch()

        self.chk_save = ToggleSwitch(); self.chk_save.setChecked(True)
        chk_save.addWidget(self.chk_save)
        
        
        vl.addLayout(chk_save)

        self.btn_load_cfg = QPushButton("📋  Load Config")
        self.btn_load_cfg.setStyleSheet(
            "QPushButton{background:#3d3d1d;font-weight:bold}"
            "QPushButton:hover{background:#4d4d2d}")
        self.btn_load_cfg.clicked.connect(self._load_config)
        vl.addWidget(self.btn_load_cfg)

        self.btn_run = QPushButton("▶  RUN")
        self.btn_run.setStyleSheet(
            "QPushButton{background:#0e639c;font-weight:bold;"
            "font-size:14px;padding:8px}"
            "QPushButton:hover{background:#1177bb}"
            "QPushButton:disabled{background:#333;color:#666}")
        self.btn_run.clicked.connect(self._run_experiment)
        self.btn_run.setEnabled(False)
        vl.addWidget(self.btn_run)

        self.btn_stop = QPushButton("⏹  STOP")
        self.btn_stop.setStyleSheet(
            "QPushButton{background:#6c1d1d;font-weight:bold}"
            "QPushButton:hover{background:#8b2222}")
        self.btn_stop.clicked.connect(self._stop_experiment)
        self.btn_stop.setEnabled(False)
        vl.addWidget(self.btn_stop)

        sr = QHBoxLayout()
        self.btn_view_seq = QPushButton("🔍 View")
        self.btn_view_seq.setToolTip("View Sequence (needs PB)")
        self.btn_view_seq.setStyleSheet(
            "QPushButton{background:#2d5f2d}QPushButton:hover{background:#3a7a3a}")
        self.btn_view_seq.clicked.connect(self._view_sequence)
        self.btn_view_seq.setEnabled(False)
        sr.addWidget(self.btn_view_seq)
        self.btn_play_seq = QPushButton("▶ Play")
        self.btn_play_seq.setToolTip("Play Sequence on PB (needs PB)")
        self.btn_play_seq.setStyleSheet(
            "QPushButton{background:#2d4a5f}QPushButton:hover{background:#3a6a7a}")
        self.btn_play_seq.clicked.connect(self._play_sequence)
        self.btn_play_seq.setEnabled(False)
        sr.addWidget(self.btn_play_seq)
        vl.addLayout(sr)
        vl.addStretch()
        return w
    
    # def _on_acquisition_mode_changed(self, idx: int):
    #     """Handle display mode change (Time/FFT)"""
    #     self.fft_mode = idx == 1
    #     if self.fft_mode:
    #         self.scope_dock.hide()
    #         self.fft_dock.show()
    #     else:
    #         self.fft_dock.hide()
    #         self.scope_dock.show()
    #     if self.display_mgr:
    #         self.display_mgr.set_fft_enabled(self.fft_mode)
    #     # Enable/disable FFT controls
    #     self.fft_tab.set_fft_mode_enabled(self.fft_mode)
    #     self.error_label.set_msg(f"Display mode: {'FFT' if self.fft_mode else 'Time'}", error=False)

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

    # ════════════════════════════════════════════════════════════════
    #  INSTRUMENTS
    # ════════════════════════════════════════════════════════════════

    def _init_instruments(self):
        try:
            self.lbl_st.setText("Initializing…")
            QApplication.processEvents()

            if self.chk_sim.isChecked():
                self.sg = SignalGenerator_sim()
                self.lbl_sg.setText("SG:  ● simulated")
            else:
                self.sg = SignalGenerator()
                self.sg.query_all()
                self.lbl_sg.setText(f"SG:  ● {self.sg.freq/1e9:.4f} GHz")
            if not self._pb_ready:
                self._do_init_pb()
            self.ao_task = AnalogOutputTask(
                dev="P6363", channels=[0, 1, 2], coil='small_confocal')
            self.lbl_ao.setText("AO:  ● ready")

            self.instruments_connected = True
            self.btn_run.setEnabled(True)
            self._update_seq_btns()
            self.lbl_st.setText("✔ Instruments initialized")
        except Exception as e:
            self.lbl_st.setText(f"❌ {e}")
            print(traceback.format_exc())

    def _init_pb_only(self):
        try:
            self.lbl_st.setText("Initializing PB…"); QApplication.processEvents()
            self._do_init_pb()
            self._update_seq_btns()
            self.lbl_st.setText("✔ PB ready — View/Play enabled")
        except Exception as e:
            self.lbl_st.setText(f"❌ PB: {e}"); print(traceback.format_exc())

    def _do_init_pb(self):
        base = {'pb':{'clk_cyc':1e3/concfg.PBclk},'scan':{},'seq':{},'mw':{}}
        self.pb = PulseBlaster(base); self.pb.configure()
        self.lbl_pb.setText("PB:  ● configured")
        self._pb_ready = True

    def _update_seq_btns(self):
        self.btn_view_seq.setEnabled(self._pb_ready)
        self.btn_play_seq.setEnabled(self._pb_ready)

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
            self.lbl_st.setText(f"✔ {name} ejected")
        except Exception as e:
            self.lbl_st.setText(f"❌ Eject {name}: {e}")

    # ════════════════════════════════════════════════════════════════
    #  TOOL LAUNCHERS (subprocess — fire and forget)
    # ════════════════════════════════════════════════════════════════

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
            self.lbl_st.setText(f"✔ Launched {script_name} (PID {proc.pid})")
        except Exception as e:
            self.lbl_st.setText(f"❌ Launch {script_name}: {e}")

    # ════════════════════════════════════════════════════════════════
    #  LOAD CONFIG
    # ════════════════════════════════════════════════════════════════

    def _load_config(self):
        cn = self.cb_cfg.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
        except Exception as e:
            self.lbl_st.setText(f"❌ Import: {e}"); return
        if not hasattr(cm,'config'):
            self.lbl_st.setText(f"❌ {cn} has no 'config'"); return
        self._loaded_config = cm.config
        self._loaded_config_name = cn
        self._update_metadata_tree(cm.config)
        self._tabs.setCurrentIndex(2)
        self.lbl_st.setText(f"✔ Loaded {cn}")

    def _ensure_config(self):
        cn = self.cb_cfg.currentText()
        if self._loaded_config is None or self._loaded_config_name != cn:
            self._load_config()
        return self._loaded_config

    # ── run ─────────────────────────────────────────────────────────────

    def _run_experiment(self):
        if not self.instruments_connected: return
        mode = self.cb_mode.currentText()
        if mode == 'Timeseries':
            self._run_timeseries()
        else:
            self._run_sweep()

    def _run_sweep(self):
        cn = self.cb_cfg.currentText(); sn = self.cb_scr.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
            sm = importlib.import_module(sn); importlib.reload(sm)
        except Exception as e:
            self.lbl_st.setText(f"❌ Import: {e}"); return

        if not hasattr(sm, 'run'):
            self.lbl_st.setText(f"❌ {sn} has no run()"); return

        # Get the ExperimentConfig object from config module
        if not hasattr(cm, 'config'):
            self.lbl_st.setText(f"❌ {cn} has no 'config' attribute"); return
        exp_config = cm.config
        self._loaded_config = exp_config
        self._loaded_config_name = cn
        self._update_metadata_tree(exp_config)

        instr = {'sg': self.sg, 'pb': self.pb, 'ao_task': self.ao_task}

        sp = fn = None
        if self.chk_save.isChecked():
            dd = Path("..") / "Saved_Data" / time.strftime("%Y-%m-%d")
            dd.mkdir(parents=True, exist_ok=True)
            fn = self.txt_f.text().strip() or "001"
            sp = dd / f"meas_{fn}"; sp.mkdir(exist_ok=True)

        # Clear plots
        self.sweep_plot.clear()
        self.sweep_plot.addLegend(offset=(10, 10))
        self.raw_curve = self.raw_plot.plot(pen=_pen(3, width=1), clear=True)
        self._curves.clear()

        # Store config ref for plotting — use ExperimentConfig directly
        self._exp_config = exp_config
        primary_name = exp_config.scan_names[0]
        self._inner_vals = exp_config.scans[primary_name].values
        self._xu = exp_config.plot.x_units
        self._Nruns = exp_config.runtime.Nruns
        self.sweep_plot.setLabel('bottom', exp_config.plot.x_label)

        # Thread
        self.experiment_thread = ExperimentThread(
            sm.run, instr, exp_config, sp, fn)
        self.experiment_thread.point_signal.connect(
            lambda *a: self._upd.emit(*a))
        self.experiment_thread.finished_signal.connect(self._on_done)
        self.experiment_thread.error_signal.connect(self._on_err)
        self.experiment_thread.start()

        self.is_running = True
        self.btn_run.setEnabled(False); self.btn_stop.setEnabled(True)
        self.btn_init.setEnabled(False)
        self.btn_view_seq.setEnabled(False); self.btn_play_seq.setEnabled(False)
        self.lbl_st.setText(f"Running {cn}…")

    # ════════════════════════════════════════════════════════════════
    #  RUN TIMESERIES
    # ════════════════════════════════════════════════════════════════

    def _run_timeseries(self):
        cn = self.cb_cfg.currentText()
        try:
            cm = importlib.import_module(cn); importlib.reload(cm)
        except Exception as e:
            self.lbl_st.setText(f"❌ Import: {e}"); return
        if not hasattr(cm, 'config'):
            self.lbl_st.setText(f"❌ {cn} has no 'config'"); return

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
            self.lbl_st.setText(f"❌ TS setup: {e}")
            print(traceback.format_exc()); return

        # Switch plots to timeseries mode
        self.sweep_plot.clear()
        self.sweep_plot.setTitle("Contrast Timeseries")
        self.sweep_plot.setLabel('bottom', 'Time', units='s')
        self.sweep_plot.setLabel('left', 'Contrast')
        self._ts_contrast_curve = self.sweep_plot.plot(
            pen=_pen(0, 220, 1.5))
        self.raw_curve = self.raw_plot.plot(
            pen=_pen(3, width=1), clear=True)

        # 30 Hz display refresh
        self._ts_display_timer = QTimer()
        self._ts_display_timer.timeout.connect(self._ts_refresh)
        self._ts_display_timer.start(33)

        # Save path
        sp = None
        if self.chk_save.isChecked():
            fn = self.txt_f.text().strip() or "001"
            dd = Path("..") / "Saved_Data" / time.strftime("%Y-%m-%d")
            dd.mkdir(parents=True, exist_ok=True)
            sp = dd / f"meas_{fn}"; sp.mkdir(exist_ok=True)
        self._ts_save_path = sp

        self._ts_exp.start()

        self.is_running = True
        self.btn_run.setEnabled(False); self.btn_stop.setEnabled(True)
        self.btn_init.setEnabled(False)
        self.btn_view_seq.setEnabled(False); self.btn_play_seq.setEnabled(False)
        self.lbl_st.setText(
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
        self.lbl_st.setText(
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
                fn = self.txt_f.text().strip() or "001"
                self._ts_exp.save(
                    self._ts_save_path / f"ts_{fn}", fmt='npz')
            self._ts_exp.teardown()

        self.is_running = False
        self.btn_run.setEnabled(True); self.btn_stop.setEnabled(False)
        self.btn_init.setEnabled(True); self._update_seq_btns()
        self.lbl_st.setText("✔ Timeseries done")

    # ════════════════════════════════════════════════════════════════
    #  VIEW SEQUENCE (PB only)
    # ════════════════════════════════════════════════════════════════

    def _view_sequence(self):
        if not self._pb_ready:
            self.lbl_st.setText("❌ PB not initialized"); return
        cfg = self._ensure_config()
        if cfg is None: return

        instr = {'sg':self.sg,'pb':self.pb,'ao_task':self.ao_task}
        from experiment_base import DiodeExperiment
        import matplotlib.pyplot as plt
        try:
            exp = DiodeExperiment(instr, cfg)
            idx = getattr(cfg.runtime,'seq_plot_indices',[0,-1])
            exp.view_sequences(indices=idx, dpi=100)
            plt.show()
            self.lbl_st.setText("✔ Sequence plots shown")
        except Exception as e:
            self.lbl_st.setText(f"❌ View: {e}"); print(traceback.format_exc())

    # ════════════════════════════════════════════════════════════════
    #  PLAY SEQUENCE (PB only)
    # ════════════════════════════════════════════════════════════════

    def _play_sequence(self):
        """Program + start PB at last seq_plot_indices value."""
        if not self._pb_ready:
            self.lbl_st.setText("❌ PB not initialized"); return
        cfg = self._ensure_config()
        if cfg is None: return

        try:
            args_names = cfg.seq.args_names
            seq_args = list(cfg.seq.args_values)
            sn = cfg.scan_names[0] if cfg.scan_names else ''
            pv = cfg.scans[sn].values if sn else np.array([0])

            indices = getattr(cfg.runtime,'seq_plot_indices',[0,-1])
            pi = indices[-1]
            if pi < 0: pi = len(pv) + pi
            val = pv[pi]

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
            self.lbl_st.setText(f"▶ PB: {cfg.seq.name} @ {sn}={val:.4g}")
        except Exception as e:
            self.lbl_st.setText(f"❌ Play: {e}"); print(traceback.format_exc())

    # ── real-time plot (main thread) ────────────────────────────────────

    @Slot(int, float, int, int, str, object, object)
    def _on_plot(self, inner_idx, val, i_run, oc_idx, oc_str, proc, raw):
        self.lbl_st.setText(
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

            if len(proc) >= 3:
                key_con = f"con_o{oc_idx}_r{i_run}"
                if key_con not in self._curves:
                    nm = f"S/R o{oc_idx}" if i_run == 0 else None
                    self._curves[key_con] = self.sweep_plot.plot(
                        pen=_pen(color_base + 1, alpha, 1.5, dash=True),
                        name=nm)
                self._curves[key_con].setData(
                    xs_all[:n], proc[2][:n, 0])

        # Raw
        if raw is not None:
            self.raw_curve.setData(np.ravel(np.array(raw)))

    @Slot(object, object)
    def _on_done(self, exp, data):
        self.is_running = False
        self.last_exp = exp  # accessible from IPython later
        self.btn_run.setEnabled(True); self.btn_stop.setEnabled(False)
        self.btn_init.setEnabled(True)
        self._update_seq_btns()

        sh = " × ".join(str(s) for s in data.shape)
        self.lbl_st.setText(f"✔ Done — data {sh}")

    @Slot(str)
    def _on_err(self, tb):
        self.is_running = False
        self.btn_run.setEnabled(True); self.btn_stop.setEnabled(False)
        self.btn_init.setEnabled(True)
        self._update_seq_btns()
        self.lbl_st.setText("❌ Error (see console)")
        print(tb)

    def _stop_experiment(self):
        if self.cb_mode.currentText() == 'Timeseries':
            self._ts_stop()
        elif self.experiment_thread:
            self.experiment_thread.request_stop()
            self.lbl_st.setText("⏹ Stopping…")

    # ── cleanup ─────────────────────────────────────────────────────────

    def closeEvent(self, ev):
        if self.is_running:
            r = QMessageBox.question(self, "Running",
                                     "Stop and quit?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if r == QMessageBox.StandardButton.No: ev.ignore(); return
            # Stop whatever is running
            if hasattr(self, '_ts_exp') and self._ts_exp and self._ts_exp.is_running:
                self._ts_stop()
            if self.experiment_thread:
                self.experiment_thread.request_stop()
                self.experiment_thread.wait(3000)
        try:
            if self.ao_task:
                self.ao_task.set_outputs_to_constant([0, 0, 0])
                time.sleep(0.2); self.ao_task.__del__()
            if self.sg and hasattr(self.sg, 'uninit'): self.sg.uninit()
            if self.pb:
                try: self.pb.stop_sequence()
                except: pass
                self.pb.closePB()
        except Exception as e:
            print(f"Cleanup: {e}")
        ev.accept()


def main():
    app = QApplication.instance() or QApplication(sys.argv)
    win = SessionManager()
    win.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
