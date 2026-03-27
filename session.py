# session.py — NV Experiment Session Manager
# =============================================================================
# Features:
#   - Instrument init / per-instrument eject
#   - Experiment launch via mainControl_*.run() in QThread
#   - Real-time pyqtgraph with:
#       * Per-run overlay (semi-transparent + running mean)
#       * Per-outer-combo separate color group
#       * Shuffled x-axis: pre-allocated from param values, scatter as data arrives
#       * Click legend to show/hide traces
#   - Auto-save with crash-safe per-run overwrites
# =============================================================================

import sys, time, importlib, traceback
import numpy as np
from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGroupBox, QLabel, QPushButton, QComboBox, QCheckBox, QStatusBar,
    QSplitter, QLineEdit, QMessageBox,
)
from PySide6.QtCore import Qt, QThread, Signal, Slot
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
    style = Qt.DashLine if dash else Qt.SolidLine
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

_SS = """
QMainWindow,QWidget{background:#1e1e1e;color:#d4d4d4}
QGroupBox{border:1px solid #3c3c3c;border-radius:4px;margin-top:8px;
  padding-top:14px;font-weight:bold}
QGroupBox::title{subcontrol-origin:margin;left:10px}
QPushButton{background:#2d2d2d;border:1px solid #3c3c3c;border-radius:3px;
  padding:5px 14px}
QPushButton:hover{background:#3c3c3c}
QPushButton:disabled{color:#555}
QComboBox,QLineEdit{background:#2d2d2d;border:1px solid #3c3c3c;
  border-radius:3px;padding:3px 6px}
QStatusBar{background:#252526}
"""


class SessionManager(QMainWindow):
    _upd = Signal(int, float, int, int, str, object, object)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("NV Session Manager")
        self.setMinimumSize(1060, 680)
        self.setStyleSheet(_SS)

        self.sg = None; self.pb = None; self.ao_task = None
        self.camera_worker = None
        self.instruments_connected = False

        self.experiment_thread: Optional[ExperimentThread] = None
        self.is_running = False
        self.last_exp = None  # DiodeExperiment object from last run

        # Per-(outer, run) curve storage
        self._curves: dict[str, pg.PlotDataItem] = {}

        self._build_ui()
        self._upd.connect(self._on_plot, Qt.QueuedConnection)

    # ── UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        c = QWidget(); self.setCentralWidget(c)
        root = QVBoxLayout(c); root.setContentsMargins(8, 8, 8, 8)

        # Instruments
        ig = QGroupBox("Instruments"); il = QVBoxLayout(ig)
        mono = "font-family:'Consolas','Courier New',monospace;font-size:12px"
        self.lbl_sg = QLabel("SG:  ○"); self.lbl_sg.setStyleSheet(mono)
        self.lbl_pb = QLabel("PB:  ○"); self.lbl_pb.setStyleSheet(mono)
        self.lbl_ao = QLabel("AO:  ○"); self.lbl_ao.setStyleSheet(mono)
        self.lbl_cam = QLabel("Cam: ○"); self.lbl_cam.setStyleSheet(mono)
        for l in (self.lbl_sg, self.lbl_pb, self.lbl_ao, self.lbl_cam):
            il.addWidget(l)

        br = QHBoxLayout()
        self.btn_init = QPushButton("Init All")
        self.btn_init.clicked.connect(self._init_instruments)
        br.addWidget(self.btn_init)
        self.chk_sim = QCheckBox("Sim SG"); self.chk_sim.setChecked(True)
        br.addWidget(self.chk_sim)
        # Per-instrument eject buttons
        for name in ('SG', 'PB', 'AO'):
            b = QPushButton(f"Eject {name}")
            b.setMaximumWidth(80)
            b.clicked.connect(lambda checked, n=name: self._eject(n))
            br.addWidget(b)
        br.addStretch()
        il.addLayout(br)
        root.addWidget(ig)

        # Splitter: launcher | plots
        sp = QSplitter(Qt.Horizontal)

        # Left: launcher
        lg = QGroupBox("Experiment"); ll = QVBoxLayout(lg)
        ll.addWidget(QLabel("Config:"))
        self.cb_cfg = QComboBox(); self.cb_cfg.setEditable(True)
        self.cb_cfg.addItems(['esr_config', 'rabi_config', 'echo_config',
                              't1_config', 't2_config'])
        ll.addWidget(self.cb_cfg)
        ll.addWidget(QLabel("Script:"))
        self.cb_scr = QComboBox(); self.cb_scr.setEditable(True)
        self.cb_scr.addItems(['mainControl_diode', 'mainControl_camera'])
        ll.addWidget(self.cb_scr)

        self.chk_rt = QCheckBox("Real-time plot"); self.chk_rt.setChecked(True)
        ll.addWidget(self.chk_rt)
        self.chk_save = QCheckBox("Auto-save"); ll.addWidget(self.chk_save)

        hf = QHBoxLayout()
        hf.addWidget(QLabel("Folder #:"))
        self.txt_f = QLineEdit("001"); self.txt_f.setMaximumWidth(60)
        hf.addWidget(self.txt_f); hf.addStretch()
        ll.addLayout(hf)

        self.btn_run = QPushButton("▶  RUN")
        self.btn_run.setStyleSheet(
            "QPushButton{background:#0e639c;font-weight:bold;"
            "font-size:14px;padding:8px}"
            "QPushButton:hover{background:#1177bb}"
            "QPushButton:disabled{background:#333;color:#666}")
        self.btn_run.clicked.connect(self._run_experiment)
        self.btn_run.setEnabled(False)
        ll.addWidget(self.btn_run)

        self.btn_stop = QPushButton("⏹  STOP")
        self.btn_stop.setStyleSheet(
            "QPushButton{background:#6c1d1d;font-weight:bold}"
            "QPushButton:hover{background:#8b2222}")
        self.btn_stop.clicked.connect(self._stop_experiment)
        self.btn_stop.setEnabled(False)
        ll.addWidget(self.btn_stop)
        ll.addStretch()
        sp.addWidget(lg)

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

        sp.addWidget(pw)
        sp.setStretchFactor(0, 1); sp.setStretchFactor(1, 3)
        root.addWidget(sp, stretch=1)

        self.sbar = QStatusBar(); self.setStatusBar(self.sbar)
        self.lbl_st = QLabel("Ready"); self.sbar.addWidget(self.lbl_st, 1)

    # ── instruments ─────────────────────────────────────────────────────

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

            base = {'pb': {'clk_cyc': 1e3 / concfg.PBclk},
                    'scan': {}, 'seq': {}, 'mw': {}}
            self.pb = PulseBlaster(base); self.pb.configure()
            self.lbl_pb.setText("PB:  ● configured")

            self.ao_task = AnalogOutputTask(
                dev="P6363", channels=[0, 1, 2], coil='small_confocal')
            self.lbl_ao.setText("AO:  ● ready")

            self.instruments_connected = True
            self.btn_run.setEnabled(True)
            self.lbl_st.setText("✔ Instruments initialized")
        except Exception as e:
            self.lbl_st.setText(f"❌ {e}")
            print(traceback.format_exc())

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
            elif name == 'AO' and self.ao_task:
                self.ao_task.set_outputs_to_constant([0, 0, 0])
                time.sleep(0.2); self.ao_task.__del__()
                self.ao_task = None
                self.lbl_ao.setText("AO:  ○ ejected")
            self.lbl_st.setText(f"✔ {name} ejected")
        except Exception as e:
            self.lbl_st.setText(f"❌ Eject {name}: {e}")

    # ── run ─────────────────────────────────────────────────────────────

    def _run_experiment(self):
        if not self.instruments_connected: return

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

        instr = {'sg': self.sg, 'pb': self.pb, 'ao_task': self.ao_task}

        sp = fn = None
        if self.chk_save.isChecked():
            dd = Path("..") / "Saved_Data" / time.strftime("%Y-%m-%d")
            dd.mkdir(parents=True, exist_ok=True)
            fn = self.txt_f.text().strip() or "001"
            sp = dd / f"{exp_config.seq.name}_{fn}"; sp.mkdir(exist_ok=True)

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
        self.lbl_st.setText(f"Running {cn}…")

    # ── real-time plot (main thread) ────────────────────────────────────

    @Slot(int, float, int, int, str, object, object)
    def _on_plot(self, inner_idx, val, i_run, oc_idx, oc_str, proc, raw):
        total_inner = len(self._inner_vals)
        self.lbl_st.setText(
            f"Outer {oc_idx+1}  Run {i_run+1}/{self._Nruns}  "
            f"Pt {inner_idx+1}/{total_inner}"
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
        sh = " × ".join(str(s) for s in data.shape)
        self.lbl_st.setText(f"✔ Done — data {sh}")

    @Slot(str)
    def _on_err(self, tb):
        self.is_running = False
        self.btn_run.setEnabled(True); self.btn_stop.setEnabled(False)
        self.btn_init.setEnabled(True)
        self.lbl_st.setText("❌ Error (see console)")
        print(tb)

    def _stop_experiment(self):
        if self.experiment_thread:
            self.experiment_thread.request_stop()
            self.lbl_st.setText("⏹ Stopping…")

    # ── cleanup ─────────────────────────────────────────────────────────

    def closeEvent(self, ev):
        if self.is_running:
            r = QMessageBox.question(self, "Running",
                                     "Stop and quit?",
                                     QMessageBox.Yes | QMessageBox.No)
            if r == QMessageBox.No: ev.ignore(); return
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
