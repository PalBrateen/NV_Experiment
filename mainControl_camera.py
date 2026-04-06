# mainControl_camera.py — v3.0 (CameraExperiment + ExperimentConfig)
"""
Thin wrapper around CameraExperiment.

Accepts ExperimentConfig with CameraConfig + FieldConfig.
Config modules must expose a `config` attribute (ExperimentConfig instance).

    from camera_esr_config import config
    exp, data = run(instruments, config)

Instruments dict must include 'camera_worker' (from Camcontrol.CameraWorker).

Phase 1: cam_levelm mode, static B field, no AO patterns.
"""

import numpy as np
import time
from pathlib import Path
from typing import Optional, Callable, TYPE_CHECKING

from camera_experiment import CameraExperiment, process_camera_data

if TYPE_CHECKING:
    from experiment_config import ExperimentConfig


def run(
    instruments: dict,
    config: 'ExperimentConfig',
    callback: Optional[Callable] = None,
    stop_check: Optional[Callable] = None,
    save_path: Optional[Path] = None,
    folder_number: Optional[str] = None,
):
    """
    Run a camera-based NV experiment.

    Args:
        instruments: dict with 'sg', 'pb', 'ao_task', 'camera_worker'
        config:      ExperimentConfig with camera= and field_cfg= populated
        callback:    callable(inner_idx, val, i_run, oc_idx, oc_str,
                             processed, last_frame)
                     processed = (mean_sig, mean_ref, contrast) or None
                     last_frame = uint16 ndarray (vsize, hsize)
        stop_check:  callable() -> bool
        save_path:   Path or None
        folder_number: str or None

    Returns:
        (exp, data_array) — exp is the CameraExperiment with all state.
        data shape: (outer_combos, Nruns, inner_pts, Nframes, vsize, hsize)
    """
    exp = CameraExperiment(instruments, config)
    try:
        exp.setup()
        data = exp.execute(
            callback=callback,
            stop_check=stop_check,
            save_path=save_path,
            folder_number=folder_number,
        )
    finally:
        exp.teardown()
    return exp, data


# ═══════════════════════════════════════════════════════════════════════
# Standalone mode
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import sys, os, psutil
    import matplotlib.pyplot as plt
    from importlib import import_module
    from SGcontrol import SignalGenerator, SignalGenerator_sim
    from PBcontrol import PulseBlaster
    from DAQcontrol import AnalogOutputTask
    from Camcontrol import CameraWorker
    import connectionConfig as concfg

    psutil.Process(os.getpid()).cpu_affinity([0, 1])

    # ── load config ──
    cfg_mod = import_module('camera_esr_config')
    config = cfg_mod.config

    trial_sg = False
    trial_cam = False

    # ── instruments ──
    sg = SignalGenerator_sim() if trial_sg else SignalGenerator()
    pb_params = {'pb': {'clk_cyc': config.pb.clk_cyc_ns},
                 'scan': {}, 'seq': {}, 'mw': {}}
    pb = PulseBlaster(pb_params); pb.configure()
    ao_task = AnalogOutputTask(dev="P6363", channels=[0, 1, 2],
                               coil='small_confocal')

    camera_worker = CameraWorker(simulate=trial_cam)
    camera_worker.init_cam()

    if not trial_sg:
        sg.enable_ntype(1)
        sg.set_amp_rf(config.mw.power)
        sg.set_freq(config.mw.freq)
        sg.setup_sg_fm('external', dev=5e5)
        sg.enable_modulation(1)

    instruments = {
        'sg': sg, 'pb': pb, 'ao_task': ao_task,
        'camera_worker': camera_worker,
    }

    try:
        exp, data = run(instruments=instruments, config=config)

        # Plot
        iv = exp.sweep.inner_values
        xu = config.plot.x_units

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5))

        for oc in range(data.shape[0]):
            for ir in range(data.shape[1]):
                proc = process_camera_data(
                    data.shape[2], data[oc, ir])
                mean_sig, mean_ref, contrast = proc

                ax1.plot(iv / xu, mean_sig, '.-', alpha=0.5,
                         label=f'S oc{oc} r{ir}')
                ax1.plot(iv / xu, mean_ref, '.--', alpha=0.5,
                         label=f'R oc{oc} r{ir}')
                ax2.plot(iv / xu, contrast, '.-', alpha=0.5,
                         label=f'oc{oc} r{ir}')

        ax1.set_xlabel(config.plot.x_label)
        ax1.set_ylabel('Intensity (pixel sum)')
        ax1.legend(fontsize=8)
        ax2.set_xlabel(config.plot.x_label)
        ax2.set_ylabel('Contrast (S/R)')

        # Show last frame
        last_frame = data[0, -1, -1, 0, :, :]
        ax3.imshow(last_frame, cmap='inferno')
        ax3.set_title('Last signal frame')

        plt.tight_layout(); plt.show()

    finally:
        camera_worker.uninit_cam()
        ao_task.set_outputs_to_constant([0, 0, 0])
        time.sleep(0.3); ao_task.__del__()
        if hasattr(sg, 'uninit'):
            sg.enable_ntype(1); sg.uninit()
        pb.stop_sequence(); pb.closePB()
        print("✔ All instruments closed")
