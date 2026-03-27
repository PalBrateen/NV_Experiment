# mainControl_diode.py — v7 (ExperimentConfig native + validation + modes)
"""
Thin wrapper around DiodeExperiment.

Accepts the new ExperimentConfig dataclass from experiment_config.py.
Config modules must expose a `config` attribute (ExperimentConfig instance).

    from esr_config import config
    exp, data = run(instruments, config)

Execution modes are controlled by config.runtime:
    RuntimeFlags(reload_pb=True,  load_all_params=False)  → Normal
    RuntimeFlags(reload_pb=False, load_all_params=False)  → Fixed-PB (ESR)
    RuntimeFlags(reload_pb=True,  load_all_params=True)   → Batch (fastest)

Validation runs automatically during setup().  To pre-check without
running, use:
    exp = DiodeExperiment(instruments, config)
    result = exp.validate(plot_errors=True)
    if not result.all_valid:
        print(result.summary())
"""

import numpy as np
import time
from pathlib import Path
from typing import Optional, Callable, TYPE_CHECKING

from experiment_base import DiodeExperiment, read_details, process_data

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
    Run a diode-based NV experiment.

    Args:
        instruments: dict with 'sg', 'pb', 'ao_task'
        config:      ExperimentConfig dataclass instance
        callback:    callable(inner_idx, val, i_run, oc_idx, oc_vals,
                             processed, raw)
        stop_check:  callable() -> bool
        save_path:   Path or None
        folder_number: str or None

    Returns:
        (exp, data_array) — exp is the DiodeExperiment with all state.
        data shape: (outer_combos, Nruns, inner_pts, daq_Nsamples * n_ch)

    Invalid cells (failed PB timing validation) are NaN-filled.
    Access exp._validation for the full ValidationResult.
    """
    exp = DiodeExperiment(instruments, config)
    try:
        exp.setup()     # validate() called automatically inside setup()
        data = exp.execute(
            callback=callback,
            stop_check=stop_check,
            save_path=save_path,
            folder_number=folder_number,
        )
    finally:
        exp.teardown()
    return exp, data


def validate_only(
    instruments: dict,
    config: 'ExperimentConfig',
    plot_errors: bool = True,
):
    """Dry-run validation without executing.

    Useful for checking parameter grids in a notebook before committing
    to a long acquisition.

    Returns:
        (exp, ValidationResult)
    """
    exp = DiodeExperiment(instruments, config)
    result = exp.validate(plot_errors=plot_errors)
    return exp, result


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
    import connectionConfig as concfg

    psutil.Process(os.getpid()).cpu_affinity([0, 1])

    # ── load config (new format) ──
    cfg_mod = import_module('rabi_config')
    config = cfg_mod.config        # ExperimentConfig instance

    trial_sg = False   # True = simulated SG

    # ── instruments ──
    sg = SignalGenerator_sim() if trial_sg else SignalGenerator()
    pb_params = {'pb': {'clk_cyc': config.pb.clk_cyc_ns},
                 'scan': {}, 'seq': {}, 'mw': {}}
    pb = PulseBlaster(pb_params); pb.configure()
    ao_task = AnalogOutputTask(dev="P6363", channels=[0, 1, 2],
                               coil='small_confocal')

    if not trial_sg:
        sg.enable_ntype(1)
        sg.set_amp_rf(config.mw.power)
        sg.set_freq(config.mw.freq)
        sg.setup_sg_fm('external', dev=5e5)
        sg.enable_modulation(1)

    instruments = {'sg': sg, 'pb': pb, 'ao_task': ao_task}

    try:
        exp, data = run(instruments=instruments, config=config)

        # Plot
        iv = exp.sweep.inner_values
        xu = config.plot.x_units
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        for oc in range(data.shape[0]):
            for ir in range(data.shape[1]):
                proc = process_data(
                    data.shape[2], exp.reads_per_cyc, data[oc, ir],
                    exp.Nsamples_cfg, exp.sequence, exp.n_channels)
                axs[0].plot(iv / xu, proc[0][:, 0], '.-',
                            alpha=0.5, label=f'oc{oc} r{ir}')
                if len(proc) > 1:
                    axs[1].plot(iv / xu, proc[2][:, 0], '.-', alpha=0.5)
        axs[0].set_xlabel(config.plot.x_label)
        axs[1].set_xlabel(config.plot.x_label)
        axs[1].set_ylabel('Contrast')
        plt.tight_layout(); plt.show()

    finally:
        ao_task.set_outputs_to_constant([0, 0, 0])
        time.sleep(0.3); ao_task.__del__()
        if hasattr(sg, 'uninit'): sg.uninit()
        pb.stop_sequence(); pb.closePB()
        print("✔ All instruments closed")
