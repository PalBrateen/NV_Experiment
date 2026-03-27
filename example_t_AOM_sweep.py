# example_t_AOM_sweep.py
"""
CASE B: Sweep t_AOM in ESR — each t_AOM value changes Nsamples,
so the AI task must be rebuilt each time.

This is a standalone script.  Run it after instruments are initialized
(e.g., from IPython after session.init_instruments(), or standalone).

What happens:
    t_AOM = 30 ms  →  Nsamples recalculated  →  setup()  →  full ESR sweep  →  teardown()
    t_AOM = 50 ms  →  Nsamples recalculated  →  setup()  →  full ESR sweep  →  teardown()
    t_AOM = 70 ms  →  Nsamples recalculated  →  setup()  →  full ESR sweep  →  teardown()

Each ESR sweep is a complete (outer × Nruns × inner) acquisition.
Results is a list of 3 data arrays (one per t_AOM), each potentially
with different daq_Nsamples dimension.
"""

import numpy as np
import time
from pathlib import Path
from spinapi import ms

# ── 1. Load config ──────────────────────────────────────────────────────
from esr_config import config     # your ExperimentConfig instance

# ── 2. Init instruments (or receive from session) ───────────────────────
from SGcontrol import SignalGenerator
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogOutputTask

pb_params = {'pb': {'clk_cyc': config.pb.clk_cyc_ns},
             'scan': {}, 'seq': {}, 'mw': {}}
sg = SignalGenerator()
pb = PulseBlaster(pb_params); pb.configure()
ao_task = AnalogOutputTask(dev="P6363", channels=[0, 1, 2],
                           coil='small_confocal')
instruments = {'sg': sg, 'pb': pb, 'ao_task': ao_task}

# ── 3. Create experiment object (ONCE) ──────────────────────────────────
from experiment_base import DiodeExperiment

exp = DiodeExperiment(instruments, config)

# ── 4. Define the config applier function ───────────────────────────────
#
# This function modifies the config IN-PLACE before each config_sweep
# iteration calls setup().  It must update everything that depends on
# the config parameter.
#
# For t_AOM in ESR:
#   - config.times.t_AOM          (the time value itself)
#   - config.seq.args_values[0]   (t_AOM is the 0th seq arg for esr_seq)
#   - config.seq.Nsamples         (depends on t_AOM via your formula)
#   - config.daq_ai.ai_samps_per_chan  (must match Nsamples)
#
# How to find the index:
#   config.seq.args_names = ['t_AOM']       ← index 0
#   config.seq.args_values = [70000000.0]   ← current value in ns

def apply_t_aom(val):
    """Called by config_sweep before each setup().  val is in ns (PB units)."""
    print(f"\n  Applying t_AOM = {val/ms:.1f} ms")

    # Update timing
    config.times.t_AOM = val

    # Update sequence args — find index by name
    idx = config.seq.args_names.index('t_AOM')
    config.seq.args_values[idx] = val

    # Recalculate Nsamples (your existing formula from esr_config.py)
    new_N = max(2, int((val - 40 * ms) / ms * 1e-3 * 10e3))
    config.seq.Nsamples = new_N
    config.daq_ai.ai_samps_per_chan = new_N

    print(f"  → Nsamples = {new_N}")


# ── 5. Run the config sweep ────────────────────────────────────────────
#
# config_sweep does this for each t_AOM value:
#   1. Calls apply_t_aom(val) to update config
#   2. Calls exp.teardown() to release old AI task
#   3. Re-reads Nsamples from updated config
#   4. Calls exp.setup() — builds new AI task, new Sweep, new data_array
#   5. Calls exp.execute() — runs the full ESR frequency sweep
#   6. Appends the data_array to results list

t_AOM_values = [30 * ms, 50 * ms, 70 * ms]   # values to sweep

try:
    results = exp.config_sweep(
        config_name='t_AOM',
        config_values=t_AOM_values,
        config_applier=apply_t_aom,
        # These kwargs go to exp.execute():
        # callback=None,        # or your real-time callback
        # stop_check=None,
        # save_path=Path("../Saved_Data/2026-03-25/t_AOM_sweep_001"),
        # folder_number="001",
    )

    # ── 6. Inspect results ──────────────────────────────────────────────
    for i, (t_val, data) in enumerate(zip(t_AOM_values, results)):
        print(f"\nt_AOM = {t_val/ms:.0f} ms: data shape = {data.shape}")
        # data.shape = (outer_combos, Nruns, inner_pts, daq_Nsamples * n_ch)
        # daq_Nsamples is DIFFERENT for each t_AOM value

finally:
    ao_task.set_outputs_to_constant([0, 0, 0])
    time.sleep(0.2); ao_task.__del__()
    pb.stop_sequence(); pb.closePB()
    print("\n✔ Done")
