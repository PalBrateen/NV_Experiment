# example_AOM_lag_rabi_sweep.py
"""
CASE A: Sweep AOM_lag in Rabi — only the PB pulse sequence changes,
DAQ shape stays the same.

AOM_lag is a timing parameter in the pulse sequence.  Changing it
shifts when the AOM fires relative to the MW pulse.  The PB needs
to be reprogrammed each scan point, but Nsamples doesn't change.

This uses the normal Sweep mechanism (not config_sweep).

For Rabi, seq.args are:
    args_names  = ['t_MW',  't_AOM', 'ro_delay', 'AOM_lag', 'MW_lag']
    args_values = [100*ns,  3*us,    300*ns,      800*ns,    150*ns  ]
                    idx 0    idx 1    idx 2        idx 3      idx 4

We want to sweep AOM_lag (index 3).

Two approaches shown below:
  A1. AOM_lag as the PRIMARY (inner) sweep — the only swept parameter
  A2. AOM_lag as an OUTER sweep with frequency as inner
"""

import numpy as np
import time
from spinapi import ns, us, ms

# ═══════════════════════════════════════════════════════════════════════
# Approach A1: AOM_lag as the ONLY swept parameter
# ═══════════════════════════════════════════════════════════════════════

def example_A1():
    """
    Simple 1D sweep of AOM_lag.

    Config file would look like:
        config = ExperimentConfig(
            scans={
                'AOM_lag': ScanAxis(
                    name='AOM_lag',
                    start=500 * ns,
                    stop=1500 * ns,
                    step=50 * ns,
                ),
            },
            seq=SequenceConfig(
                name='rabi_seq',
                args_names=['t_MW', 't_AOM', 'ro_delay', 'AOM_lag', 'MW_lag'],
                args_values=[100*ns, 3*us, 300*ns, 800*ns, 150*ns],
                Nsamples=100,
            ),
            ...
        )

    The key insight: AOM_lag is at index 3 in args_values.
    make_pb_setter(pb, param_index=3, ...) creates a closure that
    overwrites args[3] with the new AOM_lag value and reprograms PB.
    """
    from experiment_config import (
        ExperimentConfig, ScanAxis, SequenceConfig, MicrowaveConfig,
        PulseBlasterConfig, TimeConfig, DAQAIConfig, PlotConfig,
        SaveConfig, RuntimeFlags,
    )
    from connectionConfig import PBclk, laser, samp_clk, start_trig, MW

    channels = {'samp': samp_clk, 'mw': MW, 'laser': laser,
                'start': start_trig}

    config = ExperimentConfig(
        scans={
            'AOM_lag': ScanAxis(
                name='AOM_lag',
                start=500 * ns,
                stop=1500 * ns,
                step=50 * ns,
            ),
        },
        mw=MicrowaveConfig(power=8, freq=2.87e9),
        times=TimeConfig(t_AOM=3 * us, AOM_lag=800 * ns, MW_lag=150 * ns,
                         t_pi=100 * ns),
        pb=PulseBlasterConfig(
            clock_MHz=PBclk,
            channels=dict(sorted(channels.items(), key=lambda kv: kv[1])),
        ),
        seq=SequenceConfig(
            name='rabi_seq',
            args_names=['t_MW', 't_AOM', 'ro_delay', 'AOM_lag', 'MW_lag'],
            args_values=[100 * ns, 3 * us, 300 * ns, 800 * ns, 150 * ns],
            Nsamples=100,
        ),
        daq_ai=DAQAIConfig(ai_samps_per_chan=100),
        plot=PlotConfig(x_units=ns, x_label='AOM lag (ns)'),
        save_opts=SaveConfig(prefix='AOM_lag_scan'),
        runtime=RuntimeFlags(Nruns=5),
    ).use('mw', 'daq_ai')

    # ── Now: what happens inside DiodeExperiment.setup() ────────────────
    #
    # 1. setup() sees the primary scan is 'AOM_lag'
    # 2. 'rabi_seq' is NOT in FREQ_ONLY_SEQUENCES
    # 3. So it goes into the PB-reprogramming branch
    # 4. It calls make_pb_setter(pb, param_index=0, ...)
    #
    # WAIT — param_index=0 is WRONG for AOM_lag!
    #
    # The current code always uses param_index=0 because historically
    # the scanned param was always the FIRST arg.  But for AOM_lag
    # it's at index 3.
    #
    # SOLUTION: We need to tell DiodeExperiment which arg index to use.
    # The cleanest way: look up the scan name in seq.args_names.

    # This is what DiodeExperiment.setup() should do (and now does):
    #   scan_name = 'AOM_lag'
    #   param_index = config.seq.args_names.index('AOM_lag')  → 3
    #   make_pb_setter(pb, param_index=3, ...)

    print("Config created.  In a real run:")
    print("  exp, data = mainControl_diode.run(instruments, config)")
    print(f"  AOM_lag would be swept at PB arg index "
          f"{config.seq.args_names.index('AOM_lag')}")

    return config


# ═══════════════════════════════════════════════════════════════════════
# Approach A2: AOM_lag as OUTER, frequency as INNER
# ═══════════════════════════════════════════════════════════════════════

def example_A2():
    """
    2D sweep: frequency (inner, fast) × AOM_lag (outer, slow).

    For each AOM_lag value:
        PB is reprogrammed (the Rabi sequence shifts AOM timing)
        Then a full frequency sweep runs (SG changes, PB stays)
        Repeated for Nruns

    Config:
        scans = {
            'freq': ScanAxis(start=2.82*GHz, stop=2.92*GHz, step=1*MHz),
            'AOM_lag': ScanAxis(start=500*ns, stop=1500*ns, step=200*ns),
        }
        # freq is first key → inner.  AOM_lag is second → outer.

    The sweep.add() calls in setup() would be:
        sweep.add('freq', freq_values, setter=sg.set_freq)      # inner
        sweep.add('AOM_lag', lag_values, setter=pb_lag_setter)   # outer

    Where pb_lag_setter is make_pb_setter(pb, param_index=3, ...).
    """
    print("2D sweep: for each AOM_lag, run full freq sweep × Nruns")
    print("  Data shape: (n_AOM_lag_values, Nruns, n_freq_pts, daq_Nsamples)")


# ═══════════════════════════════════════════════════════════════════════
# How the setter maps to the sequence
# ═══════════════════════════════════════════════════════════════════════

def explain_pb_setter_mapping():
    """
    Visual explanation of how make_pb_setter connects config → PB.

    Your config defines:
        seq.args_names  = ['t_MW', 't_AOM', 'ro_delay', 'AOM_lag', 'MW_lag']
        seq.args_values = [100,     3000,    300,        800,       150     ]  (in ns)
                           ↑idx 0   ↑idx 1   ↑idx 2     ↑idx 3    ↑idx 4

    When you sweep 'AOM_lag', DiodeExperiment does:

        param_index = seq.args_names.index('AOM_lag')   # → 3

        pb_setter = make_pb_setter(
            pb,
            param_index=3,          # ← overwrites args[3]
            instr='diode',
            sequence='rabi_seq',
            seq_args_template=[100, 3000, 300, 800, 150],  # mutable copy
            pb_channels=channels,
        )

    When the sweep calls pb_setter(600):
        1. seq_args_template becomes [100, 3000, 300, 600, 150]
                                                      ^^^
        2. PB_program('diode', 'rabi_seq', [100,3000,300,600,150] + [channels])
           → builds new instruction list with AOM firing 600 ns after laser
        3. pb.run_sequence_for_diode(instructions)
           → PB hardware now running the updated sequence

    The next DAQ read captures data with this new timing.

    For ESR (freq sweep), the PB doesn't change at all — only SG:
        sweep.add('freq', values, setter=sg.set_freq)
        # sg.set_freq(2.87e9) sends a VISA command, PB untouched

    For Rabi (t_MW sweep at index 0):
        pb_setter = make_pb_setter(pb, param_index=0, ...)
        sweep.add('t_MW', values, setter=pb_setter)
        # pb_setter(200) → args becomes [200, 3000, 300, 800, 150]
        #                → PB reprogrammed with 200 ns MW pulse
    """
    print("See docstring for visual explanation")


if __name__ == '__main__':
    config = example_A1()
    example_A2()
    explain_pb_setter_mapping()
