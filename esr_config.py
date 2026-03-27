"""
esr_config_new.py — ESR experiment configuration.
"""

from experiment_config import (
    ExperimentConfig, ScanAxis, ScanMode, NonlinearSegment,
    MicrowaveConfig, PulseBlasterConfig, TimeConfig,
    SequenceConfig, DAQAIConfig, UHFLIConfig, PlotConfig, SaveConfig, RuntimeFlags,
    ns, us, ms, Hz, kHz, MHz, GHz,
)

try:
    from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, lia1, camera
except ImportError:
    PBclk = 500; laser = samp_clk = start_trig = MW = lia1 = camera = 0


# ── Timing ──────────────────────────────────────────────────────────────
t_AOM = 70 * ms
t_tot = 2 * t_AOM
Nsamples = int((t_AOM - 40*ms) / ms * 1e-3 * 10e3)

# ── Channel map ─────────────────────────────────────────────────────────
channels = {
    'samp':  samp_clk,
    'mw':    MW,
    'laser': laser,
    'start': start_trig,
}

# ── Build config ────────────────────────────────────────────────────────
config = ExperimentConfig(
    scans={
        'freq': ScanAxis(
            name='freq',
            start=2.82 * GHz,
            stop=2.92 * GHz,
            step=1 * MHz,
        ),
        'mw_power': ScanAxis(
            name='mw_power',
            start=-24,
            stop=-20,
            step=4,
        ),
    },

    mw=MicrowaveConfig(power=-24, freq=2.87*GHz),
    times = TimeConfig(t_AOM=t_AOM),
    pb=PulseBlasterConfig(
        clock_MHz=PBclk,
        channels=dict(sorted(channels.items(), key=lambda kv: kv[1])),
    ),
    seq=SequenceConfig(
        name='esr_seq',
        args_names=['t_AOM'],
        args_values=[t_AOM],
        Nsamples=Nsamples,
    ),
    daq_ai=DAQAIConfig(ai_samps_per_chan=Nsamples),
    plot=PlotConfig(x_units=GHz, x_label='Frequency (GHz)'),
    save_opts=SaveConfig(prefix='ESR'),
    runtime=RuntimeFlags(Nruns=2, reload_pb=True),
).use('mw', 'daq_ai', 'daq_ao')


# ── Alternative scan modes (uncomment one to use) ──────────────────────

# # Nonlinear (piecewise) sweep:
# config.scans['freq'] = ScanAxis(
#     name='freq',
#     segments=[
#         NonlinearSegment(2.82, 2.83, 10),  # coarse
#         NonlinearSegment(2.83, 2.85, 2),   # medium
#         NonlinearSegment(2.85, 2.89, 1),   # fine around resonance
#         NonlinearSegment(2.89, 2.91, 2),
#         NonlinearSegment(2.91, 2.92, 10),
#     ],
#     unit_scale=GHz,
# )

# # Explicit list of specific frequencies:
# config.scans['freq'] = ScanAxis(
#     name='freq',
#     explicit_values=[2.85*GHz, 2.87*GHz, 2.89*GHz],
# )

# ── Backward-compatible dict ────────────────────────────────────────────
params_dict = config.to_dict()
