"""
ramsey_config.py — Ramsey experiment configuration.

Demonstrates:
  - How trivially a new experiment type plugs in
  - ScanMode.LOG for logarithmic tau sweeps (T2* can span decades)
"""

from experiment_config import (
    ExperimentConfig, ScanAxis, ScanMode,
    MicrowaveConfig, TimeConfig, PulseBlasterConfig,
    SequenceConfig, DAQAIConfig, PlotConfig, SaveConfig, RuntimeFlags,
    ns, us, ms, Hz, kHz, MHz, GHz,
)

try:
    from connectionConfig import (PBclk, laser, samp_clk, start_trig,
                                  MW, camera, bx, by, bz)
except ImportError:
    PBclk = 500
    laser = samp_clk = start_trig = MW = camera = bx = by = bz = 0

# ── Timing ──────────────────────────────────────────────────────────────
t_AOM     = 3 * ms
ro_delay  = 300 * ns
AOM_lag   = 800 * ns
MW_lag    = 80 * ns
pi2_time  = 25 * ns

Nsamples = int((t_AOM - 40*ms) / ms * 1e-3 * 10e3)

channels = {
    'SampCLK':   samp_clk,
    'MW':          MW,
    'Laser':       laser,
    'StartTrig': start_trig,
    'Camera':      camera,
    'Bx': bx, 'By': by, 'Bz': bz,
}

config = ExperimentConfig(
    scans={
        'tau': ScanAxis(
            name='tau',
            start=100 * ns,
            stop=100 * us,
            Nscanpts=80,
            mode=ScanMode.LOG,    # logarithmic spacing for T2* decay
        ),
    },

    mw=MicrowaveConfig(power=8, freq=3.026*GHz),
    times=TimeConfig(t_AOM=t_AOM, t_pi=pi2_time),
    pb=PulseBlasterConfig(
        clock_MHz=PBclk,
        channels=dict(sorted(channels.items(), key=lambda kv: kv[1])),
    ),
    seq=SequenceConfig(
        name='ramsey_seq',
        args_names=['t_AOM', 'ro_delay', 'AOM_lag', 'MW_lag', 'pi2_time'],
        args_values=[t_AOM, ro_delay, AOM_lag, MW_lag, pi2_time],
        Nsamples=Nsamples,
    ),
    daq_ai=DAQAIConfig(ai_samps_per_chan=Nsamples),
    plot=PlotConfig(x_units=us, x_label='Free precession time [us]'),
    save_opts=SaveConfig(prefix='Ramsey'),
    runtime=RuntimeFlags(Nruns=10, reload_pb=True),
).use('mw', 'daq_ai',)

params_dict = config.to_dict()
