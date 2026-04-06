"""
rabi_config_new.py — Rabi experiment configuration.
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
t_AOM    = 0.1 * ms
ro_delay = 500 * ns
AOM_lag  = 800 * ns
MW_lag   = 80 * ns
Nsamples = max(2, int((t_AOM) / ms * 1e-3 * 10e3))
Nsamples = 100

# ── Channel map ─────────────────────────────────────────────────────────
channels = {
    'samp':   samp_clk,
    'mw':          MW,
    'laser':       laser,
    'start': start_trig,
    # 'Camera':      camera,
    # 'Bx': bx, 'By': by, 'Bz': bz,
}

# ── Build config ────────────────────────────────────────────────────────
config = ExperimentConfig(
    scans={
        'tau': ScanAxis(
            name='tau',
            start=10 * ns,
            stop=500 * ns,
            step=50 * ns,
        ),
        # 'mw_power': ScanAxis(
        #     name='mw_power',
        #     start=-24,
        #     stop=-20,
        #     step=4,
        # ),
        'AOM_lag': ScanAxis(
            name='AOM_lag',
            start=800,
            stop=1000,
            step=200,
        ),
        # 't_AOM': ScanAxis(
        #     name='t_AOM',
        #     start=0.1 *ms,
        #     stop=0.2 *ms,
        #     step=0.1 *ms,
        # ),
    },

    mw=MicrowaveConfig(power=8, freq=3.026*GHz),
    times=TimeConfig(t_AOM=t_AOM,),
    pb=PulseBlasterConfig(
        clock_MHz=PBclk,
        channels=dict(sorted(channels.items(), key=lambda kv: kv[1])),
    ),
    seq=SequenceConfig(
        name='rabi_seq',
        args_names=['t_AOM', 'ro_delay', 'AOM_lag', 'MW_lag'],
        args_values=[t_AOM, ro_delay, AOM_lag, MW_lag],
        Nsamples=Nsamples,
    ),
    daq_ai=DAQAIConfig(ai_samps_per_chan=Nsamples),
    plot=PlotConfig(x_units=ns, x_label='Microwave pulse length (ns)'),
    save_opts=SaveConfig(prefix='Rabi'),
    runtime=RuntimeFlags(Nruns=2, reload_pb=True, seq_plot_indices=[0, 5, -1])
).use('mw', 'daq_ai')

params_dict = config.to_dict()
