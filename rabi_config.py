"""
rabi_config.py — Rabi experiment configuration.
"""

from experiment_config import (
    ExperimentConfig, ScanAxis, ScanMode, PBpins, DAQ_INFO, MicrowaveConfig, TimeConfig,
    PulseBlasterConfig, SequenceConfig, DAQAIConfig, PlotConfig, RuntimeFlags, DAQCIConfig,
)
ns, us, ms = 1., 1e3, 1e6
Hz, kHz, MHz, GHz = 1., 1e3, 1e6, 1e9

# ── Timing ──────────────────────────────────────────────────────────────
t_AOM    = 2000 * us
ro_delay = 1000 * ns
AOM_lag  = 800 * ns
MW_lag   = 80 * ns
# Nsamples = max(2, int((t_AOM) / ms * 1e-3 * 10e3))
# Nsamples = 2 * int(1*us / us*1e-6 * 2e6)
Nsamples = 2
Ncycles = 7500

# ── Channel map ─────────────────────────────────────────────────────────
channels = [
    PBpins.samp_clk,
    PBpins.MW,
    PBpins.laser,
    PBpins.start_trig,
    PBpins.pause_trig,
    PBpins.gate,
]

# ── Build config ────────────────────────────────────────────────────────
config = ExperimentConfig(
    scans={
        'tau': ScanAxis(
            name='tau',
            start=10 * ns,
            stop=500 * ns,
            step=2 * ns,
        ),
        'mw_power': ScanAxis(
            name='mw_power',
            explicit_values=[8, -20],
        ),
        # 'AOM_lag': ScanAxis(
        #     name='AOM_lag',
        #     start=800,
        #     stop=1000,
        #     step=200,
        # ),
        # 't_AOM': ScanAxis(
        #     name='t_AOM',
        #     start=0.1 *ms,
        #     stop=0.2 *ms,
        #     step=0.1 *ms,
        # ),
    },

    mw=MicrowaveConfig(power=8, freq=2.8965*GHz),
    times=TimeConfig(t_AOM=t_AOM, ro_delay=ro_delay),
    seq=SequenceConfig(
        name='rabi_seq',
        args_names=['t_AOM', 'ro_delay', 'AOM_lag', 'MW_lag'],
        args_values=[t_AOM, ro_delay, AOM_lag, MW_lag],
        Nsamples=Nsamples,
        Ncycles=Ncycles
    ),
    # daq_ai = DAQAIConfig(
    #     voltage_ranges=[(-0.1,0.1)],
    #     sample_source='',
    #     # sample_source=DAQ_INFO['samp_clk_terminal'],
    #     sample_rate=2e6,
    #     samps_per_chan=Nsamples,
    #     start_trigger_source=DAQ_INFO['start_trig_terminal'],
    #     pause_trigger_source=DAQ_INFO['pause_trig_terminal'],
    #     # pause_trigger_source='',
    # ),
    daq_ci = DAQCIConfig(
        sample_source=DAQ_INFO['samp_clk_terminal'],
        samps_per_chan=Nsamples,
        start_trigger_source=DAQ_INFO['start_trig_terminal'],
        pause_trigger_source=DAQ_INFO['pause_trig_terminal'],
    ),
    runtime=RuntimeFlags(Nruns=2, reload_pb=True, seq_plot_indices=[-1]),

    plot=PlotConfig(x_label='MW Duration', x_label_units='s', x_units=ns),
    pb=PulseBlasterConfig(clock_MHz=ExperimentConfig.PB_CLK,
                          channels=sorted(channels),),
)

# _OPTIONAL_FIELDS = {'mw', 'times', 'daq_ai', 'daq_ao', 'daq_ci', 'uhfli', 'camera', 'field_cfg'}
# config.use('mw', 'daq_ai')

# configure detector based on input
if config.daq_ai is None and config.daq_ci:
    config.runtime.detector = 'counter'