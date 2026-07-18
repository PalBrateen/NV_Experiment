"""
esr_config.py — ESR experiment configuration.
"""

from experiment_config import (
    ExperimentConfig, ScanAxis, ScanMode, NonlinearSegment, PBpins, MicrowaveConfig,
    PulseBlasterConfig, TimeConfig, DAQ_INFO, SequenceConfig, DAQAIConfig, DAQCIConfig, 
    UHFLIConfig, PlotConfig, RuntimeFlags,
)
ns, us, ms = 1., 1e3, 1e6
Hz, kHz, MHz, GHz = 1., 1e3, 1e6, 1e9

# ── Timing ──────────────────────────────────────────────────────────────
t_AOM = 20 * ms
t_tot = 2 * t_AOM
# Nsamples = 2 * int((t_AOM/2 - 1*ms) / ms * 1e-3 * 2e6)
Nsamples = 2
Ncycles = 2


# ── PB channels used ────────────────────────────────────────────────────
channels = [
    PBpins.samp_clk,
    PBpins.MW,
    PBpins.laser,
    PBpins.start_trig,
    PBpins.pause_trig,
    PBpins.gate
]

# ── Build config ────────────────────────────────────────────────────────
config = ExperimentConfig(
    scans={
        'freq': ScanAxis(
            name='freq',
            start=2.82 * GHz,
            stop=2.92 * GHz,
            step=1 * MHz,
        ),
        # 'mw_power': ScanAxis(
        #     name='mw_power',
        #     start=-24,
        #     stop=-20,
        #     step=4,
        # ),
        # 't_AOM': ScanAxis(
        #     name='t_AOM',
        #     explicit_values=[5*ms, 25*ms, 50*ms],
        # ),
    },
    # TODO: check out multi-parameter acquisition in analog-input mode

    mw=MicrowaveConfig(power=-5, freq=2.87*GHz),
    times = TimeConfig(t_AOM=t_AOM),
    seq=SequenceConfig(
        name='esr_seq',
        args_names=['t_AOM'],
        args_values=[t_AOM],
        Nsamples=Nsamples,
        Ncycles=Ncycles,
    ),
    # daq_ai = DAQAIConfig(
    #     voltage_ranges=[(-.2,.2)],
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
    runtime=RuntimeFlags(Nruns=2, reload_pb=True),

    plot=PlotConfig(x_label='Frequency', x_label_units='Hz', x_units=GHz),
    pb=PulseBlasterConfig(clock_MHz=ExperimentConfig.PB_CLK,
                          channels=sorted(channels),),
)

# _OPTIONAL_FIELDS = {'mw', 'times', 'daq_ai', 'daq_ao', 'daq_ci', 'uhfli', 'camera', 'field_cfg'}
# config.use('mw', 'daq_ai')

# configure detector based on input
if config.daq_ai is None and config.daq_ci:
    config.runtime.detector = 'counter'


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
