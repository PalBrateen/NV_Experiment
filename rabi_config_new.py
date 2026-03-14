# rabi_config_new.py
# New simplified config format for Rabi oscillation experiment
#%% Imports
from spinapi import ns, us, ms, s
from SGcontrol import Hz, kHz, MHz, GHz
import numpy as np
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, camera

# ============================================================================
# SEQUENCE PARAMETERS
# ============================================================================

clk_cyc = 1e3/PBclk  # in ns

# Timing parameters
t_AOM = 20 * us          # AOM pulse duration
ro_delay = 1800 * ns     # Readout delay
AOM_lag = 800 * ns       # AOM lag time
MW_lag = 150 * ns        # MW lag time

PBchannels = {
    'samp': samp_clk,
    'mw': MW,
    'laser': laser,
    'start': start_trig
}

# ============================================================================
# MEASUREMENT PARAMETERS
# ============================================================================

Nsamples = 2  # Number of samples per point (signal + reference)
Nruns = 1     # Number of averaging runs

# ============================================================================
# MW PARAMETERS
# ============================================================================

MW_power = 8       # [dBm]
MW_freq = 2.87*GHz  # Fixed frequency for Rabi

# ============================================================================
# SWEEP DEFINITION (NEW FORMAT!)
# ============================================================================
# Sweep pulse_duration (MW pulse length) for Rabi oscillations

params_dict = {
    'seq': {
        'sequence': 'rabi_seq',
        't_AOM': t_AOM,
        'ro_delay': ro_delay,
        'AOM_lag': AOM_lag,
        'MW_lag': MW_lag,
        'Nsamples': Nsamples,
        'PBchannels': PBchannels,
    },

    'mw': {
        'power': MW_power,
        'freq': MW_freq,
    },

    'scan': {
        'Nruns': Nruns,
    },

    'pb': {
        'clk_cyc': clk_cyc,
        'channels': PBchannels,
    },

    # NEW: Sweep definitions
    # Sweep pulse_duration from 10 ns to 500 ns in steps of 2 ns
    'sweep': {
        'pulse_duration': np.arange(10*ns, 500*ns, 2*ns),
        # Note: PB will auto-reprogram with each pulse_duration change!
    },

    'plot': {
        'plotXaxisUnits': ns,
        'plotXaxisLabel': 'Pulse Duration (ns)',
    },

    'save': {
        'savefileprefix': "Rabi",
    },
}
