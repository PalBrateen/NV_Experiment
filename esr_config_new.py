# esr_config_new.py
# New simplified config format for ESR experiment
#%% Imports
from spinapi import ns, us, ms, s
from SGcontrol import Hz, kHz, MHz, GHz
import numpy as np
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, camera

# ============================================================================
# SEQUENCE PARAMETERS
# ============================================================================

clk_cyc = 1e3/PBclk  # in ns
t_AOM = 8000 * ms    # Duration of AOM pulse (signal acquisition time)

PBchannels = {
    'samp': samp_clk,
    'mw': MW,
    'laser': laser,
    'start': start_trig
}

# ============================================================================
# MEASUREMENT PARAMETERS
# ============================================================================

Nsamples = 8  # Number of fluorescence samples per point
Nruns = 30    # Number of averaging runs

# ============================================================================
# MW PARAMETERS
# ============================================================================

MW_power = 8  # [dBm]

# ============================================================================
# SWEEP DEFINITION (NEW FORMAT!)
# ============================================================================
# Order matters: First defined = innermost loop (fastest changing)

params_dict = {
    'seq': {
        'sequence': 'esr_seq',
        't_AOM': t_AOM,
        'Nsamples': Nsamples,
        'PBchannels': PBchannels,
    },

    'mw': {
        'power': MW_power,
        # Note: frequency is swept, not fixed here
    },

    'scan': {
        'Nruns': Nruns,
    },

    'pb': {
        'clk_cyc': clk_cyc,
        'channels': PBchannels,
    },

    # NEW: Sweep definitions
    # Order matters: First added = innermost loop
    'sweep': {
        'frequency': np.linspace(2.82*GHz, 2.92*GHz, 101),
        # Add more parameters for multi-parameter sweeps:
        # 'power': np.array([5, 8, 10]),  # Outer loop (slowest changing)
    },

    'plot': {
        'plotXaxisUnits': GHz,
        'plotXaxisLabel': 'Frequency (GHz)',
    },

    'save': {
        'savefileprefix': "ESR",
    },
}
