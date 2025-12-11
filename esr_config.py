# ESRconfig.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os, numpy as np
# from time import localtime, strftime
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, camera


clk_cyc = 1e3/PBclk #in ns
MW_power = -18         # [dBm]

scan_range = {
    'freq': {
        'range': [2.82 *GHz, 2.92 *GHz],       # [start, end] frequency in Hz
        # 'range': None,
        'step': 1 *MHz,
        # 'step': None,
        # 'Nscanpts': None,
        # 'Nscanpts': 10,
        },

    # 'mw_power': {
        # 'range': [0, 10],       # [start, end] power in dBm
        # 'step': 1,              # step in dBm
        # 'Nscanpts': 11,
        # 'values': None,
        # }
    }
shuffle_values = True

# Sequence parameters:----------------------------------------------------

t_AOM = 8000 *ms          # Duration of one-half of signal-aquisition half (camera 8000 *ms)
t_AOM = 0.08 *ms
# t_AOM = 5 *ms
t_tot = 2*t_AOM

sequenceArgs = [t_AOM]         #Sequence args

# Averaging parameters:
Nsamples = 100            # Number of FL samples to take at each frequency point
Nruns = 5                 # Number of averaging runs

reload_pb = True
load_all_params = False

PBchannels = {'samp':samp_clk, 'mw':MW, 'laser':laser, 'start':start_trig}#, 'camera': camera}

## ------------------------------------------------------------
# Can extend to other parameters using <<for param in scan_range.keys()>> structure
input_range = scan_range['freq']['range']
if input_range is not None:
    scan_start = input_range[0]
    scan_end = input_range[1] if len(input_range)>1 else scan_start

    if scan_range['freq'].get('Nscanpts', None) is None:
        scan_step = scan_range['freq']['step']
        scan_range['freq']['Nscanpts'] = round((scan_end - scan_start)/scan_step + 1)

    elif scan_range['freq'].get('step', None) is None:
        Nscanpts = scan_range['freq']['Nscanpts']
        scan_range['freq']['step'] = (scan_end - scan_start)/(Nscanpts - 1) if Nscanpts>1 else 0
    
    scan_range['freq']['values'] = list(np.linspace(scan_start, scan_end, scan_range['freq']['Nscanpts'], endpoint=True))
else:
    # For non-linear frequency sweeps, define scan_range['freq']['range'] = <list of 1st and last elements>,
    # scan_range['freq']['step'] = None,
    # scan_range['freq']['Nscanpts'] = <number of points in each segment>

    def nonlin_freq_sweep(f):
        params_dict = {}; values = []
        for i, key in enumerate(f.keys()):
            start = f[key][0]
            stop = f[key][1]
            step = f[key][2]
            # print(f"{start}, {stop}, {step}")
            n = round((stop - start)/(step/1e3))
            if i < len(f)-1:
                params_dict[i] = [n, np.linspace(start, stop, n, endpoint=False)]
                values.extend(params_dict[i][-1])
            else:
                params_dict[i] = [n, np.linspace(start, stop, n+1, endpoint=True)]
                values.extend(params_dict[i][-1])
        return params_dict, values
    
    # f = {1: [2.82, 2.85, 10], 2: [2.85, 2.86, 2], 3: [2.86, 2.88, 1], 4: [2.88, 2.89, 2], 5: [2.89, 2.92, 10]}

    f = {
        1: [2.82, 2.83, 10],
        2: [2.83, 2.85, 2],
        3: [2.85, 2.89, 1],
        # 3: [2.84, 2.90, 1],
        4: [2.89, 2.91, 2],
        5: [2.91, 2.92, 10],
    }
    params_dict, values = nonlin_freq_sweep(f)

    scan_range['freq']['range'] = [values[0], values[-1]]
    scan_range['freq']['step'] = None
    scan_range['freq']['Nscanpts'] = sum([vals[0] for vals in params_dict.values()])
    scan_range['freq']['values'] = list(np.array(values) *GHz)

# if shuffle_values:
#     np.random.shuffle(scan_range['freq']['values'])

PBchannels = dict(sorted(PBchannels.items(), key=lambda item: item[1]))

seq_times = {
    't_AOM': t_AOM,
    't_total': t_tot,
}

def update_params_dict():
    """
    Make the base parameter dictionary
    """
    
    params = {
        'scan': {
            'names': list(scan_range.keys()),
            'values': [scan_range[key]['values'] for key in scan_range.keys()],
            'Nscanpts': scan_range['freq']['Nscanpts'],
            'Nruns': Nruns,
            'reload_pb': reload_pb,
            'load_all_params': load_all_params,
        },

        'mw': {
            'power': MW_power,
            'freq': scan_range['freq']['values'],
        },
        
        'seq': {
            'Nsamples': Nsamples,
            'channels': PBchannels,
            'sequence': 'esr_seq',
            'args_name': ['t_AOM'],
            'args': sequenceArgs,

            # instructionList to be added later - only the last generated sequence will be saved
            # the variable that goes into the plot function
        },

        'pb': {
            'clk_cyc': clk_cyc,
            'channels': PBchannels,
        },

        'daq': {
            'ai': {
                'Nsamples': Nsamples,
                # add sample rate
                # add input channel numbers and voltage ranges
                # add trigger channel(s) and source(s) (multiple PB channels used)
            },
            'ao': {
                # add output channel numbers and voltage ranges
            },
        },

        'plot': {
            'plotXaxisUnits': GHz,
            'plotXaxisLabel': 'Frequency (GHz)',
        },
        'save': {
            'savefileprefix': "ESR",
        },
        'laser': {
            't_AOM': seq_times['t_AOM'],
            'power': 0,
        }
    }

    # params['seq'].update(seq_times)
    return params

params_dict = update_params_dict()

## ---------------------------------------------
# below is for YAML serialization, JSON implements this by default
# def convert_numpy_to_python(obj):
#     """Recursively convert numpy types to native Python types"""
#     if isinstance(obj, np.ndarray):
#         return obj.tolist()
#     elif isinstance(obj, np.generic):  # This catches numpy.float64, etc.
#         return obj.item()
#     elif isinstance(obj, dict):
#         return {key: convert_numpy_to_python(value) for key, value in obj.items()}
#     elif isinstance(obj, list):
#         return [convert_numpy_to_python(item) for item in obj]
#     else:
#         return obj

##---------------------------------------------
# Syntax for JSON save
# import json
# with open('params_dict.json', 'w') as f:
#     json.dump(params_dict, f, indent=4)

##---------------------------------------------
# syntax for YAML save
# import yaml
# with open('params_dict.yaml', 'w') as f:
#     yaml.dump(convert_numpy_to_python(params_dict), f)

