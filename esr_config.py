# ESRconfig.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os, numpy as np
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, lia, camera

clk_cyc = 1e3/PBclk #in ns
MW_power = -24         # [dBm]

scan = {
    'freq': {
        # 'range': [3.015 *GHz, 3.032 *GHz],
        'range': [2.9275 *GHz, 2.944 *GHz],
        # 'range': [2.74 *GHz, 3. *GHz],
        # 'range': [2.82 *GHz, 2.92 *GHz],
        # 'range': [3.02732 *GHz],
        # 'range': None,
        
        # 'range': [2.62 *GHz, 3.12 *GHz],
        # 'range': [2.995 *GHz, 3.025 *GHz],
        
        'step': 0.05 *MHz,
        # 'step': 1 *MHz,
        # 'step': None,

        # 'Nscanpts': None,
        # 'Nscanpts': 1,
        'shuffle': False,
        },

    # 'mw_power': {
        # 'range': [0, 10],       # [start, end] power in dBm
        # 'step': 1,              # step in dBm
        # 'Nscanpts': 11,
        # 'values': None,
        # }
    }

# Sequence parameters:----------------------------------------------------

# t_AOM = 8000 *ms          # Duration of one-half of signal-aquisition half (camera 8000 *ms)
# t_AOM = 0.08 *ms
# t_AOM = 2*ms
# t_AOM = 20 *us
t_AOM = 70 *ms

t_tot = 2*t_AOM

sequenceArgs = [t_AOM]         #Sequence args

# Averaging parameters:
# Nsamples = 10
Nsamples = int((t_AOM - 40*ms) /ms*1e-3 *10e3)
# Nsamples = int( 2*60 *10e3 )           # Number of FL samples to take at each frequency point
Nruns = 1                 # Number of averaging runs

# TODO: DAQ AI
# daq = {
#     'ai': {
#         'sample_rate': 10e3,    # [Sa/s]
#         'Nsamples': Nsamples,
#     }
# }
reload_pb = True
load_all_params = False

PBchannels = {'samp':samp_clk,
              'mw':MW,
              'laser':laser,
              'start':start_trig,
              #, 'camera': camera,
              'lia': lia,
            }

## ------------------------------------------------------------
# Can extend to other parameters using <<for param in scan.keys()>> structure
input_range = scan['freq']['range']
if input_range is not None:
    scan_start = input_range[0]
    scan_end = input_range[1] if len(input_range)>1 else scan_start

    if scan['freq'].get('Nscanpts', None) is None:
        scan_step = scan['freq']['step']
        scan['freq']['Nscanpts'] = round((scan_end - scan_start)/scan_step + 1)

    elif scan['freq'].get('step', None) is None:
        Nscanpts = scan['freq']['Nscanpts']
        scan['freq']['step'] = (scan_end - scan_start)/(Nscanpts - 1) if Nscanpts>1 else 0
    
    scan['freq']['values'] = list(np.linspace(scan_start, scan_end, scan['freq']['Nscanpts'], endpoint=True))
else:
    # For non-linear frequency sweeps, define scan['freq']['range'] = <list of 1st and last elements>,
    # scan['freq']['step'] = None,
    # scan['freq']['Nscanpts'] = <number of points in each segment>

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

    scan['freq']['range'] = [values[0], values[-1]]
    scan['freq']['step'] = None
    scan['freq']['Nscanpts'] = sum([vals[0] for vals in params_dict.values()])
    scan['freq']['values'] = list(np.array(values) *GHz)

if  scan['freq']['shuffle']:
    np.random.shuffle(scan['freq']['values'])

PBchannels = dict(sorted(PBchannels.items(), key=lambda item: item[1]))

seq_times = {
    't_AOM': t_AOM,
    't_total': t_tot,
}

def update_params_dict():
    """
    Make the base parameter dictionary
    """
    def remove_scan_values(dict_with_values):
        dict_without_values = {
        param: {k: v for k, v in config.items() if k != 'values'}
        for param, config in dict_with_values.items()
        }
        return dict_without_values

    params = {'scan': scan}
    add_params = {
        'names': list(scan.keys()),
        'Nruns': Nruns,
        'reload_pb': reload_pb,
        'load_all_params': load_all_params,
        }
    params['scan'].update(add_params)

    add_params = {
        'mw': {
            'power': MW_power if 'mw_power' not in params['scan']['names'] else [],
            'freq': scan['freq']['values'] if 'freq' not in params['scan']['names'] else [np.min(scan['freq']['values'])],                
        },
        
        'seq': {
            'Nsamples': Nsamples,
            'channels': PBchannels,
            'sequence': 'esr_seq',
            'args_name': ['t_AOM'],
            'args': sequenceArgs,
            't_total(s)': seq_times['t_total']/1e9

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
    params.update(add_params)
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

# FIXME: below error comes when 'mw' PBchannel key is removed which should not happen...
# Traceback (most recent call last):

#   File "D:\Brateen\NV_Experiment\mainControl_diode.py", line 474, in <module>
#     [seqArgList, instructionList] = initialize_exp(instr)

#   File "D:\Brateen\NV_Experiment\mainControl_diode.py", line 166, in initialize_exp
#     sequencecontrol.view_sequence(instr, params['seq']['sequence'], seqArgList+[params['pb']['channels']],

#   File "D:\Brateen\NV_Experiment\sequencecontrol.py", line 413, in view_sequence
#     plt.yticks(yTicks, seqArgList[-1].keys())          # Include the names of the PB channels

#   File "C:\Users\Smart\miniconda3\envs\nv-pc\lib\site-packages\matplotlib\pyplot.py", line 1719, in yticks
#     labels = ax.set_yticklabels(labels, **kwargs)

#   File "C:\Users\Smart\miniconda3\envs\nv-pc\lib\site-packages\matplotlib\axes\_base.py", line 63, in wrapper
#     return get_method(self)(*args, **kwargs)

#   File "C:\Users\Smart\miniconda3\envs\nv-pc\lib\site-packages\matplotlib\cbook\deprecation.py", line 451, in wrapper
#     return func(*args, **kwargs)

#   File "C:\Users\Smart\miniconda3\envs\nv-pc\lib\site-packages\matplotlib\axis.py", line 1796, in _set_ticklabels
#     return self.set_ticklabels(labels, minor=minor, **kwargs)

#   File "C:\Users\Smart\miniconda3\envs\nv-pc\lib\site-packages\matplotlib\axis.py", line 1717, in set_ticklabels
#     raise ValueError(

# ValueError: The number of FixedLocator locations (4), usually from a call to set_ticks, does not match the number of ticklabels (3).
# %%
