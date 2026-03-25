# Rabiconfig.py
#%%Imports
from spinapi import ns,us,ms
from SGcontrol import GHz, kHz, MHz
import numpy as np, os
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, camera, bx, by, bz

clk_cyc = 1e3/PBclk #in ns

MW_power = 8                 # Microwave power output from SRS(dBm)
MW_freq = 3.026 * GHz     # Microwave frequency (Hz)

scan = {
    'tau': {
        'range': [10 *ns, 500 *ns],
        'step': 20 *ns,
        # 'Nscanpts': None,
        'shuffle': False,
        },

    # 'mw_power': {
    #     'range': [-24, -20],       # [start, end] power in dBm
    #     'step': 4,              # step in dBm
    #     # 'Nscanpts': 11,
    #     # 'values': None,
    #     }
    }

# Sequence parameters:----------------------------------------------------

# t_AOM = 8000 *ms          # Duration of one-half of signal-aquisition half (camera 8000 *ms)
# t_AOM = 2*ms

t_AOM = 0.1 *ms                    # AOM pulse duration (optimized_initial_pulse)
# t_tot = 2*t_AOM
# ro_delay = (1800)*ns      # Readout delay (ns)
ro_delay = 500 *ns          # confocal ???

# AOM_lag = (1450)*ns     # first parameter = AOM+Preamp lag, 2nd parameter = rise/fall time of the signal as seen in PMT-Preamp-DAQ
AOM_lag = (800)*ns          # for camera
# AOM_lag = 780 *ns            # for APD in confocal

MW_lag = 80*ns

sequenceArgs = [t_AOM, ro_delay, AOM_lag, MW_lag]      #Sequence args

# Averaging parameters:
Nsamples = 2000
# Nsamples = int((t_AOM - 40*ms) /ms*1e-3 *10e3)
# Nsamples = int( 2*60 *10e3 )           # Number of FL samples to take at each frequency point
Nruns = 2                 # Number of averaging runs

# TODO: DAQ AI
# daq = {
#     'ai': {
#         'sample_rate': 10e3,    # [Sa/s]
#         'Nsamples': Nsamples,
#     }
# }
reload_pb = True
load_all_params = False

PBchannels = {
            'samp':samp_clk,
            'mw':MW,
            'laser':laser,
            'start':start_trig,
            # 'Camera':camera,
            # 'Bx':bx,
            # 'By':by,
            # 'Bz':bz
            }

## ------------------------------------------------------------
# Can extend to other parameters using <<for param in scan.keys()>> structure
input_range = scan['tau']['range']
if input_range is not None:
    scan_start = input_range[0]
    scan_end = input_range[1] if len(input_range)>1 else scan_start

    if scan['tau'].get('Nscanpts', None) is None:
        scan_step = scan['tau']['step']
        scan['tau']['Nscanpts'] = round((scan_end - scan_start)/scan_step + 1)

    elif scan['tau'].get('step', None) is None:
        Nscanpts = scan['tau']['Nscanpts']
        scan['tau']['step'] = (scan_end - scan_start)/(Nscanpts - 1) if Nscanpts>1 else 0
    
    scan['tau']['values'] = list(np.linspace(scan_start, scan_end, scan['tau']['Nscanpts'], endpoint=True))
else:
    # For non-linear frequency sweeps, define scan['tau']['range'] = <list of 1st and last elements>,
    # scan['tau']['step'] = None,
    # scan['tau']['Nscanpts'] = <number of points in each segment>

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

    f = {
        1: [2.82, 2.83, 10],
        2: [2.83, 2.85, 2],
        3: [2.85, 2.89, 1],
        # 3: [2.84, 2.90, 1],
        4: [2.89, 2.91, 2],
        5: [2.91, 2.92, 10],
    }
    params_dict, values = nonlin_freq_sweep(f)

    scan['tau']['range'] = [values[0], values[-1]]
    scan['tau']['step'] = None
    scan['tau']['Nscanpts'] = sum([vals[0] for vals in params_dict.values()])
    scan['tau']['values'] = list(np.array(values) *GHz)

if  scan['tau']['shuffle']:
    np.random.shuffle(scan['tau']['values'])

# input_range = scan['mw_power']['range']
# if input_range is not None:
#     scan_start = input_range[0]
#     scan_end = input_range[1] if len(input_range)>1 else scan_start

#     if scan['mw_power'].get('Nscanpts', None) is None:
#         scan_step = scan['mw_power']['step']
#         scan['mw_power']['Nscanpts'] = round((scan_end - scan_start)/scan_step + 1)

#     elif scan['mw_power'].get('step', None) is None:
#         Nscanpts = scan['mw_power']['Nscanpts']
#         scan['mw_power']['step'] = (scan_end - scan_start)/(Nscanpts - 1) if Nscanpts>1 else 0
    
#     scan['mw_power']['values'] = list(np.linspace(scan_start, scan_end, scan['mw_power']['Nscanpts'], endpoint=True))

PBchannels = dict(sorted(PBchannels.items(), key=lambda item: item[1]))

seq_times = {
    't_AOM': t_AOM,
    'ro_delay': ro_delay,
    'AOM_lag': AOM_lag,
    'MW_lag': MW_lag,
    # 't_total': t_tot,   # total sequence time is variable for rabi
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
            'power': MW_power if 'mw_power' not in params['scan']['names'] else [], #[np.min(scan['mw_power']['values'])],
            'freq': [MW_freq],
        },
        
        'seq': {
            'Nsamples': Nsamples,
            'channels': PBchannels,
            'sequence': 'rabi_seq',
            'args_name': ['t_AOM', 'ro_delay', 'AOM_lag', 'MW_lag'],
            'args': sequenceArgs,
            # 't_total(s)': seq_times['t_total']/1e9

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
            'plotXaxisUnits': ns,
            'plotXaxisLabel': 'Microwave pulse length (ns)',
        },
        'save': {
            'savefileprefix': "Rabi",
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


