# ESRconfig.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
# from time import localtime, strftime
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, camera

clk_cyc = 1e3/PBclk #in ns

#%% USER INPUT

# Microwave scan parameters:----------------------------------------------------
startFreq = 28.2 /10 * GHz          # Start frequency (in Hz)
endFreq = 29.2 /10 * GHz             # End frequency (in Hz)
step_size = 1 *MHz
N_scanPts = round((endFreq - startFreq)/step_size + 1)             # Number of frequency steps
# N_scanPts = 5
MW_power = 8         # Microwave power output from SRS(dBm)
SRSdisplay = 'disp2'

f = {1: [2.82, 2.85, 10], 2: [2.85, 2.86, 2], 3: [2.86, 2.88, 1], 4: [2.88, 2.89, 2], 5: [2.89, 2.92, 10]}

# f = {1: [2.82, 2.83, 10],
#       2: [2.83, 2.85, 2],
#       3: [2.85, 2.90, 1],
#       # 3: [2.84, 2.90, 1],
#       4: [2.90, 2.91, 2],
#       5: [2.91, 2.92, 10]}
params_dict = {}; scannedParam = []
for i, key in enumerate(f.keys()):
    start = f[key][0]
    stop = f[key][1]
    step = f[key][2]
    # print(f"{start}, {stop}, {step}")
    n = round((stop - start)/(step/1e3))
    if i < len(f)-1:
        params_dict[i] = [n, np.linspace(start, stop, n, endpoint=False)]
        scannedParam.extend(params_dict[i][-1])
    else:
        params_dict[i] = [n, np.linspace(start, stop, n+1, endpoint=True)]
        scannedParam.extend(params_dict[i][-1])
# N_scanpts = sum([vals[0] for vals in params_dict.values()])

# Pulse sequence parameters:----------------------------------------------------
t_AOM = 8000 *ms          # Duration of one-half of signal-aquisition half
# t_AOM = 5 *ms
t_tot = 2*t_AOM
Nsamples = 8            # Number of FL samples to take at each frequency point
Nruns = 30                    # Number of averaging runs to do
contrast_mode = 'ratio_signal_over_reference'       # Contrast mode

# Plotting options--------------------------------------------------------------
plotXaxisUnits = GHz         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Frequency (GHz)'       # Plot x axis label

# Save options------------------------------------------------------------------
# cwd = os.getcwd()
# savePath = cwd[0:cwd.rfind('\\')]+"\\Saved_Data\\"         # Path to folder where data will be saved:
saveFileName = "ESR"                           # File name for data file

# Averaging options:------------------------------------------------------------
# shotByShotNormalization = int(False)             # Option to do shot by shot contrast normalization:
# randomize = int(True)                            # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#
scannedParam = np.array(scannedParam)*1e9
# scannedParam = np.linspace(startFreq,endFreq, N_scanPts, endpoint=True)     # scannedParam = Frequency in ESR experiment
sequence = 'esr_seq'                 #Sequence string:
scanStartName = 'startFreq'         #Scan start Name
scanEndName = 'endFreq'             #Scan end Name
PBchannels = {'Samp\nCLK':samp_clk, 'MW':MW,'Laser':laser, 'Start\nTrig':start_trig, 'Camera': camera}
PBchannels = dict(sorted(PBchannels.items(), key=lambda item: item[1]))

sequenceArgs = [t_AOM]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())               #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"    #Make param file path

formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n"       #" %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList =  ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'Step_size:',step_size, 'startFreq:',scannedParam[0], 'endFreq:',scannedParam[-1], 'MW_power:',MW_power, 't_duration:',t_AOM]

def updateSequenceArgs():
    sequenceArgs = [t_AOM]
    return sequenceArgs

def updateExpParamList():
    expParamList =  ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'Step_size:',step_size, 'startFreq:',scannedParam[0], 'endFreq:',scannedParam[-1], 'MW_power:',MW_power, 't_duration:',t_AOM]   # 9 Paramters
    return expParamList