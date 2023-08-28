# ESRconfig.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
# from time import localtime, strftime
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, conv_clk, camera

clk_cyc = 1e3/PBclk #in ns

#%% USER INPUT

# Microwave scan parameters:----------------------------------------------------
startFreq = 28.0 /10 * GHz          # Start frequency (in Hz)
endFreq = 29.4 /10 * GHz             # End frequency (in Hz)
step_size = 1*MHz
N_scanPts = round((endFreq - startFreq)/step_size + 1)             # Number of frequency steps
# N_scanPts = 5
MW_power = 8         # Microwave power output from SRS(dBm)
SRSdisplay = 'disp2'

# Pulse sequence parameters:----------------------------------------------------
t_AOM = 20 *us          # Duration of one-half of signal-aquisition half
t_tot = 2*t_AOM
Nsamples = 2500             # Number of FL samples to take at each frequency point
Nruns = 1                    # Number of averaging runs to do
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

scannedParam = np.linspace(startFreq,endFreq, N_scanPts, endpoint=True)     # scannedParam = Frequency in ESR experiment
sequence = 'esr_seq'                 #Sequence string:
scanStartName = 'startFreq'         #Scan start Name
scanEndName = 'endFreq'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk,'Samp\nCLK':samp_clk, 'MW':MW,'Laser':laser, 'Start\nTrig':start_trig}#,'Camera': camera,}

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