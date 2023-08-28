# drift_config.py
#%% Imports
from spinapi import ns,us,ms
from SRScontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
# from time import localtime, strftime
from connectionConfig import PBclk, laser, samp_clk, start_trig, uW, conv_clk

clk_cyc = 1e3/PBclk #in ns

#%% USER INPUT

# Microwave scan parameters:----------------------------------------------------
freq1 = 28/10 * GHz          # MW frequency (in Hz)
freq2 = 28.7/10 * GHz
N_scanPts = 2               # Two frequency points
uW_power = 8         # Microwave power output from SRS(dBm)
SRSdisplay = 'disp2'

# Pulse sequence parameters:----------------------------------------------------
t_duration = t_tot = 100*ms                 # Duration of one-half of signal-aquisition half (in ns)
observation_time = 2 *(60)                                  # Time in seconds
Nsamples = observation_time*1e3*ms/t_duration               # Number of FL samples to take at each point
Navg = 1                    # Number of averaging runs to do
contrast_mode = 'ratio_signal_over_reference'       # Contrast mode

# Plotting options--------------------------------------------------------------
plotXaxisUnits = GHz         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Frequency (GHz)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2       # Save interval for first scan through all frequency points:
# saveSpacing_inAverages = 1      # Save interval in averaging runs:
savePath = os.getcwd()+"\\Saved_Data\\"         # Path to folder where data will be saved:
saveFileName = "drift"                           # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = int(False)             # Option to do shot by shot contrast normalization:
randomize = int(True)                            # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = [freq1, freq2]     # scannedParam = Frequency in ESR experiment
sequence = 'drift_seq'                 #Sequence string
PBchannels = {'Start\nTrig':start_trig, 'Conv\nCLK':conv_clk,
              'Samp\nCLK':samp_clk, 'uW':uW, 'Laser':laser}        #PB channels
# PBchannels = {'STARTtrig':start_trig, 'Convert\nClock':conv_clk,
              # 'Data\nSample':samp_clk, 'uW':uW,'Laser':laser}        #PB channels
sequenceArgs = [t_duration]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())               #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"    #Make param file path

formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n"       #" %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList =  ['N_scanPts:',N_scanPts, 'Navg:',Navg, 'Nsamples:',Nsamples, 'Freq1:',scannedParam[0], 'Freq2:',scannedParam[1], 'uW_power:',uW_power, 't_duration:',t_duration]       #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]              # 9 Paramters

def updateSequenceArgs():
    sequenceArgs = [t_duration]
    return sequenceArgs

def updateExpParamList():
    expParamList =  ['N_scanPts:',N_scanPts, 'Navg:',Navg, 'Nsamples:',Nsamples, 'Freq1:',scannedParam[0], 'Freq2:',scannedParam[1], 'uW_power:',uW_power, 't_duration:',t_duration]       #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]              # 9 Paramters
    return expParamList