# -*- coding: utf-8 -*-
# uW_timing_config.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
from time import localtime, strftime
from connectionConfig import *

clk_cyc = 1e3/PBclk #in ns
#%% USER INPUT

# Time gap scan
start_time = 0
# end_time = 
step_size = 100*ns
# N_scanPts = round((end_time - start_time)/step_size + 1)

# Pulse sequence parameters:----------------------------------------------------
MW_freq = 2849 /1e3 * GHz         # uW frequency corresponding to max dip in ODMR/ESR
MW_power = 8                  # uW Power in dBm

t_AOM = 30*us          # Duration of one-half of signal-aquisition
start_delay = 3*us      # Time from which the pulse sequence starts
AOM_lag = 830*ns
Nsamples = 1000             # Number of FL samples to take at each scan point
Nruns = 1                    # Number of averaging runs to do
# contrast_mode ='ratio_signal_over_reference'       # Contrast mode
t_uW = t_AOM-10*us
end_time = t_uW+50*us         # 300ns accounts for the start_trig_width
N_scanPts = round((end_time - start_time)/step_size + 1)


# Plotting options--------------------------------------------------------------
livePlotUpdate = True       # Live plot update option
plotPulseSequence = True    # Plot pulse sequence option
plotXaxisUnits = ns         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Readout-pulse Time (ns)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2     # Save interval for first scan through all freq points
# saveSpacing_inAverages = 1      # Save interval in averaging runs
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved
saveFileName = "MW_timing"        # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = int(False)     # Option to do shot by shot contrast normalization:
randomize = int(True)                    # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = np.linspace(start_time,end_time, N_scanPts, endpoint=True)     # scannedParam = 'reading_pulse' time
sequence = 'MW_timing'                 #Sequence string
scanStartName = 'start_time'         #Scan start Name
scanEndName = 'end_time'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW, 'Laser':laser,
              'Start\nTrig':start_trig}
sequenceArgs = [t_AOM, t_uW, start_delay, AOM_lag]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())               #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"    #Make param file path
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n"     # %s\t%d\n %s\t%d\n "           #Param file save settings
expParamList =  ['N_scanPts:',N_scanPts, 'Navg:',Nruns, 'Nsamples:',Nsamples, 'start_delay:', start_delay, 'start_time(after_delay):',scannedParam[0], 'end_time(after_delay):',scannedParam[-1], 'uW_freq:',MW_freq, 'uW_power:',MW_power, 't_AOM:',t_AOM, 't_wW:',t_uW, 'AOM_lag:',AOM_lag]      #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 8 parameters

def updateSequenceArgs():
    sequenceArgs = [t_AOM, t_uW, start_delay, AOM_lag]
    return sequenceArgs

def updateExpParamList():
    expParamList =  ['N_scanPts:',N_scanPts, 'Navg:',Nruns, 'Nsamples:',Nsamples, 'start_delay:', start_delay, 'start_time(after_delay):',scannedParam[0], 'end_time(after_delay):',scannedParam[-1], 'uW_freq:',MW_freq, 'uW_power:',MW_power, 't_AOM:',t_AOM, 't_wW:',t_uW, 'AOM_lag:',AOM_lag]      #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]            # 8 parameters
    return expParamList