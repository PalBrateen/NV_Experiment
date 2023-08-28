# -*- coding: utf-8 -*-
# aom_timing_config.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
from time import localtime, strftime
from connectionConfig import *

clk_cyc = 1e3/PBclk #in ns; ONE CLOCK CYCLE
# Enter the time duration > (5 clock cycles) and multiple of ONE CLOCK CYCLE
#%% USER INPUT

#%% Microwave scan parameters:----------------------------------------------------
startfreq = 28.0 /10 * GHz      # Start pulse duration (in nanoseconds)
endfreq = 29.4 /10 * GHz      # End pulse duration (in nanoseconds)
step_size = 1 * MHz
N_scanPts = round((endfreq - startfreq)/step_size + 1)
# N_scanPts = 1501              # Number of pulse length steps
MW_freq = 2.87e9
MW_power = 8          # Microwave power output from SRS(dBm)
t_duration = 150 *ns         # Duration of pi-pulse

#%% Pulse sequence parameters:----------------------------------------------------
# MW_freq = 
t_AOM = 50*us          # Duration of one-half of signal-aquisition
t_gap = 50000 *ns                # the variable parameter

start_delay = t_gap      # Time from which the pulse sequence starts
Nsamples = 500             # Number of FL samples to take at each scan point
Nruns = 1                    # Number of averaging runs to do
# contrast_mode ='ratio_signal_over_reference'       # Contrast mode

# Time gap scan
start_time = 20 *ns
end_time = 2*t_AOM + 2*t_gap #- 100 *ns      # 100 ns is the pulse width..
step_size = 1000*ns
N_scanPts = round((end_time - start_time)/step_size + 1)

# Plotting options--------------------------------------------------------------
# livePlotUpdate = True       # Live plot update option
# plotPulseSequence = True    # Plot pulse sequence option
plotXaxisUnits = ns         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Readout-pulse Time (ns)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2   # Save interval for first scan through all freq points
# saveSpacing_inAverages = 1      # Save interval in averaging runs
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved
saveFileName = "aom_timing"                           # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = False     # Option to do shot by shot contrast normalization:
randomize = True                    # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = np.linspace(start_time,end_time, N_scanPts, endpoint=True)     # scannedParam = 'reading_pulse' time
sequence = 'aom_timing'                 #Sequence string
scanStartName = 'start_time'         #Scan start Name
scanEndName = 'end_time'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW, 'Laser':laser,
              'Start\nTrig':start_trig}        #PB channels
sequenceArgs = [t_AOM, start_delay]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())       #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList =  ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay', start_delay, 'start_time (after delay):',scannedParam[0], 'end_time (after delay):',scannedParam[-1], 't_AOM:',t_AOM, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 9 parameters

def updateSequenceArgs():
    sequenceArgs = [t_AOM, start_delay]
    return sequenceArgs

def updateExpParamList():
    expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay', start_delay, 'start_time (after delay):',scannedParam[0], 'end_time (after delay):',scannedParam[-1], 't_AOM:',t_AOM, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 9 parameters
    return expParamList