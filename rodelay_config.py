# -*- coding: utf-8 -*-
# rodelayscan_config.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
from time import localtime, strftime
from connectionConfig import PBclk, laser, start_trig, samp_clk, MW, conv_clk

clk_cyc = 1e3/PBclk         #in ns; ONE CLOCK CYCLE

# Enter the time duration > (5 clock cycles) and multiple of ONE CLOCK CYCLE
#%% USER INPUT

# Pulse sequence parameters:----------------------------------------------------
t_AOM = 5*us                # Duration of AOM pulse
AOM_lag = (830)*ns         # Time lag involved in the pulsing of the AOM
# readout_delay = 350*ns
# '0.8us','1us','1.2us','2us','5us','10us','50us','100us','500us','1ms','2ms','5ms'
interval = 10 *us
# start_delay = 1*us        # Time from which the pulse sequence starts
Nsamples = 100              # Number of FL samples to take at each scan point
Nruns = 1                    # Number of averaging runs to do

# Readout Delay scan
start_rodelay = -AOM_lag
end_rodelay = 6*us
step_size = 100*ns
N_scanPts = round((end_rodelay - start_rodelay)/step_size + 1)

# Plotting options--------------------------------------------------------------
plotXaxisUnits = ns         # Plot x axis units (ns, us)
xAxisLabel = 'Readout delay (ns)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2   # Save interval for first scan through all freq points
# saveSpacing_inAverages = 1      # Save interval in averaging runs
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved
saveFileName = "rodelay"                           # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = False     # Option to do shot by shot contrast normalization:
randomize = True                    # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = np.linspace(start_rodelay, end_rodelay, N_scanPts, endpoint=True)     # readout-pulse delay
# scannedParam2 = 
# scannedParam_name = 'interval' (or 'AOM_lag')
sequence = 'rodelay'                 #Sequence string
scanStartName = 'start_delay'         #Scan start Name
scanEndName = 'end_delay'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW, 'Laser':laser,
              'Start\nTrig':start_trig}
sequenceArgs = [t_AOM, AOM_lag, interval]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())       #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n %s\t%f\n"      #" %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay:',scannedParam[0], 'end_delay:',scannedParam[-1], 't_AOM:',t_AOM, 'AOM_lag:',AOM_lag, 'Delay_Interval:',interval]      # 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 11 parameters

def updateSequenceArgs():
    sequenceArgs = [t_AOM, AOM_lag, interval]
    return sequenceArgs

def updateExpParamList():
    expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay:',scannedParam[0], 'end_delay:',scannedParam[-1], 't_AOM:',t_AOM, 'AOM_lag:',AOM_lag, 'Delay_Interval:',interval]       # 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 11 parameters
    return expParamList