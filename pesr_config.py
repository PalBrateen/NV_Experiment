# Rabiconfig.py
"""List of variables:
    * sequence: (String) Name of the sequence
    * scanStartName:
    * scanEndName:
    * sequenceargs: (List) of arguments (parameters) to be used in creating the pulse sequence
    * PBchannels: (List) of PB channels to be used in the sequence
    * expParamList: (List) of all experimental parameters used in this sequence
"""

#%%Imports
from SGcontrol import Hz, kHz, MHz, GHz
from spinapi import ns,us,ms
import numpy as np
import os
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, conv_clk
# from time import localtime, strftime

clk_cyc = 1e3/PBclk # Time resolution in ns

#-------------------------  USER INPUT  ---------------------------------------#
# See that all the experimental parameters are Python mutable object types

#%% Microwave scan parameters:----------------------------------------------------
startfreq = 28.0 /10 * GHz      # Start pulse duration (in nanoseconds)
endfreq = 29.4 /10 * GHz      # End pulse duration (in nanoseconds)
step_size = 1 * MHz
N_scanPts = round((endfreq - startfreq)/step_size + 1)
# N_scanPts = 1501              # Number of pulse length steps
MW_power = 8          # Microwave power output from SRS(dBm)
t_duration = 112 *ns         # Duration of pi-pulse

#%% Pulse sequence parameters:----------------------------------------------------
t_AOM = 5*us                    # AOM pulse duration (ns)
ro_delay = (2500)*ns      # Readout delay (ns)
# AOM_lag = (1450)*ns     # first parameter = AOM+Preamp lag, 2nd parameter = rise/fall time of the signal as seen in PMT-Preamp-DAQ
AOM_lag = (800)*ns
MW_lag = 150*ns

Nsamples = 5    # Number of signal frames to take at each scanpt
Nruns = 1                        # Number of averaging runs to do

#%% Plotting options--------------------------------------------------------------
# Contrast mode
# contrast_mode ='ratio_signal_over_reference'
# Live plot update option
# livePlotUpdate = True
# Plot pulse sequence option  - set to true to plot the pulse sequence
# plotPulseSequence = True
# Plot x axis unit multiplier (ns, us or ms)
plotXaxisUnits = GHz
# Plot x axis label
xAxisLabel = 'Microwave frequency (GHz)'

#%% Save options------------------------------------------------------------------
saveSpacing_inScanPts = 2       # Save interval for first scan through all pulse length points:
saveSpacing_inAverages = 3      # Save interval in averaging runs:
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved:
saveFileName = "PESR"      # File name for data file

#%% Averaging options:------------------------------------------------------------
# Option to do shot by shot contrast normalization:
shotByShotNormalization = int(False)
# Option to randomize order of scan points
randomize = int(True)
#------------------------- END OF USER INPUT ----------------------------------#

#%%
scannedParam = np.linspace(startfreq, endfreq, N_scanPts, endpoint=True)
sequence = 'pesr_seq'       #Sequence string
scanStartName = 'startfreq'        #Scan start Name
scanEndName = 'endfreq'            #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW, 'Laser':laser,
              'Start\nTrig':start_trig}
sequenceArgs = [t_AOM, ro_delay, AOM_lag, MW_lag, t_duration]      #Sequence args

#Make save file path
# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
#Make param file path
# paramFileName = savePath + saveFileName + dateTimeStr+'_PARAMS.txt'

#Param file save settings
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n"        # " %s\t%r\n %s\t%r\n"
#List of experimental Parameters saved in expParamList
expParamList = ['N_timePts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'startFreq:',scannedParam[0], 'endFreq:',scannedParam[-1], 'MW_power:',MW_power, 't_pi:',t_duration, 't_AOM:',t_AOM, 'ro_delay:',ro_delay, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag]        #, 'shotByShotNormalization:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]       # 11 parameters

"""-------------------------------------------------------------------------------------------------
    Note that: The following functions take the values defined above when called. However, due to their definition, the variables can be changed if called, since they are all "mutable object types" (e.g. Rabiconfig.Nsamples=100).
    If the function is called again, the new value is updated.
-------------------------------------------------------------------------------------------------"""

def updateSequenceArgs():
    sequenceArgs = [t_AOM, ro_delay, AOM_lag, MW_lag, t_duration]
    return sequenceArgs
    
def updateExpParamList():
    expParamList = ['N_timePts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'startFreq:',scannedParam[0], 'endFreq:',scannedParam[-1], 'MW_power:',MW_power, 't_pi:',t_duration, 't_AOM:',t_AOM, 'ro_delay:',ro_delay, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag]        #, 'shotByShotNormalization:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]       # 11 parameters
    return expParamList