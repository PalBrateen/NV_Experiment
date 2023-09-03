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
from spinapi import ns,us,ms
from SGcontrol import GHz, kHz, MHz
import numpy as np
import os
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, conv_clk, camera
# from time import localtime, strftime

clk_cyc = 1e3/PBclk # Time resolution in ns

#-------------------------  USER INPUT  ---------------------------------------#
# See that all the experimental parameters are Python mutable object types

#%% Microwave scan parameters:----------------------------------------------------
startPulseDuration =  10 *ns      # Start pulse duration (in nanoseconds)
endPulseDuration = 500 *ns      # End pulse duration (in nanoseconds)
step_size = 2*ns
N_scanPts = round((endPulseDuration - startPulseDuration)/step_size + 1)
# N_scanPts = 5              # Number of pulse length steps
MW_power = 8                 # Microwave power output from SRS(dBm)
MW_freq = 2870 /1e3 * GHz     # Microwave frequency (Hz)

#%% Pulse sequence parameters:----------------------------------------------------
t_AOM = 5*us                    # AOM pulse duration (ns)
ro_delay = (2500)*ns      # Readout delay (ns)
# AOM_lag = (1450)*ns     # first parameter = AOM+Preamp lag, 2nd parameter = rise/fall time of the signal as seen in PMT-Preamp-DAQ
AOM_lag = (800)*ns
MW_lag = 150*ns

Nsamples = 1    # Number of signal frames to take at each scanpt
Nruns = 1                        # Number of averaging runs to do

#%% Plotting options--------------------------------------------------------------
# Contrast mode
contrast_mode ='ratio_signal_over_reference'
# Live plot update option
# livePlotUpdate = True
# Plot pulse sequence option  - set to true to plot the pulse sequence
# plotPulseSequence = True
# Plot x axis unit multiplier (ns, us or ms)
plotXaxisUnits = ns
# Plot x axis label
xAxisLabel = 'Microwave pulse length (ns)'

#%% Save options------------------------------------------------------------------
saveSpacing_inScanPts = 2       # Save interval for first scan through all pulse length points:
saveSpacing_inAverages = 3      # Save interval in averaging runs:
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved:
saveFileName = "Rabi"      # File name for data file

#%% Averaging options:------------------------------------------------------------
# Option to do shot by shot contrast normalization:
shotByShotNormalization = int(False)
# Option to randomize order of scan points
randomize = int(True)
#------------------------- END OF USER INPUT ----------------------------------#

#%%
scannedParam = np.linspace(startPulseDuration, endPulseDuration, N_scanPts, endpoint=True)
# scannedParam = scannedParam[0:2]
# N_scanPts = len(scannedParam)
# scannedParam2 = 
sequence = 'rabi_seq'       #Sequence string
scanStartName = 'startPulseDuration'        #Scan start Name
scanEndName = 'endPulseDuration'            #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk,'Samp\nCLK':samp_clk, 'MW':MW,'Laser':laser, 'Start\nTrig':start_trig}#,'Camera': camera,}
sequenceArgs = [t_AOM, ro_delay, AOM_lag, MW_lag]      #Sequence args

#Make save file path
# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
#Make param file path
# paramFileName = savePath + saveFileName + dateTimeStr+'_PARAMS.txt'

#Param file save settings
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n" #" %s\t%r\n %s\t%r\n"
#List of experimental Parameters saved in expParamList
expParamList = ['N_timePts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'startPulseDur:',scannedParam[0], 'endPulseDur:',scannedParam[-1], 'MW_power:',MW_power, 'MW_freq:',MW_freq, 't_AOM:',t_AOM, 'ro_delay:',ro_delay, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag]     #, 'shotByShotNormalization:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]       # 11 parameters3

"""-------------------------------------------------------------------------------------------------
    Note that: The following functions take the values defined above when called. However, due to their definition, the variables can be changed if called, since they are all "mutable object types" (e.g. Rabiconfig.Nsamples=100).
    If the function is called again, the new value is updated.
-------------------------------------------------------------------------------------------------"""

def updateSequenceArgs():
    sequenceArgs = [t_AOM, ro_delay, AOM_lag, MW_lag]
    return sequenceArgs
    
def updateExpParamList():
    expParamList = ['N_timePts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'startPulseDur:',scannedParam[0], 'endPulseDur:',scannedParam[-1], 'MW_power:',MW_power, 'MW_freq:',MW_freq, 't_AOM:',t_AOM, 'ro_delay:',ro_delay, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag]     #, 'shotByShotNormalization:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]       # 11 parameters
    return expParamList