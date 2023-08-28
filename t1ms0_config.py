# T1ms0_config.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
from time import localtime, strftime
from connectionConfig import PBclk, laser, start_trig, samp_clk, conv_clk, MW

clk_cyc = 1e3/PBclk         #in ns; ONE CLOCK CYCLE

# Enter the time duration > (5 clock cycles) and multiple of ONE CLOCK CYCLE
#%% USER INPUT

# Delay/Interval scan
start_delay = 0.5 *us
end_delay = 2000 *us
step_size = 50 *us
N_scanPts = round((end_delay - start_delay)/step_size + 1)
# N_scanPts = 5

# Pulse sequence parameters:----------------------------------------------------
t_AOM = 10*us                    # AOM pulse duration (ns)
ro_delay = (2500)*ns      # Readout delay (ns)
# ro_delay = 300*ns
AOM_lag = (1450)*ns     # first parameter = AOM+Preamp lag, 2nd parameter = rise/fall time of the signal as seen in PMT-Preamp-DAQ
# AOM_lag = (800)*ns         # Time lag involved in the pulsing of the AOM
MW_lag = 150 *ns

# start_delay = 1*us        # Time from which the pulse sequence starts

Nsamples = 10000              # Number of FL samples to take at each scan point
Nruns = 1                    # Number of averaging runs to do
# contrast_mode ='ratio_signal_over_reference'       # Contrast mode

MW_power = 8                 # Microwave power output from SRS(dBm)
MW_freq = 2835 /1e3 * GHz     # Microwave frequency (Hz)
t_pi = 160 *ns

# Plotting options--------------------------------------------------------------
# livePlotUpdate = True       # Live plot update option
# plotPulseSequence = True    # Plot pulse sequence option
plotXaxisUnits = ms         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Interval (ms)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2   # Save interval for first scan through all freq points
# saveSpacing_inAverages = 1      # Save interval in averaging runs
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved
saveFileName = "T1ms0"                           # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = int(False)     # Option to do shot by shot contrast normalization:
randomize = int(True)                    # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = np.linspace(start_delay, end_delay, N_scanPts, endpoint=True)     # readout-pulse delay
scannedParam = np.array([500*round(i/500) for i in scannedParam])
# scannedParam = scannedParam[0:2]
# N_scanPts = len(scannedParam)
sequence = 'T1ms0'                 #Sequence string
scanStartName = 'start_delay'         #Scan start Name
scanEndName = 'end_delay'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'uW':MW, 'Laser':laser,
              'Start\nTrig':start_trig}
sequenceArgs = [t_AOM, AOM_lag, ro_delay, MW_lag, t_pi]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())       #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n"          #" %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay:',scannedParam[0], 'end_delay:',scannedParam[-1], 'Step_size:',step_size, 't_AOM:',t_AOM, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag, 't_pi:',t_pi, 'Readout_delay:',ro_delay]      #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 11 parameters

def updateSequenceArgs():
    sequenceArgs = [t_AOM, AOM_lag, ro_delay, MW_lag, t_pi]
    return sequenceArgs

def updateExpParamList():
    expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay:',scannedParam[0], 'end_delay:',scannedParam[-1], 'Step_size:',step_size, 't_AOM:',t_AOM, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag, 't_pi:',t_pi, 'Readout_delay:',ro_delay]      #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 11 parameters
    return expParamList