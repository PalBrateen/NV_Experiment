# T1ms1_config.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
from time import localtime, strftime
from connectionConfig import PBclk, laser, MW, start_trig, samp_clk, conv_clk

clk_cyc = 1e3/PBclk         #in ns; ONE CLOCK CYCLE

# Enter the time duration > (5 clock cycles) and multiple of ONE CLOCK CYCLE
#%% USER INPUT

# Delay/Interval scan
# Keep start delay > t_pi
start_delay = 0.5 *us
end_delay = 500 *us
step_size = 5 *us
N_scanPts = round((end_delay - start_delay)/step_size + 1)

# Pulse sequence parameters:----------------------------------------------------
t_AOM = 10 *us                # Duration of AOM pulse
AOM_lag = 830 *ns         # Time lag involved in the pulsing of the AOM
rodelay = 300 *ns
# start_delay = 1*us        # Time from which the pulse sequence starts
Nsamples = 500              # Number of FL samples to take at each scan point
Nruns = 1                    # Number of averaging runs to do
t_pi = 300 *ns              # Duration of MW pi pulse

MW_lag = 150 *ns
MW_power = 8
MW_freq =  2822 /1e3 *GHz
# contrast_mode ='ratio_signal_over_reference'       # Contrast mode

# Plotting options--------------------------------------------------------------
# livePlotUpdate = True       # Live plot update option
# plotPulseSequence = True    # Plot pulse sequence option
plotXaxisUnits = ms         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Interval (ms)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2   # Save interval for first scan through all freq points
# saveSpacing_inAverages = 1      # Save interval in averaging runs
savePath = os.getcwd()+"\\Saved_Data\\"     # Path to folder where data will be saved
saveFileName = "T1ms1"                           # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = int(False)     # Option to do shot by shot contrast normalization:
randomize = int(True)                    # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = np.linspace(start_delay, end_delay, N_scanPts, endpoint=True)     # readout-pulse delay
sequence = 'T1ms1'                 #Sequence string
scanStartName = 'start_delay'         #Scan start Name
scanEndName = 'end_delay'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW, 'Laser':laser,
              'Start\nTrig':start_trig}
sequenceArgs = [t_AOM, AOM_lag, rodelay, t_pi, MW_lag]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())       #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"
formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n"          #" %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay:',scannedParam[0], 'end_delay:',scannedParam[-1], 'Step_size:',step_size, 't_AOM:',t_AOM, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag, 'Readout_delay:',rodelay, 't_pi:',t_pi, 'MW_freq:',MW_freq, 'MW_power:',MW_power]      #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 13 parameters

def updateSequenceArgs():
    sequenceArgs = [t_AOM, AOM_lag, rodelay, t_pi, MW_lag]
    return sequenceArgs

def updateExpParamList():
    expParamList = ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'start_delay:',scannedParam[0], 'end_delay:',scannedParam[-1], 'Step_size:',step_size, 't_AOM:',t_AOM, 'AOM_lag:',AOM_lag, 'MW_lag:',MW_lag, 'Readout_delay:',rodelay, 't_pi:',t_pi, 'MW_freq:',MW_freq, 'MW_power:',MW_power]      #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]           # 13 parameters
    return expParamList