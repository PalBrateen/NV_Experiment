# ESRconfig.py
#%% Imports
from spinapi import ns,us,ms
from SGcontrol import Hz, kHz, MHz, GHz
import os
import numpy as np
# from time import localtime, strftime
from connectionConfig import PBclk, laser, samp_clk, start_trig, MW, conv_clk

clk_cyc = 1e3/PBclk #in ns

#%% USER INPUT
N_laser = 10     # Number of cycles of Laser Modulation
N_MW = 1        # Number of cycles of MW Modulation

# Microwave scan parameters:----------------------------------------------------
startFreq = 28 /10 * GHz          # Start frequency (in Hz)
endFreq = 29.4 /10 * GHz             # End frequency (in Hz)
step_size = 1*MHz
N_scanPts = round((endFreq - startFreq)/step_size + 1)             # Number of frequency steps
MW_power = 8         # Microwave power output from SRS(dBm)
SRSdisplay = 'disp2'

# Pulse sequence parameters:----------------------------------------------------
t_duration = 2.5 *us          # Duration of one-half of signal-aquisition half
# t_tot = 2*t_duration
Nsamples = N_MW             # Number of times to repeat the sequence (at each frequency point) for MW modulation
Nruns = 1                    # Number of averaging runs to do (repeating the whole scan)
# contrast_mode = 'ratio_signal_over_reference'       # Contrast mode
AOM_lag = 830*ns
MW_lag = 150*ns

# Plotting options--------------------------------------------------------------
plotXaxisUnits = GHz         # Plot x axis units (Hz, kHz, MHz or GHz)
xAxisLabel = 'Frequency (GHz)'       # Plot x axis label

# Save options------------------------------------------------------------------
# saveSpacing_inScanPts = 2       # Save interval for first scan through all frequency points:
# saveSpacing_inAverages = 1      # Save interval in averaging runs:
savePath = os.getcwd()+"\\Saved_Data\\"         # Path to folder where data will be saved:
saveFileName = "ESR_dig_mod"                           # File name for data file

# Averaging options:------------------------------------------------------------
shotByShotNormalization = int(False)             # Option to do shot by shot contrast normalization:
randomize = int(True)                            # Option to randomize order of scan points

#------------------------- END OF USER INPUT ----------------------------------#

scannedParam = np.linspace(startFreq,endFreq, N_scanPts, endpoint=True)     # scannedParam = Frequency in ESR experiment
sequence = 'esr_dig_mod_seq'                 #Sequence string:
scanStartName = 'startFreq'         #Scan start Name
scanEndName = 'endFreq'             #Scan end Name
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW,
              'Laser':laser, 'Start\nTrig':start_trig}
# PBchannels = {'STARTtrig':start_trig, 'Convert\nClock':conv_clk,
              # 'Data\nSample':samp_clk, 'MW':MW,'Laser':laser}        #PB channels
sequenceArgs = [N_laser, N_MW, t_duration, AOM_lag, MW_lag]         #Sequence args

# dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())               #Make save file path
# dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
# paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"    #Make param file path

formattingSaveString = " %s\t%d\n %s\t%d\n %s\t%d\n %s\t%d\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n %s\t%g\n"       #" %s\t%r\n %s\t%r\n"           #Param file save settings
expParamList =  ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'N_laser', N_laser, 'Step_size:',step_size, 'startFreq:',scannedParam[0], 'endFreq:',scannedParam[-1], 'MW_power:',MW_power, 't_duration:',t_duration]       #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]              # 9 Paramters

def updateSequenceArgs():
    sequenceArgs = [N_laser, N_MW, t_duration, AOM_lag, MW_lag]         #Sequence args
    return sequenceArgs

def updateExpParamList():
    expParamList =  ['N_scanPts:',N_scanPts, 'Nruns:',Nruns, 'Nsamples:',Nsamples, 'N_laser', N_laser, 'Step_size:',step_size, 'startFreq:',scannedParam[0], 'endFreq:',scannedParam[-1], 'MW_power:',MW_power, 't_duration:',t_duration]       #, 'shotByShotNorm:',shotByShotNormalization, 'randomize:',randomize] #,'saveSpacing_inScanPts:',saveSpacing_inScanPts,'saveSpacing_inAverages:',saveSpacing_inAverages]              # 9 Paramters
    return expParamList