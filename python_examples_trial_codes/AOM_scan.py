# optimizeReadoutDelay.py
"""
Optimize Readout Delay script

This script can be used to find the optimum delay between the start of the AOM pulse and the start of the DAQ pulse (see step 54 of the protocol). The script plots fluorescence emitted by an NV diamond sample as a function of the scanned delay, and saves the data as a tabulated text file (see below for saving options). 

 User inputs:
 *start_width: shortest delay in scan, in nanoseconds (must be >=5*t_min).
 *end_width:  longest delay in scan, in nanoseconds.
   *N_scanPts: number of points in the scan.
 *t_AOM: duration of AOM pulse, in ns.
 *Nsamples: number of fluorescence measurement samples to take at each delay point.
 *DAQtimeout: amount of time (in seconds) for which the DAQ will wait for the requested number of samples to become available (ie. to be acquired)
 *plotPulseSequence: if set to True, this script will generate a plot of the pulse sequence ouput by the PulseBlaster
 *savePath: path to folder where data will be saved. By default, data is saved in a folder called Saved_Data in the directory where this script is saved
 *saveFileName: file name under which to save the data. This name will later be augmented by the date and time at which the script was run.

"""
#%% Imports
from connectionConfig import *
import SRScontrol as SRSctl
import DAQcontrol as DAQctl
import PBcontrol as PBctl
import sequenceControl as seqCtl
from spinapi import ns,us,ms

from time import localtime, strftime
import matplotlib.pyplot as plt
import numpy as np
import random
import os
import sys

t_min = 1e3/PBclk
#%% User Inputs----------------------------------

# Scan parameters:----------------------------------------------------
start_width = 200  # Start pusle width (in nanoseconds), must be >=5*t_min
end_width = 20e3  # End pulse duration (in nanoseconds)
N_scanPts = 500  # Number of delay steps

# Pulse sequence parameters:----------------------------------------------------
t_readoutDelay = 2.5*us  # AOM pulse duration (in ns)
Nsamples = 250  # Number of fluorescence measurement samples to take at each delay point
DAQtimeout = 10  #DAQ timeout, in seconds

# Plotting options--------------------------------------------------------------
plotPulseSequence = True  # Plot pulse sequence option

# Save options------------------------------------------------------------------
savePath = os.getcwd()+"\\Saved_Data\\"  # Path to folder where data will be saved:
saveFileName = "AOM_scan_"  # File name for data file
# savefig = input("Save Figures (T/F)?? ")
#------------------------- END OF USER INPUT ----------------------------------#

#%% Validate Inputs
try:
    t_AOM = np.linspace(start_width,end_width,N_scanPts,endpoint=True)
    #PB channels
    PBchannels = {'AOM':AOM,'DAQ':DAQ,'STARTtrig':STARTtrig}
    #Make save file path
    dateTimeStr = strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime())
    dataFileName = savePath + saveFileName+ dateTimeStr +".txt"
    #Make param file path
    paramFileName = savePath + saveFileName+dateTimeStr+'_PARAMS'+".txt"
    #Param file save settings
    formattingSaveString = "%s\t%d\n%s\t%d\n%s\t%f\n%s\t%f\n%s\t%f\n%s\t%r\n%s\t%s\n"
    expParamList = ['N_scanPts:',N_scanPts,'Nsamples:',Nsamples,'start_width:',start_width,'end_width:',end_width,'t_readoutDelay:',t_readoutDelay,'plotPulseSequence:',plotPulseSequence,'dataFileName:',dataFileName]
    
    #Validate user input:
    if start_width<(5*t_min):
        print('Error: start_width is too short. Please set start_width>',5*t_min,'ns.')
        sys.exit()
    if start_width%t_min:
        start_width = t_min*round(start_width/t_min)
        print('Warning: start_width is not a multiple of',t_min,'ns. Rounding...\nstart_width now set to:',start_width,'ns.')
    # if end_width%t_min:
    #     end_width = t_min*round(end_width/t_min)
    #     print('Warning: end_width is not a multiple of',t_min,'ns. Rounding...\nend_width now set to:',start_width,'ns.')
    #     t_AOM = np.linspace(start_width,end_width, N_scanPts, endpoint=True) 
    stepSize = t_AOM[1]-t_AOM[0]
    if (stepSize%t_min):
        roundedStepSize = t_min*round(stepSize/t_min)
        end_width  = start_width + (N_scanPts-1)*roundedStepSize
        print('Warning: requested time step is ',stepSize,'ns, which is not an integer multiple of ',t_min,'ns. Rounding step size to the nearest multiple of ',t_min,':\nStep size is now',roundedStepSize,'.\nstart_width=',start_width,' and \nend_width=',end_width)
        t_AOM = np.linspace(start_width,end_width, N_scanPts, endpoint=True)
    
    #Configure DAQ
    DAQclosed = False
    DAQtask = DAQctl.configure_daq(Nsamples)
    
    #%% Plot Sequence
    # imgpath = os.getcwd()+"\\AOM_sweep_fig\\"
    # imageFileName = "AOM_sweep_"
    # if not (os.path.isdir(imgpath)):
    #     os.makedirs(imgpath)
    
    fluorescence = np.zeros(N_scanPts)
    if plotPulseSequence:
        # for i in range(0,N_scanPts):
        instructionList= PBctl.PB_program('AOMsweep', [t_readoutDelay,t_AOM[0]])
        [t_us,channelPulses,yTicks]=seqCtl.plot_sequence(instructionList,PBchannels)
        fig = plt.figure()
        # plt.xlim(0,end_width/1e3+5)
        for channel in channelPulses:
            plt.plot(t_us, list(channel))
        plt.yticks(yTicks, ['AOM','Start_Trig','DAQ_Trig'])
        plt.xlabel('time (us)')
        plt.ylabel('channel')
        plt.title('Pulse Sequence plot')
            # if savefig == 
            # fig.savefig(imgpath+imageFileName+str(i)+'.jpg',dpi=150,bbox_inches='tight',quality=85)
            # plt.close(fig)
        # plt.close('all')
    
    #%% Run Sequence
        #Run readout delay scan:
        for i in range (0, N_scanPts):
            instructionList= PBctl.PB_program('AOMsweep', [t_readoutDelay,t_AOM[i]])  # PB
            print('Scan point ', i+1, ' of ', N_scanPts)
            sig=DAQctl.read_daq(DAQtask,2*Nsamples,DAQtimeout)  # read DAQ
            fluorescence[i] = np.mean(sig)  #Take average of counts
    
        #Close DAQ task:
        DAQctl.close_daq_task(DAQtask)
        DAQclosed = True
    
        #Save data:
        if not (os.path.isdir(savePath)):
            os.makedirs(savePath)
            print('Warning: Save directory did not exist, creating folder named Saved_Data in the working directory. Data will be saved to this directory.')
    
        data = np.array([t_AOM,fluorescence])
        data = data.T
        dataFile = open(dataFileName, 'w')
        for item in data:
            dataFile.write("%.0f\t%f\n" % tuple(item))
        paramFile = open(paramFileName, 'w')
        paramFile.write(formattingSaveString % tuple(expParamList))
        dataFile.close()
        paramFile.close()
    
        #Plot results
        plt.figure(0)
        plt.plot(t_AOM, fluorescence)
        plt.xlabel('AOM pulse duration (ns)')
        plt.ylabel('APD Voltage (V)')
        plt.show()
except KeyboardInterrupt:
        print('User keyboard interrupt. Quitting...')
        sys.exit()
finally:
        if 'SRS' in vars():
            #Turn off SRS output
            SRSctl.disableSRS_RFOutput(SRS)
        if ('DAQtask' in vars()) and  (not DAQclosed):
            #Close DAQ task:
            DAQctl.close_daq_task(DAQtask)
            DAQclosed=True