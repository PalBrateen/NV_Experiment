# sequencecontrol.py
#%% Imports
from connectionConfig import laser, samp_clk, start_trig, conv_clk, I, Q, MW, PBclk, camera
from spinapi import ns, us, ms
import PBcontrol as PBctrl
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from importlib import import_module
# from importlib import import_module

# expCfgFile = 'esr_config'
# expCfg = import_module(expCfgFile)
conv_clk_sep = 10*us
pulse_width = 100*ns
#%% Operations
""" Create a tuple of 'PBchannel' type with entries:
    1. channelNumber
    2. startTimes
    3. pulseDurations
"""
PBchannel = namedtuple('PBchannel', ['channelNumber', 'startTimes', 'pulseDurations'])

clk_cyc = 1e3/PBclk       # Time resolution in ns

# Short pulse flags: Switch ON bits 21-23 of the Control word to enable short pulse feature
#

ONE_PERIOD = 0x200000           # 23/22/21/20 = 0010b = 2^21d = 200000x
TWO_PERIOD = 0x400000           # 23/22/21/20 = 0100 = 2^22
THREE_PERIOD = 0x600000         # 23/22/21/20 = 0110
FOUR_PERIOD = 0x800000          # 23/22/21/20 = 1000
FIVE_PERIOD = 0xA00000           # 23/22/21/20 = 1010

# plot_sequence():
# accepts:
#    'instructionList' - (List of lists) created by sequenceControl.sequenceEventCataloguer()
#    'channelMasks' - ()
# returns: [t_us,channelPulses,yTicks]

def check_params(config_file):
    expCfg = import_module(config_file)      # import the configuration file as module
    
    # Check that N_scanPts, Nsamples, Nruns are all integers and Nsamples>=1, Nruns>= 1, N_scanPts>=2
    if (not isinstance(expCfg.Nsamples, int)):
        expCfg.Nsamples = int(expCfg.Nsamples)
        print("\x10 \x1b[38;2;240;240;50mWarning\x1b[0m: Nsamples is a 'float'.")
    if  (expCfg.Nsamples<1):
        print('Error: Nsamples must be an integer >= 1.')
        sys.exit()
    if (not isinstance(expCfg.Nruns, int)) or (expCfg.Nruns<1):
        print('Error: Nruns must be an integer >= 1.')
        sys.exit()
    if (not isinstance(expCfg.N_scanPts, int)) or (expCfg.N_scanPts<2):
        print('Error: N_scanPts must be an integer >= 2.')
        sys.exit()
        
    # if len(expCfg.scannedParam) > 1:
    step_size = expCfg.scannedParam[1] - expCfg.scannedParam[0]
    if expCfg.sequence not in ['aom_timing', 'rodelay']:      # Check MW power if the seqeunce is not 'aom_timing', etc
        if expCfg.MW_power >= 9:      # Quit program if the MW power is greater/equals 9 dBm
            print("Input Microwave Power to the RF Amplifier is + "+expCfg.MW_power+" dBm. Reduce it to 8 dBm or less.")
            sys.exit()
        else:
            print("\t\x10 MW Power = "+str(expCfg.MW_power)+" dBm")
            
    if expCfg.sequence in ['aom_timing', 'rodelay']:
        if step_size < 5*clk_cyc:
            print("\x1b[1;33;41mERR: Step size = "+str(step_size)+"... Exiting..."+'\x1b[0m')
            sys.exit()
        # if step_size % clk_cyc != 0:
        #     step_size = 
        for i in range(0, expCfg.N_scanPts):
            expCfg.scannedParam[i] = round(expCfg.scannedParam[i])
            if expCfg.scannedParam[i] % clk_cyc != 0:
                # print('Rounding scannedParam to a multiple of CLK_CYC (2ns)...')
                expCfg.scannedParam[i] += 1     # actually this should be expCfg.scannedParam[i] += round(expCfg.scanned[i] % clk_cyc)
        
    if expCfg.sequence == 'MW_timing':
        if step_size < 5*clk_cyc:
            print("\x1b[1;33;41mERR: Step size = "+str(step_size)+"... Exiting..."+'\x1b[0m')
            sys.exit()
            
        for i in range(0, expCfg.N_scanPts):
            expCfg.scannedParam[i] = round(expCfg.scannedParam[i])
            if expCfg.scannedParam[i] % clk_cyc != 0:
                # print('Rounding scannedParam to a multiple of CLK_CYC (2ns)...')
                expCfg.scannedParam[i] += 1     # actually this should be expCfg.scannedParam[i] += round(expCfg.scanned[i] % clk_cyc)
    #------------------------------------------------------------
    
    #Pulse-sequence parameter checks:
    
    #Check that "IQpadding" is a multiple of clk_cyc and >5*clk_cyc:
    if expCfg.sequence in ['T2seq','XY8seq','correlSpecSeq']:
        if (expCfg.IQpadding < (5*clk_cyc)) or (expCfg.IQpadding%clk_cyc):
            print('Error: IQpadding is set to', expCfg.IQpadding,'which is either <',5*clk_cyc,'or not a multiple of',clk_cyc,'. Please edit IQpadding to ensure that it is >',5*clk_cyc,'ns and a multiple of',clk_cyc,'.')
            
    # Check t_duration in esr_seq is a multiple of (2*clk_cyc): (Why?)
    if expCfg.sequence in ['esr_seq', 'pesr_seq', 'ram_seq', 'T2_seq', 'modesr']:
        if expCfg.t_AOM%(2*clk_cyc):
            print('\x10 Warning: t_duration set to ', expCfg.t_AOM,'ns, which is not an integer multiple of ',(2*clk_cyc),'ns. Rounding t_duration to nearest multiple of ',(2*clk_cyc),'ns...')
            expCfg.t_AOM = (2*clk_cyc)*round(expCfg.t_AOM/(2*clk_cyc))
            print('t_duration now set to ', expCfg.t_AOM,'ns')
    
    #Check that t_readoutDelay and t_AOM are multiples of clk_cyc:
    if expCfg.sequence in ['rabi_seq', 'pser_seq', 'T2_seq', 'XY8seq', 'correlSpecSeq', 'T1seq', 'T1ms0', 'ram_seq']:
        if (expCfg.ro_delay%clk_cyc) or (expCfg.ro_delay < (5*clk_cyc)):
            print('Error: t_readoutDelay is set to ', expCfg.t_readoutDelay,'ns, which is not a multiple of ',clk_cyc,'ns or <',(5*clk_cyc),'ns!')
            sys.exit()
        if (expCfg.t_AOM%clk_cyc) or (expCfg.t_AOM < (5*clk_cyc)):
            print('Error: t_AOM is set to ', expCfg.t_AOM,'ns, which is not a multiple of ',clk_cyc,'ns or <',(5*clk_cyc),'ns!,')
            sys.exit()
          
    #Check that tau0 in the correlation spectroscopy sequence is an integer multiple of 2*clk_cyc:
    # if expCfg.sequence == 'correlSpecSeq':
    #     if expCfg.tau0%(2*clk_cyc):
    #         print('Error: tau0 is set to ', expCfg.tau0,'ns, which is not a multiple of ',(2*clk_cyc),'ns. Please set tau0 to an integer multiple of ',(2*clk_cyc),'ns.')
    #         sys.exit()
            
    
    # #Number of XY8 repeats check:
    # if expCfg.sequence in ['XY8seq', 'correlSpecSeq']:
    #     if expCfg.N<1 or (not isinstance(expCfg.N, int)):
    #         print('Error: number of XY8 repeats, N, must be an integer >=1.')
    #         sys.exit()
    
    
    #Pi-pulse length checks:
    # if expCfg.sequence == 'T1ms1':
    #     if expCfg.t_pi<clk_cyc or expCfg.t_pi % clk_cyc:
    #         print('Error: requested pi pulse length ',expCfg.t_pi,'ns is either <',clk_cyc,'ns or not an integer multiple of ',clk_cyc,'ns.')
    #         sys.exit()
    
    if expCfg.sequence in ['T2seq','XY8seq','correlSpecSeq']:
        # Check if the user has input a pi-pulse length which is shorter than (2*clk_cyc) or not a multiple of 2*clk_cyc:
        if expCfg.t_pi<(2*clk_cyc):
            print('Error: requested pi pulse length=',expCfg.t_pi,'ns is <',(2*clk_cyc),'ns. t_pi must be set to at least',(2*clk_cyc),'ns.')
            sys.exit()
            
        if expCfg.t_pi%(2*clk_cyc):
            print('\x1b[3;33;40m'+'Warning: t_pi set to ', expCfg.t_pi,'ns, which is not an integer multiple of ',(2*clk_cyc),'ns. Rounding t_pi to nearest multiple of ',(2*clk_cyc),'ns...')
            expCfg.t_pi = (2*clk_cyc)*round(float(expCfg.t_pi)/(2*clk_cyc))
            print('t_pi now set to ', expCfg.t_pi,'ns')
            
    # Scan step-size checks: for ESR/Rabi/T1
    # step_size = expCfg.scannedParam[1]-expCfg.scannedParam[0]
    if (expCfg.sequence in ['esr_seq', 'pesr_seq', 'modesr']):
        if ((step_size*1e6)%1):  # Round step to 1uHz if smaller; SRS freq res = 1uHz
            roundedFreqstep_size = (1e-6)*round((1e6)*step_size)
            expCfg.scannedParam[-1] = (expCfg.N_scanPts-1)*roundedFreqstep_size + expCfg.scannedParam[0] 
            print('\x1b[3;33;40m'+'Warning: Requested freq step is ',step_size,'Hz. Not integer multiple of the SRS freq resolution, 1uHz. Rounding step size to the nearest multiple of 1uHz.\nStep size is now',roundedFreqstep_size,'\n',expCfg.scanStartName,'= ',expCfg.scannedParam[0],' and \n',expCfg.scanEndName,'= ',expCfg.scannedParam[-1] + '\x1b[0m')
            expCfg.scannedParam = np.linspace(expCfg.scannedParam[0],expCfg.scannedParam[-1], expCfg.N_scanPts,endpoint= True)
            
    if expCfg.sequence in ['rabi_seq', 'T1seq', 'T1ms0', 'rodelay', 'ram_seq', 'T2_seq'] or (expCfg.sequence == 'T2seq' and expCfg.numberOfPiPulses == 1):   # why to check for numberOfPiPulses==1??
        if step_size<clk_cyc:
            print('Error: requested time step =',step_size,'ns, which is shorter than',clk_cyc,'ns. Please change N_scanPts, or ',expCfg.scanStartName,' and ',expCfg.scanEndName,' to increase time step size.')
            sys.exit()
        else:
            print('\t\x10 Step Size = '+str(step_size)+" > "+str(clk_cyc))
            
        # If requested step size is >clk_cyc but not a multiple of clk_cyc:
        if (step_size%clk_cyc):
            roundedstep_size = clk_cyc*round(step_size/clk_cyc)
            expCfg.scannedParam[-1] = (expCfg.N_scanPts-1)*roundedstep_size + expCfg.scannedParam[0] 
            print('\x1b[38;2;250;200;0mWarning: requested time step is '+str(step_size)+'ns\x1b[0m, \u2260 <int>*'+str(clk_cyc)+'ns.\nRounding step size to the nearest multiple of '+str(clk_cyc)+'...\n\x1b[38;2;250;200;0mStep size now = '+str(roundedstep_size)+'\n'+expCfg.scanStartName+' = '+str(expCfg.scannedParam[0])+'\n'+expCfg.scanEndName,'= ',str(expCfg.scannedParam[-1])+'\x1b[0m')
            expCfg.scannedParam = np.linspace(expCfg.scannedParam[0],expCfg.scannedParam[-1], expCfg.N_scanPts,endpoint= True)
            
    if expCfg.sequence =='rabi_seq':
        # Pulseblaster bug - our PulseBlaster boards do not seem to be able to output 8ns pulses. So, check if we asked for 8ns and remove this point:
        if 8 in expCfg.scannedParam:
            expCfg.scannedParam = list(expCfg.scannedParam)
            expCfg.scannedParam.remove(8)
            expCfg.N_scanPts = len(expCfg.scannedParam)
            print('\x1b[38;2;250;200;0mWarning: will not collect data at 8ns scan point \x1b[0mdue to unofficial reports of a possible issue with some PB boards whereby the instruction for outputting 8ns pulses generates 10ns pulses. Removing the 8ns scan point from the list of scan points.\x1b[0m')
    
    # if (expCfg.sequence == 'XY8seq') or (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses > 1):
    #     #Check if requested scan step is too short or not a multiple of 2*clk_cyc:
    #     if step_size<(2*clk_cyc):
    #         print('Error: requested time step is ',step_size,'ns, which is shorter than ', (2*clk_cyc),'ns. Please change N_scanPts, or ',expCfg.scanStartName,' and ',expCfg.scanEndName,' to increase time step size.')
    #         sys.exit()
            
    #     # If requested step size is >2*clk_cyc but not a multiple of clk_cyc, round to nearest multiple of (2*clk_cyc) and warn user:
    #     if (step_size%(2*clk_cyc)):
    #         roundedstep_size = (2*clk_cyc)*round(step_size/(2*clk_cyc))
    #         expCfg.scannedParam[-1] =  (expCfg.N_scanPts-1)*roundedstep_size + expCfg.scannedParam[0] 
    #         print('\x1b[3;33;40m'+'Warning: requested time step is ',step_size,'ns, which is not an integer multiple of ',(2*clk_cyc),'ns. Rounding step size to the nearest multiple of ',(2*clk_cyc),':\n Step size is now ',roundedstep_size,'\n ',expCfg.scanStartName,'= ',expCfg.scannedParam[0],' and \n',expCfg.scanEndName,'= ',expCfg.scannedParam[-1])
    #         expCfg.scannedParam = np.linspace(expCfg.scannedParam[0],expCfg.scannedParam[-1], expCfg.N_scanPts,endpoint= True)
    
    # Scan-start (minimum delay duration) checks:
    if (expCfg.sequence in ['rabi_seq', 'correlSpecSeq', 'T2_seq']) or (expCfg.sequence == 'T2seq' and expCfg.numberOfPiPulses==1):
    # Check if requested start delay/pulse length is positive and a multiple of "clk_cyc":
        if expCfg.scannedParam[0]<0:
            print('\x1b[38;2;250;2;50mERR: requested ', expCfg.scanStartName,'=', expCfg.scannedParam[0],'is <0. ', expCfg.scanStartName,' must be >=0.')
            sys.exit()
            
        if expCfg.scannedParam[0]%clk_cyc:
            print('\x1b[38;2;250;2;50mERR: ',expCfg.scanStartName,' is set to ', expCfg.scannedParam[0],', which is not a multiple of ',clk_cyc,'ns. Please set', expCfg.scanStartName,' to an integer multiple of ',clk_cyc,'ns.')
            sys.exit()
    
    #     if (expCfg.sequence == 'XY8seq') or (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses > 1):
    #     # Check if requested start delay/pulse length is a multiple of "2*clk_cyc":
    #         if expCfg.scannedParam[0]%(2*clk_cyc):
    #             print('Error: ',expCfg.scanStartName,' is set to ', expCfg.scannedParam[0],', which is not a multiple of ',(2*clk_cyc),'ns. Please set', expCfg.scanStartName,' to an integer multiple of ',(2*clk_cyc),'ns.')
    #             sys.exit()
    
    #     if expCfg.sequence == 'T1seq':
    #         if expCfg.scannedParam[0]<(expCfg.t_readoutDelay + clk_cyc*round((1*us)/clk_cyc)):
    #             print('Error: requested ',expCfg.scanStartName,' is too short.', expCfg.scanStartName,' must be >=0.')
    #             sys.exit()
    #         if expCfg.scannedParam[0]%clk_cyc:
    #             print('Error: ',expCfg.scanStartName,' is set to ', expCfg.scannedParam[0],', which is not a multiple of ',clk_cyc,'ns. Please set', expCfg.scanStartName,' to an integer multiple of',clk_cyc,'ns.')
    #             sys.exit()
    
        if expCfg.sequence in ['T2_seq','XY8seq']:
            # Check if requested start delay is shorter than 3*(5*clk_cyc):
            if expCfg.scannedParam[0]<2*(5*clk_cyc):
                print('Error: requested ',expCfg.scanStartName,' of ',expCfg.scannedParam[0],'ns is too short. For this pulse sequence, ',expCfg.scanStartName,' must be set to at least',3*(5*clk_cyc),'ns')
                sys.exit()
    
        if expCfg.sequence == 'T2seq':
            if (not isinstance(expCfg.numberOfPiPulses, int)) or (expCfg.numberOfPiPulses<1):
                print('Error: numberOfPiPulses must be a positive integer!')
                sys.exit()
            
            if expCfg.numberOfPiPulses == 1:
                if expCfg.scannedParam[0]<(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*clk_cyc)):
                    # Check if start delay is too short to allow for PB timing resolution:
                    print('Error: ',expCfg.scanStartName,' too short. For your pi_pulse length, ',expCfg.scanStartName,' must be at least', (2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*clk_cyc)),'ns.')
                    sys.exit()
                
    #         else: #no of pi pulses>1
    #             if expCfg.scannedParam[0]<(2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*clk_cyc))):
    #                 print('Error: ',expCfg.scanStartName,' too short. For your pi_pulse length, ',expCfg.scanStartName,' must be at least', (2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*clk_cyc))),'ns.')
    #                 sys.exit()
          
    #     if expCfg.sequence == 'XY8seq':
    #         if expCfg.scannedParam[0]<(2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*clk_cyc))):
    #             print('Error: ',expCfg.scanStartName,' too short. For your pi_pulse length, ',expCfg.scanStartName,' must be at least', (2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*clk_cyc))),'ns.')
    #             sys.exit()
            
    #     # Free precession time checks. The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi
    #     # or pi/2 pulse in the T2, XY8 and correlation spectroscopy sequences has to be an integer multiple of clk_cyc. The user-input 
    #     # free precession time is defined as the time between the center of subsequent pulses. Hence, for a given pi-pulse length, 
    #     # we check that the user has selected a starting free precession time which produces an edge-to-edge time that is a multiple of clk_cyc.
    #     # If not, we shift the free precession time vector by clk_cyc/2 and warn the user.
        if (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses == 1):
            if (expCfg.scannedParam[0]-(expCfg.t_pi/4))%clk_cyc:
                expCfg.scannedParam = [x+(clk_cyc/2) for x in expCfg.scannedParam]
                print('Warning: Each element of the scanned time vector has been shifted by '+str(clk_cyc/2)+'ns so that the rising-edge-to-rising-edge spacing between microwave pulses is a multiple of '+str(clk_cyc)+'ns.\
    \nDetails: The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi or pi/2 pulse in the T2 and XY8 sequences \
    has to be an integer multiple of '+str(clk_cyc)+'ns. The user-input free-precession time is defined as the time between the center of subsequent pulses.\
    Your pi pulse length,'+str(expCfg.t_pi)+'ns, produces an edge-to-edge time of'+str(expCfg.scannedParam[0]-(expCfg.t_pi/4))+'ns (at the start of the scan), which is not a multiple of '+str(clk_cyc)+'ns.\
    Hence, we shift the times by '+str(clk_cyc/2)+'ns.')	
    
    #     if (expCfg.sequence =='XY8seq') or (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses > 1):
    #         half_t_delay = expCfg.scannedParam[0]/2
    #         if (half_t_delay-(expCfg.t_pi/4))%clk_cyc:
    #             expCfg.scannedParam = [x+(clk_cyc/2) for x in expCfg.scannedParam]
    #             print('\x1b[3;33;40m'+'Warning: Each element of the scanned time vector has been shifted by ',clk_cyc/2,'ns so that the rising-edge-to-rising-edge spacing between microwave pulses is a multiple of ',clk_cyc,'ns.\
    # \nDetails: The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi or pi/2 pulse in the T2 and XY8 sequences \
    # has to be an integer multiple of ',clk_cyc,'ns. The user-input free-precession time is defined as the time between the center of subsequent pulses.\
    # Your pi pulse length,',expCfg.t_pi,'ns, produces an edge-to-edge time of', (half_t_delay-(expCfg.t_pi/4)),'ns (at the start of the scan), which is not a multiple of ',clk_cyc,'ns.\
    # Hence, we shift the times by ',clk_cyc/2,'ns.')
    
    #     if expCfg.sequence == 'correlSpecSeq':
    #         half_t_delay = expCfg.tau0/2
    #         if (half_t_delay-(expCfg.t_pi/4))%clk_cyc:
    #             expCfg.tau0 = expCfg.tau0 +(clk_cyc/2)*ns
    #             print('\x1b[3;33;40m'+'Warning: tau0 has been shifted by ',clk_cyc/2,'ns so that the rising-edge-to-rising-edge spacing between microwave pulses is a multiple of ',clk_cyc,'ns. tau0 is now set to', expCfg.tau0,'\
    # \nDetails: The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi or pi/2 pulse in the XY8 sequence \
    # has to be an integer multiple of ',clk_cyc,'ns. The user-input tau0 is defined as the time between the center of subsequent pi pulses in the XY8 sequence.\
    # For your pi pulse length,',expCfg.t_pi,'ns, your chose tau0 produces an edge-to-edge time of', half_t_delay-(expCfg.t_pi/4),'ns, which is not a multiple of ',clk_cyc,'ns.\
    # Hence, we shift the tau0 by ',clk_cyc/2,'ns.')
#-----------------------------------------------------------------------
def plot_data(x,y,xlabel,ylabel,formatting,yTicks,yTickLabels, title):
    if 'formatting' in locals():
        plt.plot(x,y, formatting)
    else:
        plt.plot(x,y)
    if 'xlabel' in locals():
        plt.xlabel(xlabel)
    if 'ylabel' in locals():
        plt.ylabel(ylabel)
    if 'yTicks' in locals() and 'yTickLabels' not in locals():
        plt.yticks(yTicks)
    elif 'yTicks' and 'yTickLabels' in locals():
        plt.yticks(yTicks, yTickLabels)
    if 'title' in locals():
        plt.title(title)
    plt.pause(0.0001)

def plot_sequence(instructions, channelMasks):
    """


    Parameters
    ----------
    instructions : list (of lists)
        DESCRIPTION.
    channelMasks : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    """
    # channelMasks = expCfg.PBchannels
    # instructions = instructionList
    scalingFactor = 0.8
    t_ns = [0, 0]
    pulses = {}
    tDone = False
    channelPulses = []
    for channelMask in channelMasks.values():
        pulses[channelMask] = [0, channelMask & instructions[0][0]]
        for i in range(0, len(instructions)):
            currentPulseLength = instructions[i][3]
            if not tDone:
                previousEdgeTime = t_ns[-1]
                nextEdgeTime = previousEdgeTime + currentPulseLength
                t_ns.append(nextEdgeTime)
                t_ns.append(nextEdgeTime)
                if i == (len(instructions)-1):
                    tDone = True
            if i == (len(instructions)-1):
                pulses[channelMask].append(channelMask & instructions[i][0])
                pulses[channelMask].append(channelMask & instructions[i][0])
            else:
                pulses[channelMask].append(channelMask & instructions[i][0])
                pulses[channelMask].append(channelMask & instructions[i+1][0])
        t_us = np.divide(t_ns, 1e3)
        channelPulses.append(list(np.add(math.log(channelMask, 2), np.multiply(list(pulses[channelMask]), scalingFactor/channelMask))))
    yTicks = np.arange(math.log(min(channelMasks.values()), 2), 1+math.log(max(channelMasks.values()), 2), 1)
    return [t_us, channelPulses, yTicks]

def sequence_event_cataloguer(allPBchannels):
    """


    Parameters
    ----------
    allPBchannels : list
        DESCRIPTION.

    Returns
    -------
    None.

    Variables:
    * allPBchannels: (list) of allPBchannels containing information on which PB channel to turn ON at what time and for what duration -> from makeSequence()
    * channel: (PBchannel) an element in the 'allPBchannels' list containing channelNumber, startTimes and endTimes
    * eventCatalog: (dictionary) keys = event times (= start/end times of a pulse)
                                 values = PB register address for the component to be pulsed
    * channelMask
    * endTimes
    * eventTime: (int) gives the start and end times of a particular component (AOM/DAQ/MW)
    * eventChannelMask: (int) stores the PB channel which changes at the 'eventTime'
    * channelBitMask: (dictionary)
"""
    # Catalogs sequence events in terms of consecutive rising edges on the allPBchannels provided. Returns a dictionary, channelBitMasks, whose keys are event (rising/falling edge) times and values are the channelBitMask which indicate which allPBchannels are on at that time.
    eventCatalog = {}  # dictionary where the keys are rising/falling edge times and the values are the channel bit masks which turn on/off at that time

    # (PBchannel) a particular entry in the 'allPBchannels' list
    for aPBchannel in allPBchannels:
        channelMask = aPBchannel.channelNumber
        endTimes = [startTime + pulseDuration for startTime, pulseDuration in zip(aPBchannel.startTimes, aPBchannel.pulseDurations)]
        eventTimes = aPBchannel.startTimes+endTimes
        # print(eventTimes)
        for eventTime in eventTimes:
            # eventTime (int) gives the start and end times of a particular component (AOM/DAQ/MW)
            # (int) stores the PB channel which goes ON at the 'eventTime'
            eventChannelMask = channelMask
            
            if eventTime in eventCatalog.keys():
                # if the eventTime is already present (i.e. the pulse duration is zero), make the eventCatalog entry corresponding to the eventTime to zero, so that the component does not get a pulse. (16^16=0, )
                eventChannelMask = eventCatalog[eventTime] ^ channelMask

                # I'm XORing instead of ORing here in case someone has a zero-length pulse in the sequence. In that case, the XOR ensures that the channel does not turn on at the pulse start/end time. If we did an OR here, it would turn on and only turn off at the next event (which would have been a rising edge), so this would have given unexpected behaviour.
            eventCatalog[eventTime] = eventChannelMask
        # print(eventCatalog)
        # input("Press any key to continue to next PBchannel")
    channelBitMasks = {}
    currentBitMask = 0
    channelBitMasks[0] = currentBitMask
    # print("Event \t CurrentBitMask \t eventCatalog[event]")
    for event in sorted(eventCatalog.keys()):
        # print(str(event)+'\t'+ str(currentBitMask) +'\t' + str(eventCatalog[event])+'\n')
        channelBitMasks[event] = currentBitMask ^ eventCatalog[event]
        currentBitMask = channelBitMasks[event]
        
    # print(channelBitMasks)
    return channelBitMasks

# the fucntion modification is completed on 17062023.. completed...
def param_err_check(instr, sequence, PBchannels, seqArgList, parameter=[0,1], Nscanpts=1):
    # Trial run to check whether the durations of all Inst < 10ns (=5*clk_cyc)
    n_error = 0; param = []; not_param = []
    print("\x10 Checking sequences for errors...")
    for i_scanpt in range (0, Nscanpts):     # scan over all the scannedParam values
        seqArgList[0] = parameter[i_scanpt]
        the_list = PBctrl.PB_program(instr, sequence, seqArgList, err_check=True)
        for i in range(0, len(the_list)):
            instructionList = the_list[i][0]    # eta chai sudhu oi parameter er sequence ta plot korar jonno...
            seq_error_count = the_list[i][1]
            inst_error_no = the_list[i][2]
            error_times = the_list[i][3]
            inst_times = the_list[i][4]
            
            # (For errors) Plot the sequence with the title format: scannedParam[i], inst_no, 
            if seq_error_count > 0:
                n_error += 1
                # -------plot and mark the region where the sequence becomes < 10ns-------
                plt.figure()
                [t_us,channelPulses,yTicks] = plot_sequence(instructionList, PBchannels)
                for channel in channelPulses:
                    plt.plot(t_us, list(channel))
                plt.xlabel('Time (us)')
                plt.ylabel('Channels')
                plt.yticks(yTicks, PBchannels.keys())
                plt.title('Err: scanParam='+str(parameter[i_scanpt])+'. Check Inst #: '+str([i for i in inst_error_no])+ ' @ ' + str([i/1e3 for i in error_times]) + 'us.\n Transitions at: ' + str([i for i in inst_times.values()]), color='r',fontsize=10)
                print('\x10 Removing \x1b[38;2;250;250;0m'+str(parameter[i_scanpt])+'\x1b[0m')
                # ---------plot done--------
                if parameter[i_scanpt] in param:
                    param.remove(parameter[i_scanpt])
                not_param.append(parameter[i_scanpt])
            else:
                if parameter[i_scanpt] not in not_param and parameter[i_scanpt] not in param:
                    param.append(parameter[i_scanpt])
    return [n_error, param]

# ------PulseBlaster Sequences------
"""Variables:
    * sequence: (String) Name of the sequence to be executed
    * args: (List) of time durations
"""

def make_sequence(instr, sequence, args):
    if instr == 'cam':
        var = {'esr_seq':make_esr_seq_camera,   'rabi_seq':make_rabi_seq_camera_FL,
               'T1ms0': make_t1_seq_camera}
    
    elif instr == 'diode':
        var = {'esr_seq':make_esr_seq,   'modesr':make_mod_esr_seq,
               'rabi_seq':make_rabi_seq,    'spin_echo':make_echo_seq_FL,
               'T2seq':makeT2Seq0,  'pesr_seq':make_pulsed_esr_seq,
               'aom_timing':make_aom_timing_seq,    'rodelay':make_rodelayscan_seq,
               'drift_seq':make_drift_analysis_sequence,    'MW_timing':make_MW_timing_seq,
               'T1ms0':make_t1_ms0_seq_MW,  'T1ms1':make_t1_ms1_seq_FL,
               'T2_seq':make_t2_seq_MW,     'ram_seq':make_ramsey_seq_MW,
               'simult_samp':check_simult_sampling,     'double_mod':make_double_mod_sequence_lcm,
               'diff_mod':make_diff_mod_sequence, 'esr_dig_mod_seq':make_dig_mod_odmr_sequence,
               'rabi_dig_mod_seq':make_dig_mod_rabi_sequence,}
        
    # print(var)
    if sequence in var.keys():
        return var[sequence](*args)
    else:
        print('Wrong sequence.. Exiting!!')
        sys.exit()

    
    
    ## elif sequence == 'T1seq':
    ##     return makeT1Seq(*args)
    # elif sequence == 'optimReadoutSeq':
    #     return makeReadoutDelaySweep(*args)
    # elif sequence == 'FLdecay':
    #     return makeFLdecaySeq(*args)
    # elif sequence == 'AOMsweep':
    #     return makeAOMsweep(*args)
    # elif sequence == 'XY8seq':
    #     return makeXY8seq(*args)
    # elif sequence == 'correlSpecSeq':
    #     return makecorrelationSpectSeq(*args)
    # print('Error: requested sequence not recognised.')
    

# All the following functions returns allPBchannels = a (list) of PBallPBchannels containing information on which PB channel to turn ON at what time and for what duration.

#------------------------------------------------------------------------------
def make_double_mod_sequence(MW_freq, laser_freq):
    MW_period = int(round((1/MW_freq)*1e9/clk_cyc)*clk_cyc)       # 1/MW_freq is in seconds
    # MW_period = 1 *us
    laser_period = int(round((1/laser_freq)*1e9/clk_cyc)*clk_cyc)
    # laser_period = 100 *ms
    
    cycles = laser_period/MW_period if laser_period > MW_period else MW_period/laser_period
    cycles = int(cycles) if ((cycles-int(cycles))%10)==0 else sys.exit("\x1b[38;2;250;40;0mError: MW and laser don't have integer cycles")      # If the laser and MW don't have integer ratio (cycles), the pattern cannot be repeated
    
    if laser_period > MW_period:
        nochannel_starts = [laser_period-MW_period/2]
        nochannel_duration = MW_period/2
        laser_starts = [0]
        MW_starts = []
        for i in range(0,cycles,1):
            # print(round((i*MW_on_time)/MW_on_time)*2*MW_on_time)
            MW_starts.append(i*MW_period)
    
    else:
        nochannel_starts = [MW_period-laser_period/2]
        nochannel_duration = laser_period/2
        MW_starts = [0]
        laser_starts = []
        for i in range(0,cycles,1):
            # print(round((i*MW_on_time)/MW_on_time)*2*MW_on_time)
            laser_starts.append(i*laser_period)
        
    AOMchannel = PBchannel(laser, laser_starts, [laser_period/2 for i in range(0,len(laser_starts))])
    MWchannel = PBchannel(MW, MW_starts, [MW_period/2 for i in range(0,len(MW_starts))])
    
    nochannel = PBchannel(0, nochannel_starts, [nochannel_duration])
    allPBchannels = [AOMchannel, MWchannel, nochannel]
    return allPBchannels

#------------------------------------------------------------------------------
# There may be cases when the ratio is not an integer but still we can program the PB. This will happen only if the main limitation of PB can be avoided -- instruction length >= 10ns.
# The function below takes the LCM of the periods of both and generates a 
def make_double_mod_sequence_lcm(MW_freq, laser_freq):
    MW_period = int(round((1/MW_freq)*1e9/clk_cyc)*clk_cyc)       # 1/MW_freq is in seconds
    # MW_period = 1 *us
    laser_period = int(round((1/laser_freq)*1e9/clk_cyc)*clk_cyc)
    # laser_period = 100 *ms
    
    common_time_period = np.lcm(MW_period, laser_period)
    print("Common time period = "+str(common_time_period)+" us.")
    MW_cycles = common_time_period/MW_period
    laser_cycles = common_time_period/laser_period
    
    if ((MW_cycles-int(MW_cycles))%10)==0 and ((laser_cycles-int(laser_cycles))%10)==0:
        MW_cycles = int(MW_cycles); laser_cycles = int(laser_cycles)
    elif ((laser_cycles-int(laser_cycles))%10)!=0:
        sys.exit("\x1b[38;2;250;40;0mError: Laser doesn't have integer cycles in common time period")
    elif ((MW_cycles-int(MW_cycles))%10)!=0:
        sys.exit("\x1b[38;2;250;40;0mError: MW and laser don't have integer cycles")
        
    MW_starts = []
    for i in range(0,MW_cycles,1):
        # print(round((i*MW_on_time)/MW_on_time)*2*MW_on_time)
        MW_starts.append(i*MW_period)

    laser_starts = []
    for i in range(0,laser_cycles,1):
        # print(round((i*MW_on_time)/MW_on_time)*2*MW_on_time)
        laser_starts.append(i*laser_period)
        
    AOMchannel = PBchannel(laser, laser_starts, [laser_period/2 for i in range(0,len(laser_starts))])
    MWchannel = PBchannel(MW, MW_starts, [MW_period/2 for i in range(0,len(MW_starts))])
    
    nochannel_starts = [(common_time_period-min(laser_period,MW_period))]
    nochannel_duration = min(laser_period, MW_period)
    nochannel = PBchannel(0, nochannel_starts, [nochannel_duration])
    allPBchannels = [AOMchannel, MWchannel, nochannel]
    return allPBchannels

#------------------------------------------------------------------------------
def make_diff_mod_sequence(MW_freq, laser_freq, diff_freq):    
    if diff_freq == 0:
        allPBchannels = make_double_mod_sequence_lcm(MW_freq, laser_freq)
    else:
        # if MW_freq == 0 or laser_freq == 0
        # 1/**_freq is in seconds; then rounding it to a multiple of clk cycle (=2ns)
        MW_period = int(round((1/MW_freq)*1e9/clk_cyc)*clk_cyc)
        laser_period = int(round((1/laser_freq)*1e9/clk_cyc)*clk_cyc) if laser_freq != 0 else MW_period
        diff_period = int(round((1/diff_freq)*1e9/clk_cyc)*clk_cyc)
        common_time_period = np.lcm(diff_period, np.lcm(MW_period, laser_period))
        
        print("Common time period = "+str(common_time_period))
        # Define the no of cycles of each signal
        MW_cycles = common_time_period/MW_period
        laser_cycles = common_time_period/laser_period if laser_freq!=0 else 1
        diff_cycles = common_time_period/diff_period
        # **_cycles should be an integer in the common_time_period. However, this should always be the case
        if ((MW_cycles-int(MW_cycles))%10)==0 and ((laser_cycles-int(laser_cycles))%10)==0:
            MW_cycles = int(MW_cycles); laser_cycles = int(laser_cycles); diff_cycles = int(diff_cycles)
        elif ((laser_cycles-int(laser_cycles))%10)!=0:
            sys.exit("\x1b[38;2;250;40;0mError: Laser doesn't have integer cycles in common time period")
        elif ((MW_cycles-int(MW_cycles))%10)!=0:
            sys.exit("\x1b[38;2;250;40;0mError: MW and laser don't have integer cycles")
            
        MW_starts = []
        for i in range(0,MW_cycles,1):
            MW_starts.append(i*MW_period)
    
        laser_starts = []
        for i in range(0,laser_cycles,1):
            laser_starts.append(i*laser_period)
        
        diff_starts = []
        for i in range(0,diff_cycles,1):
            diff_starts.append(i*diff_period)
        
        # If laser freq = 0, then the laser starts at the beginning (only once) and the duration of ON is MW_period
        AOMchannel = PBchannel(laser, laser_starts, [laser_period/2 for i in range(0,len(laser_starts))]) if laser_freq!=0 else PBchannel(laser, laser_starts, [MW_period for i in range(0,len(laser_starts))])
        
        MWchannel = PBchannel(MW, MW_starts, [MW_period/2 for i in range(0,len(MW_starts))])
        diffchannel = PBchannel(diff, diff_starts, [diff_period/2 for i in range(0,len(diff_starts))])
        
        nochannel_starts = [(common_time_period-min(laser_period,MW_period,diff_period))]
        nochannel_duration = min(laser_period,MW_period,diff_period)
        nochannel = PBchannel(0, nochannel_starts, [nochannel_duration])
        
        allPBchannels = [AOMchannel, MWchannel, diffchannel, nochannel]
    return allPBchannels

#------------------------------------------------------------------------------
def makeFLdecaySeq(t_readoutDelay, t_AOM):
    t_startTrig = clk_cyc*round(200*ns/clk_cyc)
    start_delay = clk_cyc*round(2*us/clk_cyc)  # -t_startTrig
    t_readout = clk_cyc*round(300*ns/clk_cyc)
    # PBchannel(AOM,[start_delay], [t_AOM])
    AOMchannel = PBchannel(laser, [start_delay], [t_AOM])
    # PBchannel(DAQ,[start_delay+t_AOM+t_readoutDelay],[t_readout])
    DAQchannel = PBchannel(samp_clk, [t_startTrig+t_readoutDelay], [t_readout])
    # PBchannel(STARTtrig,[start_delay+t_AOM],[14*us])          #[2*clk_cyc*round(5*us/clk_cyc)+t_startTrig]
    STARTtrigchannel = PBchannel(start_trig, [t_startTrig], [2*clk_cyc*round(6*us/clk_cyc)+t_startTrig])
    allPBchannels = [AOMchannel, DAQchannel, STARTtrigchannel]
    return allPBchannels

#------------------------------------------------------------------------------
def makeReadoutDelaySweep(t_readoutDelay, t_AOM):
    t_startTrig = clk_cyc*round(100*ns/clk_cyc)
    start_delay = clk_cyc*round(2*us/clk_cyc)  # -t_startTrig
    t_readout = clk_cyc*round(300*ns/clk_cyc)
    # PBchannel(AOM,[start_delay], [t_AOM])
    AOMchannel = PBchannel(laser, [start_delay], [t_AOM])
    # PBchannel(DAQ,[t_startTrig+t_readoutDelay],[t_readout])
    DAQchannel = PBchannel(
        samp_clk, [start_delay+t_AOM+t_readoutDelay], [t_readout])
    # [clk_cyc*round(5*us/clk_cyc)]) #[2*clk_cyc*round(5*us/clk_cyc)+t_startTrig])          #[2*clk_cyc*round(5*us/clk_cyc)+t_startTrig]
    STARTtrigchannel = PBchannel(
        start_trig, [start_delay+t_AOM], [2*clk_cyc*round(5*us/clk_cyc)+t_startTrig])
    allPBchannels = [AOMchannel, DAQchannel, STARTtrigchannel]
    return allPBchannels

#------------------------------------------------------------------------------
def makeAOMsweep(t_readoutDelay, t_AOM):
    t_startTrig = clk_cyc*round(300*ns/clk_cyc)
    start_delay = clk_cyc*round(5*us/clk_cyc)  # -t_startTrig
    t_readout = clk_cyc*round(300*ns/clk_cyc)
    # PBchannel(AOM,[start_delay], [t_AOM])
    AOMchannel = PBchannel(laser, [start_delay], [t_AOM])
    # PBchannel(DAQ,[start_delay+t_AOM+t_readoutDelay],[t_readout])
    DAQchannel = PBchannel(
        samp_clk, [start_delay+t_AOM+t_readoutDelay], [t_readout])
    # PBchannel(STARTtrig,[start_delay+t_AOM],[14*us])          #[2*clk_cyc*round(5*us/clk_cyc)+t_startTrig]
    STARTtrigchannel = PBchannel(
        start_trig, [start_delay+t_AOM], [2*clk_cyc*round(5*us/clk_cyc)+t_startTrig])
    allPBchannels = [AOMchannel, DAQchannel, STARTtrigchannel]
    return allPBchannels

#------------------------------------------------------------------------------
def make_esr_seq_camera(seq_dur):
    delay = 0.5*us
    # seq_dur = 2*seq_dur;
    # readout_width = clk_cyc*round(100*ns/clk_cyc)
    laser_start_times = [delay, delay]; laser_pulse_durations = [seq_dur,seq_dur]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel];
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.

    MW_start_times = [delay]; MW_pulse_durations = [seq_dur]
    MWchannel = PBchannel(MW, [MW_start_times[0]], [MW_pulse_durations[0]])
    allPBchannels_sig.extend([MWchannel])
    # MWchannel = PBchannel(MW, [MW_start_times[1]], [MW_pulse_durations[1]])
    # allPBchannels_ref.extend([MWchannel])
    
    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

#------------------------------------------------------------------------------
def make_esr_seq(seq_dur):
    seq_dur = 2*seq_dur;    trig_width = clk_cyc*round(100*ns/clk_cyc)
    # readout_width = clk_cyc*round(100*ns/clk_cyc)
    readout_buffer = clk_cyc*round(100*ns/clk_cyc)
    
    pd_pulse = [seq_dur/2-readout_buffer, seq_dur-readout_buffer] # seq_dur/2-readout_buffer
    # pd_pulse.extend([(seq_dur/2-readout_buffer)/2, 1.5*(seq_dur/2-readout_buffer)])
    # pd_pulse.extend([seq_dur/2-readout_buffer-1*us, seq_dur-readout_buffer-1*us])
    laser_channel = PBchannel(laser, [0], [seq_dur])
    MW_channel = PBchannel(MW, [0], [seq_dur/2]) # seq_dur/2
    # conv_clk_channel = PBchannel(conv_clk, [pd_pulse[0]-conv_clk_sep, pd_pulse[0], pd_pulse[1]-conv_clk_sep, pd_pulse[1]], [trig_width for i in range(0,4)])    
    samp_clk_channel = PBchannel(samp_clk, [pulse for pulse in pd_pulse], [trig_width for i in range(0,len(pd_pulse))]) #-conv_clk_sep-pulse_width
    start_trig_channel = PBchannel(start_trig, [0], [trig_width])
    
    allPBchannels = [laser_channel, samp_clk_channel, MW_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    return allPBchannels
#------------------------------------------------------------------------------

# def make_mod_esr_seq(t_AOM):
#     starttrig_width = clk_cyc*round(200*ns/clk_cyc)
#     readout_width = clk_cyc*round(300*ns/clk_cyc)
#     readout_buffer = clk_cyc*round(500*ns/clk_cyc)     # Readout pulse comes at 500ns before the end of t_AOM
#     AOMchannel = PBchannel(AOM, [0], [2*t_AOM])
#     MWchannel = PBchannel(MW, [0, 2*t_AOM], [t_AOM, t_AOM])
#     DAQchannel = PBchannel(DAQ, [t_AOM-readout_buffer, 2*t_AOM-readout_buffer, 3*t_AOM-readout_buffer], [readout_width, readout_width, readout_width])
#     STARTtrigchannel = PBchannel(STARTtrig, [0], [starttrig_width])
#     allPBchannels = [AOMchannel, DAQchannel, MWchannel, STARTtrigchannel]
#     return allPBchannels
#------------------------------------------------------------------------------
def make_mod_esr_seq(t_AOM):
    starttrig_width = clk_cyc*round(200*ns/clk_cyc)
    readout_width = clk_cyc*round(300*ns/clk_cyc)
    readout_buffer = clk_cyc*round(500*ns/clk_cyc)     # Readout pulse comes at 500ns before the end of t_AOM
    laser_ch = PBchannel(laser, [0], [3*t_AOM])
    MW_ch = PBchannel(MW, [0, 2*t_AOM], [t_AOM, 2*t_AOM])
    samp_clk_ch = PBchannel(samp_clk, [t_AOM-readout_buffer, 2*t_AOM-readout_buffer, 3*t_AOM-readout_buffer, 4*t_AOM-readout_buffer], [readout_width, readout_width, readout_width, readout_width])
    start_trig_ch = PBchannel(start_trig, [0], [starttrig_width])
    allPBchannels = [laser_ch, samp_clk_ch, MW_ch, start_trig_ch]
    return allPBchannels

#------------------------------------------------------------------------------
def make_dig_mod_odmr_sequence(N_laser, N_MW, t_AOM, AOM_lag, MW_lag):
    t_MW = 0
    RO_width = 300*ns; conv_clk_width = 20*ns
    pulse_width = conv_clk_width
    sig_ref_buffer = 1000*ns
    t_total = (t_AOM + t_MW + RO_width + 2*pulse_width + sig_ref_buffer)
    
    a = [[0+i*t_total, (t_AOM+t_MW+i*t_total)] for i in range(0,N_laser,1)]
    laser_starts = [a[i][j] for i in range(0,len(a),1) for j in range(0,len(a[0]),1)]
    a = [[N_laser*t_total+i*t_total, (N_laser*t_total+t_AOM+t_MW+i*t_total)] for i in range(0,N_laser,1)]
    laser_starts += [a[i][j] for i in range(0,len(a),1) for j in range(0,len(a[0]),1)]
    
    a = [[t_AOM, RO_width] for i in range(0,2*N_laser)]
    laser_dur = [a[i][j] for i in range(0,len(a)) for j in range(0,len(a[0]))]
    
    MW_starts = list(np.zeros(N_laser)) + [N_laser*t_total+AOM_lag-MW_lag+i*t_total for i in range(0,N_laser,1)]
    MW_dur = list(np.zeros(N_laser)) + [t_AOM+t_MW for i in range(0,N_laser)]
    
    # samp_clk_starts = [0+(t_AOM)+t_MW+i*t_total for i in range(0,N_laser)]
    # samp_clk_starts += [N_laser*t_total+(t_AOM)+t_MW+i*t_total for i in range(0,N_laser)]
    # samp_clk_dur = [RO_width for i in range(0,2*N_laser,1)]
    # samp_clk_dur = [a[i][j] for i in range(0, len(a)) for j in range(0, len(a[0]))]
    
    a = [[0+(t_AOM+AOM_lag)+t_MW+RO_width+i*t_total, 0+(t_AOM+AOM_lag)+t_MW+RO_width+pulse_width+sig_ref_buffer+i*t_total] for i in range(0,N_laser,1)]
    samp_clk_starts = [a[i][j] for i in range(0,len(a)) for j in range(0,len(a[0]))]
    a = [[N_laser*t_total+(t_AOM+AOM_lag)+t_MW+RO_width+i*t_total, N_laser*t_total+(t_AOM+AOM_lag)+t_MW+RO_width+pulse_width+sig_ref_buffer+i*t_total] for i in range(0,N_laser,1)]
    samp_clk_starts += [a[i][j] for i in range(0,len(a)) for j in range(0,len(a[0]))]
    samp_clk_dur = [pulse_width for i in range(0,2*2*N_laser)]
    
    laser_ch = PBchannel(laser, laser_starts, laser_dur)
    MW_ch = PBchannel(MW, MW_starts, MW_dur)
    start_trig_ch = PBchannel(start_trig, [0], [200*ns])
    samp_clk_ch = PBchannel(samp_clk, samp_clk_starts, samp_clk_dur)
    # conv_clk_ch = PBchannel(conv_clk, conv_clk_starts, conv_clk_dur)
    
    allPBchannels = [laser_ch, MW_ch, start_trig_ch, samp_clk_ch]
    return allPBchannels
    
#------------------------------------------------------------------------------
# def make_aom_check_seq(t_readout_pulse):
#     AOM_delay = 800*ns
#     start_delay = 1*us          # AOM first pulse
#     start_trig_width = 300*ns
#     daq_readout_width = 300*ns
#     t_AOM = 2.5*us
#     t_gap = 250*ns
    
#     AOM_readout_pulse = start_delay + t_AOM + AOM_delay + t_gap         # AOM 2nd pulse
#     # total_duration = start_delay + t_AOM + t_gap + t_AOM + AOM_delay
    
#     AOMchannel = PBchannel(AOM,[start_delay, AOM_readout_pulse],[t_AOM, t_AOM])
#     STARTtrigchannel = PBchannel(STARTtrig, [AOM_readout_pulse-t_gap], [start_trig_width])
#     DAQchannel = PBchannel(DAQ, [AOM_readout_pulse+AOM_delay+t_readout_pulse], [daq_readout_width])
    
#     allPBchannels = [AOMchannel, STARTtrigchannel, DAQchannel]
#     return allPBchannels
#------------------------------------------------------------------------------
# def make_aom_timing_seq(reading_pulse, t_AOM, start_delay):
#     start_trig_width = 100*ns       # The time of the pulse rising edge matters
#     daq_readout_width = 100*ns      # The time of the pulse rising edge matters
#     t_gap = 20*us
    
#     AOM_second_pulse = start_delay + t_AOM + t_gap         # AOM 2nd pulse
#     # total_duration = start_delay + t_AOM + t_gap + t_AOM + 2*us
#     # if reading_pulse > t_AOM:
#     #     reading_pulse = t_AOM
    
#     laserchannel = PBchannel(laser,[start_delay, AOM_second_pulse],[t_AOM, t_AOM])
#     # laserchannel = PBchannel(laser,[start_delay],[])
#     # start_trig_channel = PBchannel(start_trig, [start_delay], [start_trig_width])
#     start_trig_channel = PBchannel(start_trig, [0], [start_trig_width])
#     # samp_clk_channel = PBchannel(samp_clk, [start_delay+reading_pulse], [daq_readout_width])
    # samp_clk_channel = PBchannel(samp_clk, [0+reading_pulse], [daq_readout_width])
    
#     allPBchannels = [laserchannel, start_trig_channel, samp_clk_channel]
#     return allPBchannels

def make_aom_timing_seq(reading_pulse, t_AOM, t_gap):       # for Rabi...
    allPBchannels = []
    
    t_drive = t_gap
    
    # laser_channel = PBchannel(laser, [t_drive],[])
    laser_channel = PBchannel(laser, [t_drive, 2*t_drive+t_AOM], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    # # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    
    # if t_MW <= 5*clk_cyc and t_MW > 0:
    #     # (PBchannel) for MW pulse
    #     MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag], [5*clk_cyc])
    #     # SHORT PULSE duration
    #     shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
    #     shortPulseChannel = PBchannel(shortpulseFLAG, [0*us+AOM_lag-MW_lag], [5*clk_cyc])     # (PBchannel) for SHORT MW pulse
    #     allPBchannels = [shortPulseChannel, MWchannel]  # Short pulse feature
    # else:
    #     MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag], [t_MW])
    #     allPBchannels = [MWchannel]
    
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    # conv_clk_channel = PBchannel(conv_clk, [0], [])
    samp_clk_channel = PBchannel(samp_clk, [0+reading_pulse], [pulse_width])
    
    allPBchannels.extend([laser_channel, samp_clk_channel, start_trig_channel])
    # allPBchannels.extend([conv_clk_channel])
    # print(allPBchannels)
    return allPBchannels
    
# def make_aom_timing_seq(reading_pulse, t_AOM, t_gap):       # for PESR...
#     # (t_AOM, t_ro_delay, t_pi):
#     """
#     Make pulse sequence for Pulsed ODMR

#     Parameters
#     ----------
#     t_MW : float
#         Width of MW pulse.
#     t_AOM : float
#         AOM 'ON' time or AOM width.
#     t_ro_delay : float
#         Length of (AOM) readout pulse.
#     AOM_lag: float
        
#     MW_lag: float
        
#     Returns
#     -------
#     allPBchannels : list
#         A list of PBchannel object types.

#     """
#     t_pi = t_gap
#     t_flip = t_pi
   
#     laser_channel = PBchannel(laser, [t_flip, 2*t_flip+t_AOM], [t_AOM,t_AOM])
#     start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
#     # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
#     MWchannel = PBchannel(MW, [AOM_lag-MW_lag], [t_flip])
#     # conv_clk_channel = PBchannel(conv_clk,[apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
#     # samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width for i in range(0,2)])
    
#     samp_clk_channel = PBchannel(samp_clk,[apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,2)])
#     allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]
#------------------------------------------------------------------------------
def make_MW_timing_seq(reading_pulse, t_AOM, t_MW, start_delay, AOM_lag):
    # start_trig_width = 300*ns       # The time of the pulse rising edge matters
    # daq_readout_width = 300*ns      # The time of the pulse rising edge matters
    # AOM_to_MW = 0*us
    first_half = start_delay + AOM_lag + t_MW + 10*us        # 5*us: for the readout-pulse time, as it is being scanned
    if (reading_pulse + start_delay) > first_half:
        first_half = reading_pulse + start_delay
    laser_channel = PBchannel(laser, [0],[2*first_half])
    start_trig_channel = PBchannel(start_trig, [0],[pulse_width])
    samp_clk_channel = PBchannel(samp_clk, [start_delay+reading_pulse, first_half + start_delay + reading_pulse],[pulse_width, pulse_width])
    MWchannel = PBchannel(MW, [start_delay+2*us], [t_MW])
    allPBchannels = [laser_channel, start_trig_channel, samp_clk_channel, MWchannel]
    return allPBchannels
#------------------------------------------------------------------------------c
def make_rabi_seq_camera_FL(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag):
    """
    Make pulse sequence for Rabi oscillations

    Parameters
    ----------
    t_MW : float
        Width of MW pulse.
    t_AOM : float
        AOM 'ON' time or AOM width.
    t_ro_delay : float
        Length of (AOM) readout pulse.
    AOM_lag: float
        
    MW_lag: float
        
    Returns
    -------
    allPBchannels : list
        A list of PBchannel object types.

    """
    if t_MW < 10:
        t_drive = 10+0*ns
    else:
        t_drive = t_MW + 0*ns
    
    laser_start_times = [t_drive, t_drive]; laser_pulse_durations = [t_AOM, t_AOM]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel];
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    
    if t_MW <= 5*clk_cyc and t_MW > 0:
        MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [5*clk_cyc for i in range(0,2)]
        # SHORT PULSE duration
        shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
        shortPulseChannel = PBchannel(shortpulseFLAG, [MW_start_times[0]], [MW_pulse_durations[0]])
        allPBchannels_sig.extend([shortPulseChannel])
        shortPulseChannel = PBchannel(shortpulseFLAG, [MW_start_times[1]], [MW_pulse_durations[1]])
        # allPBchannels_ref.extend([shortPulseChannel])
    else:
        MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [t_MW for i in range(0,2)]
    MWchannel = PBchannel(MW, [MW_start_times[0]], [MW_pulse_durations[0]])
    allPBchannels_sig.extend([MWchannel])
    MWchannel = PBchannel(MW, [MW_start_times[1]], [MW_pulse_durations[1]])
    # allPBchannels_ref.extend([MWchannel])
    
    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

def make_rabi_seq_camera_MW(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag):
    """
    Make pulse sequence for Rabi oscillations

    Parameters
    ----------
    t_MW : float
        Width of MW pulse.
    t_AOM : float
        AOM 'ON' time or AOM width.
    t_ro_delay : float
        Length of (AOM) readout pulse.
    AOM_lag: float
        
    MW_lag: float
        
    Returns
    -------
    allPBchannels : list
        A list of PBchannel object types.

    """
    if t_MW < 10:
        t_drive = 10 +0*ns
    else:
        t_drive = t_MW + 0*ns
    
    laser_start_times = [t_drive, t_drive]; laser_pulse_durations = [t_AOM, t_AOM]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel];
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    
    if t_MW <= 5*clk_cyc and t_MW > 0:
        MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [5*clk_cyc for i in range(0,2)]
        # SHORT PULSE duration
        shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
        shortPulseChannel = PBchannel(shortpulseFLAG, [MW_start_times[0]], [MW_pulse_durations[0]])
        allPBchannels_sig.extend([shortPulseChannel])
        shortPulseChannel = PBchannel(shortpulseFLAG, [MW_start_times[1]], [MW_pulse_durations[1]])
        # allPBchannels_ref.extend([shortPulseChannel])
    else:
        MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [t_MW for i in range(0,2)]
    MWchannel = PBchannel(MW, [MW_start_times[0]], [MW_pulse_durations[0]])
    allPBchannels_sig.extend([MWchannel])
    MWchannel = PBchannel(MW, [MW_start_times[1]], [MW_pulse_durations[1]])
    # allPBchannels_ref.extend([MWchannel])
    
    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

#------------------------------------------------------------------------------
def make_rabi_seq(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag):
    """
    Make pulse sequence for Rabi oscillations

    Parameters
    ----------
    t_MW : float
        Width of MW pulse.
    t_AOM : float
        AOM 'ON' time or AOM width.
    t_ro_delay : float
        Length of (AOM) readout pulse.
    AOM_lag: float
        
    MW_lag: float
        
    Returns
    -------
    allPBchannels : list
        A list of PBchannel object types.

    """
    if t_MW < 10:
        t_drive = 10+0*ns
    else:
        t_drive = t_MW + 0*ns

    apd_pulse = [t_drive+AOM_lag+t_ro_delay, 2*t_drive+t_AOM+AOM_lag+t_ro_delay]
    
    # laser_channel = PBchannel(laser, [t_drive],[])
    laser_channel = PBchannel(laser, [t_drive, 2*t_drive+t_AOM], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    allPBchannels = []
    if t_MW <= 5*clk_cyc and t_MW > 0:
        # (PBchannel) for MW pulse
        MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag], [5*clk_cyc])
        # SHORT PULSE duration
        shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
        shortPulseChannel = PBchannel(shortpulseFLAG, [0*us+AOM_lag-MW_lag], [5*clk_cyc])     # (PBchannel) for SHORT MW pulse
        allPBchannels = [shortPulseChannel, MWchannel]  # Short pulse feature
    else:
        # MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag], [t_MW])
        MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag, 0*us+AOM_lag-MW_lag+t_drive+t_AOM], [t_MW,t_MW])
        allPBchannels = [MWchannel]
    
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    # conv_clk_channel = PBchannel(conv_clk, [0], [])
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,len(apd_pulse))])
    allPBchannels.extend([laser_channel, samp_clk_channel, start_trig_channel])
    # allPBchannels.extend([conv_clk_channel])
    # print(allPBchannels)
    return allPBchannels

#------------------------------------------------------------------------------
def make_dig_mod_rabi_sequence(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag, N_laser, N_MW):
    RO_width = 300*ns; conv_clk_width = 20*ns
    pulse_width = conv_clk_width
    sig_ref_buffer = 500*ns
    if t_MW < 10:
        t_MW = 10 +0*ns
    else:
        t_MW = t_MW +0*ns
    
    t_total = (t_AOM + t_MW + RO_width + 2*pulse_width + sig_ref_buffer)
    
    a = [[0+i*t_total, (t_AOM+t_MW+i*t_total)] for i in range(0,N_laser,1)]
    laser_starts = [a[i][j] for i in range(0,len(a),1) for j in range(0,len(a[0]),1)]
    a = [[N_laser*t_total+i*t_total, (N_laser*t_total+t_AOM+t_MW+i*t_total)] for i in range(0,N_laser,1)]
    laser_starts += [a[i][j] for i in range(0,len(a),1) for j in range(0,len(a[0]),1)]
    
    a = [[t_AOM, RO_width] for i in range(0,2*N_laser)]
    laser_dur = [a[i][j] for i in range(0,len(a)) for j in range(0,len(a[0]))]
    
    MW_starts = list(np.zeros(N_laser)) + [N_laser*t_total+AOM_lag-MW_lag+t_AOM+i*t_total for i in range(0,N_laser,1)]
    MW_dur = list(np.zeros(N_laser)) + [t_MW for i in range(0,N_laser)]
    
    # samp_clk_starts = [0+(t_AOM)+t_MW+i*t_total for i in range(0,N_laser)]
    # samp_clk_starts += [N_laser*t_total+(t_AOM)+t_MW+i*t_total for i in range(0,N_laser)]
    # samp_clk_dur = [RO_width for i in range(0,2*N_laser,1)]
    # samp_clk_dur = [a[i][j] for i in range(0, len(a)) for j in range(0, len(a[0]))]
    
    a = [[0+(t_AOM+AOM_lag)+t_MW+RO_width+i*t_total, 0+(t_AOM+AOM_lag)+t_MW+RO_width+pulse_width+sig_ref_buffer+i*t_total] for i in range(0,N_laser,1)]
    samp_clk_starts = [a[i][j] for i in range(0,len(a)) for j in range(0,len(a[0]))]
    a = [[N_laser*t_total+(t_AOM+AOM_lag)+t_MW+RO_width+i*t_total, N_laser*t_total+(t_AOM+AOM_lag)+t_MW+RO_width+pulse_width+sig_ref_buffer+i*t_total] for i in range(0,N_laser,1)]
    samp_clk_starts += [a[i][j] for i in range(0,len(a)) for j in range(0,len(a[0]))]
    samp_clk_dur = [pulse_width for i in range(0,2*2*N_laser)]
    
    laser_ch = PBchannel(laser, laser_starts, laser_dur)
    MW_ch = PBchannel(MW, MW_starts, MW_dur)
    start_trig_ch = PBchannel(start_trig, [0], [200*ns])
    samp_clk_ch = PBchannel(samp_clk, samp_clk_starts, samp_clk_dur)
    # conv_clk_ch = PBchannel(conv_clk, conv_clk_starts, conv_clk_dur)
    
    allPBchannels = [laser_ch, MW_ch, start_trig_ch, samp_clk_ch]
    return allPBchannels

def make_pulsed_esr_seq(t_AOM, t_ro_delay, AOM_lag, MW_lag, t_pi):
    """
    Make pulse sequence for Pulsed ODMR

    Parameters
    ----------
    t_MW : float
        Width of MW pulse.
    t_AOM : float
        AOM 'ON' time or AOM width.
    t_ro_delay : float
        Length of (AOM) readout pulse.
    AOM_lag: float
        
    MW_lag: float
        
    Returns
    -------
    allPBchannels : list
        A list of PBchannel object types.

    """
    t_flip = t_pi
    apd_pulse = [t_flip+AOM_lag+t_ro_delay, 2*t_flip+t_AOM+AOM_lag+t_ro_delay]
    laser_channel = PBchannel(laser, [t_flip, 2*t_flip+t_AOM], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    MWchannel = PBchannel(MW, [AOM_lag-MW_lag], [t_flip])
    # conv_clk_channel = PBchannel(conv_clk,[apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    # samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width for i in range(0,2)])
    
    samp_clk_channel = PBchannel(samp_clk,[apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,2)])
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    # print(allPBchannels)
    return allPBchannels

#------------------------------------------------------------------------------
def make_echo_seq_MW(delay_2nd_half, delay_1st_half, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi):
    """ Spin-echo seq. 
    Ref = FL w/o MW"""
    t_delay = delay_1st_half + delay_2nd_half + 2*t_pi# + 5*us     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
    apd_pulse = [t_delay+AOM_lag+ro_delay, 2*t_delay+t_AOM+AOM_lag+ro_delay]

    laser_channel = PBchannel(laser, [t_delay, 2*t_delay+t_AOM], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag, 0*us+AOM_lag+t_pi/2+delay_1st_half-MW_lag, 0*us+AOM_lag+(t_delay-0*us)-MW_lag-t_pi/2], [t_pi/2, t_pi, t_pi/2])
    # The third MW pulse 
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width, pulse_width])
    conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel, conv_clk_channel]
    return allPBchannels

#------------------------------------------------------------------------------
def make_echo_seq_FL(delay_2nd_half, delay_1st_half, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi):
    """ Spin-echo seq. 
    Ref = Saturated FL"""
    seq_start_delay = 1*us      # = t_AOM_init_pulse
    # daq_ref_pulse = seq_start_delay + 5*us#dur_AOM + AOM_lag - (1000)*ns      For ref pulse ahead of signal pulse
    t_delay = delay_1st_half + delay_2nd_half + 2*t_pi# + 5*us     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
    t_AOM_readout_pulse = seq_start_delay + t_AOM + t_delay
    daq_sig_pulse = t_AOM_readout_pulse + AOM_lag + ro_delay
    MW_pulse_start = seq_start_delay + t_AOM + AOM_lag# + 5*us
    daq_ref_pulse = t_AOM_readout_pulse + 7*us
    apd_pulse = [daq_sig_pulse, daq_ref_pulse]
    
    laser_channel = PBchannel(laser, [seq_start_delay, t_AOM_readout_pulse], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [seq_start_delay], [pulse_width])
    MWchannel = PBchannel(MW, [MW_pulse_start-MW_lag, MW_pulse_start+t_pi/2+delay_1st_half-MW_lag, MW_pulse_start+(t_delay-0*us)-MW_lag-t_pi/2], [t_pi/2, t_pi, t_pi/2])
    
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width for i in range(0,2)])
    
    conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel, conv_clk_channel]
    return allPBchannels

#------------------------------------------------------------------------------
def make_ramsey_seq_MW(t_precession, t_AOM, ro_delay, AOM_lag, MW_lag, t_piby2):
    t_precession = t_precession + 2*t_piby2     # increase the actual laser off time to accomodate the non-zero width of the pi/2 pulses so that the actual precession time is as defined by the param variable of mainControl
    apd_pulse = [t_precession+AOM_lag+ro_delay, 2*t_precession+t_AOM+AOM_lag+ro_delay]
    laser_channel = PBchannel(laser, [t_precession, 2*t_precession+t_AOM], [t_AOM, t_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    MWchannel = PBchannel(MW, [AOM_lag-MW_lag, t_precession+AOM_lag-MW_lag-t_piby2], [t_piby2, t_piby2])
    
    # samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width,apd_pulse[1]-pulse_width], [pulse_width for i in range(0,2)])
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,2)])
    
    allPBchannels = [MWchannel, laser_channel,samp_clk_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    return allPBchannels

def make_ramsey_seq_FL(t_delay, t_AOM, ro_delay, AOM_lag, MW_lag, t_piby2):
    """ Spin-echo seq. 
    Ref = Saturated FL"""
    seq_start_delay = 1*us      # = t_AOM_init_pulse
    # daq_ref_pulse = seq_start_delay + 5*us#dur_AOM + AOM_lag - (1000)*ns      For ref pulse ahead of signal pulse
    t_delay = t_delay + 2*t_piby2     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
    t_AOM_readout_pulse = seq_start_delay + t_AOM + t_delay
    daq_sig_pulse = t_AOM_readout_pulse + AOM_lag + ro_delay
    MW_pulse_start = seq_start_delay + t_AOM + AOM_lag
    daq_ref_pulse = t_AOM_readout_pulse + 7*us
    apd_pulse = [daq_sig_pulse, daq_ref_pulse]
    
    laser_channel = PBchannel(laser, [seq_start_delay, t_AOM_readout_pulse], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [seq_start_delay], [pulse_width])
    MWchannel = PBchannel(MW, [MW_pulse_start-MW_lag,  MW_pulse_start+t_delay-MW_lag-t_piby2], [t_piby2 for i in range(0,2)])
    
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width for i in range(0,2)])
    
    conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel, conv_clk_channel]
    return allPBchannels
#------------------------------------------------------------------------------
# def make_rabi_seq(t_MW, t_AOM, t_ro_delay):
#     """
#     Make pulse sequence for Rabi oscillations

#     Parameters
#     ----------
#     t_MW : float
#         Width of MW pulse.
#     t_AOM : float
#         AOM 'ON' time or AOM width.
#     t_ro_delay : float
#         Length of (AOM) readout pulse.

#     Returns
#     -------
#     allPBchannels : list
#         A list of PBchannel object types.

#     """
#     AOM_lag = (860+300)*ns
#     MW_lag = 250*ns
#     t_drive = 1*us
#     starttrig_width = 300*ns
#     daq_readout_width = 300*ns

#     AOMchannel = PBchannel(AOM, [t_drive, 2*t_drive+t_AOM], [t_AOM,t_AOM])
#     STARTtrigchannel = PBchannel(STARTtrig, [AOM_lag], [starttrig_width])
#     DAQchannel = PBchannel(DAQ, [t_drive+AOM_lag+t_ro_delay, 2*t_drive+t_AOM+AOM_lag+t_ro_delay], [daq_readout_width, daq_readout_width])
    
#     # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
#     if t_MW <= 5*clk_cyc and t_MW > 0:
#         # (PBchannel) for MW pulse
#         MWchannel = PBchannel(MW, [AOM_lag-MW_lag+500], [5*clk_cyc])
#         # SHORT PULSE duration
#         shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
#         shortPulseChannel = PBchannel(shortpulseFLAG, [AOM_lag-MW_lag+500], [5*clk_cyc])     # (PBchannel) for SHORT MW pulse
#         allPBchannels = [shortPulseChannel, MWchannel]  # Short pulse feature
#     else:
#         MWchannel = PBchannel(MW, [AOM_lag-MW_lag+500], [t_MW])
#         allPBchannels = [MWchannel]
    
#     allPBchannels.extend([AOMchannel, DAQchannel, STARTtrigchannel])
#     # print(allPBchannels)
#     return allPBchannels
#------------------------------------------------------------------------------
# def make_rabi_seq(t_MW, t_AOM, t_readoutDelay):       # ORIGINAL
#     start_delay = t_readoutDelay
#     t_startTrig = clk_cyc*round(100*ns/clk_cyc)    # startTrigger time
#     t_readout_width = clk_cyc*round(100*ns/clk_cyc)      #
#     # AOM_delay = (800+300+500)*ns ... Put in t_readoutDelay in rabi_config.py
#     MWtoAOM_delay = clk_cyc*round(200*ns/clk_cyc)
#     firstHalfDuration = start_delay + t_MW + MWtoAOM_delay + t_AOM

#     AOMchannel = PBchannel(AOM, [firstHalfDuration-t_AOM, 2*firstHalfDuration-t_AOM], [t_AOM, t_AOM])
#     DAQchannel = PBchannel(DAQ, [firstHalfDuration-t_AOM+t_readoutDelay,2*firstHalfDuration-t_AOM+t_readoutDelay], [t_readout_width, t_readout_width])
#     STARTtrigchannel = PBchannel(STARTtrig, [0], [t_startTrig])
    
#     # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
#     if t_MW <= 5*clk_cyc and t_MW > 0:
#         # (PBchannel) for MW pulse
#         MWchannel = PBchannel(MW, [start_delay], [5*clk_cyc])
#         # SHORT PULSE duration
#         shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
#         shortPulseChannel = PBchannel(shortpulseFLAG, [start_delay], [5*clk_cyc])     # (PBchannel) for SHORT MW pulse
#         allPBchannels = [shortPulseChannel, MWchannel]  # Short pulse feature
#     else:
#         MWchannel = PBchannel(MW, [start_delay], [t_MW])
#         allPBchannels = [MWchannel]
    
#     allPBchannels.extend([AOMchannel, DAQchannel, STARTtrigchannel])
#     # print(allPBchannels)
#     return allPBchannels
#------------------------------------------------------------------------------
# allPBchannels: (List) of PBallPBchannels
    #        Each (PBchannel) consists of 1.channelNumber 2.starttime 3.Pulse duration for a component (MW, AOM, DAQ)
    #   E.g.:
    #    PBchannel(MW,[100],[500]) means channelNumber 32 will be turned ON at 100ns (starttime) for a duration of 500ns.
    #    PBchannel(AOM,[10,200],[2000,350]) means channelNumber 16 will be turned ON at 10ns for 2us and again turned ON at 200ns for 350ns.

def make_rodelayscan_seq(rodelay, dur_AOM, AOM_lag, delay):
    seq_start_delay = 1*us      # = t_AOM_init_pulse
    
    daq_ref_pulse = seq_start_delay + dur_AOM + AOM_lag # - (600)*ns
    t_AOM_readout_pulse = seq_start_delay + dur_AOM + delay
    daq_sig_pulse = t_AOM_readout_pulse + AOM_lag + rodelay
    # tot_time = 
    
    laser_channel = PBchannel(laser, [seq_start_delay, t_AOM_readout_pulse], [dur_AOM, dur_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    samp_clk_channel = PBchannel(samp_clk, [daq_ref_pulse, daq_sig_pulse], [pulse_width, pulse_width])
    
    allPBchannels = [laser_channel, start_trig_channel, samp_clk_channel]
    return allPBchannels

def make_t1_seq_camera(relax_time, dur_AOM, AOM_lag, rodelay, MW_lag, t_pi):
    # laser_start_times = [relax_time, 0]; laser_pulse_durations = [dur_AOM, relax_time + dur_AOM]
    laser_start_times = [relax_time, relax_time]; laser_pulse_durations = [dur_AOM, dur_AOM]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel]
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]

    MW_start_times = [AOM_lag-MW_lag]; MW_durations = [t_pi]
    MW_channel = PBchannel(MW, [MW_start_times[0]], [MW_durations[0]])
    allPBchannels_sig.extend([MW_channel])

    # MW_start_times = [dur_AOM+relax_time+AOM_lag-MW_lag]; MW_durations = [t_pi]
    # MW_channel = PBchannel(MW, [MW_start_times[0]], [MW_durations[0]])
    # allPBchannels_ref.append([MW_channel])

    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

def make_t1_ms0_seq_FL(relax_time, dur_AOM, AOM_lag, rodelay, MW_lag, t_pi):
    """
    Generate the pulse sequence to study T1 relaxation from m_s=0 state. Sequence implmented as signal and reference with MW and Saturated FL. 

    Parameters
    ----------
    delay : float
        The 'wait' interval between the initialization pulse and the readout pulse in the sequence.
    dur_AOM : float
        The ON time of the AOM initialization/readout pulse.
    AOM_lag : float
        The time lag involved in the AOM pulse from PB reaching the AOM and the actual start of the light pulse
    ro_delay : float
        The time to wait before the reading the FL signal after the AOM readout pulse occurs.

    Returns
    -------
    allPBchannels : list of PBchannels
        Contains information of all the PB channels involved in the sequence with their start times and pulse durations.

    """
    seq_start_delay = 1*us      # = t_AOM_init_pulse
    
    daq_ref_pulse = seq_start_delay + 8*us#dur_AOM + AOM_lag - (1000)*ns
    t_AOM_readout_pulse = seq_start_delay + dur_AOM + relax_time
    daq_sig_pulse = t_AOM_readout_pulse + AOM_lag + rodelay
    apd_pulse = [daq_ref_pulse, daq_sig_pulse]
    # tot_time = 
    
    laser_channel = PBchannel(laser, [seq_start_delay, t_AOM_readout_pulse], [dur_AOM])
    start_trig_channel = PBchannel(start_trig, [seq_start_delay], [pulse_width])
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width, pulse_width])
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [laser_channel, samp_clk_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    return allPBchannels

def make_t1_ms0_seq_MW(relax_time, dur_AOM, AOM_lag, rodelay, MW_lag, t_pi):
    """
    Generate the pulse sequence to study T1 relaxation from m_s=+-1 state. Sequence implemented as signal and reference with/without MW 

    Parameters
    ----------
    relax_time : float
        The 'wait' interval between the initialization pulse and the readout pulse in the sequence.
    dur_AOM : float
        The ON time of the AOM initialization/readout pulse.
    AOM_lag : float
        The time lag involved in the AOM pulse from PB reaching the AOM and the actual start of the light pulse
    rodelay : float
        The time to wait before the reading the FL signal after the AOM readout pulse occurs.

    Returns
    -------
    allPBchannels : list of PBchannels
        Contains information of all the PB channels involved in the sequence with their start times and pulse durations.

    """
    relax_time = relax_time + t_pi     # increase the actual laser off time to accomodate the non-zero width of the pi pulses so that the actual precession time is as defined by the param variable of mainControl
    apd_pulse = [relax_time+AOM_lag+rodelay, 2*relax_time+dur_AOM+AOM_lag+rodelay]
    
    laser_channel = PBchannel(laser, [relax_time, 2*relax_time+dur_AOM], [dur_AOM,dur_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    # MWchannel = PBchannel(MW, [AOM_lag-MW_lag, relax_time+dur_AOM+AOM_lag-MW_lag], [t_pi, t_pi])
    MWchannel = PBchannel(MW, [AOM_lag-MW_lag], [t_pi])
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width, pulse_width])
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0, 4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    
    return allPBchannels

def make_t1_ms1_seq_FL(relax_time, dur_AOM, AOM_lag, rodelay, t_pi, MW_lag):
    """
    Generate the pulse sequence to study T1 relaxation from m_s=+-1 state. Sequence implmented as signal and reference with MW and Saturated FL. 

    Parameters
    ----------
    delay : float
        The 'wait' interval between the initialization pulse and the readout pulse in the sequence.
    dur_AOM : float
        The ON time of the AOM initialization/readout pulse.
    AOM_lag : float
        The time lag involved in the AOM pulse from PB reaching the AOM and the actual start of the light pulse
    ro_delay : float
        The time to wait before the reading the FL signal after the AOM readout pulse occurs.

    Returns
    -------
    allPBchannels : list of PBchannels
        Contains information of all the PB channels involved in the sequence with their start times and pulse durations.

    """
    seq_start_delay = 1*us      # = t_AOM_init_pulse
    
    daq_ref_pulse = seq_start_delay + 5*us#dur_AOM + AOM_lag - (1000)*ns
    t_AOM_readout_pulse = seq_start_delay + dur_AOM + relax_time
    daq_sig_pulse = t_AOM_readout_pulse + AOM_lag + rodelay
    apd_pulse = [daq_ref_pulse, daq_sig_pulse]
    # tot_time = 
    
    laser_channel = PBchannel(laser, [seq_start_delay, t_AOM_readout_pulse], [dur_AOM, dur_AOM])
    start_trig_channel = PBchannel(start_trig, [seq_start_delay], [pulse_width])
    MWchannel = PBchannel(MW, [seq_start_delay+dur_AOM+AOM_lag-MW_lag], [t_pi])
    
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width, pulse_width])
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    return allPBchannels



def make_t2_seq_FL(t_delay, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi):
    """ Spin-echo seq. 
    Ref = Saturated FL"""
    seq_start_delay = 1*us      # = t_AOM_init_pulse
    daq_ref_pulse = seq_start_delay + 5*us#dur_AOM + AOM_lag - (1000)*ns      For ref pulse ahead of signal pulse
    t_delay = t_delay + 2*t_pi     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
    t_AOM_readout_pulse = seq_start_delay + t_AOM + t_delay
    daq_sig_pulse = t_AOM_readout_pulse + AOM_lag + ro_delay
    MW_pulse_start = seq_start_delay + t_AOM + AOM_lag
    # daq_ref_pulse = t_AOM_readout_pulse + 7*us
    apd_pulse = [daq_sig_pulse, daq_ref_pulse]
    
    laser_channel = PBchannel(laser, [seq_start_delay, t_AOM_readout_pulse], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [seq_start_delay], [pulse_width])
    MWchannel = PBchannel(MW, [MW_pulse_start-MW_lag, MW_pulse_start+t_delay/2-MW_lag-t_pi/2, MW_pulse_start+t_delay-MW_lag-t_pi/2], [t_pi/2, t_pi, t_pi/2])
    
    # samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width for i in range(0,2)])    
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,2)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]
    # allPBchannels.extend([conv_clk_channel])
    return allPBchannels

def make_t2_seq_MW(t_delay, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi):
    """ T2 seq. 
    Ref = FL w/o MW"""
    t_delay = t_delay + 2*t_pi     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
    apd_pulse = [t_delay+AOM_lag+ro_delay, 2*t_delay+t_AOM+AOM_lag+ro_delay]
    
    laser_channel = PBchannel(laser, [t_delay, 2*t_delay+t_AOM], [t_AOM,t_AOM])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    MWchannel = PBchannel(MW, [AOM_lag-MW_lag, AOM_lag+t_delay/2-MW_lag-t_pi/2, t_delay+AOM_lag-MW_lag-t_pi/2], [t_pi/2, t_pi, t_pi/2])
    # samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width, pulse_width])
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,2)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]
     # allPBchannels.extend([conv_clk_channel])
    return allPBchannels
    
# def make_t1_ms0_seq(delay, t_AOM, AOM_lag):
#     starttrig_width = 300*ns
#     daq_readout_width = 300*ns
#     seq_start_delay = 1*us
#     daq_ref_pulse = seq_start_delay + t_AOM + AOM_lag - 500*ns - daq_readout_width
#     daq_sig_pulse = seq_start_delay + t_AOM + AOM_lag + delay
#     # tot_time = 
    
#     AOMchannel = PBchannel(AOM, [seq_start_delay], [t_AOM])
#     STARTtrigchannel = PBchannel(STARTtrig, [seq_start_delay], [starttrig_width])
#     DAQchannel = PBchannel(DAQ, [daq_ref_pulse,daq_sig_pulse], [daq_readout_width,daq_readout_width])
    
#     allPBchannels = [AOMchannel, STARTtrigchannel, DAQchannel]
#     return allPBchannels

def make_t1_ms1_seq0(t_delay,t_AOM,t_readoutDelay,t_pi):
    t_startTrig = clk_cyc*round(300*ns/clk_cyc)
    t_readout = clk_cyc*round(300*ns/clk_cyc)
    MWtoAOM_delay =clk_cyc*round(1*us/clk_cyc)
    AOMstartTime1 = t_delay
    firstHalfDuration=AOMstartTime1+t_AOM
    AOMstartTime2 =  firstHalfDuration+AOMstartTime1 
    AOMchannel = PBchannel(laser,[AOMstartTime1,AOMstartTime2],[t_AOM,t_AOM])
    MWchannel = PBchannel(MW,[firstHalfDuration +t_readoutDelay + clk_cyc*round(1*us/clk_cyc)],[t_pi])
    DAQchannel = PBchannel(samp_clk,[AOMstartTime1+t_readoutDelay, AOMstartTime2+t_readoutDelay],[t_readout,t_readout])
    STARTtrigchannel = PBchannel(start_trig,[0],[t_startTrig])
    allPBchannels = [AOMchannel, DAQchannel, MWchannel, STARTtrigchannel]
    return allPBchannels
# -----------------------------------------------------------------------------

def check_simult_sampling():
    interval = 1*us
    samp_clk_edge = interval + 0.5*interval
    conv_clk_edge = samp_clk_edge+0.1*interval
    laser_channel = PBchannel(laser, [2*interval],[1*interval])
    samp_clk_channel = PBchannel(samp_clk, [samp_clk_edge], [pulse_width])
    conv_clk_channel = PBchannel(conv_clk, [conv_clk_edge, conv_clk_edge+2*us], [pulse_width for i in range(0,2)])
    allPBchannels = [laser_channel, samp_clk_channel, conv_clk_channel]
    return allPBchannels

def make_drift_analysis_sequence(duration):
    
    laser_channel = PBchannel(laser, [0], [])
    MW_channel = PBchannel(MW, [0], [])
    start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
    samp_clk_channel = PBchannel(samp_clk, [duration-pulse_width], [pulse_width])
    conv_clk_channel = PBchannel(conv_clk, [duration, duration+2*us-pulse_width], [pulse_width for i in range(0,2)])
    
    allPBchannels = [laser_channel, MW_channel, samp_clk_channel, conv_clk_channel, start_trig_channel]
    return allPBchannels

def makeT2Seq0(t_delay,t_AOM,t_readoutDelay,t_pi,IQpadding, numberOfPiPulses):
    t_piby2=t_pi/2
    t_startTrig = clk_cyc*round(300*ns/clk_cyc)
    t_readout = clk_cyc*round(300*ns/clk_cyc)
    MWtoAOM_delay =clk_cyc*round(1*us/clk_cyc)
    start_delay = (clk_cyc*round(1*us/clk_cyc) + t_readoutDelay)
    #Make pulses for signal half of the sequence:
    [MWstartTimes1,MWdurations,IstartTimes1,Idurations,QstartTimes1,Qdurations]= makeCPMGpulses(start_delay,numberOfPiPulses,t_delay,t_pi, t_piby2,IQpadding)
    CPMGduration = MWstartTimes1[-1]+t_piby2-start_delay
    AOMstartTime1 = start_delay+CPMGduration +MWtoAOM_delay
    DAQstartTime1 = AOMstartTime1+t_readoutDelay
    firstHalfDuration = AOMstartTime1+t_AOM
    #Make pulses for background half of the sequence:
    MWstartTimes2 =[x+firstHalfDuration for x in MWstartTimes1]
    #IstartTimes2 = [x+firstHalfDuration for x in IstartTimes1]
    QstartTimes2nd = QstartTimes1[:-1]
    QstartTimes2 = [x+firstHalfDuration for x in QstartTimes2nd]
    AOMstartTime2 = firstHalfDuration + AOMstartTime1
    DAQstartTime2 = firstHalfDuration + DAQstartTime1
    #Make full MW, I and Q pulse lists
    MWstartTimes=MWstartTimes1 +MWstartTimes2
    # Make allPBchannels:
    laser_channel 		 = PBchannel(laser,[AOMstartTime1,AOMstartTime2],[t_AOM,t_AOM])
    samp_clk_channel 		 = PBchannel(samp_clk,[DAQstartTime1,DAQstartTime2],[t_readout,t_readout])
    MWchannel  		 = PBchannel(MW,MWstartTimes1 +MWstartTimes2,MWdurations+MWdurations)
    Ichannel   		 = PBchannel(I,IstartTimes1,Idurations)
    Qchannel   		 = PBchannel(Q,QstartTimes1+QstartTimes2,Qdurations+Qdurations[:-1])
    start_trig_channel = PBchannel(start_trig,[0],[t_startTrig])
    allPBchannels=[laser_channel, samp_clk_channel, MWchannel, Ichannel, Qchannel, start_trig_channel]
    return allPBchannels
#
#------------------------------------------------------------------------------
# def makeXY8seq(t_delay,t_AOM,t_readoutDelay,t_pi,IQpadding, numberOfRepeats):
# 	t_piby2=t_pi/2
# 	t_startTrig = clk_cyc*round(300*ns/clk_cyc)
# 	t_readout = clk_cyc*round(300*ns/clk_cyc)
# 	MWtoAOM_delay =clk_cyc*round(1*us/clk_cyc)
# 	start_delay = (clk_cyc*round(1*us/clk_cyc) + t_readoutDelay)
# 	#Make pulses for signal half of the sequence:
# 	[MWstartTimes1,MWdurations,IstartTimes1,Idurations,QstartTimes1,Qdurations]= makeXY8pulses(start_delay,numberOfRepeats,t_delay,t_pi, t_piby2,IQpadding)
# 	XY8duration = MWstartTimes1[-1]+t_pi/2-start_delay
# 	AOMstartTime1 = start_delay+XY8duration +MWtoAOM_delay
# 	DAQstartTime1 = AOMstartTime1+t_readoutDelay
# 	firstHalfDuration = AOMstartTime1+t_AOM
# 	#Make pulses for background half of the sequence:
# 	MWstartTimes2 =[x+firstHalfDuration for x in MWstartTimes1]
# 	QstartTimes2nd = QstartTimes1[:-1]
# 	QstartTimes2 = [x+firstHalfDuration for x in QstartTimes2nd]
# 	AOMstartTime2 = firstHalfDuration + AOMstartTime1
# 	DAQstartTime2 = firstHalfDuration + DAQstartTime1
# 	#Make allPBchannels:
# 	AOMchannel 		 = PBchannel(AOM,[AOMstartTime1,AOMstartTime2],[t_AOM,t_AOM])
# 	DAQchannel 		 = PBchannel(DAQ,[DAQstartTime1,DAQstartTime2],[t_readout,t_readout])
# 	MWchannel  		 = PBchannel(MW,MWstartTimes1 +MWstartTimes2,MWdurations+MWdurations)
# 	Ichannel   		 = PBchannel(I,IstartTimes1,Idurations)
# 	Qchannel   		 = PBchannel(Q,QstartTimes1+QstartTimes2,Qdurations+Qdurations[:-1])
# 	STARTtrigchannel = PBchannel(STARTtrig,[0],[t_startTrig])
# 	allPBchannels=[AOMchannel,DAQchannel, MWchannel, Ichannel, Qchannel, STARTtrigchannel]
# 	return allPBchannels
#
#------------------------------------------------------------------------------
# def makecorrelationSpectSeq(t_delay_betweenXY8seqs,t_delay, t_AOM,t_readoutDelay,t_pi,IQpadding,numberOfRepeats):
# 	t_piby2=t_pi/2
# 	t_startTrig = clk_cyc*round(300*ns/clk_cyc)
# 	t_readout = clk_cyc*round(300*ns/clk_cyc)
# 	MWtoAOM_delay =clk_cyc*round(1*us/clk_cyc)
# 	start_delay = clk_cyc*round(2*us/clk_cyc)
# 	#Make pulses for first XY8 in the first half of the sequence (I only pulses in second XY8 of second half, so we will take the I times here and shift them in time):
# 	[MWstartTimes1a,MWdurations1a,IstartTimes,Idurations,QstartTimes1a,Qdurations1a]= makeXY8pulses(start_delay,numberOfRepeats,t_delay,t_pi, t_piby2,IQpadding)
# 	#Make pulses for second XY8 in the first half of the sequence:
# 	firstXY8duration = MWstartTimes1a[-1]+t_piby2
# 	MWstartTimes1b = [x+firstXY8duration + t_delay_betweenXY8seqs for x in MWstartTimes1a]
# 	QstartTimes1b = [x +firstXY8duration + t_delay_betweenXY8seqs for x in QstartTimes1a]
# 	Qdurations1b = Qdurations1a
# 	#Make AOM pulse and DAQ pulse for signal half
# 	XY8duration = MWstartTimes1b[-1]+t_pi/2-start_delay
# 	AOMstartTime1 = start_delay+XY8duration +MWtoAOM_delay
# 	DAQstartTime1 = AOMstartTime1+t_readoutDelay
# 	firstHalfDuration = AOMstartTime1+t_AOM
#
# 	#Make pulses for first XY8 in the second half of the sequence (no I's on this half):
# 	MWstartTimes2 = [x+firstHalfDuration for x in MWstartTimes1a+MWstartTimes1b]
# 	QstartTimes2 = [x+firstHalfDuration for x in QstartTimes1a+QstartTimes1b[:-1]]
# 	AOMstartTime2 = firstHalfDuration + AOMstartTime1
# 	DAQstartTime2 = firstHalfDuration + DAQstartTime1
# 	IstartTimes = [x +firstHalfDuration+firstXY8duration+t_delay_betweenXY8seqs for x in IstartTimes]
#
# 	#concatenate pulse times:
# 	MWstartTimes = MWstartTimes1a+MWstartTimes1b+MWstartTimes2
# 	MWdurations = MWdurations1a*4
# 	QstartTimes = QstartTimes1a+QstartTimes1b+QstartTimes2
# 	Qdurations = Qdurations1a+Qdurations1b+Qdurations1a+Qdurations1b[:-1]
#
# 	#Make allPBchannels:
# 	AOMchannel 		 = PBchannel(AOM,[AOMstartTime1,AOMstartTime2],[t_AOM,t_AOM])
# 	DAQchannel 		 = PBchannel(DAQ,[DAQstartTime1,DAQstartTime2],[t_readout,t_readout])
# 	MWchannel  		 = PBchannel(MW,MWstartTimes,MWdurations)
# 	Ichannel   		 = PBchannel(I,IstartTimes,Idurations)
# 	Qchannel   		 = PBchannel(Q,QstartTimes,Qdurations)
# 	STARTtrigchannel = PBchannel(STARTtrig,[0],[t_startTrig])
# 	allPBchannels=[AOMchannel,DAQchannel, MWchannel, Ichannel, Qchannel, STARTtrigchannel]
# 	return allPBchannels
#
#------------------------------------------------------------------------------
def makeCPMGpulses(start_delay,numberOfPiPulses,t_delay,t_pi, t_piby2,IQpadding):
    t_piby4 = t_piby2/2;
    if numberOfPiPulses == 1:
        MWstartTimes = [start_delay, start_delay +t_delay-t_piby4, start_delay +2*t_delay]
        MWdurations =  [t_piby2, t_pi, t_piby2]
    else:
         half_t_delay = t_delay/2
         #Start off the sequence by adding the initial pi/2 and first pi pulse
         MWstartTimes =[start_delay, start_delay +half_t_delay-t_piby4]
         MWdurations =[t_piby2, t_pi]
         #Add remaining pi pulses:
         for i in range(1,numberOfPiPulses):
             currentEdgeTime = MWstartTimes[-1]+t_delay
             MWstartTimes.append(currentEdgeTime)
             MWdurations.append(t_pi)
         #Append the final pi/2 pulse:
         MWstartTimes.append(MWstartTimes[-1]+half_t_delay+t_piby4)
         MWdurations.append(t_piby2)
    #Make the I and Q channel pulses:
    #Q is ON during pi(y) pulses and the final pi/2(-x), but not the first pi/2(x) pulse
    QstartTimes = [x-IQpadding for x in MWstartTimes[1:]]
    Qdurations=[x +2*IQpadding for x in MWdurations[1:]]
 	#I is only on during the final pi/2(-x) pulse:
    IstartTimes =[x-IQpadding for x in [MWstartTimes[-1]]]
    Idurations =[x +2*IQpadding for x in [MWdurations[-1]]]
    return [MWstartTimes,MWdurations,IstartTimes,Idurations,QstartTimes,Qdurations]
#
#------------------------------------------------------------------------------
# def makeXY8pulses(start_delay,numberOfRepeats,t_delay,t_pi, t_piby2,IQpadding):
# 	t_piby4=t_piby2/2
# 	#Start off the sequence by adding the initial pi/2:
# 	half_t_delay = t_delay/2
# 	MWstartTimes =[start_delay]
# 	MWdurations =[t_piby2]
# 	QstartTimes=[]
# 	Qdurations =[]
# 	#Add remaining pi pulses:
# 	firstPiPulseDone = False
# 	for i in range(0,numberOfRepeats):
# 		#Make next 8 pi pulses:
# 		next8piPulseStartTimes =[]
# 		next8piPulseDurations=[]
# 		#Add the first pulse in the set of 8
# 		currentEdgeTime=0
# 		if not firstPiPulseDone:
# 			currentEdgeTime = MWstartTimes[-1]+half_t_delay-t_piby4
# 			firstPiPulseDone=True
# 		else:
# 			currentEdgeTime = MWstartTimes[-1]+t_delay
# 		next8piPulseStartTimes.append(currentEdgeTime)
# 		next8piPulseDurations.append(t_pi)
# 		for j in range (1,8):
# 			newEdgeTime = next8piPulseStartTimes[-1]+t_delay
# 			next8piPulseStartTimes.append(newEdgeTime)
# 			next8piPulseDurations.append(t_pi)
# 		# Make next 8 Q start times (Q is only on for pulses 1,3,4,6 of the xy8 pi pulses, for a 0-indexed sequence):
# 		next8QstartTimes = list(next8piPulseStartTimes[i] for i in [1,3,4,6])
# 		next8Qdurations = list(next8piPulseDurations[i] for i in [1,3,4,6])
# 		# Append next 8 pi pulses and Q pulses to start times lists:
# 		MWstartTimes.extend(next8piPulseStartTimes)
# 		MWdurations.extend(next8piPulseDurations)
# 		QstartTimes.extend(next8QstartTimes)
# 		Qdurations.extend(next8Qdurations)
#
# 	#Append the final pi/2 pulse:
# 	MWstartTimes.append(MWstartTimes[-1]+half_t_delay+t_piby4)
# 	MWdurations.append(t_piby2)
# 	#Append the final pi/2 pulse to Q channel since (in the signal bin) Q is on for this pulse as it is a -x pulse.
# 	QstartTimes.append(MWstartTimes[-1])
# 	Qdurations.append(MWdurations[-1])
# 	#Pad the Q channel pulses:
# 	QstartTimes = [x-IQpadding for x in QstartTimes]
# 	Qdurations=[x +2*IQpadding for x in Qdurations]
# 	#Make I channel pulses. I is only on during the final pi/2(-x) pulse:
# 	IstartTimes =[x-IQpadding for x in [MWstartTimes[-1]]]
# 	Idurations =[x +2*IQpadding for x in [MWdurations[-1]]]
# 	return [MWstartTimes,MWdurations,IstartTimes,Idurations,QstartTimes,Qdurations]

# %%
