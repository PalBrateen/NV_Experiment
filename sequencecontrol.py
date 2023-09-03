# sequencecontrol.py
#%% Imports
from connectionConfig import laser, samp_clk, start_trig, conv_clk, I, Q, MW, PBclk, camera
from spinapi import ns, us, ms
import PBcontrol as PBctrl
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from importlib import import_module
import control_camera_sequences as cam_seq; import control_daq_sequences as daq_seq
# from importlib import import_module

# expCfgFile = 'esr_config'
# expCfg = import_module(expCfgFile)

#%% Operations
clk_cyc = 1e3/PBclk       # Time resolution in ns

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
        var = {'esr_seq':cam_seq.make_esr_seq_camera,   'rabi_seq':cam_seq.make_rabi_seq_camera_MW,
               'pesr_seq': cam_seq.make_pulsed_esr_seq_camera_FL,  'T1ms0': cam_seq.make_t1_seq_camera}
    
    elif instr == 'diode':
        var = {'esr_seq':daq_seq.make_esr_seq,   'modesr':daq_seq.make_mod_esr_seq,
               'rabi_seq':daq_seq.make_rabi_seq,    'spin_echo':daq_seq.make_echo_seq_FL,
               'T2seq':daq_seq.makeT2Seq0,  'pesr_seq':daq_seq.make_pulsed_esr_seq,
               'aom_timing':daq_seq.make_aom_timing_seq,    'rodelay':daq_seq.make_rodelayscan_seq,
               'drift_seq':daq_seq.make_drift_analysis_sequence,    'MW_timing':daq_seq.make_MW_timing_seq,
               'T1ms0':daq_seq.make_t1_ms0_seq_MW,  'T1ms1':daq_seq.make_t1_ms1_seq_FL,
               'T2_seq':daq_seq.make_t2_seq_MW,     'ram_seq':daq_seq.make_ramsey_seq_MW,
               'simult_samp':daq_seq.check_simult_sampling,     'double_mod':daq_seq.make_double_mod_sequence_lcm,
               'diff_mod':daq_seq.make_diff_mod_sequence, 'esr_dig_mod_seq':daq_seq.make_dig_mod_odmr_sequence,
               'rabi_dig_mod_seq':daq_seq.make_dig_mod_rabi_sequence,}
        
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

# %%
