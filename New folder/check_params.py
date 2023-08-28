# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 23:10:27 2022

@author: PC
"""

from importlib import import_module
from PBcontrol import PBclk
import sys
import numpy as np

def check_params(config_file):
    expCfg = import_module(config_file)      # import the configuration file as module
    clk_cyc = 1e3/PBclk
    
    # Check that N_scanPts, Nsamples, Navg are all integers and Nsamples>=1, Navg>= 1, N_scanPts>=2
    if (not isinstance(expCfg.Nsamples, int)):
        expCfg.Nsamples = int(expCfg.Nsamples)
        print("\x10 \x1b[38;2;240;240;50mWarning\x1b[0m: Nsamples is a 'float'.")
    if  (expCfg.Nsamples<1):
        print('Error: Nsamples must be an integer >= 1.')
        sys.exit()
    if (not isinstance(expCfg.Navg, int)) or (expCfg.Navg<1):
        print('Error: Navg must be an integer >= 1.')
        sys.exit()
    if (not isinstance(expCfg.N_scanPts, int)) or (expCfg.N_scanPts<2):
        print('Error: N_scanPts must be an integer >= 2.')
        sys.exit()
        
    step_size = expCfg.scannedParam[1] - expCfg.scannedParam[0]
    if expCfg.sequence not in ['aom_timing', 'T1ms0', 'rodelay']:      # Check uW power if the seqeunce is not 'aom_timing', 'T1ms0', etc
        if expCfg.uW_power >= 9:      # Quit program if the uW power is greater/equals 9 dBm
            print("Input Microwave Power to the RF Amplifier is + "+expCfg.uW_power+" dBm. Reduce it to 8 dBm or less.")
            sys.exit()
        else:
            print("\t\x10 MW Power = "+str(expCfg.uW_power)+" dBm")
            
    if expCfg.sequence in ['aom_timing', 'T1ms0', 'rodelay']:
        if step_size <= 5*clk_cyc:
            print("\x1b[1;33;41mERR: Step size = "+str(step_size)+"... Exiting..."+'\x1b[0m')
            sys.exit()
        # if step_size % clk_cyc != 0:
        #     step_size = 
        for i in range(0, expCfg.N_scanPts):
            expCfg.scannedParam[i] = round(expCfg.scannedParam[i])
            if expCfg.scannedParam[i] % clk_cyc != 0:
                print('Rounding scannedParam to a multiple of CLK_CYC (2ns)...')
                expCfg.scannedParam[i] += 1     # actually this should be expCfg.scannedParam[i] += round(expCfg.scanned[i] % clk_cyc)
        
    if expCfg.sequence == 'uW_timing':
        if step_size <= 5*clk_cyc:
            print("\x1b[1;33;41mERR: Step size = "+str(step_size)+"... Exiting..."+'\x1b[0m')
            sys.exit()
            
        for i in range(0, expCfg.N_scanPts):
            expCfg.scannedParam[i] = round(expCfg.scannedParam[i])
            if expCfg.scannedParam[i] % clk_cyc != 0:
                print('Rounding scannedParam to a multiple of CLK_CYC (2ns)...')
                expCfg.scannedParam[i] += 1     # actually this should be expCfg.scannedParam[i] += round(expCfg.scanned[i] % clk_cyc)
    #------------------------------------------------------------
    
    #Pulse-sequence parameter checks:
    
    #Check that "IQpadding" is a multiple of clk_cyc and >5*clk_cyc:
    if expCfg.sequence in ['T2seq','XY8seq','correlSpecSeq']:
        if (expCfg.IQpadding < (5*clk_cyc)) or (expCfg.IQpadding%clk_cyc):
            print('Error: IQpadding is set to', expCfg.IQpadding,'which is either <',5*clk_cyc,'or not a multiple of',clk_cyc,'. Please edit IQpadding to ensure that it is >',5*clk_cyc,'ns and a multiple of',clk_cyc,'.')
            
    # Check t_duration in esr_seq is a multiple of (2*clk_cyc): (Why?)
    if expCfg.sequence in ['esr_seq', 'pesr_seq', 'ram_seq', 'T2_seq', 'modesr']:
        if expCfg.t_duration%(2*clk_cyc):
            print('\x10 Warning: t_duration set to ', expCfg.t_duration,'ns, which is not an integer multiple of ',(2*clk_cyc),'ns. Rounding t_duration to nearest multiple of ',(2*clk_cyc),'ns...')
            expCfg.t_duration = (2*clk_cyc)*round(expCfg.t_duration/(2*clk_cyc))
            print('t_duration now set to ', expCfg.t_duration,'ns')
    
    #Check that t_readoutDelay and t_AOM are multiples of clk_cyc:
    if expCfg.sequence in ['rabi_seq', 'pser_seq', 'T2_seq', 'XY8seq', 'correlSpecSeq', 'T1seq', 'ram_seq']:
        if (expCfg.readout_delay%clk_cyc) or (expCfg.readout_delay < (5*clk_cyc)):
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
    if expCfg.sequence == 'T1ms1':
        if not expCfg.t_pi % clk_cyc: # or expCfg.t_pi<clk_cyc:
            print('Error: requested pi pulse length '+str(expCfg.t_pi)+'ns \u2260 n('+str(clk_cyc)+') ns.')
            sys.exit()
    
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
            
    if expCfg.sequence in ['rabi_seq', 'T1seq', 'T1ms0', 'T1ms1', 'rodelay', 'ram_seq', 'T2_seq'] or (expCfg.sequence == 'T2seq' and expCfg.numberOfPiPulses == 1):   # why to check for numberOfPiPulses==1??
        if step_size<clk_cyc:
            print('Error: requested time step =',step_size,'ns, < ',clk_cyc,'ns. Please change N_scanPts, or ',expCfg.scanStartName,' and ',expCfg.scanEndName,' to increase time step size.')
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
