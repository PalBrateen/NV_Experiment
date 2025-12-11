# sequencecontrol.py
from connectionConfig import laser, samp_clk, start_trig, I, Q, MW, PBclk, camera
import sys, math, numpy as np, matplotlib.pyplot as plt, control_camera_sequences as cam_seq,\
    control_daq_sequences as daq_seq
# from PBcontrol import PulseBlaster
import PBcontrol

# expCfgFile = 'esr_config'
# expCfg = import_module(expCfgFile)

# self.clk_cyc = 1e3/PBclk       # Time resolution in ns
# pb = PBctrl.PulseBlaster()

# plot_sequence():
# accepts:
#    'instructionList' - (List of lists) created by sequenceControl.sequenceEventCataloguer()
#    'channelMasks' - ()
# returns: [t_us,channelPulses,yTicks]
class sequencecontrol:

    def __init__(self, parameter_dict) -> None:
        self.parameter_dict: dict = parameter_dict
        self.clk_cyc: float = self.parameter_dict['pb']['clk_cyc']
        # self._pbctrl: PulseBlaster = PulseBlaster(self.parameter_dict)  # Will be injected
        # self.daq_sequences = daq_seq.DAQsequences(self.parameter_dict['pb']['channels'])
    
    # def set_pbcontrol(self, pbctrl_obj):
    #     self._pbctrl = pbctrl_obj
        
    def return_params(self):
        return self.parameter_dict
    
    def check_params(self):
        """
        Docstring for check_params
        
        :param self: Description
        """
        
        # Check that N_scanPts, Nsamples, Nruns are all integers and Nsamples>=1, Nruns>= 1, N_scanPts>=2
        if (not isinstance(self.parameter_dict['seq']['Nsamples'], int)):
            self.parameter_dict['seq']['Nsamples'] = int(self.parameter_dict['seq']['Nsamples'])
            print("⚠ \x1b[38;2;240;240;50mWarning\x1b[0m: Nsamples is a 'float'.")
        
        if  (self.parameter_dict['seq']['Nsamples']<1):
            print('❌ Error: Nsamples must be an integer >= 1.')
            sys.exit()

        if (not isinstance(self.parameter_dict['scan']['Nruns'], int)) or\
            (self.parameter_dict['scan']['Nruns']<1):
            print('❌ Error: Nruns must be an integer >= 1.')
            sys.exit()
        
        Nscanpts = self.parameter_dict['scan']['Nscanpts']
        if (not isinstance(Nscanpts, int)) or (Nscanpts < 2):
            print('❌ Error: N_scanPts must be an integer >= 2.')
            sys.exit()
        
        scan = self.parameter_dict['scan']['values']

        # if len(scan) > 1:
        # step_size = expCfg.scannedParam[1] - expCfg.scannedParam[0]

        if self.parameter_dict['seq']['sequence'] not in ['aom_timing', 'T1ms0_train']:      # Check MW power if the seqeunce is not 'aom_timing', etc
            if self.parameter_dict['mw']['power'] >= 9:      # Quit program if the MW power is greater/equals 9 dBm
                print("❌ Input Microwave Power to the RF Amplifier is + "+self.parameter_dict['mw']['power']+" dBm. Reduce it to 8 dBm or less.")
                sys.exit()
            # else:
            #     # print("\t MW Power = "+str(expCfg.MW_power)+" dBm")
            #     None
        
        # # small step size checks:
        # if self.parameter_dict['seq']['sequence'] in ['aom_timing', 'rodelay']:
        #     if step_size < 1*self.clk_cyc:
        #         print("\x1b[1;33;41mERR: Step size = "+str(step_size)+"... Exiting..."+'\x1b[0m')
        #         sys.exit()
        #     # if step_size % self.clk_cyc != 0:
        #     #     step_size = 
        #     for i in range(0, expCfg.N_scanPts):
        #         expCfg.scannedParam[i] = round(expCfg.scannedParam[i])
        #         if expCfg.scannedParam[i] % self.clk_cyc != 0:
        #             # print('Rounding scannedParam to a multiple of CLK_CYC (2ns)...')
        #             expCfg.scannedParam[i] += 1     # actually this should be expCfg.scannedParam[i] += round(expCfg.scanned[i] % self.clk_cyc)
            
        # if expCfg.sequence == 'MW_timing':
        #     if step_size < 5*self.clk_cyc:
        #         print("\x1b[1;33;41mERR: Step size = "+str(step_size)+"... Exiting..."+'\x1b[0m')
        #         sys.exit()
                
        #     for i in range(0, expCfg.N_scanPts):
        #         expCfg.scannedParam[i] = round(expCfg.scannedParam[i])
        #         if expCfg.scannedParam[i] % self.clk_cyc != 0:
        #             # print('Rounding scannedParam to a multiple of CLK_CYC (2ns)...')
        #             expCfg.scannedParam[i] += 1     # actually this should be expCfg.scannedParam[i] += round(expCfg.scanned[i] % self.clk_cyc)
        #------------------------------------------------------------
        
        #Pulse-sequence parameter checks:
        
        # #Check that "IQpadding" is a multiple of self.clk_cyc and >5*self.clk_cyc:
        # if expCfg.sequence in ['T2seq','XY8seq','correlSpecSeq']:
        #     if (expCfg.IQpadding < (5*self.clk_cyc)) or (expCfg.IQpadding%self.clk_cyc):
        #         print('Error: IQpadding is set to', expCfg.IQpadding,'which is either <',5*self.clk_cyc,'or not a multiple of',self.clk_cyc,'. Please edit IQpadding to ensure that it is >',5*self.clk_cyc,'ns and a multiple of',self.clk_cyc,'.')
                
        # # Check t_duration in esr_seq is a multiple of (2*self.clk_cyc): (Why?)
        # if expCfg.sequence in ['esr_seq', 'pesr_seq', 'ram_seq', 'T2_seq', 'modesr']:
        #     if expCfg.t_AOM%(2*self.clk_cyc):
        #         print(' Warning: t_duration set to ', expCfg.t_AOM,'ns, which is not an integer multiple of ',(2*self.clk_cyc),'ns. Rounding t_duration to nearest multiple of ',(2*self.clk_cyc),'ns...')
        #         expCfg.t_AOM = (2*self.clk_cyc)*round(expCfg.t_AOM/(2*self.clk_cyc))
        #         print('t_duration now set to ', expCfg.t_AOM,'ns')
        
        # #Check that t_readoutDelay and t_AOM are multiples of self.clk_cyc:
        # if expCfg.sequence in ['rabi_seq', 'pser_seq', 'T2_seq', 'XY8seq', 'correlSpecSeq', 'T1seq', 'T1ms0', 'ram_seq']:
        #     if (expCfg.ro_delay%self.clk_cyc) or (expCfg.ro_delay < (5*self.clk_cyc)):
        #         print('Error: t_readoutDelay is set to ', expCfg.t_readoutDelay,'ns, which is not a multiple of ',self.clk_cyc,'ns or <',(5*self.clk_cyc),'ns!')
        #         sys.exit()
        #     if (expCfg.t_AOM%self.clk_cyc) or (expCfg.t_AOM < (5*self.clk_cyc)):
        #         print('Error: t_AOM is set to ', expCfg.t_AOM,'ns, which is not a multiple of ',self.clk_cyc,'ns or <',(5*self.clk_cyc),'ns!,')
        #         sys.exit()
            
        #Check that tau0 in the correlation spectroscopy sequence is an integer multiple of 2*self.clk_cyc:
        # if expCfg.sequence == 'correlSpecSeq':
        #     if expCfg.tau0%(2*self.clk_cyc):
        #         print('Error: tau0 is set to ', expCfg.tau0,'ns, which is not a multiple of ',(2*self.clk_cyc),'ns. Please set tau0 to an integer multiple of ',(2*self.clk_cyc),'ns.')
        #         sys.exit()
                
        
        # #Number of XY8 repeats check:
        # if expCfg.sequence in ['XY8seq', 'correlSpecSeq']:
        #     if expCfg.N<1 or (not isinstance(expCfg.N, int)):
        #         print('Error: number of XY8 repeats, N, must be an integer >=1.')
        #         sys.exit()
        
        
        #Pi-pulse length checks:
        # if expCfg.sequence == 'T1ms1':
        #     if expCfg.t_pi<self.clk_cyc or expCfg.t_pi % self.clk_cyc:
        #         print('Error: requested pi pulse length ',expCfg.t_pi,'ns is either <',self.clk_cyc,'ns or not an integer multiple of ',self.clk_cyc,'ns.')
        #         sys.exit()
        
        # if expCfg.sequence in ['T2seq','XY8seq','correlSpecSeq']:
        #     # Check if the user has input a pi-pulse length which is shorter than (2*self.clk_cyc) or not a multiple of "2*self.clk_cyc":
        #     if expCfg.t_pi<(2*self.clk_cyc):
        #         print('Error: requested pi pulse length=',expCfg.t_pi,'ns is <',(2*self.clk_cyc),'ns. t_pi must be set to at least',(2*self.clk_cyc),'ns.')
        #         sys.exit()
                
        #     if expCfg.t_pi%(2*self.clk_cyc):
        #         print('\x1b[3;33;40m'+'Warning: t_pi set to ', expCfg.t_pi,'ns, which is not an integer multiple of ',(2*self.clk_cyc),'ns. Rounding t_pi to nearest multiple of ',(2*self.clk_cyc),'ns...')
        #         expCfg.t_pi = (2*self.clk_cyc)*round(float(expCfg.t_pi)/(2*self.clk_cyc))
        #         print('t_pi now set to ', expCfg.t_pi,'ns')
                
        # # Scan step-size checks: for ESR/Rabi/T1
        # step_size = expCfg.scannedParam[1]-expCfg.scannedParam[0]
        # if (expCfg.sequence in ['esr_seq', 'pesr_seq', 'modesr']):
        #     if ((step_size*1e6)%1):  # Round step to 1uHz if smaller; SRS freq res = 1uHz
        #         roundedFreqstep_size = (1e-6)*round((1e6)*step_size)
        #         expCfg.scannedParam[-1] = (expCfg.N_scanPts-1)*roundedFreqstep_size + expCfg.scannedParam[0] 
        #         print('\x1b[3;33;40m'+'Warning: Requested freq step is ',step_size,'Hz. Not integer multiple of the SRS freq resolution, 1uHz. Rounding step size to the nearest multiple of 1uHz.\nStep size is now',roundedFreqstep_size,'\n',expCfg.scanStartName,'= ',expCfg.scannedParam[0],' and \n',expCfg.scanEndName,'= ',expCfg.scannedParam[-1] + '\x1b[0m')
        #         expCfg.scannedParam = np.linspace(expCfg.scannedParam[0],expCfg.scannedParam[-1], expCfg.N_scanPts,endpoint= True)
                
        # if expCfg.sequence in ['rabi_seq', 'T1seq', 'T1ms0', 'rodelay', 'ram_seq', 'T2_seq'] or (expCfg.sequence == 'T2seq' and expCfg.numberOfPiPulses == 1):   # why to check for numberOfPiPulses==1??
        #     if step_size<self.clk_cyc:
        #         print('Error: requested time step =',step_size,'ns, which is shorter than',self.clk_cyc,'ns. Please change N_scanPts, or ',expCfg.scanStartName,' and ',expCfg.scanEndName,' to increase time step size.')
        #         sys.exit()
        #     else:
        #         # print('\t Step Size = '+str(step_size)+" > "+str(self.clk_cyc))
        #         None
                
        #     # If requested step size is >self.clk_cyc but not a multiple of self.clk_cyc:
        #     if (step_size%self.clk_cyc):
        #         roundedstep_size = self.clk_cyc*round(step_size/self.clk_cyc)
        #         expCfg.scannedParam[-1] = (expCfg.N_scanPts-1)*roundedstep_size + expCfg.scannedParam[0] 
        #         print('\x1b[38;2;250;200;0mWarning: requested time step is '+str(step_size)+'ns\x1b[0m, \u2260 <int>*'+str(self.clk_cyc)+'ns.\nRounding step size to the nearest multiple of '+str(self.clk_cyc)+'...\n\x1b[38;2;250;200;0mStep size now = '+str(roundedstep_size)+'\n'+expCfg.scanStartName+' = '+str(expCfg.scannedParam[0])+'\n'+expCfg.scanEndName,'= ',str(expCfg.scannedParam[-1])+'\x1b[0m')
        #         expCfg.scannedParam = np.linspace(expCfg.scannedParam[0],expCfg.scannedParam[-1], expCfg.N_scanPts,endpoint= True)
                
        if self.parameter_dict['seq']['sequence'] =='rabi_seq':
            # Pulseblaster bug - our PulseBlaster boards do not seem to be able to output 8ns pulses. So, check if we asked for 8ns and remove this point:
            if 8 in self.parameter_dict['scan']['values']:
                # expCfg.scannedParam = list(expCfg.scannedParam)
                self.parameter_dict['scan']['values'].remove(8)
                self.parameter_dict['scan']['Nscanpts'] = len(self.parameter_dict['scan']['values'])
                print('⚠ \x1b[38;2;250;200;0mWarning: will not collect data at 8ns scan point \x1b[0mdue to unofficial reports of a possible issue with some PB boards whereby the instruction for outputting 8ns pulses generates 10ns pulses. Removing the 8ns scan point from the list of scan points.\x1b[0m')
        
        # if (expCfg.sequence == 'XY8seq') or (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses > 1):
        #     #Check if requested scan step is too short or not a multiple of 2*self.clk_cyc:
        #     if step_size<(2*self.clk_cyc):
        #         print('Error: requested time step is ',step_size,'ns, which is shorter than ', (2*self.clk_cyc),'ns. Please change N_scanPts, or ',expCfg.scanStartName,' and ',expCfg.scanEndName,' to increase time step size.')
        #         sys.exit()
                
        #     # If requested step size is >2*self.clk_cyc but not a multiple of self.clk_cyc, round to nearest multiple of (2*self.clk_cyc) and warn user:
        #     if (step_size%(2*self.clk_cyc)):
        #         roundedstep_size = (2*self.clk_cyc)*round(step_size/(2*self.clk_cyc))
        #         expCfg.scannedParam[-1] =  (expCfg.N_scanPts-1)*roundedstep_size + expCfg.scannedParam[0] 
        #         print('\x1b[3;33;40m'+'Warning: requested time step is ',step_size,'ns, which is not an integer multiple of ',(2*self.clk_cyc),'ns. Rounding step size to the nearest multiple of ',(2*self.clk_cyc),':\n Step size is now ',roundedstep_size,'\n ',expCfg.scanStartName,'= ',expCfg.scannedParam[0],' and \n',expCfg.scanEndName,'= ',expCfg.scannedParam[-1])
        #         expCfg.scannedParam = np.linspace(expCfg.scannedParam[0],expCfg.scannedParam[-1], expCfg.N_scanPts,endpoint= True)
        
        # # Scan-start (minimum delay duration) checks:
        # if (expCfg.sequence in ['rabi_seq', 'correlSpecSeq', 'T2_seq']) or (expCfg.sequence == 'T2seq' and expCfg.numberOfPiPulses==1):
        # # Check if requested start delay/pulse length is positive and a multiple of "self.clk_cyc":
        #     if expCfg.scannedParam[0]<0:
        #         print('\x1b[38;2;250;2;50mERR: requested ', expCfg.scanStartName,'=', expCfg.scannedParam[0],'is <0. ', expCfg.scanStartName,' must be >=0.')
        #         sys.exit()
                
        #     if expCfg.scannedParam[0]%self.clk_cyc:
        #         print('\x1b[38;2;250;2;50mERR: ',expCfg.scanStartName,' is set to ', expCfg.scannedParam[0],', which is not a multiple of ',self.clk_cyc,'ns. Please set', expCfg.scanStartName,' to an integer multiple of ',self.clk_cyc,'ns.')
        #         sys.exit()
        
        #     if (expCfg.sequence == 'XY8seq') or (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses > 1):
        #     # Check if requested start delay/pulse length is a multiple of "2*self.clk_cyc":
        #         if expCfg.scannedParam[0]%(2*self.clk_cyc):
        #             print('Error: ',expCfg.scanStartName,' is set to ', expCfg.scannedParam[0],', which is not a multiple of ',(2*self.clk_cyc),'ns. Please set', expCfg.scanStartName,' to an integer multiple of ',(2*self.clk_cyc),'ns.')
        #             sys.exit()
        
        #     if expCfg.sequence == 'T1seq':
        #         if expCfg.scannedParam[0]<(expCfg.t_readoutDelay + self.clk_cyc*round((1*us)/self.clk_cyc)):
        #             print('Error: requested ',expCfg.scanStartName,' is too short.', expCfg.scanStartName,' must be >=0.')
        #             sys.exit()
        #         if expCfg.scannedParam[0]%self.clk_cyc:
        #             print('Error: ',expCfg.scanStartName,' is set to ', expCfg.scannedParam[0],', which is not a multiple of ',self.clk_cyc,'ns. Please set', expCfg.scanStartName,' to an integer multiple of',self.clk_cyc,'ns.')
        #             sys.exit()
        
            # if expCfg.sequence in ['T2_seq','XY8seq']:
            #     # Check if requested start delay is shorter than 3*(5*self.clk_cyc):
            #     if expCfg.scannedParam[0]<2*(5*self.clk_cyc):
            #         print('Error: requested ',expCfg.scanStartName,' of ',expCfg.scannedParam[0],'ns is too short. For this pulse sequence, ',expCfg.scanStartName,' must be set to at least',3*(5*self.clk_cyc),'ns')
            #         sys.exit()
        
            # if expCfg.sequence == 'T2seq':
            #     if (not isinstance(expCfg.numberOfPiPulses, int)) or (expCfg.numberOfPiPulses<1):
            #         print('Error: numberOfPiPulses must be a positive integer!')
            #         sys.exit()
                
            #     if expCfg.numberOfPiPulses == 1:
            #         if expCfg.scannedParam[0]<(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*self.clk_cyc)):
            #             # Check if start delay is too short to allow for PB timing resolution:
            #             print('Error: ',expCfg.scanStartName,' too short. For your pi_pulse length, ',expCfg.scanStartName,' must be at least', (2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*self.clk_cyc)),'ns.')
            #             sys.exit()
                    
        #         else: #no of pi pulses>1
        #             if expCfg.scannedParam[0]<(2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*self.clk_cyc))):
        #                 print('Error: ',expCfg.scanStartName,' too short. For your pi_pulse length, ',expCfg.scanStartName,' must be at least', (2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*self.clk_cyc))),'ns.')
        #                 sys.exit()
            
        #     if expCfg.sequence == 'XY8seq':
        #         if expCfg.scannedParam[0]<(2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*self.clk_cyc))):
        #             print('Error: ',expCfg.scanStartName,' too short. For your pi_pulse length, ',expCfg.scanStartName,' must be at least', (2*(2*expCfg.IQpadding + (3/4)*expCfg.t_pi + (5*self.clk_cyc))),'ns.')
        #             sys.exit()
                
        #     # Free precession time checks. The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi
        #     # or pi/2 pulse in the T2, XY8 and correlation spectroscopy sequences has to be an integer multiple of self.clk_cyc. The user-input 
        #     # free precession time is defined as the time between the center of subsequent pulses. Hence, for a given pi-pulse length, 
        #     # we check that the user has selected a starting free precession time which produces an edge-to-edge time that is a multiple of self.clk_cyc.
        #     # If not, we shift the free precession time vector by self.clk_cyc/2 and warn the user.
        #     if (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses == 1):
        #         if (expCfg.scannedParam[0]-(expCfg.t_pi/4))%self.clk_cyc:
        #             expCfg.scannedParam = [x+(self.clk_cyc/2) for x in expCfg.scannedParam]
        #             print('Warning: Each element of the scanned time vector has been shifted by '+str(self.clk_cyc/2)+'ns so that the rising-edge-to-rising-edge spacing between microwave pulses is a multiple of '+str(self.clk_cyc)+'ns.\
        # \nDetails: The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi or pi/2 pulse in the T2 and XY8 sequences \
        # has to be an integer multiple of '+str(self.clk_cyc)+'ns. The user-input free-precession time is defined as the time between the center of subsequent pulses.\
        # Your pi pulse length,'+str(expCfg.t_pi)+'ns, produces an edge-to-edge time of'+str(expCfg.scannedParam[0]-(expCfg.t_pi/4))+'ns (at the start of the scan), which is not a multiple of '+str(self.clk_cyc)+'ns.\
        # Hence, we shift the times by '+str(self.clk_cyc/2)+'ns.')	
        
        #     if (expCfg.sequence =='XY8seq') or (expCfg.sequence=='T2seq' and expCfg.numberOfPiPulses > 1):
        #         half_t_delay = expCfg.scannedParam[0]/2
        #         if (half_t_delay-(expCfg.t_pi/4))%self.clk_cyc:
        #             expCfg.scannedParam = [x+(self.clk_cyc/2) for x in expCfg.scannedParam]
        #             print('\x1b[3;33;40m'+'Warning: Each element of the scanned time vector has been shifted by ',self.clk_cyc/2,'ns so that the rising-edge-to-rising-edge spacing between microwave pulses is a multiple of ',self.clk_cyc,'ns.\
        # \nDetails: The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi or pi/2 pulse in the T2 and XY8 sequences \
        # has to be an integer multiple of ',self.clk_cyc,'ns. The user-input free-precession time is defined as the time between the center of subsequent pulses.\
        # Your pi pulse length,',expCfg.t_pi,'ns, produces an edge-to-edge time of', (half_t_delay-(expCfg.t_pi/4)),'ns (at the start of the scan), which is not a multiple of ',self.clk_cyc,'ns.\
        # Hence, we shift the times by ',self.clk_cyc/2,'ns.')
        
        #     if expCfg.sequence == 'correlSpecSeq':
        #         half_t_delay = expCfg.tau0/2
        #         if (half_t_delay-(expCfg.t_pi/4))%self.clk_cyc:
        #             expCfg.tau0 = expCfg.tau0 +(self.clk_cyc/2)*ns
        #             print('\x1b[3;33;40m'+'Warning: tau0 has been shifted by ',self.clk_cyc/2,'ns so that the rising-edge-to-rising-edge spacing between microwave pulses is a multiple of ',self.clk_cyc,'ns. tau0 is now set to', expCfg.tau0,'\
        # \nDetails: The spacing between the rising edge of a pi or pi/2 pulse and the rising edge of the subsequent pi or pi/2 pulse in the XY8 sequence \
        # has to be an integer multiple of ',self.clk_cyc,'ns. The user-input tau0 is defined as the time between the center of subsequent pi pulses in the XY8 sequence.\
        # For your pi pulse length,',expCfg.t_pi,'ns, your chose tau0 produces an edge-to-edge time of', half_t_delay-(expCfg.t_pi/4),'ns, which is not a multiple of ',self.clk_cyc,'ns.\
        # Hence, we shift the tau0 by ',self.clk_cyc/2,'ns.')
    #-----------------------------------------------------------------------

    # ------PulseBlaster Sequences------
    """Variables:
        * sequence: (String) Name of the sequence to be executed
        * args: (List) of time durations
    """
    @staticmethod
    def make_sequence(instr, sequence, args):
        # coming from PBcontrol.py > PulseBlaster.PB_program__()
        var = {}
        if instr == 'cam' or instr == 'cam_levelm' or instr == 'cam_syncm':
            # the incoming 'instr' can only be 'cam' or 'cam_levelm' to simplify the function
            var = {'esr_seq':cam_seq.make_esr_seq_camera,   'rabi_seq':cam_seq.make_rabi_seq_camera_FL,
                'pesr_seq': cam_seq.make_pulsed_esr_seq_camera_FL,  'T1ms0': cam_seq.make_t1_seq_camera,
                'T2_seq': cam_seq.make_t2_seq}
            
        elif instr == 'cam_level1':
            var = {'esr_seq':cam_seq.make_esr_seq_camera_level_trigger,   'rabi_seq':cam_seq.make_rabi_seq_camera_level_trigger,
                'pesr_seq': cam_seq.make_pulsed_esr_seq_camera_level_trigger,  'T1ms0': cam_seq.make_t1_seq_camera_level_trigger}
            
        elif instr == 'diode':
            var = {'esr_seq':daq_seq.make_esr_seq,   'modesr':daq_seq.make_mod_esr_seq,
                'rabi_seq':daq_seq.make_rabi_seq,    'spin_echo':daq_seq.make_echo_seq_FL,
                'T2seq':daq_seq.makeT2Seq0,  'pesr_seq':daq_seq.make_pulsed_esr_seq,
                'aom_timing':daq_seq.make_aom_timing_seq,    'rodelay':daq_seq.make_opt_readout_time_sequence,
                'drift_seq':daq_seq.make_drift_analysis_sequence,    'MW_timing':daq_seq.make_MW_timing_seq,
                'T1ms0':daq_seq.make_t1_ms0_seq_MW2,  'T1ms1':daq_seq.make_t1_ms1_seq_FL,
                # 'T1ms0_train': daq_seq.make_t1ms0_pulse_train,
                'T2_seq':daq_seq.make_t2_seq_MW,     'ram_seq':daq_seq.make_ramsey_seq_MW,
                'simult_samp':daq_seq.check_simult_sampling,     'double_mod':daq_seq.make_double_mod_sequence_lcm,
                'diff_mod':daq_seq.make_diff_mod_sequence, 'esr_dig_mod_seq':daq_seq.make_dig_mod_odmr_sequence,
                'rabi_dig_mod_seq':daq_seq.make_dig_mod_rabi_sequence,
                'rabi_contrast_seq':daq_seq.make_rabi_contrast_sequence,
                }
            
        # print(var)
        if sequence in var.keys():
            # instead of directly assigning the name of the function to the dictionary, return the name when the function is called
            # at the caller function, receive it only when required
            # this change can me some functions in sequececontrol and PulseBlaster @staticmethod
            # self.parameter_dict['seq']['seqctrl_name'] = 
            # print(args)
            return var[sequence].__name__, var[sequence](*args)
        else:
            print('❌ Wrong sequence.. Exiting!!')
            sys.exit()
    
    @staticmethod
    def plot_data(x,y, xlabel, ylabel, formatting, yTicks, yTickLabels, title):
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

    @staticmethod
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
        t_ns = [0, 0]; t_us = []
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

    @staticmethod
    def view_sequence(instr, sequence, seqArgList, only_plot=False, parameter=[0], seq_no_plot:list=[0], plot_dpi:int=100):
        # the_list=[]
        for i in range(0, len(seq_no_plot)):
            if only_plot == False:
                if  sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                    seqArgList[0] = parameter[seq_no_plot[i]]
            # plt.figure(num=f"{sequence} sequence plot {time.strftime(' [%H%M%S]', time.localtime())}", dpi = plot_dpi)
            plt.figure()
            _, the_list = PBcontrol.PulseBlaster.PB_program(instr, sequence, seqArgList)  # last element of seqArgList is PBchannels
            for j in range(0, len(the_list)):
                instructionList = the_list[j][0]
                inst_times = the_list[j][4]
                # instructionList=pbctrl.PB_program_camera(sequence,seqArgList)[0][0]
                # inst_times = pbctrl.PB_program_camera(sequence,seqArgList)[0][4]
                # plt.subplot(1,2,(j+1))
                [t_us,channelPulses,yTicks] = sequencecontrol.plot_sequence(instructionList, seqArgList[-1])
                for channel in channelPulses:
                    plt.plot(t_us, list(channel))
                plt.yticks(yTicks, seqArgList[-1].keys())          # Include the names of the PB channels
                plt.xlabel('Time (us)')
                plt.ylabel('Channel')
                plt.title(f"{sequence} Pulse Seq plot. Param val @ {str(seq_no_plot[i])}: {str(parameter[seq_no_plot[i]])} ns\n\
                          Transitions at: {str([vals for vals in inst_times.values()])}", fontsize=12)
        # Way to display the total no of instructions in the plot, since it will vary with different scanPts:
            # print("Total %d inst" % len(inst_times.keys()))

    @staticmethod
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
        * channel: (PBchannel) an element in the 'allPBchannels' list containing channel_number, start_times and end_times
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
            channelMask = aPBchannel.channel_number
            if channelMask < 0:
                continue  # skip this PBchannel if channel_number is negative
            end_times = [startTime + pulseDuration for startTime, pulseDuration in zip(aPBchannel.start_times, aPBchannel.pulse_durations)]
            eventTimes = aPBchannel.start_times+end_times
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
    @staticmethod
    def param_err_check(instr, sequence, seqArgList, parameter=[0,1], Nscanpts=1):
        # Trial run over all parameters to check whether the durations of all Inst < 10ns (=5*self.clk_cyc)
        n_error = 0; param = []; not_param = []
        print("🔃 Checking sequences for errors...")
        for i_scanpt in range (0, Nscanpts):     # scan over all the scannedParam values
            seqArgList[0] = parameter[i_scanpt]
            _, the_list = PBcontrol.PulseBlaster.PB_program(instr, sequence, seqArgList[0:-1], err_check=True)  # last element of seqArgList is PBchannels
            for i in range(0, len(the_list)):
                instructionList = the_list[i][0]    # eta chai sudhu oi parameter er sequence ta plot korar jonno...
                seq_error_count = the_list[i][1]
                inst_error_no = the_list[i][2]
                error_times = the_list[i][3]
                inst_times = the_list[i][4]
                
                # (For errors) Plot the sequence with the title format: scannedParam[i], inst_no, 
                if seq_error_count > 0:
                    n_error += 1
                    # # -------plot and mark the region where the sequence becomes < 10ns-------
                    # plt.figure()
                    # [t_us,channelPulses,yTicks] = plot_sequence(instructionList, PBchannels)
                    # for channel in channelPulses:
                    #     plt.plot(t_us, list(channel))
                    # plt.xlabel('Time (us)')
                    # plt.ylabel('Channels')
                    # # print(yticks)
                    # # print(PBchannels.keys())
                    # # plt.yticks(yTicks, PBchannels.keys())
                    # plt.title('Err: scanParam='+str(parameter[i_scanpt])+'. Check Inst #: '+str([i for i in inst_error_no])+ ' @ ' + str([i/1e3 for i in error_times]) + 'us.\n Transitions at: ' + str([i for i in inst_times.values()]), color='r',fontsize=10)
                    # print(' Removing \x1b[38;2;250;250;0m'+str(parameter[i_scanpt])+'\x1b[0m')
                    # # ---------plot done--------
                    if parameter[i_scanpt] in param:
                        param.remove(parameter[i_scanpt])
                    not_param.append(parameter[i_scanpt])
                else:
                    if parameter[i_scanpt] not in not_param and parameter[i_scanpt] not in param:
                        param.append(parameter[i_scanpt])
        return [n_error, param]
    
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

