# control_daq_sequences.py

from connectionConfig import laser, samp_clk, start_trig, I, Q, MW, PBclk, camera, bx, by, bz
from spinapi import ns, us, ms
from collections import namedtuple; import sys; import numpy as np

PBchannel = namedtuple('PBchannel', ['channelNumber', 'startTimes', 'pulseDurations'])
conv_clk_sep = 5*us
pulse_width = 100*ns

clk_cyc = 1e3/PBclk       # Time resolution in ns
ONE_PERIOD = 0x200000           # 23/22/21/20 = 0010b = 2^21d = 200000x
TWO_PERIOD = 0x400000           # 23/22/21/20 = 0100 = 2^22
THREE_PERIOD = 0x600000         # 23/22/21/20 = 0110
FOUR_PERIOD = 0x800000          # 23/22/21/20 = 1000
FIVE_PERIOD = 0xA00000           # 23/22/21/20 = 1010

def make_esr_seq(seq_dur):
    seq_dur = 2*seq_dur;    trig_width = clk_cyc*round(100*ns/clk_cyc)
    # readout_width = clk_cyc*round(100*ns/clk_cyc)
    readout_buffer = clk_cyc*round(100*ns/clk_cyc)
    # pd_pulse refers to the readout pulse timings; in the case of multi-channel acquisition, it refers to the last channel readout pulse...
    pd_pulse = [seq_dur/2-readout_buffer, seq_dur-readout_buffer] # seq_dur/2-readout_buffer
    # pd_pulse.extend([(seq_dur/2-readout_buffer)/2, 1.5*(seq_dur/2-readout_buffer)])
    # pd_pulse.extend([seq_dur/2-readout_buffer-1*us, seq_dur-readout_buffer-1*us])
    laser_channel = PBchannel(laser, [0], [seq_dur])
    MW_channel = PBchannel(MW, [0], [seq_dur/2]) # seq_dur/2
    start_trig_channel = PBchannel(start_trig, [0], [trig_width])

    samp_clk_channel = PBchannel(samp_clk, [(pulse-conv_clk_sep-pulse_width) for pulse in pd_pulse], [trig_width for i in range(0,len(pd_pulse))]) #-conv_clk_sep-pulse_width
    # conv_clk_channel = PBchannel(conv_clk, [pd_pulse[0]-conv_clk_sep, pd_pulse[0], pd_pulse[1]-conv_clk_sep, pd_pulse[1]], [trig_width for i in range(0,4)])
    
    bx_channel = PBchannel(bx, [0], [seq_dur])
    by_channel = PBchannel(by, [0], [seq_dur])
    bz_channel = PBchannel(bz, [0], [seq_dur])
    allPBchannels = [laser_channel, samp_clk_channel, MW_channel, start_trig_channel, bx_channel, by_channel, bz_channel]
    # allPBchannels.extend([conv_clk_channel])
    return allPBchannels
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
    seq_dur = (2*t_drive + 2*t_AOM) if (2*t_drive + 2*t_AOM) > (2*t_drive+t_AOM+AOM_lag+t_ro_delay+pulse_width) else (2*t_drive+t_AOM+AOM_lag+t_ro_delay+pulse_width)
    
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
        MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag], [t_MW])
        allPBchannels = [MWchannel]
    
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    # conv_clk_channel = PBchannel(conv_clk, [0], [])
    samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0], apd_pulse[1]], [pulse_width for i in range(0,len(apd_pulse))])
    bx_channel = PBchannel(bx, [0], [seq_dur])
    by_channel = PBchannel(by, [0], [seq_dur])
    bz_channel = PBchannel(bz, [0], [seq_dur])
    allPBchannels.extend([laser_channel, samp_clk_channel, start_trig_channel, bx_channel, by_channel, bz_channel])
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
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]#, conv_clk_channel]
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
    
    # conv_clk_channel = PBchannel(conv_clk, [apd_pulse[0], apd_pulse[0]+conv_clk_sep, apd_pulse[1], apd_pulse[1]+conv_clk_sep], [pulse_width for i in range(0,4)])
    
    allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel]#, conv_clk_channel]
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

