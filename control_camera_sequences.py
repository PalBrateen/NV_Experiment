# control_camera_sequences.py

from connectionConfig import laser, samp_clk, start_trig, I, Q, MW, PBclk, camera
from spinapi import ns, us, ms
from collections import namedtuple

""" Create a tuple of 'PBchannel' type with entries:
    1. channelNumber
    2. startTimes
    3. pulseDurations
"""
PBchannel = namedtuple('PBchannel', ['channelNumber', 'startTimes', 'pulseDurations'])
conv_clk_sep = 10*us
pulse_width = 100*ns
camera_delay = (87.7+0)*us

clk_cyc = 1e3/PBclk       # Time resolution in ns
# Short pulse flags: Switch ON bits 21-23 of the Control word to enable short pulse feature
ONE_PERIOD = 0x200000           # 23/22/21/20 = 0010b = 2^21d = 200000x
TWO_PERIOD = 0x400000           # 23/22/21/20 = 0100 = 2^22
THREE_PERIOD = 0x600000         # 23/22/21/20 = 0110
FOUR_PERIOD = 0x800000          # 23/22/21/20 = 1000
FIVE_PERIOD = 0xA00000           # 23/22/21/20 = 1010

def make_esr_seq_camera_level_trigger(seq_dur):
    exp_time = 1.05*ms
    wait_time = 1.05*ms
    camera_channel = PBchannel(camera, [0, wait_time+exp_time], [exp_time, exp_time])
    
    laser_start_times = [0]; laser_pulse_durations = [2*(wait_time+exp_time)]
    laser_channel = PBchannel(laser, [i for i in laser_start_times], [i for i in laser_pulse_durations])
    
    MW_start_times = [0,2*(wait_time+exp_time)-1*us]; MW_pulse_durations = [1.5*exp_time, 1.5*us]
    MW_channel = PBchannel(MW, [i for i in MW_start_times], [i for i in MW_pulse_durations])

    allPBchannels = [camera_channel, laser_channel, MW_channel]
    return allPBchannels

def make_esr_seq_camera(seq_dur):
    delay = 500*ns
    seq_dur -= 500*ns;
    # readout_width = clk_cyc*round(100*ns/clk_cyc)
    laser_start_times = [delay, delay]; laser_pulse_durations = [seq_dur,seq_dur]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel];
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.

    MW_start_times = [delay,delay]; MW_pulse_durations = [seq_dur,seq_dur]
    MWchannel = PBchannel(MW, [MW_start_times[0]], [MW_pulse_durations[0]])
    allPBchannels_sig.extend([MWchannel])
    MWchannel = PBchannel(MW, [MW_start_times[1]], [MW_pulse_durations[1]])
    # allPBchannels_ref.extend([MWchannel])
    
    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

def make_rabi_seq_camera_level_trigger(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag):
    allPBchannels = []

    if t_MW < 10:
        t_drive = 10+0*ns
    else:
        t_drive = t_MW + 0*ns
    exp_time = 1.05*ms
    wait_time = 3*ms
    camera_channel = PBchannel(camera, [0, wait_time+exp_time], [exp_time, exp_time])
    
    laser_offset = camera_delay-t_AOM-AOM_lag
    laser_start_times = [laser_offset, laser_offset+t_AOM+t_MW, (wait_time+exp_time)+laser_offset, (wait_time+exp_time)+laser_offset+t_AOM+t_MW]; laser_pulse_durations = [t_AOM, t_ro_delay, t_AOM, t_ro_delay]

    laser_channel = PBchannel(laser, [i for i in laser_start_times], [i for i in laser_pulse_durations])
    # laser_channel = PBchannel(laser, [],[])
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    MW_offset = camera_delay-MW_lag
    if t_MW <= 5*clk_cyc and t_MW > 0:
        MW_start_times = [MW_offset, 2*(wait_time+exp_time)-t_MW]; MW_pulse_durations = [5*clk_cyc, 5*clk_cyc]
        # SHORT PULSE duration
        shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
        shortPulseChannel = PBchannel(shortpulseFLAG, [i for i in MW_start_times], [i for i in MW_pulse_durations])
        allPBchannels = [shortPulseChannel]
    else:
        MW_start_times = [MW_offset, 2*(wait_time+exp_time)-t_MW]; MW_pulse_durations = [t_MW, t_MW]
    MW_channel = PBchannel(MW, [i for i in MW_start_times], [i for i in MW_pulse_durations])
    
    allPBchannels.extend([camera_channel, laser_channel, MW_channel])
    return allPBchannels

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

def make_t2_seq(t_delay, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi):
    """ T2 seq. 
    Ref = FL w/o MW"""
    t_delay = t_delay + 2*t_pi     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl

    laser_start_times = [t_delay, t_delay]; laser_pulse_durations = [t_AOM, t_AOM]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel]
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]

    MW_start_times = [[AOM_lag-MW_lag, AOM_lag+t_delay/2-MW_lag-t_pi/2, t_delay+AOM_lag-MW_lag-t_pi/2]]
    MW_pulse_durations = [[t_pi/2, t_pi, t_pi/2]]
    MWchannel = PBchannel(MW, MW_start_times[0], MW_pulse_durations[0])
    allPBchannels_sig.extend([MWchannel])
    
    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels


# bhul achhe kichhu ekta...
# def make_echo_seq_MW(delay_2nd_half, delay_1st_half, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi):
#     """ Spin-echo seq. 
#     """

#     laser_start_times = [t_drive, t_drive]; laser_pulse_durations = [t_AOM, t_AOM]
#     laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
#     allPBchannels_sig = [laser_channel];
#     laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
#     allPBchannels_ref = [laser_channel]

#     t_delay = delay_1st_half + delay_2nd_half + 2*t_pi# + 5*us     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
#     apd_pulse = [t_delay+AOM_lag+ro_delay, 2*t_delay+t_AOM+AOM_lag+ro_delay]

#     laser_channel = PBchannel(laser, [t_delay, 2*t_delay+t_AOM], [t_AOM,t_AOM])
#     start_trig_channel = PBchannel(start_trig, [0], [pulse_width])
#     MWchannel = PBchannel(MW, [0*us+AOM_lag-MW_lag, 0*us+AOM_lag+t_pi/2+delay_1st_half-MW_lag, 0*us+AOM_lag+(t_delay-0*us)-MW_lag-t_pi/2], [t_pi/2, t_pi, t_pi/2])
#     # The third MW pulse 
#     samp_clk_channel = PBchannel(samp_clk, [apd_pulse[0]-pulse_width, apd_pulse[1]-pulse_width], [pulse_width, pulse_width])
    
#     allPBchannels = [MWchannel, laser_channel, samp_clk_channel, start_trig_channel, conv_clk_channel]
#     return allPBchannels

def make_pulsed_esr_seq_camera_level_trigger(t_AOM, t_ro_delay, AOM_lag, MW_lag, t_pi):
    exp_time = 1.05*ms
    camera_channel = PBchannel(camera, [0, 2*exp_time,], [exp_time, exp_time])
    
    laser_offset = camera_delay-t_AOM-AOM_lag
    laser_start_times = [laser_offset, laser_offset+t_AOM+t_pi, 2*exp_time+laser_offset, 2*exp_time+laser_offset+t_AOM+t_pi]; laser_pulse_durations = [t_AOM, t_ro_delay, t_AOM, t_ro_delay]
    laser_channel = PBchannel(laser, [i for i in laser_start_times], [i for i in laser_pulse_durations])
    
    MW_offset = camera_delay-MW_lag
    MW_start_times = [MW_offset, 4*exp_time-t_pi]; MW_pulse_durations = [t_pi, t_pi]
    MW_channel = PBchannel(MW, [MW_start_times], [MW_pulse_durations])

    allPBchannels = [camera_channel, laser_channel, MW_channel]
    return allPBchannels

def make_pulsed_esr_seq_camera_FL(t_AOM, t_ro_delay, AOM_lag, MW_lag, t_pi):
    """
    Make pulse sequence for Rabi oscillations

    Parameters
    ----------
    t_AOM : float
        AOM 'ON' time or AOM width.
    t_ro_delay : float
        Length of (AOM) readout pulse.
    AOM_lag: float
        
    MW_lag: float
        
    t_pi : float
        Width of MW pi pulse.
    Returns
    -------
    allPBchannels : list
        A list of PBchannel object types.

    """
    t_flip = t_pi
    
    laser_start_times = [t_flip, t_flip]; laser_pulse_durations = [t_AOM, t_AOM]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel];
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]
    # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    
    MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [t_flip for i in range(0,2)]
    MWchannel = PBchannel(MW, [MW_start_times[0]], [MW_pulse_durations[0]])
    allPBchannels_sig.extend([MWchannel])
    MWchannel = PBchannel(MW, [MW_start_times[1]], [MW_pulse_durations[1]])
    # allPBchannels_ref.extend([MWchannel])
    
    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

def make_t1_seq_camera_level_trigger(relax_time, t_AOM, AOM_lag, rodelay, MW_lag, t_pi):
    # relax_time += t_pi
    exp_time = 1.05*ms
    camera_channel = PBchannel(camera, [0, 2*exp_time,], [exp_time, exp_time])
    
    laser_offset = camera_delay-t_AOM-AOM_lag
    laser_start_times = [laser_offset, laser_offset+t_AOM+relax_time, 2*exp_time+laser_offset, 2*exp_time+laser_offset+t_AOM+relax_time]; laser_pulse_durations = [t_AOM, rodelay, t_AOM, rodelay]
    laser_channel = PBchannel(laser, [i for i in laser_start_times], [i for i in laser_pulse_durations])
    
    MW_offset = camera_delay-MW_lag
    MW_start_times = [MW_offset,4*exp_time-t_pi]; MW_pulse_durations = [t_pi, t_pi]
    MW_channel = PBchannel(MW, [i for i in MW_start_times], [i for i in MW_pulse_durations])

    allPBchannels = [camera_channel, laser_channel, MW_channel]
    return allPBchannels

def make_t1_seq_camera(relax_time, dur_AOM, AOM_lag, rodelay, MW_lag, t_pi):
    # laser_start_times = [relax_time, 0]; laser_pulse_durations = [dur_AOM, relax_time + dur_AOM]
    relax_time += t_pi
    laser_start_times = [relax_time, relax_time]; laser_pulse_durations = [dur_AOM, dur_AOM]
    laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
    allPBchannels_sig = [laser_channel]
    laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
    allPBchannels_ref = [laser_channel]

    MW_start_times = [AOM_lag-MW_lag, AOM_lag-MW_lag]; MW_durations = [t_pi, t_pi]
    MW_channel = PBchannel(MW, [MW_start_times[0]], [MW_durations[0]])
    allPBchannels_sig.extend([MW_channel])
    MW_channel = PBchannel(MW, [MW_start_times[1]], [MW_durations[1]])
    # allPBchannels_ref.extend([MW_channel])

    allPBchannels = [allPBchannels_sig, allPBchannels_ref]
    return allPBchannels

# def make_rabi_seq_camera_MW(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag):
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
#     AOM_lag: float
        
#     MW_lag: float
        
#     Returns
#     -------
#     allPBchannels : list
#         A list of PBchannel object types.

#     """
#     if t_MW < 10:
#         t_drive = 10 +0*ns
#     else:
#         t_drive = t_MW + 0*ns
    
#     laser_start_times = [t_drive, t_drive]; laser_pulse_durations = [t_AOM, t_AOM]
#     laser_channel = PBchannel(laser, [laser_start_times[0]], [laser_pulse_durations[0]])
#     allPBchannels_sig = [laser_channel];
#     laser_channel = PBchannel(laser, [laser_start_times[1]], [laser_pulse_durations[1]])
#     allPBchannels_ref = [laser_channel]
#     # if the MW pulse duration is less than 5*clk_cyc, use the SHORT PULSE feature. Otherwise, use the original time duration.
    
#     if t_MW <= 5*clk_cyc and t_MW > 0:
#         MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [5*clk_cyc for i in range(0,2)]
#         # SHORT PULSE duration
#         shortpulseFLAG = int((t_MW/2)*ONE_PERIOD)
#         shortPulseChannel = PBchannel(shortpulseFLAG, [MW_start_times[0]], [MW_pulse_durations[0]])
#         allPBchannels_sig.extend([shortPulseChannel])
#         shortPulseChannel = PBchannel(shortpulseFLAG, [MW_start_times[1]], [MW_pulse_durations[1]])
#         # allPBchannels_ref.extend([shortPulseChannel])
#     else:
#         MW_start_times = [0*us+AOM_lag-MW_lag for i in range(0,2)]; MW_pulse_durations = [t_MW for i in range(0,2)]
#     MWchannel = PBchannel(MW, [MW_start_times[0]], [MW_pulse_durations[0]])
#     allPBchannels_sig.extend([MWchannel])
#     MWchannel = PBchannel(MW, [MW_start_times[1]], [MW_pulse_durations[1]])
#     # allPBchannels_ref.extend([MWchannel])
    
#     allPBchannels = [allPBchannels_sig]#, allPBchannels_ref]
#     return allPBchannels
