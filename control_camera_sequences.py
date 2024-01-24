# control_camera_sequences.py

from connectionConfig import laser, samp_clk, start_trig, conv_clk, I, Q, MW, PBclk, camera
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

clk_cyc = 1e3/PBclk       # Time resolution in ns
# Short pulse flags: Switch ON bits 21-23 of the Control word to enable short pulse feature
ONE_PERIOD = 0x200000           # 23/22/21/20 = 0010b = 2^21d = 200000x
TWO_PERIOD = 0x400000           # 23/22/21/20 = 0100 = 2^22
THREE_PERIOD = 0x600000         # 23/22/21/20 = 0110
FOUR_PERIOD = 0x800000          # 23/22/21/20 = 1000
FIVE_PERIOD = 0xA00000           # 23/22/21/20 = 1010


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
