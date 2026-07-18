# control_counter_sequences.py

import sys, numpy as np
from experiment_config import PB_CLK, PBpins
from collections import namedtuple
ns, us, ms = 1., 1e3, 1e6

PBchannel = namedtuple('PBchannel', ['channel_number', 'start_times', 'pulse_durations'])

def compile_sequence(pbchannels: list[PBpins], channel_timings: dict) -> list[PBchannel]:
    """Compact sequence maker."""
    allPBchannels = []
    active_pbchannels = {m.name.lower(): m.value for m in pbchannels}       # type: ignore
    # Build channels
    for channel, (starts, durations) in channel_timings.items():
        if channel in active_pbchannels:
            allPBchannels.append(PBchannel(active_pbchannels[channel], starts, durations))
    
    return allPBchannels

conv_clk_sep = 5*us
pulse_width = 100*ns

clk_cyc = 1e3/PB_CLK            # Time resolution in ns
ONE_PERIOD = 0x200000           # 23/22/21/20 = 0010b = 2^21d = 200000x
TWO_PERIOD = 0x400000           # 23/22/21/20 = 0100 = 2^22
THREE_PERIOD = 0x600000         # 23/22/21/20 = 0110
FOUR_PERIOD = 0x800000          # 23/22/21/20 = 1000
FIVE_PERIOD = 0xA00000           # 23/22/21/20 = 1010

# TODO: DEBUG: return channel_timings variable and save to param file..

def make_esr_seq(seq_dur, pb_channels: list):
    gate2clock = 1 *us   # wait time after gate OFF for sample clock
    offset_time = 10 *us

    channel_timings = {
        'start_trig':   ([0],
                         [pulse_width]),
        'laser':        ([0],
                         [seq_dur + pulse_width + gate2clock]),
        'mw':           ([0],
                         [seq_dur/2]),
        
        'pause_trig':   ([offset_time, seq_dur/2 + offset_time],
                         [seq_dur/2-offset_time]*2),
        
        # counter input special
        'gate':   ([offset_time-60*ns, seq_dur/2 + offset_time - 60*ns],
                         [seq_dur/2-offset_time]*2),
        'samp_clk':     ([seq_dur/2 + gate2clock, seq_dur + gate2clock],
                         [pulse_width]*2),
    }

    return compile_sequence(pbchannels=pb_channels, channel_timings=channel_timings)

def make_rabi_seq(t_MW, t_AOM, t_ro_delay, AOM_lag, MW_lag, pb_channels):
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
    gate2clock = 1*us
    if t_MW < 10:
        t_drive = 10+0*ns
    else:
        t_drive = t_MW + 0*ns

    # apd_pulse = [t_drive+AOM_lag+t_ro_delay, 2*t_drive+t_AOM+AOM_lag+t_ro_delay]
    apd_pulse = [t_drive+AOM_lag, 2*t_drive+t_AOM+AOM_lag]
    # seq_dur = (2*t_drive + 2*t_AOM) if (2*t_drive + 2*t_AOM) > (2*t_drive+t_AOM+AOM_lag+t_ro_delay+pulse_width) else (2*t_drive+t_AOM+AOM_lag+t_ro_delay+pulse_width)
    
    channel_timings = {
        'start_trig':       ([0],
                             [pulse_width]),
        'laser':            ([t_drive, 2*t_drive+t_AOM],
                             [t_AOM,t_AOM]),
        'mw':               ([0*us+AOM_lag-MW_lag],
                             [t_MW]),
        'pause_trig':       ([apd_pulse[0], apd_pulse[1]],
                             [t_ro_delay]*2),
        # 'samp_clk':       ([apd_pulse[0], apd_pulse[1]], [pulse_width]*2),
        # 'samp_clk':         ([0.2*ms, 1.7*ms], [seq_dur/2-0.2*ms]*2),
        # 'lia':              ([0], [seq_dur]),
        # 'bx':               ([0], [seq_dur]),
        # 'by':               ([0], [seq_dur]),
        # 'bz':               ([0], [seq_dur]),
        # for counter
        'gate':       ([apd_pulse[0] - 60*ns, apd_pulse[1] - 60*ns],
                             [t_ro_delay + 120*ns]*2),
        'samp_clk':     ([apd_pulse[0] + t_ro_delay + gate2clock, apd_pulse[1] + t_ro_delay + gate2clock],
                         [pulse_width]*2),
    }
    return compile_sequence(pbchannels=pb_channels, channel_timings=channel_timings)

def make_pulsed_esr_seq(t_AOM, t_ro_delay, AOM_lag, MW_lag, t_pi, pb_channels):
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

    t_flip = 2*us + t_pi
    gate2clock = 1*us
    apd_pulse = [t_flip+AOM_lag+t_ro_delay, 2*t_flip+t_AOM+AOM_lag+t_ro_delay]
    
    channel_timings = {
        'start_trig':       ([0],
                             [pulse_width]),
        'laser':            ([t_flip, 2*t_flip + t_AOM],
                             [t_AOM, t_AOM]),
        'mw':               ([AOM_lag - MW_lag],
                             [t_pi]),
        'pause_trig':       ([apd_pulse[0], apd_pulse[1]],
                             [t_ro_delay]*2),
        # 'samp_clk':       ( [t_flip+AOM_lag, 2*t_flip+t_AOM+AOM_lag], [t_AOM, t_AOM]),
        # 'samp_clk':         ([0.2*ms, 1.7*ms], [seq_dur/2-0.2*ms]*2),
        # 'lia':              ([0], [seq_dur]),
        # 'bx':               ([0], [seq_dur]),
        # 'by':               ([0], [seq_dur]),
        # 'bz':               ([0], [seq_dur]),
        # for counter
        'gate':       ([apd_pulse[0] - 60*ns, apd_pulse[1] - 60*ns],
                             [t_ro_delay + 120*ns]*2),
        'samp_clk':     ([apd_pulse[0] + t_ro_delay + gate2clock, apd_pulse[1] + t_ro_delay + gate2clock],
                         [pulse_width]*2),
    }
    return compile_sequence(pbchannels=pb_channels, channel_timings=channel_timings)

#------------------------------------------------------------------------------
def make_echo_seq_MW(delay_2nd_half, delay_1st_half, t_AOM, ro_delay, AOM_lag, MW_lag, t_pi, pb_channels):
    """ Spin-echo seq. 
    Ref = FL w/o MW"""
        
    gate2clock = 1*us
    t_delay = delay_1st_half + delay_2nd_half + 2*t_pi# + 5*us     # increase the actual laser off time to accomodate the non-zero width of the pi/2 and pi pulses so that the actual precession time is as defined by the param variable of mainControl
    apd_pulse = [t_delay+AOM_lag+ro_delay, 2*t_delay+t_AOM+AOM_lag+ro_delay]
    
    channel_timings = {
        'start_trig':       ([0],
                             [pulse_width]),
        'laser':            ([t_delay, 2*t_delay+t_AOM],
                             [t_AOM,t_AOM]),
        'mw':               ([0*us+AOM_lag-MW_lag, 0*us+AOM_lag+t_pi/2+delay_1st_half-MW_lag, 0*us+AOM_lag+(t_delay-0*us)-MW_lag-t_pi/2],
                             [t_pi/2, t_pi, t_pi/2]),
        'pause_trig':       ([apd_pulse[0], apd_pulse[1]],
                             [ro_delay]*2),
        # 'samp_clk':       ( [t_flip+AOM_lag, 2*t_flip+t_AOM+AOM_lag], [t_AOM, t_AOM]),
        # 'samp_clk':         ([0.2*ms, 1.7*ms], [seq_dur/2-0.2*ms]*2),
        # 'lia':              ([0], [seq_dur]),
        # 'bx':               ([0], [seq_dur]),
        # 'by':               ([0], [seq_dur]),
        # 'bz':               ([0], [seq_dur]),
        # for counter
        'gate':       ([apd_pulse[0] - 60*ns, apd_pulse[1] - 60*ns],
                             [ro_delay + 120*ns]*2),
        'samp_clk':     ([apd_pulse[0] + ro_delay + gate2clock, apd_pulse[1] + ro_delay + gate2clock],
                         [pulse_width]*2),
    }
    return compile_sequence(pbchannels=pb_channels, channel_timings=channel_timings)
