# sequencecontrol.py
"""
PulseBlaster sequence construction and visualization.

Three public functions:

    make_sequence(instr, sequence, args)
        Dispatch to the correct sequence builder in control_daq_sequences
        or control_camera_sequences.  Returns (func_name, allPBchannels).

    event_cataloguer(pb_channels)
        Convert a list of PBchannel namedtuples into a time-sorted dict
        of {event_time: combined_bitmask}.  This is the input to
        PBcontrol.create_PBinstruction().

    plot_sequence(instructions, channel_map)
        Render a PB instruction list as a multi-channel timing diagram.
        Returns (t_us, channel_traces, y_ticks) for matplotlib.

    view_sequence(instr, sequence, seq_args, ...)
        Quick visualization: compile a PB program and plot it.

All validation (parameter bounds, Cartesian-product PB timing checks)
lives in experiment_base.DiodeExperiment.validate().
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

import control_daq_sequences as daq_seq
import control_camera_sequences as cam_seq
import PBcontrol


# ═══════════════════════════════════════════════════════════════════════
#  Sequence dispatch tables
# ═══════════════════════════════════════════════════════════════════════

# Maps (instr_type, sequence_name) → builder function.
# Each builder takes (*args) and returns a list of PBchannel namedtuples.

_DIODE_SEQUENCES = {
    'esr_seq':              daq_seq.make_esr_seq,
    'modesr':               daq_seq.make_mod_esr_seq,
    'rabi_seq':             daq_seq.make_rabi_seq,
    'spin_echo':            daq_seq.make_echo_seq_FL,
    'T2seq':                daq_seq.makeT2Seq0,
    'pesr_seq':             daq_seq.make_pulsed_esr_seq,
    'aom_timing':           daq_seq.make_aom_timing_seq,
    'rodelay':              daq_seq.make_opt_readout_time_sequence,
    'drift_seq':            daq_seq.make_drift_analysis_sequence,
    'MW_timing':            daq_seq.make_MW_timing_seq,
    'T1ms0':                daq_seq.make_t1_ms0_seq_MW2,
    'T1ms1':                daq_seq.make_t1_ms1_seq_FL,
    'T2_seq':               daq_seq.make_t2_seq_MW,
    'ram_seq':              daq_seq.make_ramsey_seq_MW,
    'simult_samp':          daq_seq.check_simult_sampling,
    'double_mod':           daq_seq.make_double_mod_sequence_lcm,
    'diff_mod':             daq_seq.make_diff_mod_sequence,
    'esr_dig_mod_seq':      daq_seq.make_dig_mod_odmr_sequence,
    'rabi_dig_mod_seq':     daq_seq.make_dig_mod_rabi_sequence,
    'rabi_contrast_seq':    daq_seq.make_rabi_contrast_sequence,
    't1ms0_train':          daq_seq.make_t1ms0_train_sequence,

}

_CAM_SEQUENCES = {
    'esr_seq':  cam_seq.make_esr_seq_camera,
    'rabi_seq': cam_seq.make_rabi_seq_camera_FL,
    'pesr_seq': cam_seq.make_pulsed_esr_seq_camera_FL,
    'T1ms0':    cam_seq.make_t1_seq_camera,
    'T2_seq':   cam_seq.make_t2_seq,
}

_CAM_LEVEL1_SEQUENCES = {
    'esr_seq':  cam_seq.make_esr_seq_camera_level_trigger,
    'rabi_seq': cam_seq.make_rabi_seq_camera_level_trigger,
    'pesr_seq': cam_seq.make_pulsed_esr_seq_camera_level_trigger,
    'T1ms0':    cam_seq.make_t1_seq_camera_level_trigger,
}

_DISPATCH = {
    'diode':      _DIODE_SEQUENCES,
    'cam':        _CAM_SEQUENCES,
    'cam_levelm': _CAM_SEQUENCES,
    'cam_syncm':  _CAM_SEQUENCES,
    'cam_level1': _CAM_LEVEL1_SEQUENCES,
}


# ═══════════════════════════════════════════════════════════════════════
#  make_sequence
# ═══════════════════════════════════════════════════════════════════════

def make_sequence(instr: str, sequence: str, args: list):
    """Dispatch to the correct sequence builder.

    Parameters
    ----------
    instr : str
        Instrument mode: 'diode', 'cam', 'cam_levelm', 'cam_syncm',
        'cam_level1'.
    sequence : str
        Sequence name (e.g. 'esr_seq', 'rabi_seq', 'ram_seq').
    args : list
        Positional arguments forwarded to the builder function.
        Typically: [timing_params..., pb_channels_dict].

    Returns
    -------
    (func_name, pb_channels) : tuple[str, list[PBchannel]]
        func_name is the builder's __name__ (for logging/display).
        pb_channels is whatever the builder returns — typically a list
        of PBchannel namedtuples (one set for diode, possibly nested
        for camera signal/reference halves).

    Raises
    ------
    ValueError
        If instr or sequence is not recognized.
    """
    table = _DISPATCH.get(instr)
    if table is None:
        raise ValueError(
            f"Unknown instrument type '{instr}'. "
            f"Valid: {list(_DISPATCH.keys())}")

    builder = table.get(sequence)
    if builder is None:
        raise ValueError(
            f"Unknown sequence '{sequence}' for instr='{instr}'. "
            f"Valid: {list(table.keys())}")

    return builder.__name__, builder(*args)


# ═══════════════════════════════════════════════════════════════════════
#  event_cataloguer
# ═══════════════════════════════════════════════════════════════════════

def event_cataloguer(pb_channels) -> dict[float, int]:
    """Convert PBchannel pulse definitions to a time-sorted bitmask dict.

    Each PBchannel namedtuple specifies which hardware channel turns on
    at which times and for how long.  This function merges all channels
    into a single timeline of composite bitmasks — the format needed by
    PBcontrol.create_PBinstruction().

    Parameters
    ----------
    pb_channels : list[PBchannel]
        Each PBchannel has: channel_number (bitmask), start_times (list),
        pulse_durations (list).

    Returns
    -------
    dict[float, int]
        Keys: event times in ns (sorted).
        Values: combined channel bitmask active after that event.
        Always starts with {0: 0} (all off at t=0).

    Notes
    -----
    Uses XOR (not OR) when merging overlapping edges. This correctly
    handles zero-duration pulses: a start edge followed immediately by
    an end edge at the same time cancels out, so the channel never
    actually turns on. With OR, it would turn on and only turn off at
    the next unrelated event.
    """
    # Phase 1: collect all edges into an event catalog.
    # For each event time, XOR the contributing channel mask.
    edge_catalog = {}

    for ch in pb_channels:
        mask = ch.channel_number
        if mask < 0:
            continue  # channels with negative mask are disabled

        for t_start, duration in zip(ch.start_times, ch.pulse_durations):
            t_end = t_start + duration
            for t in (t_start, t_end):
                if t in edge_catalog:
                    # XOR: toggling the same channel twice at the same
                    # time cancels it (handles zero-length pulses).
                    edge_catalog[t] ^= mask
                else:
                    edge_catalog[t] = mask

    # Phase 2: walk events in time order, accumulating the running
    # bitmask via XOR at each edge.
    bitmasks = {0: 0}
    running = 0
    for t in sorted(edge_catalog):
        running ^= edge_catalog[t]
        bitmasks[t] = running

    return bitmasks


# ═══════════════════════════════════════════════════════════════════════
#  plot_sequence
# ═══════════════════════════════════════════════════════════════════════

def plot_sequence(instructions: list, channel_map: dict):
    """Render PB instructions as a multi-channel timing diagram.

    Parameters
    ----------
    instructions : list[list]
        PB instruction list from create_PBinstruction().
        Each entry: [bitmask, inst_type, inst_data, duration_ns].
    channel_map : dict[str, int]
        Channel names → bitmask values.
        E.g. {'laser': 8, 'mw': 4, 'samp': 2, 'start': 16}.

    Returns
    -------
    (t_us, channel_traces, y_ticks) : tuple
        t_us: ndarray of time points in microseconds (step-function edges).
        channel_traces: list of lists, one per channel, with y-values
            offset by channel bit index for stacked display.
        y_ticks: ndarray of y-tick positions for labelling channels.
    """
    n_inst = len(instructions)
    if n_inst == 0:
        return np.array([0]), [], np.array([0])

    SCALE = 0.8  # vertical scaling per channel

    # Build time axis (pairs of points per instruction for step function)
    t_ns = [0.0, 0.0]
    for inst in instructions:
        prev = t_ns[-1]
        nxt = prev + inst[3]   # inst[3] = duration_ns
        t_ns.append(nxt)
        t_ns.append(nxt)
    t_us = np.array(t_ns) / 1e3

    # Build per-channel traces
    channel_traces = []
    for ch_name, ch_mask in channel_map.items():
        bit_index = int(math.log2(ch_mask))
        trace = [0.0, float(bool(ch_mask & instructions[0][0]))]

        for i in range(n_inst):
            current_on = float(bool(ch_mask & instructions[i][0]))
            if i < n_inst - 1:
                next_on = float(bool(ch_mask & instructions[i + 1][0]))
            else:
                next_on = current_on  # last instruction: hold state
            trace.append(current_on)
            trace.append(next_on)

        # Offset and scale for stacked display
        scaled = [bit_index + v * SCALE for v in trace]
        channel_traces.append(scaled)

    # Y-ticks at each channel's bit index
    bit_indices = [int(math.log2(m)) for m in channel_map.values()]
    y_ticks = np.arange(min(bit_indices), max(bit_indices) + 1, 1)

    return t_us, channel_traces, y_ticks


# ═══════════════════════════════════════════════════════════════════════
#  view_sequence  —  quick compile-and-plot utility
# ═══════════════════════════════════════════════════════════════════════

def view_sequence(instr: str, sequence: str, seq_args: list,
                  param_values=None, indices=None,
                  dpi: int = 100):
    """Compile and plot PB sequences at one or more parameter values.

    Parameters
    ----------
    instr : str
        Instrument type ('diode', 'cam', etc.).
    sequence : str
        Sequence name.
    seq_args : list
        Full sequence args list (last element = pb_channels dict).
        For PB-timing sweeps, seq_args[0] is overwritten per value.
    param_values : array-like, optional
        Sweep values to plot.  If None, plots the sequence as-is.
    indices : list[int], optional
        Which indices into param_values to plot (default: [0]).
    dpi : int
        Figure DPI.

    Notes
    -----
    For freq-only sequences (ESR, etc.), the PB program doesn't change
    with parameter value, so param_values is ignored.
    """
    freq_only = sequence in [
        'esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']

    if param_values is None:
        param_values = [seq_args[0]]
        indices = [0]
    if indices is None:
        indices = [0]

    pb_channels = seq_args[-1]  # last arg is always the channel map

    for idx in indices:
        if not freq_only:
            seq_args[0] = param_values[idx]

        _, the_list = PBcontrol.PulseBlaster.PB_program(
            instr, sequence, seq_args)

        for sub in the_list:
            inst_list = sub[0]
            inst_times = sub[4]

            fig, ax = plt.subplots(dpi=dpi)
            t_us, traces, y_ticks = plot_sequence(inst_list, pb_channels)

            for trace in traces:
                ax.plot(t_us, trace)

            ax.set_yticks(y_ticks, list(pb_channels.keys()))
            ax.set_xlabel('Time (µs)')
            ax.set_ylabel('Channel')
            ax.set_title(
                f"{sequence}  |  param[{idx}] = "
                f"{param_values[idx]:.4g} ns\n"
                f"Transitions: {list(inst_times.values())}",
                fontsize=10)


# ═══════════════════════════════════════════════════════════════════════
#  Backward-compatible class wrapper
# ═══════════════════════════════════════════════════════════════════════
#
# PBcontrol.py calls sequencecontrol.sequencecontrol.make_sequence(...)
# and sequencecontrol.sequencecontrol.sequence_event_cataloguer(...).
#
# Rather than touching PBcontrol.py, we expose a thin class that
# delegates to the module-level functions above.  The class holds no
# state and needs no __init__ args (the old one only stored
# parameter_dict for check_params, which is gone).

class sequencecontrol:
    """Thin class wrapper for PBcontrol.py compatibility.

    PBcontrol imports this as:
        import sequencecontrol
        sequencecontrol.sequencecontrol.make_sequence(...)
        sequencecontrol.sequencecontrol.sequence_event_cataloguer(...)

    All methods delegate to module-level functions.
    """

    def __init__(self, parameter_dict=None):
        # parameter_dict accepted for backward compat but unused.
        pass

    @staticmethod
    def make_sequence(instr, sequence, args):
        return make_sequence(instr, sequence, args)

    @staticmethod
    def sequence_event_cataloguer(pb_channels):
        return event_cataloguer(pb_channels)

    @staticmethod
    def plot_sequence(instructions, channel_map):
        return plot_sequence(instructions, channel_map)

    @staticmethod
    def view_sequence(instr, sequence, seq_args, only_plot=False,
                      parameter=None, seq_no_plot=None, plot_dpi=100):
        """Legacy-compatible signature that maps to view_sequence()."""
        if parameter is None:
            parameter = [0]
        if seq_no_plot is None:
            seq_no_plot = [0]
        # Map old args to new
        param_values = parameter
        indices = seq_no_plot
        view_sequence(instr, sequence, seq_args,
                      param_values=param_values,
                      indices=indices, dpi=plot_dpi)
