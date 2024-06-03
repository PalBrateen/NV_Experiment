# sequenceControl.py
from connectionConfig import AOM, DAQ, STARTtrig, I, Q, uW, PBclk
from spinapi import ns, us, ms
import sys, math
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from importlib import import_module
expCfgFile = 'esr_config'; expCfg = import_module(expCfgFile)
#%% Operations
""" Create a tuple of 'PBchannel' type with entries:
    1. channelNumber
    2. startTimes
    3. pulseDurations
"""
PBchannel = namedtuple('PBchannel', ['channelNumber', 'startTimes', 'pulseDurations'])
clk_cyc = 1e3/PBclk       # Time resolution in ns
# Short pulse flags: Switch ON bits 21-23 of the Control word to enable short pulse feature
ONE_PERIOD = 0x200000           # 23/22/21/20 = 0010b = 2^21d = 200000x
TWO_PERIOD = 0x400000           # 23/22/21/20 = 0100 = 2^22
THREE_PERIOD = 0x600000         # 23/22/21/20 = 0110
FOUR_PERIOD = 0x800000          # 23/22/21/20 = 1000
FIVE_PERIOD = 0xA00000           # 23/22/21/20 = 1010

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
    
# plot_sequence():
# accepts:
#    'instructionList' - (List of lists) created by sequenceControl.sequenceEventCataloguer()
#    'channelMasks' - ()
# returns: [t_us,channelPulses,yTicks]
def plot_sequence(instructions, channelMasks):
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
    # Catalogs sequence events in terms of consecutive rising edges on the allPBchannels provided. Returns a dictionary, channelBitMasks, whose keys are event (rising/falling edge) times and values are the channelBitMask which indicate which allPBchannels are on at that time.
    eventCatalog = {}  # dictionary where the keys are rising/falling edge times and the values are the channel bit masks which turn on/off at that time

    # (PBchannel) a particular entry in the 'allPBchannels' list
    for aPBchannel in allPBchannels:
        channelMask = aPBchannel.channelNumber
        endTimes = [startTime + pulseDuration for startTime, pulseDuration in zip(aPBchannel.startTimes, aPBchannel.pulseDurations)]
        eventTimes = aPBchannel.startTimes+endTimes
        # print(eventTimes)
        for eventTime in eventTimes:
            # eventTime (int) gives the start and end times of a particular component (AOM/DAQ/uW)
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

# ------PulseBlaster Sequences------
"""Variables:
    * sequence: (String) Name of the sequence to be executed
    * args: (List) of time durations
"""
def make_sequence(sequence, args):
    if sequence == 'esr_seq':
        return makeESRseq(*args)
    elif sequence == 'modesr':
        return make_mod_esr_seq(*args)
    elif sequence == 'rabi_seq':
        return make_rabi_seq(*args)
    ...
    ...
    ...
    
#PBcontrol.py
from ctypes import *
from spinapi import *
import numpy as np
import sequenceControl as seqctrl
from connectionConfig import *
import sys

start = [0]
def errorCatcher(statusVar):
    if statusVar<0:
        print ('\x1b[1;37;41m'+'Error:'+ pb_get_error())
        sys.exit()
def configurePB():
    pb_set_debug(1)
    status = pb_init()
    errorCatcher(status)
    pb_core_clock(PBclk)
    return 0
# Define pb_inst_pbonly in a new way. It is already defined in the spinapi.py file. This method accepts inputs that are easily understandable and then converts it into the types that the function in spinapi.py file understands.
def pb_inst_pbonly(flags,inst,inst_data,length):
    flags = c_uint(flags)
    inst = c_int(inst)
    inst_data = c_int(inst_data)
    length = c_double(length)
    return spinapi.pb_inst_pbonly(flags, inst, inst_data, length)

def PB_program(sequence,sequenceArgs):
    allPBchannels=seqctrl.make_sequence(sequence, sequenceArgs)       # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
    channelBitMasks = seqctrl.sequence_event_cataloguer(allPBchannels)
    instructionList = create_PBinstruction(channelBitMasks)
    return instructionList

# Let channelBitMasks={0:0, 100:1, 250:3, 300:2, 350:34, 1000:2}
def create_PBinstruction(channelBitMasks):
    eventTimes = list(channelBitMasks.keys())       # (List) giving the (start/end) times = [0, 100, 250, 300, 350, 1000] (6 items)
    numEvents = len(eventTimes)          # (int) No. of events (i.e. no. of times the PB has to be turned ON/OFF) = 6
    numIntervals = numEvents-1           # (int) no. of instructions = 5; actually this should be numIntervals
    eventDurations = list(np.zeros(numIntervals))    # (list) with (numEvents-1) zeros = [0,0,0,0,0]
    seq_error_count = 0
    inst_error_no = []
    error_times = []
    inst_times = {}
    for i in range(0,numIntervals):              # <<<<something wrong here.....>>>>
        #if i == numIntervals-1:
        eventDurations[i] = eventTimes[i+1] - eventTimes[i]
        inst_times[i] = eventTimes[i]/1e3
        if eventDurations[i]<10:
            inst_error_no.extend([i])
            error_times.extend([eventTimes[i]])
            # print("Instruction at "+str(i)+" = "+'\x1b[1;37;41m'+str(eventDurations[i])+" < "+str(10)+"ns.")
            seq_error_count += 1
            # eventDurations[i] = 10
            # sys.exit()
        # else:
        eventDurations[i] = eventTimes[i+1]-eventTimes[i]
    instructionList = []
    bitMasks = list(channelBitMasks.values())      # (List) = [0, 1, 3, 2, 34, 2]
    # Put start = [0] in the main program so that it is not limited to this function namespace.. Changed...
    # print('bitMasks \t Event Duration \n')
    # print('---------\t-------\n')
    for i in range(0,numIntervals):
        if i == (numIntervals-1):
            instructionList.extend([[bitMasks[i], Inst.BRANCH, start[0], eventDurations[i]]])
        else:
            instructionList.extend([[bitMasks[i], Inst.CONTINUE, 0, eventDurations[i]]])
        # print(str(bitMasks[i])+'\t'+str(eventDurations[i])+'\n')
    # if seq_error_count > 0:
    #     print("Total instructions #: ",numIntervals)
    return [instructionList, seq_error_count, inst_error_no, error_times, inst_times]
    
def run_sequence(instructionList):
    configurePB()
    status = pb_start_programming(PULSE_PROGRAM)
    errorCatcher(status)
    startDone = False
    # start = [0]                         # (List) with one element; start[0]=0
    # pb_inst_pbonly(int flags, int inst, int inst_data, double time_period)
    for i in range(0, len(instructionList)):
        if startDone:
            status = pb_inst_pbonly(instructionList[i][0],instructionList[i][1],instructionList[i][2],instructionList[i][3])
            errorCatcher(status)
        else:
            start[0]= pb_inst_pbonly(instructionList[0][0],instructionList[0][1],instructionList[0][2],instructionList[0][3])
            errorCatcher(start[0])
            startDone = True
    status = pb_stop_programming()
    errorCatcher(status)
    status = pb_start()
    errorCatcher(status)
    status = pb_close()
    errorCatcher(status)
    # print('\x1b[38;2;250;0;0m\x10 PB STARTED\x1b[0m')
    # return instructionList
