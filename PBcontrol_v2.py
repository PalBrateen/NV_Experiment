#PBcontrol.py

from ctypes import *
from spinapi import *
import numpy as np
import sequencecontrol as seqctrl
import connectionConfig as concfg
import sys
# the time unit is in ns.. as followed by the Pulse Blaster
# start = [0]
clk_cyc = (1/concfg.PBclk)*1e3      # in ns

def errorCatcher(statusVar):
    """
    Display the error occuring due to programming the PB board.

    Parameters
    ----------
    statusVar : int
        Value indicates the occurence of a PB error
        0 implies no error.
        -1 implies an error.
        
    Returns
    -------
    None.

    """
    if statusVar<0:
        print ('\x1b[1;37;41m'+'Hi.. Error:'+ pb_get_error())
        # pb_init();
        pb_stop(); pb_close()
        sys.exit()

def configurePB():
    """
    Initialize the PB board for programming.

    Returns
    -------
    int (0)
        Indicates the PB board is successfully initialized for programming. 0 is according to the PB manual.
        -ve indicates failure in initialization.

    """
    pb_set_debug(1)
    status = pb_init()
    errorCatcher(status)
    pb_core_clock(concfg.PBclk)
    return 0

# Define pb_inst_pbonly in a new way. It is already defined in the spinapi.py file. This method accepts inputs that are easily understandable and then converts it into the types that the function in spinapi.py file understands.
def pb_inst_pbonly(flags,inst,inst_data,length):
    """   

    Parameters
    ----------
    flags : TYPE
        DESCRIPTION.
    inst : TYPE
        DESCRIPTION.
    inst_data : TYPE
        DESCRIPTION.
    length : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    flags = c_uint(flags)
    inst = c_int(inst)
    inst_data = c_int(inst_data)
    length = c_double(length)
    return spinapi.pb_inst_pbonly(flags, inst, inst_data, length)

def stop_sequence():
    # pb_init()
    status = pb_stop();     errorCatcher(status);
    pb_close()

def PB_program(instr, sequence, sequenceArgs, err_check=False):
    if instr == 'cam':
        return PB_program_camera(sequence, sequenceArgs, err_check)
    elif instr == 'diode':
        return PB_program_diode(sequence, sequenceArgs, err_check)

def PB_program_camera(sequence, sequenceArgs, err_check=False):
    List = []
    allPBchannels = seqctrl.make_sequence('cam',sequence, sequenceArgs)       # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
    # print(allPBchannels)
    for i in range(0,2):    # 2: one for signal sequence, other for reference sequence.. duto 'allPBchannels' alada kore produce kora hochhe.. tai eta..
        channelBitMasks = seqctrl.sequence_event_cataloguer(allPBchannels[i])
        List.append(create_PBinstruction(channelBitMasks, err_check))
    # print(List)
    return List

def PB_program_diode(sequence, sequenceArgs, err_check=False):
    """ Make the Pulse Sequence instruction to be sent to the PB.
    Parameters
    ----------
    sequence : String
        Name of the sequence to be executed.
    sequenceArgs : ??
        DESCRIPTION.

    Returns
    -------
    instructionList : List
        DESCRIPTION.    """
    List = []
    allPBchannels = seqctrl.make_sequence('diode',sequence, sequenceArgs)       # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
    # print(allPBchannels)
    channelBitMasks = seqctrl.sequence_event_cataloguer(allPBchannels)
    List.append(create_PBinstruction(channelBitMasks, err_check))
    # print(List)
    return List

# Let channelBitMasks={0:0, 100:1, 250:3, 300:2, 350:34, 1000:2}
def create_PBinstruction(channelBitMasks, err_check):
    """     Creates the list of instructions sent to the PB consisting of 4 arguments for each line of instruction.
    Parameters
    ----------
    channelBitMasks : TYPE
        DESCRIPTION.

    Returns
    -------
    instructionList : List
        DESCRIPTION.    """
        
    eventTimes = list(channelBitMasks.keys())       # (List) giving the (start/end) times = [0, 100, 250, 300, 350, 1000] (6 items)
    numEvents = len(eventTimes)          # (int) No. of events (i.e. no. of times the PB has to be turned ON/OFF) = 6
    numIntervals = numEvents-1           # (int) no. of instructions = 5; actually this should be numIntervals
    eventDurations = list(np.zeros(numIntervals))    # (list) with (numEvents-1) zeros = [0,0,0,0,0]
    
    seq_error_count = 0; inst_error_no = []; error_times = []
    inst_times = {}
    for i in range(0,numIntervals):              # <<<<something wrong here.....>>>>
        #if i == numIntervals-1:
        eventDurations[i] = eventTimes[i+1] - eventTimes[i]
        inst_times[i] = eventTimes[i]/1e3
        if eventDurations[i]<10 and err_check == True:
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
    # niche instructionList ta banano hochhe...
    for i in range(0,numIntervals):
        if i == (numIntervals-1):
            # jodi list er last instruction hoye, then use Inst.BRANCH
            # instructionList.extend([[bitMasks[i], Inst.BRANCH, start[0], eventDurations[i]]]) # changed on 23062023... start[0] = 0 generally... hence the below statement...
            instructionList.extend([[bitMasks[i], Inst.BRANCH, 0, eventDurations[i]]])

        else:
            # jodi list er last instruction na hoye, then use Inst.CONTINUE
            instructionList.extend([[bitMasks[i], Inst.CONTINUE, 0, eventDurations[i]]])
        # print(str(bitMasks[i])+'\t'+str(eventDurations[i])+'\n')
    # if seq_error_count > 0:
    #     print("Total instructions #: ",numIntervals)
    return [instructionList, seq_error_count, inst_error_no, error_times, inst_times]
    
def run_sequence_for_diode(instructionList):
    """    Program and run the PB sequence for camera. Requires the PB board to be connected to the PC.

    Parameters
    ----------
    instructionList : List
        DESCRIPTION.

    Returns
    -------
    None.    """
    # configurePB()
    pb_reset()
    status = pb_start_programming(PULSE_PROGRAM);   errorCatcher(status)
    instructionList = instructionList[0]    # there is only one element in 'instructionList'.. hence assigning this as instructionList instead of signal_instruction or reference_instruction..
    started = False
    start = [0]                                    # (List) with one element; start[0]=0
    # pb_inst_pbonly(int flags, int inst, int inst_data, double time_period)
    # print(instructionList)
    for i in range(0, len(instructionList)):
        if started:
            status = pb_inst_pbonly(instructionList[i][0],instructionList[i][1],instructionList[i][2],instructionList[i][3]);
            # print(instructionList[i][0],instructionList[i][1],instructionList[i][2],instructionList[i][3])
            errorCatcher(status)
        else:
            start[0] = pb_inst_pbonly(instructionList[0][0],instructionList[0][1],instructionList[0][2],instructionList[0][3])
            # print(instructionList[0][0],instructionList[0][1],instructionList[0][2],instructionList[0][3])
            errorCatcher(start[0])
            started = True
    status = pb_stop_programming();     errorCatcher(status)
    status = pb_start();    errorCatcher(status)
    # status = pb_close();                errorCatcher(status)

def run_sequence_for_camera(instructionList,t_exposure,t_seq_total,N_total):
    # t_seq_total should be a list: t_seq_total = [t_seq_tot_sig, t_seq_tot_ref]
    # print(instructionList)
    # instructionList               => len(instructionList)=2; contains PB instructions for 'signal' and 'reference' sequences
    # instructionList[i]            => i=0: signal, i=1: reference
    # instructionList[0][i]         => i=0: 1st instruction of the signal part, i=1: 2nd instruction...
    # isntructionList[0][0][i]   => i=0: PBchannel number, i=1: Inst, i=2: Inst data, i=3: time (delay)

    # configurePB()
    pb_reset()
    status = pb_start_programming(PULSE_PROGRAM);    errorCatcher(status)

    start = [0]                                    # (List) with one element; start[0]=0
    #-----------------------------korlam.. 15/07/2023-----------------------------
    t_cam_response = 38.96 *us
    
    # below are two lists, with 2 elements each, one for signal and other for reference...
    # the operation should use np.floor() and not np.rint().. think..
    if N_total == []:
        N_trigger = [int(np.floor(t_cam_response/element)) for element in t_seq_total]
        N_remaining = [int(np.floor((t_exposure - t_cam_response)/element)) for element in t_seq_total]
        N_total = [int(np.floor(t_exposure/element)) for element in t_seq_total]
        # need to check if (N_trigger + N_remaining) <= N_total) ??

        # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
        t_buffer0 = (t_exposure - t_cam_response) - N_remaining[0]*t_seq_total[0]
        t_buffer1 = t_cam_response - N_trigger[0]*t_seq_total[0]
        t_buffer2 = (t_exposure - t_cam_response) - N_remaining[1]*t_seq_total[1]
        t_buffer3 = t_cam_response - N_trigger[1]*t_seq_total[1]
        t_buffer = [t_buffer0, t_buffer1, t_buffer2, t_buffer3]
        t_buffer = [clk_cyc*round(time/clk_cyc) if time>10*ns else 10.0*ns for time in t_buffer]
        
    else:
        N_trigger = [int(np.floor(t_cam_response/element)) for element in t_seq_total]
        N_remaining = [int(N_total[i] - N_trigger[i]) for i in range(0,len(N_total))]
        t_exposure = [(t_cam_response+N_remaining[i]*t_seq_total[i]) for i in range(0,len(N_total))]    # this is in ns..
        # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
        t_buffer1 = t_cam_response - N_trigger[0]*t_seq_total[0]
        t_buffer3 = t_cam_response - N_trigger[1]*t_seq_total[1]
        t_buffer0 = 0
        t_buffer2 = 0

        t_buffer1 = t_buffer1 if t_buffer1>=10*ns else 10*ns
        t_buffer3 = t_buffer3 if t_buffer3>=10*ns else 10*ns

        t_buffer = [t_buffer0, clk_cyc*round(t_buffer1/clk_cyc), t_buffer2, clk_cyc*round(t_buffer3/clk_cyc)]

    N = [N_trigger, N_remaining, N_total]

    # print(N)
    # print(t_buffer)
    # print(t_exposure)

    #-----------------------------korlam.. 15/07/2023-----------------------------

    signal_instruction = instructionList[0]
    # print(signal_instruction)
    reference_instruction = instructionList[1]
    # print(reference_instruction)
    
    # follow notebook for the details of the sequence...
    #------------------------------------------------------------------------------------
    if N_trigger[1] >= 1:
        started = False
        # 1st: sequence starts with camera HIGH and the 'reference' sequence (instructionList[1]) running (assume, as a part of previous cycle)... so the commands include a XOR with concfg.camera...
        for i in range(0, len(reference_instruction)-1):
            # print(reference_instruction[i])
            if started:
                # the instructions other than the first and last contain Inst.CONTINUE from their source...
                status = pb_inst_pbonly(reference_instruction[i][0]^concfg.camera, reference_instruction[i][1],reference_instruction[i][2], reference_instruction[i][3])
                # print(reference_instruction[i][0]^concfg.camera, reference_instruction[i][1],reference_instruction[i][2], reference_instruction[i][3])
                errorCatcher(status)
            else:
                # the first instruction contains Inst.LOOP with repetitions N_trigger(_reference)[1], bypassing what is contained in instructionList[0][1] (Inst) & instructionList[0][2] (Inst_data)
                start[0] = pb_inst_pbonly(reference_instruction[0][0]^concfg.camera, Inst.LOOP, N_trigger[1],reference_instruction[0][3])
                # print(reference_instruction[0][0]^concfg.camera, Inst.LOOP, N_trigger[1],reference_instruction[0][3])
                errorCatcher(start[0])
                started = True
        # the last instruction.. the stack of instructions needs to be repeated.. put Inst.END_LOOP with the start point (start[0])
        status = pb_inst_pbonly(reference_instruction[i+1][0]^concfg.camera, Inst.END_LOOP, start[0], reference_instruction[i+1][3]);
        # print(reference_instruction[i+1][0]^concfg.camera, Inst.END_LOOP, start[0], reference_instruction[i+1][3]);
        errorCatcher(status)
        # wait for the buffer time (t_buffer[3])... use Inst.CONTINUE as there are more instructions...
        status = pb_inst_pbonly(0^concfg.camera, Inst.CONTINUE, 0, t_buffer[3]) if t_buffer[3]>=10 else 0
        # print(0, Inst.CONTINUE, 0, t_buffer[3])
        errorCatcher(status)
    else:
        start[0] = pb_inst_pbonly(concfg.camera, Inst.CONTINUE, 0, t_cam_response); errorCatcher(start[0])
    #------------------------------------------------------------------------------------
    # 2nd: sequence has 'no' camera, with the sequence (instructionList[0]) running...
    started = False
    # the buffer time is at the beginning in this case...
    status = pb_inst_pbonly(0, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0]>=10 else 0
    # print(0, Inst.CONTINUE, 0, t_buffer[0])
    errorCatcher(status)
    for i in range(0, len(signal_instruction)-1):
        # print(signal_instruction[i])
        if started:
            # the instructions other than the first and the last...
            status = pb_inst_pbonly(signal_instruction[i][0], signal_instruction[i][1], signal_instruction[i][2],signal_instruction[i][3])
            # print(signal_instruction[i][0], signal_instruction[i][1], signal_instruction[i][2],signal_instruction[i][3])
            errorCatcher(status)
        else:
            # first instruction with Inst.LOOP and repetitions N_remaining(_signal)[0]
            start.append(pb_inst_pbonly(signal_instruction[0][0], Inst.LOOP, N_remaining[0], signal_instruction[0][3]))
            # print(signal_instruction[0][0], Inst.LOOP, N_remaining[0], signal_instruction[0][3])
            errorCatcher(start[1])
            started = True
    # the last instruction.. with Inst.END_LOOP and start point (start[1])...
    status = pb_inst_pbonly(signal_instruction[i+1][0], Inst.END_LOOP, start[1], signal_instruction[i+1][3]);
    # print(signal_instruction[i+1][0], Inst.END_LOOP, start[1], signal_instruction[i+1][3])
    errorCatcher(status)
    #------------------------------------------------------------------------------------
    # 3rd: sequence has 'camera' has again (XORed with concfg.camera)... with the signal sequence (instructionList[0]) running...
    if N_trigger[0] >= 1:
        started = False
        for i in range(0, len(signal_instruction)-1):
            # print(signal_instruction[i])
            if started:
                # middle instructions...
                status = pb_inst_pbonly(signal_instruction[i][0]^concfg.camera, signal_instruction[i][1], signal_instruction[i][2], signal_instruction[i][3])
                # print(signal_instruction[i][0]^concfg.camera, signal_instruction[i][1], signal_instruction[i][2], signal_instruction[i][3])
                errorCatcher(status)
            else:
                # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                start.append(pb_inst_pbonly(signal_instruction[0][0]^concfg.camera, Inst.LOOP, N_trigger[0], signal_instruction[0][3]))
                # print(signal_instruction[0][0]^concfg.camera, Inst.LOOP, N_trigger[0], signal_instruction[0][3])
                errorCatcher(start[2])
                started = True
                
        # the last instruction.. with Inst.END_LOOP and start point (start[2])...
        status = pb_inst_pbonly(signal_instruction[i+1][0]^concfg.camera, Inst.END_LOOP, start[2], signal_instruction[i+1][3]); 
        # print(signal_instruction[i+1][0]^concfg.camera, Inst.END_LOOP, start[2], signal_instruction[i+1][3])
        errorCatcher(status)
        # the buffer comes here.. 
        status = pb_inst_pbonly(0^concfg.camera, Inst.CONTINUE, 0, t_buffer[1]) if t_buffer[1]>=10 else 0
        # print(0^concfg.camera, Inst.CONTINUE, 0, t_buffer[1])
        errorCatcher(status)
    else:
        start.append(pb_inst_pbonly(concfg.camera, Inst.CONTINUE, 0, t_cam_response)); errorCatcher(start[2])
    #------------------------------------------------------------------------------------
    # 4th: sequence runs w/o 'camera'... with the 'reference' sequence (instructionList[1])...
    started = False
    # buffer comes here...
    status = pb_inst_pbonly(0, Inst.CONTINUE, 0, t_buffer[2]) if t_buffer[2]>=10 else 0
    # print(0, Inst.CONTINUE, 0, t_buffer[2])
    errorCatcher(status)
    for i in range(0, len(reference_instruction)-1):
        # print(reference_instruction[i])
        if started:
            # instructions in the middle..
            status = pb_inst_pbonly(reference_instruction[i][0], reference_instruction[i][1], reference_instruction[i][2],reference_instruction[i][3])
            # print(reference_instruction[i][0], reference_instruction[i][1], reference_instruction[i][2],reference_instruction[i][3])
            errorCatcher(status)
        else:
            # first instruction... Inst.LOOP and repetition N_remaining(_reference)[1]...
            start.append(pb_inst_pbonly(reference_instruction[0][0], Inst.LOOP, N_remaining[1], reference_instruction[0][3]))
            # print(reference_instruction[0][0], Inst.LOOP, N_remaining[1], reference_instruction[0][3])
            errorCatcher(start[3])
            started = True
    # last instruction being written.. with Inst.END_LOOP and start point (start[3])...
    status = pb_inst_pbonly(reference_instruction[i+1][0], Inst.END_LOOP, start[3], reference_instruction[i+1][3]);
    # print(reference_instruction[i+1][0], Inst.END_LOOP, start[3], reference_instruction[i+1][3])
    errorCatcher(status)
    # a last instruction with Inst.BRANCH of length 10 ns... returns the control to start[0] (the first instruction)
    status = pb_inst_pbonly(reference_instruction[i+1][0], Inst.BRANCH, start[0], 10*ns);
    # print(reference_instruction[i+1][0], Inst.BRANCH, start[0], 10*ns)
    errorCatcher(status)
    # 10 ns is very small compared to the time-scales involed (response time and esp. exposure time) in the camera program...
    #---------------------------------------------DONE-------------------------------------------

    status = pb_stop_programming();     errorCatcher(status)
    status = pb_start();                errorCatcher(status)
    # status = pb_close();                errorCatcher(status)
    # print('\x1b[38;2;250;0;0m\x10 PB STARTED\x1b[0m')
    # return instructionList
    return [t_exposure, N]
        


    # # ekhan theke 'signal' frame er part ta shuru hochhe...
    # start_trigger = pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.CONTINUE, 0, 38.96*us)     # camera trigger with pulse width of 38.96 us
    # for i in range(0, len(instructionList[0])-1):
    #     if started:
    #         status = pb_inst_pbonly(instructionList[0][i][0],instructionList[0][i][1],instructionList[0][i][2],instructionList[0][i][3])
    #         errorCatcher(status)
    #     else:
    #         start[0] = pb_inst_pbonly(instructionList[0][0][0], Inst.LOOP, N_repeat,instructionList[0][0][3])    # start korar somoy mention N_repeat bar loop ta cholbe...
    #         errorCatcher(start[0])
    #         started = True
    # loop_end = pb_inst_pbonly(instructionList[0][i+1][0], Inst.END_LOOP, start[0], instructionList[0][i+1][3])
    # # ebar final trigger -- signal frame ta sesh kore reference frame ta shuru hobe...
    # # sathe sathe laser ta o ON kore dite hobe...
    # final_trigger = pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.CONTINUE, 0, 38.96*us)
    # # 'signal' frame er part ta sesh holo...

    # # ebar 'reference' frame er part ta shuru hochhe.......... ki kore hobe ei part ta???
    # started = False
    # for i in range(0, len(instructionList[1])-1):
    #     if started:
    #         status = pb_inst_pbonly(instructionList[1][i][0],instructionList[1][i][1],instructionList[1][i][2],instructionList[1][i][3])
    #         errorCatcher(status)
    #     else:
    #         start.append(pb_inst_pbonly(instructionList[1][0][0], Inst.LOOP, N_repeat,instructionList[1][0][3]))    # start korar somoy mention N_repeat bar loop ta cholbe...
    #         errorCatcher(start[0])
    #         started = True
    # loop_end = pb_inst_pbonly(instructionList[1][i+1][0], Inst.END_LOOP, start[1], instructionList[1][i+1][3])
    # final_trigger = pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.BRANCH, start[0], 1*us)
    # # 'reference΄' frame er part ta sesh holo...
    # # ebar sequence end kora hobe.. 