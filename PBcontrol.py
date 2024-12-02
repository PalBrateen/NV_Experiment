#PBcontrol.py
#%%
from ctypes import *
from spinapi import *
import numpy as np, sequencecontrol as seqctrl, connectionConfig as concfg, sys, time

class PulseBlaster():
    
    def __init__(self) -> None:
        # the time unit is in ns.. as followed by the Pulse Blaster
        self.pbclk = concfg.PBclk
        self.clk_cyc = (1 / self.pbclk) * 1e3  # in ns
        self.pb_status = "closed"
        # self.configure()

        # the state of the PulseBlaster is as follows, state can be enquired after initialization
        # error         - if error occurs, replace the sys.exit() with this
        # init          - after pb_init()
        # closed        - after pb_close()
        # configured    - after pb_core_clock()
        # programmed    - after pb_start_programming() pb_stop_programming()
        # running       - after pb_start()
        # stopped       - after pb_stop()

    def errorCatcher(self, status):
        """
        Display the error occuring due to programming the PB board.

        Parameters
        ----------
        statusVar : int
            Value indicates the occurence of a PB error.
            0 implies no error.
            -1 implies an error.
            
        Returns
        -------
        None.

        """
        if status < 0:
            print('\x1b[1;37;41m' + 'Error:' + pb_get_error() + '\x1b[0m')
            # pb_init();
            pb_stop()
            pb_close()
            # sys.exit()
            self.pb_status = "error"
            sys.exit()
        else:
            return 0


    def configure(self):
        """
        Initialize the PB board for programming.

        Returns
        -------
        int (0)
            Indicates the PB board is successfully initialized for programming. 0 is according to the PB manual.
            -ve indicates failure in initialization.

        """
        pb_set_debug(1)
        # status = pb_init()
        # while status < 0:
        #     print(pb_get_error() + '\x1b[0m')
        #     time.sleep(0.5)
        #     pb_close()
        #     status = pb_init()
        if self.pb_status == "closed":
            status = pb_init(); self.errorCatcher(status)
            self.pb_status = "init"
            # print("PB Init'd..")
            print('\x10 PB: \x1b[38;2;250;250;0mv' + pb_get_version() + '\x1b[0m')
        else:
            print("Already init'd")
            status = 0
        
        pb_core_clock(self.pbclk)
        self.pb_status = "configured"
        return self.pb_status

    def closePB(self):
        status = pb_close()
        self.errorCatcher(status)
        self.pb_status = "closed"

    # Define pb_inst_pbonly in a new way. It is already defined in the spinapi.py file. This method accepts inputs that are easily understandable and then converts it into the types that the function in spinapi.py file understands.
    def pb_inst_pbonly(self, flags, inst, inst_data, length):
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

    def start_sequence(self):
        # pb_init()
        # if self.pb_status
        status = pb_start()
        self.errorCatcher(status)
        self.pb_status = "running"
        # pb_close()
        # self.pb_status = "closed"
        return self.pb_status

    def stop_sequence(self):
        # pb_init()
        # if self.pb_status
        status = pb_stop()
        self.errorCatcher(status)
        self.pb_status = "stopped"
        # pb_close()
        # self.pb_status = "closed"
        return self.pb_status


    def PB_program(self, instr, sequence, sequenceArgs, err_check=False):
        self.seqctrl = seqctrl.sequencecontrol()
        if instr == 'cam':
            return self.PB_program_camera(sequence, sequenceArgs, err_check)
        elif instr == 'cam_level1':
            return self.PB_program_camera_level_trigger_1(sequence, sequenceArgs, err_check)
        elif 'levelm' in instr:
            return self.PB_program_camera_level_trigger_many(sequence, sequenceArgs, err_check)
        elif instr == 'diode':
            return self.PB_program_diode(sequence, sequenceArgs, err_check)
        elif 'timeseries' in instr:
            return self.PB_program_diode(sequence, sequenceArgs, err_check)      # dummy return statement, required for program flow


    def PB_program_camera(self, sequence, sequenceArgs, err_check=False):
        List = []
        allPBchannels = self.seqctrl.make_sequence('cam', sequence, sequenceArgs)
        # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
        # print(allPBchannels)
        for i in range(0, len(allPBchannels)):
            # 2: one for signal sequence, other for reference sequence.. duto 'allPBchannels' alada kore produce kora hochhe.. tai eta..
            channelBitMasks = self.seqctrl.sequence_event_cataloguer(allPBchannels[i])
            List.append(self.create_PBinstruction(channelBitMasks, err_check))
        # print(List)
        return List


    def PB_program_camera_level_trigger_many(self, sequence, sequenceArgs, err_check=False):
        List = []
        allPBchannels = self.seqctrl.make_sequence('cam_levelm', sequence, sequenceArgs)
        # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
        # print(allPBchannels)
        for i in range(0, len(allPBchannels)):
            # 2: one for signal sequence, other for reference sequence.. duto 'allPBchannels' alada kore produce kora hochhe.. tai eta..
            channelBitMasks = self.seqctrl.sequence_event_cataloguer(allPBchannels[i])
            List.append(self.create_PBinstruction(channelBitMasks, err_check))
        # print(List)
        return List


    def PB_program_camera_level_trigger_1(self, sequence, sequenceArgs, err_check=False):
        List = []
        allPBchannels = self.seqctrl.make_sequence('cam_level1', sequence, sequenceArgs)
        # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
        # print(allPBchannels)
        channelBitMasks = self.seqctrl.sequence_event_cataloguer(allPBchannels)
        List.append(self.create_PBinstruction(channelBitMasks, err_check))
        # print(List)
        return List


    def PB_program_diode(self, sequence, sequenceArgs, err_check=False):
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
        allPBchannels = self.seqctrl.make_sequence('diode', sequence, sequenceArgs)
        # (List) of (PBchannels) containing information on which PB channel to turn ON at what time and for what duration.
        # print(allPBchannels)
        channelBitMasks = self.seqctrl.sequence_event_cataloguer(allPBchannels)
        List.append(self.create_PBinstruction(channelBitMasks, err_check))
        # print(List)
        return List


    # Let channelBitMasks={0:0, 100:1, 250:3, 300:2, 350:34, 1000:2}
    def create_PBinstruction(self, channelBitMasks, err_check):
        """     Creates the list of instructions sent to the PB consisting of 4 arguments for each line of instruction.
        Parameters
        ----------
        channelBitMasks : TYPE
            DESCRIPTION.

        Returns
        -------
        instructionList : List
            DESCRIPTION.    """

        eventTimes = list(channelBitMasks.keys())  # (List) giving the (start/end) times = [0, 100, 250, 300, 350, 1000] (6 items)
        numEvents = len(eventTimes)  # (int) No. of events (i.e. no. of times the PB has to be turned ON/OFF) = 6
        numIntervals = numEvents - 1  # (int) no. of instructions = 5; actually this should be numIntervals
        eventDurations = list(np.zeros(numIntervals))  # (list) with (numEvents-1) zeros = [0,0,0,0,0]

        seq_error_count = 0;
        inst_error_no = [];
        error_times = []
        inst_times = {}
        for i in range(0, numIntervals):  # <<<<something wrong here.....>>>>
            #if i == numIntervals-1:
            eventDurations[i] = eventTimes[i + 1] - eventTimes[i]
            inst_times[i] = eventTimes[i] / 1e3
            if eventDurations[i] < 10 and err_check == True:
                inst_error_no.extend([i])
                error_times.extend([eventTimes[i]])
                # print("Instruction at "+str(i)+" = "+'\x1b[1;37;41m'+str(eventDurations[i])+" < "+str(10)+"ns.")
                seq_error_count += 1
                # eventDurations[i] = 10
                # sys.exit()
            # else:
            eventDurations[i] = eventTimes[i + 1] - eventTimes[i]
        instructionList = []
        bitMasks = list(channelBitMasks.values())  # (List) = [0, 1, 3, 2, 34, 2]
        # Put start = [0] in the main program so that it is not limited to this function namespace.. Changed...
        # print('bitMasks \t Event Duration \n')
        # print('---------\t-------\n')
        # niche instructionList ta banano hochhe...
        for i in range(0, numIntervals):
            if i == (numIntervals - 1):
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


    # def run_sequence(instr, instructionList,t_exposure,t_seq_total,N_total):
    # def run_sequence(instr, instructionList,*sequenceparams): # --------better implementation---------
    #     if instr == 'cam':
    #         return run_sequence_for_camera(instructionList,t_exposure,t_seq_total,N_total)
    #     elif instr == 'cam_level':
    #         return run_sequence_for_camera_level_trigger(instructionList)
    #     elif instr == 'diode':
    #         return self.run_sequence_for_diode(instructionList)

    def run_sequence_for_diode(self, instructionList):
        """    Program and run the PB sequence for camera. Requires the PB board to be connected to the PC.

        Parameters
        ----------
        instructionList : List
            DESCRIPTION.

        Returns
        -------
        None.    """
        # configure()
        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);
        self.errorCatcher(status)
        # there is only one element in 'instructionList'.. hence assigning this as instructionList instead of signal_instruction or reference_instruction..
        instructionList = instructionList[0]
        started = False
        start = [0]  # (List) with one element; start[0]=0
        # pb_inst_pbonly(int flags, int inst, int inst_data, double time_period)
        # print(instructionList)
        for i in range(0, len(instructionList)):
            if started:
                status = pb_inst_pbonly(instructionList[i][0], instructionList[i][1], instructionList[i][2], instructionList[i][3]);
                # print(instructionList[i][0],instructionList[i][1],instructionList[i][2],instructionList[i][3])
                self.errorCatcher(status)
            else:
                start[0] = pb_inst_pbonly(instructionList[0][0], instructionList[0][1], instructionList[0][2], instructionList[0][3])
                # print(instructionList[0][0],instructionList[0][1],instructionList[0][2],instructionList[0][3])
                self.errorCatcher(start[0])
                started = True
        status = pb_stop_programming();
        self.errorCatcher(status)
        status = pb_start();
        self.errorCatcher(status)
        self.pb_status = "running"
        # status = pb_close();                self.errorCatcher(status)


    def run_sequence_for_camera_level_trigger_1(self, instructionList):
        self.run_sequence_for_diode(instructionList)


    def run_sequence_for_camera(self, instructionList, t_exposure, t_seq_total, N_total):
        # t_seq_total should be a list: t_seq_total = [t_seq_tot_sig, t_seq_tot_ref]
        # print(instructionList)
        # instructionList               => len(instructionList)=2; contains PB instructions for 'signal' and 'reference' sequences
        # instructionList[i]            => i=0: signal, i=1: reference
        # instructionList[0][i]         => i=0: 1st instruction of the signal part, i=1: 2nd instruction...
        # isntructionList[0][0][i]   => i=0: PBchannel number, i=1: Inst, i=2: Inst data, i=3: time (delay)

        # configure()

        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);
        self.errorCatcher(status)

        start = [0]  # (List) with one element; start[0]=0
        #-----------------------------korlam.. 15/07/2023-----------------------------
        t_cam_response = 38.96 * us
        t_jitter = 0 * us
        # below are two lists, with 2 elements each, one for signal and other for reference...
        # the operation should use np.floor() and not np.rint().. think..
        if N_total == []:
            N_trigger = [int(np.floor((t_cam_response - t_jitter) / element)) for element in t_seq_total]
            N_remaining = [int(np.floor((t_exposure - t_cam_response - t_jitter) / element)) for element in t_seq_total]
            N_total = [int(np.floor(t_exposure / element)) for element in t_seq_total]
            # print(N_total)
            # need to check if (N_trigger + N_remaining) <= N_total) ??

            # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period

            t_buffer0 = (t_exposure - t_cam_response) - N_remaining[0] * t_seq_total[0]
            t_buffer1 = t_cam_response - N_trigger[0] * t_seq_total[0]
            t_buffer2 = (t_exposure - t_cam_response) - N_remaining[1] * t_seq_total[1]
            t_buffer3 = t_cam_response - N_trigger[1] * t_seq_total[1]
            t_buffer = [t_buffer0, t_buffer1, t_buffer2, t_buffer3]
            t_buffer = [self.clk_cyc * round(time / self.clk_cyc) if time > 10 * ns else 10.0 * ns for time in t_buffer]

        else:
            # change t_exposure keeping N_total constant
            N_trigger = [int(np.floor(t_cam_response / element)) for element in t_seq_total]
            N_remaining = [int(N_total[i] - N_trigger[i]) for i in range(0, len(N_total))]
            t_exposure = [(t_cam_response + N_remaining[i] * t_seq_total[i]) for i in
                        range(0, len(N_total))]  # this is in ns..
            # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
            t_buffer1 = t_cam_response - N_trigger[0] * t_seq_total[0]
            t_buffer3 = t_cam_response - N_trigger[1] * t_seq_total[1]
            t_buffer0 = 0
            t_buffer2 = 0

            t_buffer1 = t_buffer1 if t_buffer1 >= 10 * ns else 10 * ns
            t_buffer3 = t_buffer3 if t_buffer3 >= 10 * ns else 10 * ns

            t_buffer = [t_buffer0, self.clk_cyc * round(t_buffer1 / self.clk_cyc), t_buffer2, self.clk_cyc * round(t_buffer3 / self.clk_cyc)]

        N = [N_trigger, N_remaining, N_total]

        print(N)
        print(t_buffer)
        print(t_exposure)

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
            for i in range(0, len(reference_instruction) - 1):
                # print(reference_instruction[i])
                if started:
                    # the instructions other than the first and last contain Inst.CONTINUE from their source...
                    status = pb_inst_pbonly(reference_instruction[i][0] ^ concfg.camera, reference_instruction[i][1],
                                            reference_instruction[i][2], reference_instruction[i][3])
                    # print(reference_instruction[i][0]^concfg.camera, reference_instruction[i][1],reference_instruction[i][2], reference_instruction[i][3])
                    self.errorCatcher(status)
                else:
                    # the first instruction contains Inst.LOOP with repetitions N_trigger(_reference)[1], bypassing what is contained in instructionList[0][1] (Inst) & instructionList[0][2] (Inst_data)
                    start[0] = pb_inst_pbonly(reference_instruction[0][0] ^ concfg.camera, Inst.LOOP, N_trigger[1],
                                            reference_instruction[0][3])
                    # print(reference_instruction[0][0]^concfg.camera, Inst.LOOP, N_trigger[1],reference_instruction[0][3])
                    self.errorCatcher(start[0])
                    started = True
            # the last instruction.. the stack of instructions needs to be repeated.. put Inst.END_LOOP with the start point (start[0])
            status = pb_inst_pbonly(reference_instruction[i + 1][0] ^ concfg.camera, Inst.END_LOOP, start[0],
                                    reference_instruction[i + 1][3]);
            # print(reference_instruction[i+1][0]^concfg.camera, Inst.END_LOOP, start[0], reference_instruction[i+1][3]);
            self.errorCatcher(status)
            # wait for the buffer time (t_buffer[3])... use Inst.CONTINUE as there are more instructions...
            status = pb_inst_pbonly(0 ^ concfg.camera, Inst.CONTINUE, 0, t_buffer[3]) if t_buffer[3] >= 10 else 0
            # print(0, Inst.CONTINUE, 0, t_buffer[3])
            self.errorCatcher(status)
        else:
            start[0] = pb_inst_pbonly(concfg.camera, Inst.CONTINUE, 0, t_cam_response);
            self.errorCatcher(start[0])
        #------------------------------------------------------------------------------------
        # 2nd: sequence has 'no' camera, with the sequence (signal_instruction) running...
        started = False
        # the buffer time is at the beginning in this case...
        status = pb_inst_pbonly(0, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else 0
        # print(0, Inst.CONTINUE, 0, t_buffer[0])
        self.errorCatcher(status)
        for i in range(0, len(signal_instruction) - 1):
            # print(signal_instruction[i])
            if started:
                # the instructions other than the first and the last...
                status = pb_inst_pbonly(signal_instruction[i][0], signal_instruction[i][1], signal_instruction[i][2],
                                        signal_instruction[i][3])
                # print(signal_instruction[i][0], signal_instruction[i][1], signal_instruction[i][2],signal_instruction[i][3])
                self.errorCatcher(status)
            else:
                # first instruction with Inst.LOOP and repetitions N_remaining(_signal)[0]
                start.append(pb_inst_pbonly(signal_instruction[0][0], Inst.LOOP, N_remaining[0], signal_instruction[0][3]))
                # print(signal_instruction[0][0], Inst.LOOP, N_remaining[0], signal_instruction[0][3])
                self.errorCatcher(start[1])
                started = True
        # the last instruction.. with Inst.END_LOOP and start point (start[1])...
        status = pb_inst_pbonly(signal_instruction[i + 1][0], Inst.END_LOOP, start[1], signal_instruction[i + 1][3]);
        # print(signal_instruction[i+1][0], Inst.END_LOOP, start[1], signal_instruction[i+1][3])
        self.errorCatcher(status)
        #------------------------------------------------------------------------------------
        # 3rd: sequence has 'camera' has again (XORed with concfg.camera)... with the signal sequence (signal_instruction) running...
        if N_trigger[0] >= 1:
            started = False
            for i in range(0, len(signal_instruction) - 1):
                # print(signal_instruction[i])
                if started:
                    # middle instructions...
                    status = pb_inst_pbonly(signal_instruction[i][0] ^ concfg.camera, signal_instruction[i][1],
                                            signal_instruction[i][2], signal_instruction[i][3])
                    # print(signal_instruction[i][0]^concfg.camera, signal_instruction[i][1], signal_instruction[i][2], signal_instruction[i][3])
                    self.errorCatcher(status)
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    start.append(pb_inst_pbonly(signal_instruction[0][0] ^ concfg.camera, Inst.LOOP, N_trigger[0],
                                                signal_instruction[0][3]))
                    # print(signal_instruction[0][0]^concfg.camera, Inst.LOOP, N_trigger[0], signal_instruction[0][3])
                    self.errorCatcher(start[2])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            status = pb_inst_pbonly(signal_instruction[i + 1][0] ^ concfg.camera, Inst.END_LOOP, start[2],
                                    signal_instruction[i + 1][3]);
            # print(signal_instruction[i+1][0]^concfg.camera, Inst.END_LOOP, start[2], signal_instruction[i+1][3])
            self.errorCatcher(status)
            # the buffer comes here.. 
            status = pb_inst_pbonly(0 ^ concfg.camera, Inst.CONTINUE, 0, t_buffer[1]) if t_buffer[1] >= 10 else 0
            # print(0^concfg.camera, Inst.CONTINUE, 0, t_buffer[1])
            self.errorCatcher(status)
        else:
            start.append(pb_inst_pbonly(concfg.camera, Inst.CONTINUE, 0, t_cam_response));
            self.errorCatcher(start[2])
        #------------------------------------------------------------------------------------
        # 4th: sequence runs w/o 'camera'... with the 'reference' sequence (instructionList[1])...
        started = False
        # buffer comes here...
        status = pb_inst_pbonly(0, Inst.CONTINUE, 0, t_buffer[2]) if t_buffer[2] >= 10 else 0
        # print(0, Inst.CONTINUE, 0, t_buffer[2])
        self.errorCatcher(status)
        for i in range(0, len(reference_instruction) - 1):
            # print(reference_instruction[i])
            if started:
                # instructions in the middle..
                status = pb_inst_pbonly(reference_instruction[i][0], reference_instruction[i][1],
                                        reference_instruction[i][2], reference_instruction[i][3])
                # print(reference_instruction[i][0], reference_instruction[i][1], reference_instruction[i][2],reference_instruction[i][3])
                self.errorCatcher(status)
            else:
                # first instruction... Inst.LOOP and repetition N_remaining(_reference)[1]...
                start.append(
                    pb_inst_pbonly(reference_instruction[0][0], Inst.LOOP, N_remaining[1], reference_instruction[0][3]))
                # print(reference_instruction[0][0], Inst.LOOP, N_remaining[1], reference_instruction[0][3])
                self.errorCatcher(start[3])
                started = True
        # last instruction being written.. with Inst.END_LOOP and start point (start[3])...
        status = pb_inst_pbonly(reference_instruction[i + 1][0], Inst.END_LOOP, start[3], reference_instruction[i + 1][3]);
        # print(reference_instruction[i+1][0], Inst.END_LOOP, start[3], reference_instruction[i+1][3])
        self.errorCatcher(status)
        # a last instruction with Inst.BRANCH of length 10 ns... returns the control to start[0] (the first instruction)
        status = pb_inst_pbonly(reference_instruction[i + 1][0], Inst.BRANCH, start[0], 10 * ns);
        # print(reference_instruction[i+1][0], Inst.BRANCH, start[0], 10*ns)
        self.errorCatcher(status)
        # 10 ns is very small compared to the time-scales involed (response time and esp. exposure time) in the camera program...
        #---------------------------------------------DONE-------------------------------------------

        status = pb_stop_programming();
        self.errorCatcher(status)
        status = pb_start();
        self.errorCatcher(status)
        self.pb_status = "running"
        # status = pb_close();                self.errorCatcher(status)
        # print('\x1b[38;2;250;0;0m\x10 PB STARTED\x1b[0m')
        # return instructionList
        return [t_exposure, N]

        # # ekhan theke 'signal' frame er part ta shuru hochhe...
        # start_trigger = pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.CONTINUE, 0, 38.96*us)     # camera trigger with pulse width of 38.96 us
        # for i in range(0, len(instructionList[0])-1):
        #     if started:
        #         status = pb_inst_pbonly(instructionList[0][i][0],instructionList[0][i][1],instructionList[0][i][2],instructionList[0][i][3])
        #         self.errorCatcher(status)
        #     else:
        #         start[0] = pb_inst_pbonly(instructionList[0][0][0], Inst.LOOP, N_repeat,instructionList[0][0][3])    # start korar somoy mention N_repeat bar loop ta cholbe...
        #         self.errorCatcher(start[0])
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
        #         self.errorCatcher(status)
        #     else:
        #         start.append(pb_inst_pbonly(instructionList[1][0][0], Inst.LOOP, N_repeat,instructionList[1][0][3]))    # start korar somoy mention N_repeat bar loop ta cholbe...
        #         self.errorCatcher(start[0])
        #         started = True
        # loop_end = pb_inst_pbonly(instructionList[1][i+1][0], Inst.END_LOOP, start[1], instructionList[1][i+1][3])
        # final_trigger = pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.BRANCH, start[0], 1*us)
        # # 'reference΄' frame er part ta sesh holo...
        # # ebar sequence end kora hobe..
        


    def run_sequence_for_camera_level_trigger_many_dc(self, instructionList, t_exposure, t_align, t_seq_total, N_total):
        # t_seq_total should be a list: t_seq_total = [t_seq_tot_sig, t_seq_tot_ref]
        # print(instructionList)
        # instructionList               => len(instructionList)=2; PB instructions for 'signal' and 'reference' sequences
        # instructionList[i]            => i=0: signal, i=1: reference
        # instructionList[0][i]         => i=0: 1st instruction of the signal part, i=1: 2nd instruction...
        # isntructionList[0][0][i]   => i=0: PBchannel number, i=1: Inst, i=2: Inst data, i=3: time (delay)
        # t_exposure = 100 *ms
        # configure()
        # print("Loading PB...")
        # print(t_exposure)
        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);
        self.errorCatcher(status)

        # -----------------------------korlam.. 24/05/2023-----------------------------
        t_cam_response = 87.7 * us
        t_jitter = 0 * us
        # below are two lists, with 2 elements each, one for signal and other for reference...
        # the operation should use np.floor() and not np.rint().. think..
        if N_total == []:
            # N_trigger = [int(np.floor((t_cam_response-t_jitter)/element)) for element in t_seq_total]
            # N_total = [int(np.floor((t_exposure - t_cam_response - t_jitter)/element)) for element in t_seq_total]
            # print(t_exposure)
            N_total = [int(np.floor(t_exposure / element)) for element in t_seq_total]
            # the total time is not t_exposure-t_cam_response but only t_exposure, the t_exposure need to be adjusted to
            # t_exposure + 87.7us
            
            # print(N_total)
            # need to check if (N_trigger + N_remaining) <= N_total) ??

            # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
            t_buffer = [t_exposure - N_t * t for (t, N_t) in zip(t_seq_total, N_total)]
            # t_exposure is just for reference.. Ask camera for actual value.. this causes problems in the odmr sequence logic
            # t_exposure = [(t_cam_response + N_total[i] * t_seq_total[i] + t_buffer[i]) for i in range(0, len(N_total))]  # this is in ns..
            t_exposure = [t_exposure, t_exposure]
            # print(t_exposure)
            # print(N_total)

        # -------------------------------------------------------------------------------------------
        # -------------------------revisit this below------------------------------------------------
        # -----this is the case when the sequence time varies by orders of magnitude and the exposure time needs to be modified to keep the exposure (number of photons collected) the same.....----------------
        # ------------------this can also be used to do experiments with some fixed number of repetitions of 'normal' sequences-----------------------
        else:
            # change t_exposure keeping N_total constant
            # this assignment is just for user to refer... The pulse sequence loop is controlled only by 'N'
            t_exposure = [(t_cam_response + N_total[i] * t_seq_total[i]) for i in range(0, len(N_total))]  # this is in ns..
            # print('this also runs!!!')
            # print(t_exposure)
            t_buffer = [0, 10*ns]

        N = N_total

        # print('N: ', N)
        # print(f't_seq_total = {t_seq_total}')
        # print(f't_buffer: {t_buffer}')
        # print('t_exposure: ', t_exposure)
        # print(f't_align = {t_align}')

        #-----------------------------korlam.. 24/05/2023-----------------------------

        signal_instruction = instructionList[0]
        reference_instruction = instructionList[1]
        # print(signal_instruction)
        # print(reference_instruction)

        # follow notebook for the details of the sequence...
        #------------------------------------------------------------------------------------

        # part 1: aligning the propeller
        DAQ_AO_CH = concfg.bx ^ concfg.by
        # DAQ_AO_test = concfg.by ^ concfg.bz

        pbchannels = concfg.bz #^ concfg.laser# ^ concfg.MW
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_align/2);
        self.errorCatcher(status)
        pbchannels = concfg.bx ^ concfg.by #^ concfg.laser# ^ concfg.MW
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_align/2);
        self.errorCatcher(status)

        # pbchannels = concfg.bz #^ concfg.laser# ^ concfg.MW
        # status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_align/2);
        # self.errorCatcher(status)
        
        # print(pbchannels, Inst.CONTINUE, 0, t_align)
        start = [0]  # store the 'start' points of the sequence: 1st instruction and loop start points

        pbchannels = 0 #^ concfg.laser
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, 10*ms);
        self.errorCatcher(status)

        # part 2: camera exposure edge and subsequent delay
        # do we need laser here?
        pbchannels = concfg.camera #^ concfg.laser #^ concfg.MW 
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_cam_response);
        self.errorCatcher(status)
        # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 3: the SIGNAL sequence loop while camera exposure
        instruction = signal_instruction
        if N[0] >= 1:  # use the else case for ODMR!!!!!!!!!
            started = False
            for i in range(0, len(instruction) - 1):
                # print(f'instruction[{i}]= {instruction[i]}')
                if started:
                    # middle instructions...
                    pbchannels = instruction[i][0] ^ concfg.camera
                    status = pb_inst_pbonly(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3]);
                    self.errorCatcher(status)
                    # print(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3])
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    # start[1]
                    pbchannels = instruction[0][0] ^ concfg.camera
                    start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, N[0], instruction[0][3]));
                    self.errorCatcher(start[1])
                    # print(pbchannels, Inst.LOOP, N[0], instruction[0][3])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            pbchannels = instruction[i + 1][0] ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[1], instruction[i + 1][3]);
            self.errorCatcher(status)
            # print(pbchannels, Inst.END_LOOP, start[1], instruction[i + 1][3])

            # the buffer comes here..
            pbchannels = 0 ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else 0
            self.errorCatcher(status)
            # print(pbchannels, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else None

        else:  # -------------check this out for normal ODMR!!!!!!--------------
            pbchannels = concfg.camera ^ concfg.laser ^ concfg.MW
            # start[1]
            # there is no buffer in this case
            # Exposure and laser happens for the same time
            start.append(pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_exposure[0]));
            self.errorCatcher(start[1])
            # print(pbchannels, Inst.CONTINUE, 0, t_exposure[0])

        # # part 4: aligning the propeller:
        # pbchannels = concfg.bz #^ concfg.laser
        # status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_align/2);
        # self.errorCatcher(status)
        # # print(pbchannels, Inst.CONTINUE, 0, t_align)
        # pbchannels = concfg.bx ^ concfg.by #^ concfg.laser
        # status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_align/2);
        # self.errorCatcher(status)

        # # part 5: camera exposure edge and subsequent delay
        # # do we need laser here?
        # pbchannels = concfg.camera #^ concfg.laser
        # status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_cam_response);
        # self.errorCatcher(status)
        # # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 4: gap between frames
        pbchannels = 0 #^ concfg.laser
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, 8*ms)
        self.errorCatcher(status)

        # part 5: camera exposure edge and subsequent delay
        # do we need laser here?
        pbchannels = concfg.camera #^ concfg.laser
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_cam_response);
        self.errorCatcher(status)
        # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 6: the REFERENCE sequence loop while camera exposure
        instruction = reference_instruction
        if N[1] >= 1:  # use the else case for ODMR!!!!!!!!!
            started = False
            for i in range(0, len(instruction) - 1):
                # print(f'instruction[{i}]= {instruction[i]}')
                if started:
                    # middle instructions...
                    pbchannels = instruction[i][0] ^ concfg.camera
                    status = pb_inst_pbonly(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3]);
                    self.errorCatcher(status)
                    # print(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3])
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    # start[2]
                    pbchannels = instruction[0][0] ^ concfg.camera
                    start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, N[1], instruction[0][3]));
                    self.errorCatcher(start[2])
                    # print(pbchannels, Inst.LOOP, N[1], instruction[0][3])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            pbchannels = instruction[i + 1][0] ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[2], instruction[i + 1][3]);
            self.errorCatcher(status)
            # print(pbchannels, Inst.END_LOOP, start[2], instruction[i + 1][3])

            # the buffer comes here..
            # the last instruction with Inst.BRANCH... return the control to start[0] (the first instruction)
            pbchannels = 0 ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.BRANCH, start[0], t_buffer[1]) if t_buffer[1] >= 10 else 0;
            self.errorCatcher(status)
            # print(pbchannels, Inst.BRANCH, start[0], t_buffer[1]) if t_buffer[1] >= 10 else 0

        else:  # check this out for normal ODMR!!!!!!
            pbchannels = concfg.camera ^ concfg.laser
            # start[2]
            start.append(pb_inst_pbonly(pbchannels, Inst.BRANCH, start[0], t_exposure[1]));
            self.errorCatcher(start[2])
            # print(pbchannels, Inst.BRANCH, start[0], t_exposure[1])
        #---------------------------------------------DONE-------------------------------------------

        status = pb_stop_programming();
        self.errorCatcher(status)
        status = pb_start();
        self.errorCatcher(status)
        self.pb_status = "running"
        # status = pb_close();                self.errorCatcher(status)
        # print('\x1b[38;2;250;0;0m\x10 PB STARTED\x1b[0m')
        # return instructionList
        return [t_exposure, N]

    def run_sequence_for_camera_level_trigger_many_ac_level_trigger(self, instructionList, t_exposure, t_align, t_seq_total, N_total):
        # t_seq_total should be a list: t_seq_total = [t_seq_tot_sig, t_seq_tot_ref]
        # print(instructionList)
        # instructionList               => len(instructionList)=2; PB instructions for 'signal' and 'reference' sequences
        # instructionList[i]            => i=0: signal, i=1: reference
        # instructionList[0][i]         => i=0: 1st instruction of the signal part, i=1: 2nd instruction...
        # isntructionList[0][0][i]   => i=0: PBchannel number, i=1: Inst, i=2: Inst data, i=3: time (delay)
        # t_exposure = 100 *ms
        # configure()
        # print("Loading PB...")
        # print(t_exposure)
        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);
        self.errorCatcher(status)

        # -----------------------------korlam.. 24/05/2023-----------------------------
        t_cam_response = 87.7 * us
        t_jitter = 0 * us
        # below are two lists, with 2 elements each, one for signal and other for reference...
        # the operation should use np.floor() and not np.rint().. think..
        if N_total == []:
            # N_trigger = [int(np.floor((t_cam_response-t_jitter)/element)) for element in t_seq_total]
            # N_total = [int(np.floor((t_exposure - t_cam_response - t_jitter)/element)) for element in t_seq_total]
            # print(t_exposure)
            N_total = [int(np.floor(t_exposure / element)) for element in t_seq_total]
            # the total time is not t_exposure-t_cam_response but only t_exposure, the t_exposure need to be adjusted to
            # t_exposure + 87.7us
            
            # print(N_total)
            # need to check if (N_trigger + N_remaining) <= N_total) ??

            # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
            t_buffer = [t_exposure - N_t * t for (t, N_t) in zip(t_seq_total, N_total)]
            # t_exposure is just for reference.. Ask camera for actual value.. this causes problems in the odmr sequence logic
            # t_exposure = [(t_cam_response + N_total[i] * t_seq_total[i] + t_buffer[i]) for i in range(0, len(N_total))]  # this is in ns..
            t_exposure = [t_exposure, t_exposure]
            # print(t_exposure)
            # print(N_total)

        # -------------------------------------------------------------------------------------------
        # -------------------------revisit this below------------------------------------------------
        # -----this is the case when the sequence time varies by orders of magnitude and the exposure time needs to be modified to keep the exposure (number of photons collected) the same.....----------------
        # ------------------this can also be used to do experiments with some fixed number of repetitions of 'normal' sequences-----------------------
        else:
            # change t_exposure keeping N_total constant
            # this assignment is just for user to refer... The pulse sequence loop is controlled only by 'N'
            t_exposure = [(t_cam_response + N_total[i] * t_seq_total[i]) for i in range(0, len(N_total))]  # this is in ns..
            # print('this also runs!!!')
            # print(t_exposure)
            t_buffer = [0, 10*ns]

        N = N_total

        # print('N: ', N)
        # print(f't_seq_total = {t_seq_total}')
        # print(f't_buffer: {t_buffer}')
        # print('t_exposure: ', t_exposure)
        # print(f't_align = {t_align}')

        #-----------------------------korlam.. 24/05/2023-----------------------------

        signal_instruction = instructionList[0]
        reference_instruction = instructionList[1]
        # print(signal_instruction)
        # print(reference_instruction)

        # follow notebook for the details of the sequence...
        #------------------------------------------------------------------------------------

        # part 1: aligning the propeller
        b_channels = concfg.bx ^ concfg.by ^ concfg.bz ^ concfg.start_trig
        status = pb_inst_pbonly(b_channels, Inst.CONTINUE, 0, t_align);
        self.errorCatcher(status)
        # print(b_channels, Inst.CONTINUE, 0, t_align)

        start = [0]  # store the 'start' points of the sequence: 1st instruction and loop start points

        # part 2: camera exposure edge and subsequent delay
        # do we need laser here?
        pbchannels = concfg.camera #^ concfg.laser #^ concfg.MW
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_cam_response);
        self.errorCatcher(status)
        # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 3: the SIGNAL sequence loop while camera exposure
        instruction = signal_instruction
        if N[0] >= 1:  # use the else case for ODMR!!!!!!!!!
            started = False
            for i in range(0, len(instruction) - 1):
                # print(f'instruction[{i}]= {instruction[i]}')
                if started:
                    # middle instructions...
                    pbchannels = instruction[i][0] ^ concfg.camera
                    status = pb_inst_pbonly(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3]);
                    self.errorCatcher(status)
                    # print(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3])
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    # start[1]
                    pbchannels = instruction[0][0] ^ concfg.camera
                    start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, N[0], instruction[0][3]));
                    self.errorCatcher(start[1])
                    # print(pbchannels, Inst.LOOP, N[0], instruction[0][3])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            pbchannels = instruction[i + 1][0] ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[1], instruction[i + 1][3]);
            self.errorCatcher(status)
            # print(pbchannels, Inst.END_LOOP, start[1], instruction[i + 1][3])

            # the buffer comes here..
            pbchannels = 0 ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else 0
            self.errorCatcher(status)
            # print(pbchannels, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else None

        else:  # -------------check this out for normal ODMR!!!!!!--------------
            pbchannels = concfg.camera ^ concfg.laser ^ concfg.MW
            # start[1]
            # there is no buffer in this case
            # Exposure and laser happens for the same time
            start.append(pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_exposure[0]));
            self.errorCatcher(start[1])
            # print(pbchannels, Inst.CONTINUE, 0, t_exposure[0])

        # part 4: gap between frames
        pbchannels = 0 #^ concfg.laser
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, 8*ms)
        # print(pbchannels, Inst.CONTINUE, 0, 8*ms)
        self.errorCatcher(status)

        # part 5: camera exposure edge and subsequent delay
        # do we need laser here?
        pbchannels = concfg.camera #^ concfg.laser
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_cam_response);
        self.errorCatcher(status)
        # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 6: the REFERENCE sequence loop while camera exposure
        instruction = reference_instruction
        if N[1] >= 1:  # use the else case for ODMR!!!!!!!!!
            started = False
            for i in range(0, len(instruction) - 1):
                # print(f'instruction[{i}]= {instruction[i]}')
                if started:
                    # middle instructions...
                    pbchannels = instruction[i][0] ^ concfg.camera
                    status = pb_inst_pbonly(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3]);
                    self.errorCatcher(status)
                    # print(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3])
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    # start[2]
                    pbchannels = instruction[0][0] ^ concfg.camera
                    start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, N[1], instruction[0][3]));
                    self.errorCatcher(start[2])
                    # print(pbchannels, Inst.LOOP, N[1], instruction[0][3])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            pbchannels = instruction[i + 1][0] ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[2], instruction[i + 1][3]);
            self.errorCatcher(status)
            # print(pbchannels, Inst.END_LOOP, start[2], instruction[i + 1][3])

            # the buffer comes here..
            # the last instruction with Inst.BRANCH... return the control to start[0] (the first instruction)
            pbchannels = 0 ^ concfg.camera
            status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, start[0], t_buffer[1]) if t_buffer[1] >= 10 else 0;
            self.errorCatcher(status)
            # print(pbchannels, Inst.BRANCH, start[0], t_buffer[1]) if t_buffer[1] >= 10 else 0

        else:  # check this out for normal ODMR!!!!!!
            pbchannels = concfg.camera ^ concfg.laser
            # start[2]
            start.append(pb_inst_pbonly(pbchannels, Inst.CONTINUE, start[0], t_exposure[1]));
            self.errorCatcher(start[2])
            # print(pbchannels, Inst.BRANCH, start[0], t_exposure[1])
        pb_inst_pbonly(0, Inst.CONTINUE, 0, 14*ns)
        pb_inst_pbonly(0, Inst.STOP, 0, 10*ns)
        #---------------------------------------------DONE-------------------------------------------

        status = pb_stop_programming();
        self.errorCatcher(status)
        status = pb_start();
        self.errorCatcher(status)
        self.pb_status = "running"
        # status = pb_close();                self.errorCatcher(status)
        # print('\x1b[38;2;250;0;0m\x10 PB STARTED\x1b[0m')
        # return instructionList
        return [t_exposure, N]

    def run_sequence_for_camera_level_trigger_many_ac(self, instructionList, t_exposure, t_align, t_seq_total, N_total, Nsamples):
        # t_seq_total should be a list: t_seq_total = [t_seq_tot_sig, t_seq_tot_ref]
        # print(instructionList)
        # instructionList               => len(instructionList)=2; PB instructions for 'signal' and 'reference' sequences
        # instructionList[i]            => i=0: signal, i=1: reference
        # instructionList[0][i]         => i=0: 1st instruction of the signal part, i=1: 2nd instruction...
        # isntructionList[0][0][i]   => i=0: PBchannel number, i=1: Inst, i=2: Inst data, i=3: time (delay)
        # t_exposure = 100 *ms
        # configure()
        # print("Loading PB...")
        # print(t_exposure)
        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);
        self.errorCatcher(status)

        # -----------------------------korlam.. 24/05/2023-----------------------------
        t_cam_response = 38.96 * us
        t_jitter = 0 * us
        # below are two lists, with 2 elements each, one for signal and other for reference...
        # the operation should use np.floor() and not np.rint().. think..
        if N_total == []:
            # N_trigger = [int(np.floor((t_cam_response-t_jitter)/element)) for element in t_seq_total]
            # N_total = [int(np.floor((t_exposure - t_cam_response - t_jitter)/element)) for element in t_seq_total]
            # print(t_exposure)
            N_total = [int(np.floor(t_exposure / element)) for element in t_seq_total]
            # the total time is not t_exposure-t_cam_response but only t_exposure, the t_exposure need to be adjusted to
            # t_exposure + 87.7us
            
            # print(N_total)
            # need to check if (N_trigger + N_remaining) <= N_total) ??

            # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
            t_buffer = [t_exposure - N_t * t for (t, N_t) in zip(t_seq_total, N_total)]
            # t_exposure is just for reference.. Ask camera for actual value.. this causes problems in the odmr sequence logic
            # t_exposure = [(t_cam_response + N_total[i] * t_seq_total[i] + t_buffer[i]) for i in range(0, len(N_total))]  # this is in ns..
            t_exposure = [t_exposure, t_exposure]
            # print(t_exposure)
            # print(N_total)

        # -------------------------------------------------------------------------------------------
        # -------------------------revisit this below------------------------------------------------
        # -----this is the case when the sequence time varies by orders of magnitude and the exposure time needs to be modified to keep the exposure (number of photons collected) the same.....----------------
        # ------------------this can also be used to do experiments with some fixed number of repetitions of 'normal' sequences-----------------------
        else:
            # change t_exposure keeping N_total constant
            # this assignment is just for user to refer... The pulse sequence loop is controlled only by 'N'
            t_exposure = [(t_cam_response + N_total[i] * t_seq_total[i]) for i in range(0, len(N_total))]  # this is in ns..
            # print('this also runs!!!')
            # print(t_exposure)
            t_buffer = [0, 10*ns]

        N = N_total

        # print('N: ', N)
        # print(f't_seq_total = {t_seq_total}')
        # print(f't_buffer: {t_buffer}')
        # print('t_exposure: ', t_exposure)
        # print(f't_align = {t_align}')

        #-----------------------------korlam.. 24/05/2023-----------------------------

        signal_instruction = instructionList[0]
        reference_instruction = instructionList[1]
        # print(signal_instruction)
        # print(reference_instruction)

        # follow notebook for the details of the sequence...
        #------------------------------------------------------------------------------------

        # part 1: aligning the propeller
        b_channels = concfg.bx ^ concfg.by ^ concfg.bz ^ concfg.start_trig
        b_channels_only = concfg.bx ^ concfg.by ^ concfg.bz
        status = pb_inst_pbonly(b_channels, Inst.CONTINUE, 0, t_align);
        self.errorCatcher(status)
        # print(b_channels, Inst.CONTINUE, 0, t_align)
        # consider the inductance of the coil and adjust the timescale
        status = pb_inst_pbonly(b_channels_only, Inst.CONTINUE, 0, 5*ms);
        self.errorCatcher(status)
        start = []  # store the 'start' points of the sequence: 1st instruction and loop start point

        # part 2: camera exposure edge and subsequent delay
        # do we need laser here?
        pbchannels = concfg.camera ^ concfg.laser ^ b_channels_only#^ concfg.MW
        start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, Nsamples, t_cam_response))
        self.errorCatcher(start[0])
        # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 3: the SIGNAL sequence loop while camera exposure
        instruction = signal_instruction
        if N[0] >= 1:  # use the else case for ODMR!!!!!!!!!
            started = False
            for i in range(0, len(instruction) - 1):
                # print(f'instruction[{i}]= {instruction[i]}')
                if started:
                    # middle instructions...
                    pbchannels = instruction[i][0] ^ b_channels_only
                    status = pb_inst_pbonly(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3]);
                    self.errorCatcher(status)
                    # print(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3])
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    # start[1]
                    pbchannels = instruction[0][0] ^ b_channels_only
                    start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, N[0], instruction[0][3]));
                    self.errorCatcher(start[1])
                    # print(pbchannels, Inst.LOOP, N[0], instruction[0][3])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            pbchannels = instruction[i + 1][0] ^ b_channels_only
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[1], instruction[i + 1][3]);
            self.errorCatcher(status)
            # print(pbchannels, Inst.END_LOOP, start[1], instruction[i + 1][3])

            # the buffer comes here..
            pbchannels = 0 ^ b_channels_only
            status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else 0
            self.errorCatcher(status)
            # print(pbchannels, Inst.CONTINUE, 0, t_buffer[0]) if t_buffer[0] >= 10 else None

        else:  # -------------check this out for normal ODMR!!!!!!--------------
            pbchannels = concfg.laser ^ concfg.MW ^ b_channels_only
            # start[1]
            # there is no buffer in this case
            # Exposure and laser happens for the same time
            start.append(pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_exposure[0]));
            self.errorCatcher(start[1])
            # print(pbchannels, Inst.CONTINUE, 0, t_exposure[0])

        # part 5: camera exposure edge and subsequent delay
        # do we need laser here?
        pbchannels = concfg.camera ^ concfg.laser ^ b_channels_only
        status = pb_inst_pbonly(pbchannels, Inst.CONTINUE, 0, t_cam_response);
        self.errorCatcher(status)
        # print(pbchannels, Inst.CONTINUE, 0, t_cam_response)

        # part 6: the REFERENCE sequence loop while camera exposure
        instruction = reference_instruction
        if N[1] >= 1:  # use the else case for ODMR!!!!!!!!!
            started = False
            for i in range(0, len(instruction) - 1):
                # print(f'instruction[{i}]= {instruction[i]}')
                if started:
                    # middle instructions...
                    pbchannels = instruction[i][0] ^ b_channels_only
                    status = pb_inst_pbonly(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3]);
                    self.errorCatcher(status)
                    # print(pbchannels, instruction[i][1], instruction[i][2], instruction[i][3])
                else:
                    # first instruction.. with Inst.LOOP and repetitions N_trigger(_signal)[0]
                    # start[2]
                    pbchannels = instruction[0][0] ^ b_channels_only
                    start.append(pb_inst_pbonly(pbchannels, Inst.LOOP, N[1], instruction[0][3]));
                    self.errorCatcher(start[2])
                    # print(pbchannels, Inst.LOOP, N[1], instruction[0][3])
                    started = True

            # the last instruction.. with Inst.END_LOOP and start point (start[2])...
            pbchannels = instruction[i + 1][0] ^ b_channels_only
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[2], instruction[i + 1][3]);
            self.errorCatcher(status)
            # print(pbchannels, Inst.END_LOOP, start[2], instruction[i + 1][3])

            # the buffer comes here..
            # the last instruction with Inst.BRANCH... return the control to start[0] (the first instruction)
            pbchannels = 0 ^ b_channels_only
            status = pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[0], t_buffer[1]) if t_buffer[1] >= 10 else 0;
            self.errorCatcher(status)
            # print(pbchannels, Inst.BRANCH, start[0], t_buffer[1]) if t_buffer[1] >= 10 else 0

        else:  # check this out for normal ODMR!!!!!!
            pbchannels = concfg.laser ^ b_channels_only
            # start[2]
            start.append(pb_inst_pbonly(pbchannels, Inst.END_LOOP, start[0], t_exposure[1]));
            self.errorCatcher(start[2])
            # print(pbchannels, Inst.BRANCH, start[0], t_exposure[1])
        pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.CONTINUE, 0, 14*ns)
        pb_inst_pbonly(concfg.laser, Inst.STOP, 0, 10*ns)
        #---------------------------------------------DONE-------------------------------------------

        status = pb_stop_programming();
        self.errorCatcher(status)
        # status = pb_start();
        # self.errorCatcher(status)
        # self.pb_status = "running"
        # status = pb_close();                self.errorCatcher(status)
        # print('\x1b[38;2;250;0;0m\x10 PB STARTED\x1b[0m')
        # return instructionList
        return [t_exposure, N]

    def run_only_daq(self, t_align):      # incoming t_align in ms
        ### Run this function to keep the manipulation ON during scanpt change
        # DAQ_AO_CH = concfg.bx ^ concfg.by
        if t_align < 50*ms:
            instructionList = [[concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, t_align/2],
                            [concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align/2]]
        else:
            # t_align is already in ns
            duty = 0.07
            x = int(duty*t_align/(1-duty)/10)*10       # in ms
            # print(x)
            instructionList = [[concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, x],
                            [concfg.bz, Inst.CONTINUE, 0, (t_align/2 - x)],
                            [concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align/2]]
        self.run_sequence_for_diode([instructionList])

    def custom_trigger(self, t_exposure, t_align):
        # instructionList, t_exposure, t_align, t_seq_total, N_total
        # t_exposure = 5 *ms
        # t_field = 125 *ms
        # t_cam_response = 38.96 *us
        
        print(f"t_exposure = {t_exposure}")
        print(f"t_align = {t_align}")
        t_field = t_align/2
        trigger_width = 1 *ms
        t_total = 2000 *ms
        start = []

        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);    self.errorCatcher(status)

        number_of_frames = int(int(t_field/t_exposure)/2)
        print(number_of_frames)
        # first the Bz part
        start.append(pb_inst_pbonly(concfg.camera ^ concfg.bz ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width));    self.errorCatcher(start[0])
        status = pb_inst_pbonly(concfg.bz ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.camera ^ concfg.bz ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.bz ^ concfg.MW ^ concfg.laser, Inst.END_LOOP, start[0], (t_exposure-trigger_width));    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.camera ^ concfg.bz ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.bz ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        # next the Bx part
        start.append(pb_inst_pbonly(concfg.camera ^ concfg.bx ^ concfg.MW ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width));    self.errorCatcher(start[1])
        status = pb_inst_pbonly(concfg.bx ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.camera ^ concfg.bx ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.bx ^ concfg.laser, Inst.END_LOOP, start[1], (t_exposure-trigger_width));    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.camera ^ concfg.bx ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.bx ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        # remaining part without field
        number_of_frames = int(int(t_total/t_exposure)/2)
        print(number_of_frames)
        start.append(pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width));    self.errorCatcher(start[2])
        status = pb_inst_pbonly(0 ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.camera ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        status = pb_inst_pbonly(0 ^ concfg.MW ^ concfg.laser, Inst.END_LOOP, start[2], (t_exposure-trigger_width));    self.errorCatcher(status)

        # status = pb_inst_pbonly(0 ^ concfg.laser, Inst.BRANCH, start[2], 100 *ns);    self.errorCatcher(status)

        status = pb_stop_programming();    self.errorCatcher(status)
        status = pb_start();    self.errorCatcher(status)
        self.pb_status = "running"

    def custom_trigger_rot_field(self, t_exposure, t_align, t_meas, t_align_extended):
        # t_exposure is the exposure time returned by camera
        # t_align is the extended rotating field pattern time determined by the pattern, t_align_rot_extended in main program
        # t_meas is the measurement time ~100 ms
        # t_exposure = 5 *ms
        # t_field = 125 *ms
        # t_cam_response = 38.96 *us
        print(f"Target measurement time = {t_meas} ns")
        print(f"t_exposure = {t_exposure} ns")
        print(f"t_align = {t_align} ns")

        t_field = t_align       # the field channels are kept open for the whole time of t_align
        trigger_width = 100 *us
        start = []
        b_channels = concfg.bx ^ concfg.by ^ concfg.bz ^ concfg.start_trig      # open/close all PB channels for rotating field aignment

        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);    self.errorCatcher(status)

        number_of_frames = int(np.ceil(np.ceil(t_field/t_exposure)/2))       # (number of frames)/2 during the field ON time
        print(f"{number_of_frames*2} in field ON time")
        # print(concfg.camera ^ b_channels ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width)
        start.append(pb_inst_pbonly(concfg.camera ^ b_channels ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width));    self.errorCatcher(start[0])
        status = pb_inst_pbonly(b_channels ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        # print(b_channels ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width))
        status = pb_inst_pbonly(concfg.camera ^ b_channels ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        # print(concfg.camera ^ b_channels ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width)
        status = pb_inst_pbonly(b_channels ^ concfg.MW ^ concfg.laser, Inst.END_LOOP, start[0], (t_exposure-trigger_width));    self.errorCatcher(status)
        # print(b_channels ^ concfg.MW ^ concfg.laser, Inst.END_LOOP, start[0], (t_exposure-trigger_width))
        # status = pb_inst_pbonly(concfg.camera ^ b_channels ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        # status = pb_inst_pbonly(b_channels ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)

        # start.append(pb_inst_pbonly(concfg.camera ^ b_channels ^ concfg.MW ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width));    self.errorCatcher(start[1])
        # status = pb_inst_pbonly(b_channels ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        # status = pb_inst_pbonly(concfg.camera ^ b_channels ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        # status = pb_inst_pbonly(b_channels ^ concfg.laser, Inst.END_LOOP, start[1], (t_exposure-trigger_width));    self.errorCatcher(status)
        # status = pb_inst_pbonly(concfg.camera ^ b_channels ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        # status = pb_inst_pbonly(b_channels ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)

        # number_of_frames = int(np.ceil(np.ceil((t_meas-t_field)/t_exposure)/2))     # (number of frames)/2 during field OFF time
        # print(f"{number_of_frames*2} in field OFF time ideally")
        # if (t_meas - t_align_extended) > 0:

        t_meas_modified = t_meas-(t_align_extended-t_align)
        print(f"t_align extended = {t_align_extended} ns")
        print(f"Modified t_meas = {t_meas_modified}")
        number_of_frames = int(np.ceil(np.ceil(t_meas_modified/t_exposure)/2))
        print(f"{number_of_frames*2} in field OFF time")
        start.append(pb_inst_pbonly(concfg.camera ^ concfg.laser, Inst.CONTINUE, number_of_frames, trigger_width));    self.errorCatcher(start[1])
        status = pb_inst_pbonly(0 ^ concfg.laser, Inst.CONTINUE, 0, (t_exposure-trigger_width));    self.errorCatcher(status)
        status = pb_inst_pbonly(concfg.camera ^ concfg.MW ^ concfg.laser, Inst.CONTINUE, 0, trigger_width);    self.errorCatcher(status)
        status = pb_inst_pbonly(0 ^ concfg.MW ^ concfg.laser, Inst.BRANCH, start[1], (t_exposure-trigger_width));    self.errorCatcher(status)

        # status = pb_inst_pbonly(0 ^ concfg.laser, Inst.BRANCH, start[2], 100 *ns);    self.errorCatcher(status)

        status = pb_stop_programming();    self.errorCatcher(status)
        status = pb_start();    self.errorCatcher(status)
        self.pb_status = "running"

    def custom_trigger_rot_field_continuous(self, t_exposure, t_align, t_meas, t_align_extended):
        # t_exposure is the exposure time returned by camera
        # t_align is the extended rotating field pattern time determined by the pattern, t_align_rot_extended in main program
        # t_meas is the total measurement time ~100 ms
        # t_exposure = 5 *ms
        # t_field = 125 *ms
        # t_cam_response = 38.96 *us
        print(f"Target measurement time = {t_meas} ns")
        print(f"t_exposure = {t_exposure} ns")
        print(f"t_align_rot = {t_align} ns")
        t_field = t_align       # the field channels are kept open for the whole time of t_align
        t_meas_modified = t_meas-(t_align_extended-t_align)
        print(f"t_align extended = {t_align_extended} ns")
        print(f"Modified t_meas = {t_meas_modified}")

        trigger_width = 1 *ms
        start = []
        b_channels = concfg.bx ^ concfg.by ^ concfg.bz ^ concfg.start_trig      # open/close all PB channels for rotating field aignment
        b_channels_only = concfg.bx ^ concfg.by ^ concfg.bz
        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);    self.errorCatcher(status)

        number_of_frames = int(np.ceil(np.ceil(t_field/t_exposure)/2))       # (number of frames)/2 during the field ON time
        print(f"{number_of_frames*2} in field ON time")

        # print(b_channels ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width)
        start.append(pb_inst_pbonly(b_channels ^ concfg.laser, Inst.CONTINUE, 0, t_field));    self.errorCatcher(start[0])
        status = pb_inst_pbonly(b_channels_only ^ concfg.laser, Inst.BRANCH, 0, t_meas_modified);    self.errorCatcher(status)

        status = pb_stop_programming();    self.errorCatcher(status)
        status = pb_start();    self.errorCatcher(status)
        self.pb_status = "running"

    def custom_trigger_rot_field_continuous_uncorrected(self, t_exposure, t_align, t_meas, t_align_extended):
        # t_exposure is the exposure time returned by camera
        # t_align is the extended rotating field pattern time determined by the pattern, t_align_rot_extended in main program
        # t_meas is the total measurement time ~100 ms
        # t_exposure = 5 *ms
        # t_field = 125 *ms
        # t_cam_response = 38.96 *us
        print(f"Target measurement time = {t_meas} ns")
        print(f"t_exposure = {t_exposure} ns")
        print(f"t_align_rot = {t_align} ns")
        t_field = t_align       # the field channels are kept open for the whole time of t_align
        t_meas_modified = t_meas-(t_align_extended-t_align)
        print(f"t_align extended = {t_align_extended} ns")
        print(f"Modified t_meas = {t_meas_modified}")

        if t_meas_modified<0:
            n_frames_rot_alignment = np.floor(t_align/t_exposure)      # UNITS???? t_exposure is seconds, and t_align_rot is in seconds (calculated from dt in seconds)
            t_align_extended = n_frames_rot_alignment*t_exposure       # extended/modified t_align_rot to fit even number of frames (signals and references)
            t_meas_modified = t_meas+(t_align-t_align_extended)
            print(f"t_align extended adjusted = {t_align_extended} ns")
            print(f"Modified t_meas = {t_meas_modified}")

        trigger_width = 1 *ms
        start = []
        b_channels = concfg.bx ^ concfg.by ^ concfg.bz ^ concfg.start_trig      # open/close all PB channels for rotating field aignment
        b_channels_only = concfg.bx ^ concfg.by ^ concfg.bz
        pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);    self.errorCatcher(status)

        number_of_frames = int(np.ceil(np.ceil(t_field/t_exposure)/2))       # (number of frames)/2 during the field ON time
        print(f"{number_of_frames*2} in field ON time")

        # print(b_channels ^ concfg.laser, Inst.LOOP, number_of_frames, trigger_width)
        start.append(pb_inst_pbonly(b_channels ^ concfg.laser, Inst.CONTINUE, 0, t_field));    self.errorCatcher(start[0])
        status = pb_inst_pbonly(b_channels_only ^ concfg.laser, Inst.BRANCH, 0, t_meas_modified);    self.errorCatcher(status)

        status = pb_stop_programming();    self.errorCatcher(status)
        status = pb_start();    self.errorCatcher(status)
        self.pb_status = "running"

    def focus_adjustment_sequence(self, t_exposure, t_align, Nsamples, focus_time: float):
        """PB sequence during the focus adjustment.
        Designs and loads the sequence but does not start the sequence.
        
        Parameters
        ----------
            t_exposure : float
                frame exposure time / sampling time in [ns]
            t_align : float
                alignment time / rotating field sequence time in [ns]
            Nsamples : int
                total number of samples to acquire, signal+ref
        """

        self.min_focus_time = focus_time       # minimum focus time before restarting again if required
        self.sequence_time = t_align+2*t_exposure       # total sequence time; 2 exposure frames: signal + ref

        start = []
        b_channels = concfg.bx ^ concfg.by ^ concfg.bz ^ concfg.start_trig      # open/close all PB channels for rotating field alignment
        # b_channels_only = concfg.bx ^ concfg.by ^ concfg.bz
        # try:
        # pb_reset()
        status = pb_start_programming(PULSE_PROGRAM);    self.errorCatcher(status)

        self.n_repeat = int(np.ceil(self.min_focus_time/self.sequence_time))       # number of repetitions of the sequence before looping
        self.min_focus_time = self.n_repeat*self.sequence_time

        start.append(pb_inst_pbonly(b_channels ^ concfg.laser, Inst.LOOP, self.n_repeat, t_align));    self.errorCatcher(start[0])
        status = pb_inst_pbonly(concfg.laser, Inst.END_LOOP, 0, 5*ms + Nsamples*t_exposure);    self.errorCatcher(status)

        status = pb_stop_programming();    self.errorCatcher(status)
        self.pb_status = "programmed"
        # except Exception as e:
