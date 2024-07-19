# %% Initialization and Definition
# reset
import PBcontrol_v2 as pbctrl, sequencecontrol as seqctrl, SGcontrol as sgctrl, \
    connectionConfig as conCfg, matplotlib.pyplot as plt, numpy as np, sys, time, dialog
from spinapi import ns, us, ms, Inst
from os.path import isdir, isfile;
from importlib import import_module
import dcam
# from matlab import engine as eng
# %matplotlib qt5

global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, SG, expCfg, filenumber, hdcamcon, instr, data_raw

# def main():
print("\x10 \x1b[0mImports Successful...")
plt.rcParams.update({'figure.max_open_warning': 0})  # No warnings on opening mult fig windows

instr = 'cam_levelm'  # optons: cam, cam_level1, cam_levelm
expCfgFile = 'esr' + '_config'
N_total = [2000, 2000]  # a 2-element list (sig and ref) for total number of repetitions.. if empty then N_total is allowed to change for each scanpt..
N_total = []

trial_run = ['n']  # 1st = SG
# all times here are in seconds
t_exposure = 50e-3
t_meas = t_exposure;
t_align = 12e-3
align_field = [50, 0, 0]

# ------------- Plotting details----------------------------
seq_no_plot = [-1]
voltage_unit = 1  # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100  # The dpi of the displayed pulse sequence plot
plotPulseSequence = False
# t_exposure = 0.1     # in seconds.. eta dorkar o nei.. user will define it..
livePlotUpdate = False


def initialize_instr(sequence):
    global hdcamcon
    try:
        pbctrl.pb_close()
        pbctrl.configurePB()
        print(
            '\x10 PB: \x1b[38;2;250;250;0mv' + pbctrl.pb_get_version() + '\x1b[0m')  # Display the PB board version using pb_get_version()
    except:
        print("Error Initializing PB !!")
    # finally:
    #     sys.exit()
    if trial_run[0] == 'n' and sequence not in ['aom_timing', 'rodelay']:
        # Do not initialize SG if it is a trial run or the sequence is present in the list ['aom_timing', 'rodelay']
        SG = sgctrl.init_sg(conCfg.serialaddr,conCfg.model_name)
        if SG != '':
            sgctrl.enable_SG_op();
            print("SG Output Enabled...")
            sgctrl.set_SG_amp(expCfg.MW_power)
            sgctrl.setup_SG_pulse_mod();
            print("SG Ext Pulse Mod Enabled...")
    else:
        SG = None
    return [SG]


def close_all(SG):
    """End the measurement
    Closes all the instruments.
    """
    if (trial_run[0] == 'n') and (SG is not None):
        SG.write('freq2.87ghz');
        # SG.write('LCAL')
        sgctrl.disable_SG_op()
        sgctrl.uninit_sg()

    # pbctrl.pb_init();
    closed = True
    pbctrl.pb_stop();
    pbctrl.pb_close();
    print("Pulse Blaster closed...\x1b[0m")
    return True


def initialize_exp(instr, expCfg):
    """Initialize the experimental parameters and check whether the sequence can be sent to the PB...
    Also, start the initial PB sequence.
    Returns: 
    param_save_format - contains the save format of the parameters in the parmater file
    savePath          -
    seqArgList        - 
    expParamList      -
    Nscanpts          -
    param             -
    instructionList   - 
    Include trial run check so that the error plots are only displayed when it is not a trial run
    """
    param = expCfg.scannedParam
    n_error = 0
    param_save_format = expCfg.formattingSaveString
    sequenceArgs = expCfg.updateSequenceArgs()  # Variables used in the pulse sequence
    expParamList = expCfg.updateExpParamList()  # List of experimental parameters

    # ------------------------------------------------------------------
    # est_time = expCfg.t_tot * expCfg.Nsamples * expCfg.N_scanPts ??
    # check the scanned parameters for errors, if errors are found, remove those params
    if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr',
                               'drift_seq']:  # for sequences except ESR, set MW frequency
        seqArgList = [param[seq_no_plot[-1]]]
        seqArgList.extend(sequenceArgs)  # Make a dummy seqArgList just to create it
        [n_error, param] = seqctrl.param_err_check(instr, expCfg.sequence, expCfg.PBchannels, seqArgList,
                                                   expCfg.scannedParam, expCfg.N_scanPts)
        if n_error > 0:
            print('\x1b[1;37;41m' + 'Err: Check Sequences...\x1b[0m')
            print('\x1b[38;2;250;0;0m' + str(n_error) + '\x1b[0m parameters removed...')
            # Close Error plots??
            close_plots = dialog.yesno_box('Close Plots', 'Close the Error Plots?')
            if close_plots == 'yes':
                plt.close('all')
            print("\x10 Sequences checked... \x1b[38;2;100;250;0mErrors removed...\x1b[0m")
            # est_time = (2*expCfg.t_AOM*Nscanpts + sum(param))*expCfg.Nsamples
        else:
            print('\x1b[38;2;100;250;0m----No Errors----\x1b[0m')
    else:  # for ESR exp, set MW freq to start pt of the scan -> there is no scannedParam
        seqArgList = sequenceArgs

    # define 'Nscanpts' as the length of the 'param' variable...
    Nscanpts = len(param)
    if Nscanpts > 0:
        print("\x10 \x1b[38;2;250;100;10m%d\x1b[0m scan pts" % Nscanpts)
    else:
        print("\x1b[38;2;200;200;10mSubtle errors...\x1b[0m")

    # Plotting the pulse sequences -------------------------------------------
    if plotPulseSequence:
        print(expCfg.PBchannels)
        view_sequence(expCfg.sequence, expCfg.PBchannels, seqArgList, False, param, seq_no_plot, seq_plot_dpi)

    # Start the initial sequence now ------------
    print("Starting Initial Sequence...")
    if trial_run[0] == 'n':
        instructionList = [start_initial_PB_seq()]
        pbctrl.run_sequence_for_diode(instructionList)
        print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
        if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq', 'aom_timing',
                                   'rodelay'] and trial_run[0] == 'n':
            sgctrl.set_SG_freq(expCfg.MW_freq)
    else:
        print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    # stop the initial PB sequence after initializing the parameters
    status = pbctrl.pb_stop()
    pbctrl.errorCatcher(status)

    return [seqArgList, expParamList, Nscanpts, param, instructionList]


def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[conCfg.laser, Inst.CONTINUE, 0, 500 * ms],
                       [conCfg.laser, Inst.BRANCH, 0, 500 * ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...

    return instructionList


# def calc_contrast(signal, reference, op):
# if op == '+-':
#     contrast = signal - reference
# elif op == '-+':
#     contrast = -signal + reference
# elif op == 's/r':
#     contrast = signal/reference
# elif op == 'r/s':
#     contrast = reference/signal
# return contrast

# the view_sequence needs to be modified to (sub)plot both the signal and reference sequence for camera
def view_sequence(sequence, PBchannels, seqArgList, only_plot=False, parameter=[0], seq_no_plot=[0], plot_dpi=100):
    """Plot the pulse sequence for visualization. Cross-check with an oscilloscope."""
    for i in range(0, len(seq_no_plot)):
        if only_plot == False:
            if sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                # expCfg.t_AOM = t_AOM*parameter[seq_no_plot[i]]/parameter[0]
                # sequenceArgs = expCfg.updateSequenceArgs()
                # seqArgList = [parameter[seq_no_plot[i]]]
                # seqArgList.extend(sequenceArgs)
                # expCfg.t_AOM = t_AOM
                seqArgList[0] = parameter[seq_no_plot[i]]
        plt.figure(dpi=plot_dpi)
        the_list = pbctrl.PB_program(instr, sequence, seqArgList)
        for j in range(0, len(the_list)):
            instructionList = the_list[j][0]
            inst_times = the_list[j][4]

            if instr == 'cam' or instr == 'cam_levelm':
                plt.subplot(1, 2, (j + 1))
            [t_us, channelPulses, yTicks] = seqctrl.plot_sequence(instructionList, PBchannels)
            for channel in channelPulses:
                plt.plot(t_us, list(channel))
            plt.yticks(yTicks, PBchannels.keys())  # Include the names of the PB channels
            plt.xlabel('Time (us)')
            # plt.ylabel('Channel')
            plt.title(sequence + ' Plot. Param @ ' + str(seq_no_plot[i]) + ': ' + str(
                parameter[seq_no_plot[i]]) + 'ns\nTransitions [us]: ' + str([vals for vals in inst_times.values()]),
                      fontsize=10)
            plt.show()
    # Way to display the total no of instructions in the plot, since it will vary with different scanPts:
    # print("Total %d inst" % len(inst_times.keys()))

# ----------------------------------------------------------------------------

def pulse_sequence_loop(*timings, t_seq_total, parameters, sequence, seqArgList):#, trigger_event, condition):
    """ Control the scanpoint loop and pulse sequence for experiment as well as camera trigger
    
    """
    Nscanpts = len(parameters)
    inst_set_time_list = []
    # time.sleep(0.02)
    print("starting loop..")
    
    for i_scanpt in range(0, Nscanpts):
        t_seq_total_i = t_seq_total[i_scanpt]
        param = parameters[i_scanpt]
        print("\x10 ", i_scanpt + 1, ' / ', Nscanpts, ': ', param)
        
        # TODO: thread this part as well for parallel operation..
        inst_set_start_time = time.perf_counter()
        instructionList = []
        if sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
            if trial_run[0] == 'n':
                sgctrl.set_SG_freq(param)
        else:
            seqArgList[0] = param
        # inst_set_time_list.append(time.perf_counter()-inst_set_start_time)
        inst_set_start_time = time.perf_counter()
        the_list = pbctrl.PB_program(instr, sequence, seqArgList)
        # print('the_list: ', the_list)
        for i in range(0, len(the_list)):
            instructionList.append(the_list[i][0])
        # print('instructionList: ', instructionList)

        if trial_run[0] == 'n':
            # print('Starting PB sequence...')
            if instr == 'cam':
                pbctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total_i, N_total)
            elif instr == 'cam_level1':
                pbctrl.run_sequence_for_camera_level_trigger_1(instructionList)
            elif instr == 'cam_levelm':
                t_align = timings[2]
                pbctrl.run_sequence_for_camera_level_trigger_many(instructionList, t_exposure, t_align, t_seq_total_i, N_total)
            
            inst_set_time_list.append(time.perf_counter()-inst_set_start_time)
                
    return [i_scanpt, inst_set_time_list]


if __name__ == '__main__':
    # global data_raw

    expCfg = import_module(expCfgFile)
    # expCfg.N_scanPts = len(expCfg.scannedParam)
    t_AOM = expCfg.t_AOM
    clk_cyc = 1e3 / conCfg.PBclk  # One clk cycle of PB = inverse of the clk freq of PB
    print('\x10 \x1b[38;2;250;250;0mRunning ' + expCfg.saveFileName + ' sequence\x1b[0m')
    [SG] = initialize_instr(expCfg.sequence)

    seqctrl.check_params(expCfgFile)
    print("\x10 Parameter checks completed...")

    [seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr,expCfg)
    if (expCfg.sequence == 'esr_seq'):
        t_manip = np.zeros(len(param))
    elif (expCfg.sequence == 'pesr_seq'):
        t_manip = np.ones(len(param)) * expCfg.t_duration
    else:
        # t_manip = [param for ]
        t_manip = np.array(param)

    # below 2 lines are for T1 with varying 'init' duration such that the duty cycle remains same...
    # t_seq_total_sig = [t_AOM*param[i]/param[0] for i in range(0,len(param))] + t_manip
    # t_seq_total_ref = [t_AOM*param[i]/param[0] for i in range(0,len(param))] + t_manip

    t_seq_total_sig = expCfg.t_AOM + t_manip
    t_seq_total_ref = expCfg.t_AOM + t_manip  # for other sequences..
    # t_seq_total_ref = np.zeros(len(t_manip),dtype='float64')        # for T1 seq with laser ALWAYS ON
    t_seq_total = np.transpose(np.array([t_seq_total_sig, t_seq_total_ref]))

    # -------------------------- 20062023-------------------------
    
    closed = False
    display_parameters = dialog.yesno_box(expCfg.saveFileName + ' Params', 'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' % (str(conCfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples, param[0] / expCfg.plotXaxisUnits, param[-1] / expCfg.plotXaxisUnits))

    if display_parameters == 'yes':
        if trial_run[0] == 'n':
            try:
                # continue_run = 'yes'
                # save_flag = False
                start_time = time.perf_counter()    # entire measurement time
                PB_wait_time = []
                
                print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)    # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
                
                for i_run in range(0, expCfg.Nruns):     # TODO: add Nruns loop as well
                    print("Run: ", i_run + 1, ' / ', expCfg.Nruns)
                    
                    acq_start_time = time.perf_counter()
                    # time.sleep(5)
                    print("going to loop...")
                    [i_scanpt, inst_set_time] = pulse_sequence_loop(t_exposure * 1e9, N_total, t_align * 1e9, t_seq_total=t_seq_total, parameters=param, sequence=expCfg.sequence, seqArgList=seqArgList)
                    
                    # [i_scanpt_inst, inst_set_time_list] = scan_thread.result
                    acq_time = time.perf_counter() - acq_start_time
                    
                # all runs complete.. process for display...
                end_time = time.perf_counter()
                exec_time = end_time - start_time  # entire measurement time
                # Close all after full acquisition
                closed = close_all(SG)                    

            except KeyboardInterrupt:
                end_time = time.perf_counter()
                exec_time = end_time - start_time  # entire measurement time
                print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
                # Plot the last run data upto the point where it was interrupted...
            
            finally:
                if not closed:
                    closed = close_all(SG)
            
            instructionList = []
            if expCfg.sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                None
            else:
                seqArgList[0] = param[seq_no_plot[-1]]
            print(seqArgList)
            the_list = pbctrl.PB_program(instr, expCfg.sequence, seqArgList)
            # print(the_list)
            # print(t_seq_total)
            for i in range(0, len(the_list)):
                instructionList.append(the_list[i][0])
            if instr == 'cam':
                pbctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total[seq_no_plot[-1]], N_total)
            elif instr == 'cam_level':
                pbctrl.run_sequence_for_camera_level_trigger_1(instructionList)
            elif instr == 'cam_levelm':
                pbctrl.run_sequence_for_camera_level_trigger_many(instructionList, t_exposure * 1e9, t_align * 1e9, t_seq_total[seq_no_plot[-1]], N_total)
            # pbctrl.pb_close(); print("\x10 Pulse Blaster closed...\x1b[0m")
    else:
        print("Measurement aborted...")
        if not closed:
            closed = close_all(SG)
    plt.plot(inst_set_time)

    # %%
