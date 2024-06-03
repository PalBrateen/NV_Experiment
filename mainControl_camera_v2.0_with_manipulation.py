# %% Initialization and Definition
# reset
import Camcontrol_v2_0 as camctrl, PBcontrol_v2 as PBctrl, sequencecontrol as seqctrl, SGcontrol as SGctrl, \
    connectionConfig as conCfg, matplotlib.pyplot as plt, numpy as np, sys, time, dialog, cv2, os, skimage as sk, \
    tifffile as tifff, DAQcontrol as daqctrl
from spinapi import ns, us, ms, Inst
from os.path import isdir, isfile;
from importlib import import_module
import dcam

# %matplotlib qt5

# from matlab import engine as eng

# def main():
print("\x10 \x1b[0mImports Successful...")
plt.rcParams.update({'figure.max_open_warning': 0})  # No warnings on opening mult fig windows
global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, SG, expCfg, filenumber, hdcamcon

instr = 'cam_levelm'  # optons: cam, cam_level1, cam_levelm
expCfgFile = 'esr' + '_config'
N_total = [2000,
           2000]  # a 2-element list (sig and ref) for total number of repetitions.. if empty then N_total is allowed to change for each scanpt..
N_total = []

trial_run = ['n', 'n']  # 1st = SG, 2nd = camera
# all times here are in seconds
t_exposure = 1e-3 if trial_run[1] == 'y' else None
t_meas = t_exposure;
t_align = 500e-3
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
        PBctrl.pb_close()
        PBctrl.configurePB()
        print(
            '\x10 PB: \x1b[38;2;250;250;0mv' + PBctrl.pb_get_version() + '\x1b[0m')  # Display the PB board version using pb_get_version()
    except:
        print("Error Initializing PB !!")
    # finally:
    #     sys.exit()
    if trial_run[0] == 'n' and sequence not in ['aom_timing',
                                                'rodelay']:  # Do not initialize SG if it is a trial run or the sequence is present in the list ['aom_timing', 'rodelay']
        SG = SGctrl.init_sg(conCfg.serialaddr,
                            conCfg.model_name)  # Initialize SG using RS-232 address and the model name
        if SG != '':
            SGctrl.enable_SG_op();
            print("SG Output Enabled...")
            SGctrl.set_SG_amp(expCfg.MW_power)
            SGctrl.setup_SG_pulse_mod();
            print("SG Ext Pulse Mod Enabled...")
    else:
        SG = None
    hdcamcon = camctrl.init_cam() if trial_run[1] == 'n' else None
    return [SG, hdcamcon]


def close_all(SG, hdcamcon, ao_task):  # arguments are 'SG' and 'hdcamcon'
    """End the measurement
    Closes all the instruments.
    """
    if (trial_run[0] == 'n') and (SG is not None):
        SG.write('freq2.87ghz');
        SG.write('lcal')
        SGctrl.disable_SG_op()
        SGctrl.uninit_sg()

    if (trial_run[1] == 'n') and (hdcamcon is not None):
        camctrl.uninit_cam();
        print("\x1b[38;2;10;250;50mCamera Closed...")

    if (trial_run[1] == 'n') and (hdcamcon is not None):
        daqctrl.close_daq_task(ao_task)
        print('AO stopped...')

    # PBctrl.pb_init(); 
    PBctrl.pb_stop();
    PBctrl.pb_close();
    print("Pulse Blaster closed...\x1b[0m")
    return True


def initialize_exp(instr, expCfg):
    # Initialize the experimental parameters and check whether the sequence can be sent to the PB...
    # Also, start the initial PB sequence
    # Return variables: 
    # param_save_format - contains the save format of the parameters in the parmater file
    # savePath          -
    # seqArgList        - 
    # expParamList      -
    # Nscanpts          -
    # param             -
    # instructionList   - 
    ###### Include trial run check so that the error plots are only displayed when it is not a trial run

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
    if trial_run[0] == 'n' or trial_run[1] == 'n':
        instructionList = [start_initial_PB_seq()]
        PBctrl.run_sequence_for_diode(instructionList)
        print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
        if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq', 'aom_timing',
                                   'rodelay'] and trial_run[0] == 'n':
            SGctrl.set_SG_freq(expCfg.MW_freq)
    else:
        print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    # Create the data save folder...
    cwd = os.getcwd()  # current working directory - generally 'NV_Experiment' folder. Data folder is created a level up
    # rfind() gives the last occurence of '\\'
    savePath = cwd[0:cwd.rfind('\\')] + "\\Saved_Data\\" + time.strftime("%Y-%m-%d", time.localtime()) + '\\'

    if not (isdir(savePath)):
        os.makedirs(savePath)
    print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')

    # viewing live image and adjust the exposure time
    if trial_run[1] == 'n':
        camctrl.live_frames(timeout_ms=1000)

    # now select the ROI
    roi = camctrl.select_roi() if trial_run[1] == 'n' else None
    # [716, 880, 300, 300]

    # stop the initial PB sequence after initializing the parameters
    status = PBctrl.pb_stop()
    PBctrl.errorCatcher(status)

    # something wrong in the return statement, if else both returns the same parameters!!!!!!!!!!!
    return [roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] if trial_run[
                                                                                                                 1] == 'n' else [
        roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList]


def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[conCfg.laser, Inst.CONTINUE, 0, 500 * ms],
                       [conCfg.laser, Inst.BRANCH, 0, 500 * ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...

    return instructionList


def read_save_details(roi, N_scanpts, Nsamples_expCfg):
    """defines the 'frame_per_cyc' to be captured, and the write formats of the 'scannedParam' and the 'data'
    Returns: 
    frames_per_cyc            - define the number of frames captured in one cycle of the sequence, 1 signal, 1 ref, so 2. 
    scannedparam_write_format - write format of the scannedparams in the data file
    """
    hsize = roi[2];
    vsize = roi[3];
    frames_per_cyc = [2]  # signal and reference (2)
    Nsamples = frames_per_cyc[0] * Nsamples_expCfg

    scannedparam_write_format = "%g\t" * N_scanpts
    scannedparam_write_format = scannedparam_write_format[0:-1] + "\n"

    return [frames_per_cyc, Nsamples, scannedparam_write_format]


def extra_param_save_details(scan_time_list, exec_time):
    extra_params_format = ' %s\t%0.2f\n %s\t%0.2f\n %s\t%0.4g\n %s\t%0.4g\n %s\t%g\n %s\t%d\n %s\t%d\n %s\t%d\n %s\t%d\n'  # %s\t%f\n
    # scan_time_list in ms; exec_time in seconds
    param_list = ['Max_scan_time(us):', max(scan_time_list) * 1e6, 'Min_scan_time(us):', min(scan_time_list) * 1e6,
                  'Total_scan_time(s):', np.sum(scan_time_list), 'Total_run_time(s):', exec_time, 'Step:',
                  (expCfg.scannedParam[1] - expCfg.scannedParam[0]), 'X0:', roi[0], 'Y0:', roi[1], 'W:', roi[2], 'H:',
                  roi[3]]
    # param_list[1] = i_scanpt+1..... Ekhane ki hbe?? Kon parameter save korbo??
    # ekhane ekta if kore, jodi 'i'=1 hoe, thle first param_list ta return korbe, nhle porer param_list ta return korbe.
    # Eta jodi kora hoe, thle runs>1 hole ba multi-param scan hole, notun param er sathe tar details save kora jabe...
    # Config file e ekta variable lagbe jeta dekhabe je kon variable ta scan hoechhe, other than scannedParam
    param_list = tuple(param_list)
    return [extra_params_format, param_list]


def prepare_for_saving(savePath):
    global file_number
    if 'file_number' in vars():  # Returns a dict of all local variables
        print("Previous file number: \x1b[38;2;250;150;0m" + file_number + '\x1b[0m')

    while True:
        file_number = input("File name:" + expCfg.saveFileName + "_camera_#: ")
        datafilename = savePath + expCfg.saveFileName + "_camera_" + file_number + ".tiff"
        if (isfile(datafilename)):
            print('\x1b[38;2;250;250;0mFile exists. Retry...\x1b[0m')
            # continue
        else:
            break

    paramfilename = savePath + expCfg.saveFileName + "_camera_" + "params_" + file_number + ".txt"

    return [paramfilename, datafilename, file_number]


def save_data(datafilename, data_raw, scannedparam_write_format, param):
    print("\x10 Saving to file.....")

    tifff.imwrite(datafilename, data_raw, imagej=True)
    # sk.io.imsave(datafilename, data_raw, plugin='tifffile')

    print("\x10 Data saved to\x1b[38;2;100;250;50m %s_%s\x1b[0m !!!" % (expCfg.saveFileName, file_number))
    return True


def save_parameters(paramfilename, param_save_format, param_list):
    paramFile = open(paramfilename, 'a')
    paramFile.write(param_save_format % tuple(param_list))
    paramFile.close()


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
        the_list = PBctrl.PB_program(instr, sequence, seqArgList)
        for j in range(0, len(the_list)):
            instructionList = the_list[j][0]
            inst_times = the_list[j][4]

            if instr == 'cam':
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


def acquire_data(*timings, t_seq_total_i, roi, Nsamples, parameter, sequence, seqArgList):
    # use starred expression here.. very difficult to generalize otherwise.....
    # setup next scan iteration (e.g. for ESR experiment, change microwave frequency; for T2 experiment, reprogram pulseblaster with new delay)
    global trial_run
    t_exposure = timings[0]
    N_total = timings[1]
    instructionList = []
    if sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
        if trial_run[0] == 'n':
            SGctrl.set_SG_freq(parameter)
    else:
        seqArgList[0] = parameter
    the_list = PBctrl.PB_program(instr, sequence, seqArgList)
    print('the_list: ', the_list)
    for i in range(0, len(the_list)):
        instructionList.append(the_list[i][0])
    print('instructionList: ', instructionList)

    timeout_ms = 2000  # revisit...

    if trial_run[1] == 'n':

        print('Starting PB sequence...')
        scan_start_time = time.perf_counter()  # time in seconds
        if instr == 'cam':
            PBctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total_i, N_total)
        elif instr == 'cam_level1':
            PBctrl.run_sequence_for_camera_level_trigger_1(instructionList)
        elif instr == 'cam_levelm':
            t_align = timings[2]
            PBctrl.run_sequence_for_camera_level_trigger_many(instructionList, t_exposure, t_align, t_seq_total_i,
                                                              N_total)
        # print("PB req start time: ", (time.perf_counter()-scan_start_time))

        hsize = roi[2]
        vsize = roi[3]
        # frames=np.zeros((vsize,hsize,Nsamples),dtype=np.uint16)
        frames = np.zeros((Nsamples, vsize, hsize), dtype='uint16')
        t1 = time.perf_counter()
        print('Starting wait_capevent_frameready...')
        result = hdcamcon.wait_capevent_frameready(timeout_ms)
        PB_wait_time.append(time.perf_counter() - t1)

        timeout_happened = 0
        if (result) is not True:
            # frame does not come
            if result != camctrl.dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                print('-NG: Dcam.wait_event() failed with error {}'.format(result))
                # break
            # TIMEOUT error happens
            timeout_happened += 1
            if timeout_happened == 1:
                print('Waiting for a frame to arrive.', end='')
                if hdcamcon.get_propertyvalue(
                        camctrl.dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == camctrl.dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
                    print(' Check your trigger source.', end='')
                else:
                    print(' Check <timeout_ms>.', end='')
                print(' Press Ctrl+C to abort.')
            else:
                print('.')
                if timeout_happened > 5:
                    timeout_happened = 0
            # continue

        for j in range(0, Nsamples):
            print(j, end='')
            if j > 0:
                result = hdcamcon.wait_capevent_frameready(timeout_ms)
            # wait_capevent_frameready() succeeded
            frame = hdcamcon.get_lastframedata()
            cv2.imshow("Image", frame)
            cv2.waitKey(1)
            if frame is not False:
                print(' frame elo...')
                frames[j, :, :] = frame
            else:
                print("No frame...")
        scan_end_time = time.perf_counter()
        print('Single scan point end...')
    else:
        if instr == 'cam':
            PBctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total_i, N_total)
        elif instr == 'cam_level':
            PBctrl.run_sequence_for_camera_level_trigger_1(instructionList)
        elif instr == 'cam_levelm':
            t_align = timings[3]
            PBctrl.run_sequence_for_camera_level_trigger_many(instructionList, t_exposure, t_align, t_seq_total_i,
                                                              N_total)
    status = PBctrl.pb_stop();
    PBctrl.errorCatcher(status);

    scan_time = (scan_end_time - scan_start_time)  # in seconds
    # if ith_scan_pt==0:
    #     time.sleep(100)
    return [frames, scan_time, instructionList, PB_wait_time]


def process_data(i_max, roi, n_frames, data_raw):
    """Process data_raw for plotting. Find the mean signal, reference and the contrast
    Note: This is done after the whole scan is performed."""
    hsize = roi[2];
    vsize = roi[3];
    mean_sig = np.zeros(i_max);
    mean_ref = np.zeros(i_max);
    contrast = np.zeros(i_max);
    for i_scanpt in range(0, i_max):
        # signal_frames = np.array([data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize] for j in range(0,n_frames,2)])
        # reference_frames = np.array([data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize] for j in range(1,n_frames,2)])
        # calculate the sum of pixels for each signal and ref frame and form an array
        # sum_of_signal_frames = np.array([np.sum(data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize]) for j in range(0,n_frames,2)])
        # sum_of_reference_frames = np.array([np.sum(data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize]) for j in range(1,n_frames,2)])

        signal_frames = data_raw[i_scanpt, 0::2, :, :];
        reference_frames = data_raw[i_scanpt, 1::2, :, :];

        # mean_sig[i] = np.mean(sum_of_signal_frames)
        # # print(mean_sig[i])
        # mean_ref[i] = np.mean(sum_of_reference_frames)
        # # print(mean_ref[i])
        mean_sig[i_scanpt] = np.sum(np.mean(signal_frames, 0))
        mean_ref[i_scanpt] = np.sum(np.mean(reference_frames, 0))

        contrast[i_scanpt] = mean_sig[i_scanpt] / mean_ref[i_scanpt]
    return [mean_sig, mean_ref, contrast]


def plot_data(i_max, param, processed_data, x_label, x_unit, roi, n_frames, live=False):
    """Plot the signal, reference and the contrast"""
    # n_frames = total frames in one cycle of sequence..
    [mean_sig, mean_ref, contrast] = processed_data
    xValues = param[0:i_max]
    if live:
        plt.plot([x / x_unit for x in xValues], contrast[0:i_max], 'b--')
        plt.show()
        plt.pause(0.0001)
    else:
        # plt.figure("Signal, Reference & Contrast Plots: %d points" %(i_max))
        plt.figure()
        plt.subplot(121)
        plt.plot([x / x_unit for x in xValues], mean_sig[0:i_max], '.-', [x / x_unit for x in xValues],
                 mean_ref[0:i_max], '.-')
        plt.legend(['Avg Sig', 'Avg Ref'])
        plt.xlabel(x_label)

        plt.subplot(122)
        plt.plot([x / x_unit for x in xValues], contrast[0:i_max], '.-b')
        plt.xlabel(x_label);
        plt.ylabel('Contrast')
        plt.show()


# ----------------------------------------------------------------------------

# if __name__ == '__main__':
# instructionList = main()

expCfg = import_module(expCfgFile)
# expCfg.N_scanPts = len(expCfg.scannedParam)
t_AOM = expCfg.t_AOM
clk_cyc = 1e3 / conCfg.PBclk  # One clk cycle of PB = inverse of the clk freq of PB
print('\x10 \x1b[38;2;250;250;0mRunning ' + expCfg.saveFileName + ' sequence\x1b[0m')
[SG, hdcamcon] = initialize_instr(expCfg.sequence)

seqctrl.check_params(expCfgFile)
print("\x10 Parameter checks completed...")

[roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr,
                                                                                                                expCfg)
if trial_run[1] == 'n':
    hsize = roi[2];
    vsize = roi[3];

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
if trial_run[1] == 'n':
    [frames_per_cyc, Nsamples, scannedparam_write_format] = read_save_details(roi, Nscanpts, expCfg.Nsamples)

# -------------------------- 20062023-------------------------

# expCfg.Nsamples = 1 mane holo at each scan pt 1 frame (like N APD readouts for averaging). So multiply by frames_per_cyc[0]
# thle Nsamples = 2 hbe... at each scan pt 2 frames.. then go to next scan pt..
# Nsamples = ekta scanpt e total kotogulo signal (ba reference) frames...
# Nsamples = frames_per_cyc[0]*expCfg.Nsamples        # feed this as the 'n_frames' in capture() of Camcontrol.py... eta read_save_details() theke ashbe ebar...
save_flag = False
if trial_run[1] == 'n':  # same as if hdcamcon is not None...
    # check this part and simplify..
    # just coverting everything to new library for now..----------------------

    camctrl.configure_camera(instr)  # configure the camera for triggered acquisition

    # the ROI is already set in initialize_experiment()
    # camctrl.set_roi(roi,True)   # set the ROI for the acquisition

    t_exposure = hdcamcon.get_propertyvalue(
        propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME)  # mind the value... in seconds !!!

# status = PBctrl.pb_stop(); PBctrl.errorCatcher(status) already stopped...
if trial_run[1] == 'n':
    ao_task = daqctrl.config_ao()
    output_aovoltage = daqctrl.coil_calibration(output_field=align_field)

closed = False
display_parameters = dialog.yesno_box(expCfg.saveFileName + ' Params',
                                      'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' % (
                                          str(conCfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples,
                                          param[0] / expCfg.plotXaxisUnits, param[-1] / expCfg.plotXaxisUnits))
if display_parameters == 'yes':
    if trial_run[1] == 'n':
        try:
            # continue_run = 'yes'
            # save_flag = False

            camctrl.configure_camera(instr)
            camctrl.query_cam_values()
            t_exposure = hdcamcon.setget_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME, val=0.001)
            PB_wait_time = []
            # set buffer and start capture sequence..
            while not camctrl.startingcapture(size_buffer=Nsamples, sequence=True):
                print("Retrying..")

            # start DC field
            daqctrl.start_ao(ao_task, data=output_aovoltage)
            time.sleep(t_align * 3)

            print(
                "\x10 Camera configured for %d frames..." % expCfg.Nsamples)  # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
            print('\x10 %d frames in each cycles...' % frames_per_cyc[0])
            # from above, the total number of frames to acquire is thus frames_per_cyc[0]*expCfg.Nsamples = Nsamples

            print("------Acquiring %dx%d------" % tuple([roi[2], roi[3]]))
            for i_run in range(0, expCfg.Nruns):
                # data_raw = np.zeros((vsize*Nscanpts,hsize*Nsamples))   # note the dimension...
                data_raw = np.zeros((Nscanpts, Nsamples, vsize, hsize), dtype='uint16')
                scan_time_list = []  # time (in seconds) for each scannedParam
                # print("\x10",i_run+1,'/',expCfg.Nruns)
                start_time = time.perf_counter()  # in seconds
                for i_scanpt in range(0, Nscanpts):
                    # time.sleep(5)
                    # if i_scanpt>0:
                    #     break
                    print("\x10 ", i_scanpt + 1, ' / ', Nscanpts, ': ', param[i_scanpt])
                    # print("Exp [ms]: ",(cam["t_exposure"].value)*1e3)
                    # expCfg.t_AOM = t_AOM*param[i_scanpt]/param[0]
                    # print("t_AOM [ns] ",expCfg.t_AOM)
                    # sequenceArgs = expCfg.updateSequenceArgs()
                    # seqArgList = [param[i_scanpt]]
                    # seqArgList.extend(sequenceArgs)
                    # print(seqArgList)
                    [frames, scan_time, instructionList, PB_wait_time] = acquire_data(t_exposure * 1e9, N_total,
                                                                                      t_align * 1e9,
                                                                                      t_seq_total_i=t_seq_total[
                                                                                          i_scanpt], roi=roi,
                                                                                      Nsamples=Nsamples,
                                                                                      parameter=param[i_scanpt],
                                                                                      sequence=expCfg.sequence,
                                                                                      seqArgList=seqArgList)
                    # data_raw[i_scanpt*vsize:(i_scanpt+1)*vsize,:] = np.concatenate(tuple(frames[:,:,j] for j in range(0,Nsamples)), axis=1)       # 2d-array... ekhane for j in range(0,frames_per_cyc[0]) chhilo... otake range(0,Nsamples) kore dilam.. think.. concat korte gele j should run over all the frames.. jodi total 4 frames thake (2 sig, 2 ref), then j=0,1,2,3...
                    data_raw[i_scanpt, :, :, :] = frames
                    scan_time_list.append(scan_time)  # in seconds

                end_time = time.perf_counter()  # in seconds
                exec_time = end_time - start_time  # in seconds
                processed_data = process_data(i_scanpt, roi, Nsamples,
                                              data_raw)  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples,
                          live=False)
                print("\x10 Execution time = %.2fs" % exec_time)  # in seconds
                print("\x10 Scan time = %.2fs" % np.sum(scan_time_list))  # in seconds; scan_time_list in seconds
                # Close all if this is the last run
                if i_run + 1 == expCfg.Nruns and not closed:
                    closed = close_all(SG, hdcamcon, ao_task)

                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Data file saving
                    if save_flag == False:  # Ask for file number iff no save was performed
                        [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
                    save_flag = save_data(datafilename, data_raw, scannedparam_write_format, param)
                else:
                    print("\x10 \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")
                if (i_run + 1) < expCfg.Nruns:
                    continue_run = dialog.yesno_box('Continue', "Continue Run?")
                    if continue_run == 'yes':
                        continue
                    else:
                        sys.exit("Run Interrupted...")

        except KeyboardInterrupt:
            end_time = time.perf_counter()
            exec_time = end_time - start_time  # in seconds
            print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
            # Plot the last run data upto the point where it was interrupted...
            processed_data = process_data(i_scanpt, roi, Nsamples,
                                          data_raw)  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
            plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples,
                      live=False)
            # Then ask whether to save it...
            savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
            if savefile_yn == 'yes':
                # Ask for filename only if there was no save operation
                if save_flag == False:
                    [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
                save_flag = save_data(datafilename, data_raw, scannedparam_write_format, param)
        # nicher ei 'else' ta kno dewa hoyechhilo??? 18-06-2023
        # else:
        #     if (i_run+1)==expCfg.Nruns: # Save data of the final run here
        #         save_data(len(conCfg.input_terminals), datafilename, data_raw, expCfg.Nsamples, reads_per_cyc)
        finally:
            if not closed:
                closed = close_all(SG, hdcamcon, ao_task)

            print(
                "\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))

            print("\x10 shape(frames) = \x1b[38;2;250;150;50m%d x %d x %d\x1b[0m" % np.shape(frames))
            # Save parameters (only if there was one save operation)...
            if save_flag:
                expParamList[1] = i_scanpt + 1  # expParamList[1] -> value of N_scanPts
                expParamList[3] = i_run + 1  # expParamList[3] -> value of Nruns
                save_parameters(paramfilename, param_save_format, expParamList)

                [next_params_format, expParamList] = extra_param_save_details(scan_time_list, exec_time)
                save_parameters(paramfilename, next_params_format, expParamList)

            # scan_time_list in seconds; exec_time in seconds
            # ----------------------------------------------------------------------

            # matlab_engine = eng.start_matlab()

            # if expCfg.sequence == 'esr':
            #     analysis = eng.ESR_raw_analysis()
    elif trial_run[1] == 'y':
        instructionList = []
        if expCfg.sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
            None
        else:
            seqArgList[0] = param[seq_no_plot[-1]]
        print(seqArgList)
        the_list = PBctrl.PB_program(instr, expCfg.sequence, seqArgList)
        # print(the_list)
        # print(t_seq_total)
        for i in range(0, len(the_list)):
            instructionList.append(the_list[i][0])
        if instr == 'cam':
            PBctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total[seq_no_plot[-1]], N_total)
        elif instr == 'cam_level':
            PBctrl.run_sequence_for_camera_level_trigger_1(instructionList)
        elif instr == 'cam_levelm':
            PBctrl.run_sequence_for_camera_level_trigger_many(instructionList, t_exposure * 1e9, t_align * 1e9,
                                                              t_seq_total[seq_no_plot[-1]], N_total)
        # PBctrl.pb_close(); print("\x10 Pulse Blaster closed...\x1b[0m")
else:
    print("Measurement aborted...")
    if not closed:
        closed = close_all(SG, hdcamcon, ao_task=None)

    # ther lower statement were done as the None values were not assigned previously
    # check before discarding...
    # if expCfg.sequence not in ['aom_timing','rodelay']:
    #     if trial_run[0] == 'n' and trial_run[1] == 'n':
    #         close_all(cam, SG)
    #     elif trial_run[0] == 'y' and trial_run[1] == 'n':
    #         close_all(cam)
    #     elif trial_run[1] == 'y' and trial_run[0] == 'n':
    #         close_all(SG)
    # sys.exit()

# return instructionList

# def camera_triggers(changing_N, t_exposure, t_seq_total):
#     """Should N (repetitions) vary in one exposure time or remain constant??"""
#     # take the t_exposure as a list and then return as a list after changing it (if required)...
#     # also, take the t_seq_total as the whole list and do the full calculation at one go..
#     # then supply the changed values to the run_sequence_for_camera()..
#     t_cam_response = 38.96 *us
#     if changing_N:
#         # here, the N changes and the t_exposure remains the same
#         # exposure time ta ke adjust korte hbe depending on the  


#         # below are two lists, with 2 elements each, one for signal and other for reference...
#         # the operation should use np.floor() and not np.rint().. think..
#         N_trigger = [int(np.floor(t_cam_response/element)) for element in t_seq_total]
#         N_remaining = [int(np.floor((t_exposure - t_cam_response)/element)) for element in t_seq_total]
#         N_total = [int(np.floor(t_exposure/element)) for element in t_seq_total]
#         # need to check if (N_trigger + N_remaining) <= N_total) ??
#         print(N_trigger)
#         print(N_remaining)
#         print(N_total)
#         # defining the buffer times (in nanoseconds) -  time when there will be no sequence running - quiet period
#         t_buffer1 = (t_exposure - t_cam_response) - N_remaining[0]*t_seq_total[0]
#         t_buffer2 = t_cam_response - N_trigger[0]*t_seq_total[0]
#         t_buffer3 = (t_exposure - t_cam_response) - N_remaining[1]*t_seq_total[1]
#         t_buffer4 = t_cam_response - N_trigger[1]*t_seq_total[1]
#         t_buffer = [t_buffer1, t_buffer2, t_buffer3, t_buffer4]
#     else:
#         N_total = [int(np.floor(t_exposure/element)) for element in t_seq_total]

#     N = [N_trigger, N_remaining, N_total]

#     return [N, t_exposure, t_buffer]


# %%
