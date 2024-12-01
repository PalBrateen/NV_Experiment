# mainControl_camera -v2.0
# Initialization and Definition
# reset
import Camcontrol as camctrl, PBcontrol_v2 as pbctrl, sequencecontrol as seqctrl, SGcontrol as sgctrl, \
    connectionConfig as concfg, matplotlib.pyplot as plt, numpy as np, sys, time, dialog, cv2, os, skimage as sk, \
    tifffile as tifff, DAQcontrol as daqctrl, threading, logging
from spinapi import ns, us, ms, Inst
from os.path import isdir, isfile;
from importlib import import_module
import dcam,  warnings, matplotlib
# from matlab import engine as eng
# %matplotlib qt5

global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, SG, expCfg, f_number, hdcamcon, instr, data#, data_raw_time

# def main():
# print("\x10 \x1b[0mImports Successful...")
plt.rcParams.update({'figure.max_open_warning': 0})  # No warnings on opening mult fig windows

instr = 'cam_levelm_trigger_ao'
# options: cam, cam_level1, cam_levelm, cam_timeseries, cam_timeseries_trigger_ao, cam_levelm_trigger_ao
expCfgFile = 'esr' + '_config'
N_total = [2, 2]  # a 2-element list (sig and ref) for total number of repetitions.. if empty then N_total is allowed to change for each scanpt..
# N_total = [] #if expCfgFile == 'esr_config' else [12354,12354]

trial_run = ['n', 'n']  # 1st = SG, 2nd = camera
# all times here are in seconds
fps = 996.3
rot_field_amp = 50      # field amplitude in [gauss]
rot_field_freq = 29     # field frequency [Hz]

direction = 'z'
rot_angle = 45      # rotation angle in degrees
t_meas = 50        # measurement time [ms] for cam_timeseries_trigger_ao only

t_exposure = 3     # exposure [ms]; set to 5 ms for 'timeseries' measurements
t_align_dc = 250         # [ms]; set to 250 ms for 'timeseries' mmesurements
align_field = [ 50 , 0, +50 ]          # [bx] G
# align_field = [ 0 , 0 , 0 ]
# test_field = [15, 15]        # [by, bz] G
# test_field = [0, 0]


# cam_timeseries_trigger_ao: take the timeseries measurement at 5ms exposure time for some time interval after performing the 3x(Bz, Bx) + fractional rotating field alignment
# cam_levelm_trigger_ao: perform ODMR with triggered rotating field for diffusion control with level triggered camnera acquisition
t_exposure /= 1e3       # [s]
# ------------- Plotting details----------------------------
seq_no_plot = [-1]
voltage_unit = 1  # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100  # The dpi of the displayed pulse sequence plot
plotPulseSequence = False
livePlotUpdate = False
ao_task = None

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
        SG = sgctrl.init_sg(concfg.serialaddr,concfg.model_name)
        if SG != '':
            sgctrl.enable_SG_op();
            print("SG Output Enabled...")
            sgctrl.set_SG_amp(expCfg.MW_power)
            SG.write("freq2.87ghz")
            sgctrl.setup_SG_pulse_mod();
            print("SG Ext Pulse Mod Enabled...")
    else:
        SG = None
    hdcamcon = camctrl.init_cam() if trial_run[1] == 'n' else None
    return [SG, hdcamcon]


def close_all(SG, hdcamcon, ao_task):
    """End the measurement
    Closes all the instruments.
    """
    if (trial_run[1] == 'n') and (hdcamcon is not None):
        # stops capture, releases buffer, closes camera and uninitializes DCAM-API
        camctrl.uninit_cam()
        # print("\x1b[38;2;10;250;50mCamera Closed...")

    if (ao_task is not None):
        daqctrl.set_outputs_to_constant(ao_task, [0,0,0])
        daqctrl.close_daq_task(ao_task)
        print('AO stopped...')
    
    if (trial_run[0] == 'n') and (SG is not None):
        SG.write('freq2.87ghz')
        # sgctrl.disable_SG_op()
        sgctrl.uninit_sg()

    # pbctrl.pb_init();
    if t_align_dc < 50:
        pbctrl.run_only_daq(250 *ms)        # let t_align_dc=250
    else:
        pbctrl.run_only_daq(t_align_dc *ms)
    closed = True
    # pbctrl.pb_stop();
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

    def is_folder_empty(folder_path):
        return len(os.listdir(folder_path)) == 0
    
    # ------------------------------------------------------------------
    # est_time = expCfg.t_tot * expCfg.Nsamples * expCfg.N_scanPts ??
    # check the scanned parameters for errors, if errors are found, remove those params
    if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr',
                               'drift_seq', 'timeseries_seq']:  # for sequences except ESR, set MW frequency
        seqArgList = [param[seq_no_plot[-1]]]
        seqArgList.extend(sequenceArgs)  # Make a dummy seqArgList just to create it
        [n_error, param] = seqctrl.param_err_check(instr, expCfg.sequence, expCfg.PBchannels, seqArgList,
                                                   expCfg.scannedParam, expCfg.N_scanPts)
        if n_error > 0:
            print('\x1b[1;37;41m' + 'Err: Check Sequences...\x1b[0m')
            print('\x1b[38;2;250;0;0m' + str(n_error) + '\x1b[0m parameters removed...')
            # Close Error plots??
            # close_plots = dialog.yesno_box('Close Plots', 'Close the Error Plots?')
            # if close_plots == 'yes':
            #     plt.close('all')
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
        # print(expCfg.PBchannels)
        view_sequence(expCfg.sequence, expCfg.PBchannels, seqArgList, False, param, seq_no_plot, seq_plot_dpi)

    # Start the initial sequence now ------------
    print("Starting Initial Sequence...")
    if trial_run[0] == 'n' or trial_run[1] == 'n':
        # instructionList = [start_initial_PB_seq()]
        # pbctrl.run_sequence_for_diode(instructionList)
        if t_align_dc < 50:
            instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, t_align_dc*ms/2],
                            [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align_dc*ms/2]]
        else:
            # t_align_dc is in ms
            duty = 0.07
            x = np.ceil(duty*t_align_dc/(1-duty)/10)*10       # in ms
            instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, x],
                            [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_dc*ms/2 - x)],
                            [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align_dc*ms/2]]
        pbctrl.run_sequence_for_diode([instructionList])
                
        print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
        if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq', 'aom_timing', 'rodelay', 'timeseries_seq'] and trial_run[0] == 'n':
            sgctrl.set_SG_freq(expCfg.MW_freq)
    else:
        print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    # Create the data save folder...
    print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')
    cwd = os.getcwd()  # current working directory - generally 'NV_Experiment' folder. Data folder is created a level up
    # rfind() gives the last occurence of '\\'
    savePath = cwd[0:cwd.rfind('\\')] + "\\Saved_Data\\" + time.strftime("%Y-%m-%d", time.localtime()) + '\\'
    if not (isdir(savePath)):
        os.makedirs(savePath)
    if expCfg.Nruns > 1:
        global f_number
        f_number = input("Enter folder number:: "+expCfg.sequence+"_#: ")
        # while isdir(savePath + expCfg.sequence + '_' + f_number + '\\'):
        #     f_number = input("Folder exists.. Re-enter folder number:: "+expCfg.sequence+"_#: ")
        #     # while not is_folder_empty(savePath + expCfg.sequence + '_' + f_number):
        #     #     f_number = input("Folder exists.. Re-enter folder number:: "+expCfg.sequence+"_#: ")

        while True:
            full_path = savePath + expCfg.sequence + '_' + f_number + '\\'
            
            if not isdir(full_path):
                os.makedirs(full_path)
                break
            elif not os.listdir(full_path):  # Check if the directory is empty
                break
            else:
                f_number = input("Folder exists and is not empty. Re-enter folder number:: "+expCfg.sequence+"_#: ")
        
        savePath = full_path
        if not (isdir(savePath)):
            os.makedirs(savePath)
    
    roi = [1228, 1344, 88, 52]
    roi = []
    
    # viewing live image and adjust the exposure time
    if trial_run[1] == 'n' and hdcamcon is not None:
        pbchannels = concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.bz
        # manage the duty cycle of MW
        instructionList = [[pbchannels ^ concfg.MW, Inst.CONTINUE, 0, 40 * ms],
                       [pbchannels, Inst.BRANCH, 0, 600 * ms]]
        pbctrl.run_sequence_for_diode([instructionList])
        ao_task = daqctrl.config_ao(dev='P6363')
        
        from_liveframes = camctrl.live_frames(ao_task, roi=roi, exposure=t_exposure, t_align=t_align_dc, field=align_field, running=False)
        pbctrl.run_only_daq(t_align_dc *ms)
    else:
        from_liveframes = [None, None, None, None, None, None]
    # now select the ROI
    if roi == [] and trial_run[1] == 'n':
        roi = camctrl.select_roi(data=from_liveframes[-1], roi=roi) if (trial_run[1] == 'n' and hdcamcon is not None) else None

    # stop the initial PB sequence after initializing the parameters
    # status = pbctrl.pb_stop()
    # pbctrl.errorCatcher(status)
    # instead of stopping PB, run the manipulation fields
    pbctrl.run_only_daq(t_align_dc *ms)
    
    # something wrong in the return statement, if and else both returns the same parameters!!!!!!!!!!!
    return [from_liveframes, roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] #if trial_run[1] == 'n' else [roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList]


def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[concfg.laser ^ concfg.MW, Inst.CONTINUE, 0, 500 * ms],
                       [concfg.laser, Inst.BRANCH, 0, 500 * ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...

    return instructionList


def read_save_details(roi, N_scanpts, Nsamples_expCfg):
    """defines the 'frame_per_cyc' to be captured, and the write formats of the 'scannedParam' and the 'data'
    Returns:
    frames_per_cyc            - define the number of frames captured in one cycle of the sequence, 1 signal, 1 ref, so 2. 
    scannedparam_write_format - write format of the scannedparams in the data file
    """
    frames_per_cyc = [2]  # signal and reference (2)
    Nsamples = frames_per_cyc[0] * Nsamples_expCfg

    # scannedparam_write_format = "%g\t" * N_scanpts
    # scannedparam_write_format = scannedparam_write_format[0:-1] + "\n"

    return [frames_per_cyc, Nsamples]


def extra_param_save_details(scan_time_list, exec_time, roi, exposure):
    extra_params_format = ' %s\t%0.2f\n %s\t%0.2f\n %s\t%0.4g\n %s\t%0.4g\n %s\t%g\n %s\t%d\n %s\t%d\n %s\t%d\n %s\t%d\n %s\t%f\n'  # %s\t%f\n
    # scan_time_list in ms; exec_time in seconds
    param_list = ['Max_scan_time(us):', max(scan_time_list) * 1e6, 'Min_scan_time(us):', min(scan_time_list) * 1e6,
                  'Total_scan_time(s):', np.sum(scan_time_list), 'Total_run_time(s):', exec_time, 'Step:',
                  (expCfg.scannedParam[1] - expCfg.scannedParam[0]), 'X0:', roi[0], 'Y0:', roi[1], 'W:', roi[2], 'H:',
                  roi[3], 'exposure (s):', exposure]
    # param_list[1] = i_scanpt+1..... Ekhane ki hbe?? Kon parameter save korbo??
    # ekhane ekta if kore, jodi 'i'=1 hoe, thle first param_list ta return korbe, nhle porer param_list ta return korbe.
    # Eta jodi kora hoe, thle runs>1 hole ba multi-param scan hole, notun param er sathe tar details save kora jabe...
    # Config file e ekta variable lagbe jeta dekhabe je kon variable ta scan hoechhe, other than scannedParam
    param_list = tuple(param_list)
    return [extra_params_format, param_list]


def prepare_for_saving(savePath):
    global f_number
    # if 'f_number' in vars():  # Returns a dict of all local variables
    #     print("Previous file number: \x1b[38;2;250;150;0m" + f_number + '\x1b[0m')

    while True:

        f_number = input("File name:" + expCfg.saveFileName + "_camera_#: ")
        datafilename = savePath + expCfg.saveFileName + "_camera_" + f_number + ".tiff"
        if (isfile(datafilename)):
            print('\x1b[38;2;250;250;0mFile exists. Retry...\x1b[0m')
            # continue
        else:
            break

    paramfilename = savePath + expCfg.saveFileName + "_camera_" + "params_" + f_number + ".txt"

    return [paramfilename, datafilename, f_number]


def save_data(datafilename, data, f_number):
    print("\x10 Saving to file.....")

    tifff.imwrite(datafilename, data, ome=True)
    # sk.io.imsave(datafilename, data, plugin='tifffile')

    print("\x10 Data saved to\x1b[38;2;100;250;50m %s_%s\x1b[0m !!!" % (expCfg.saveFileName, f_number))
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
        plt.figure(num=f"{sequence} sequence plot", dpi = plot_dpi)
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


def acquire_data(trigger_event, condition, roi, Nsamples, Nscanpts, i_run):
    # use starred expression here.. very difficult to generalize otherwise.....
    # setup next scan iteration (e.g. for ESR experiment, change microwave frequency; for T2 experiment, reprogram pulseblaster with new delay)

    scan_time_list = []  # time (in seconds) for each scannedParam
    global trial_run, data#, data_raw_time

    timeout_ms = 100000  # revisit...

    if trial_run[1] == 'n':

        for i_scanpt_cam in range(0, Nscanpts):

            trigger_event.wait()        # waiting for trigger to be set...

            # if i_scanpt_cam > 0:
            #     break
            
            # print('Starting wait_capevent_frameready...')
            t1 = time.perf_counter()
            # result = False 
            result = hdcamcon.wait_capevent_frameready(timeout_ms)

            timeout_happened = 0
            while result is not True:
            # if (result) is not True:
                # frame does not come
                if result != camctrl.dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                    print('-NG: Dcam.wait_event() failed with error {}'.format(result))
                    break

                # TIMEOUT error happens
                timeout_happened += 1
                if timeout_happened == 1:
                    print('Waiting for a frame to arrive.', end='')
                    if hdcamcon.get_propertyvalue(camctrl.dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == camctrl.dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
                        print(' Check your trigger source.', end='')
                    else:
                        print(' Check <timeout_ms>.', end='')
                    print(' Press Ctrl+C to abort.')
                else:
                    print('.')
                    if timeout_happened > 5:
                        timeout_happened = 0

                result = hdcamcon.wait_capevent_frameready(timeout_ms)
                # data_raw_time[i_run,i_scanpt_cam,0] = time.perf_counter()
                # continue
            if i_scanpt_cam>0:
                discard_frame = hdcamcon.get_lastframedata()
                result = hdcamcon.wait_capevent_frameready(timeout_ms)
            
            # 
            for i_sample in range(0, Nsamples):
                # wait_capevent_frameready() succeeded
                if i_sample+1 == Nsamples:
                    # Last frame for a scanpt:: stop PB and take out the last frame
                    status = pbctrl.pb_stop();    pbctrl.errorCatcher(status);
                    # pbctrl.run_only_daq(t_align_dc *ms)
                    frame = hdcamcon.get_lastframedata()
                    # print(i_sample, end='')
                    if frame is not False:
                        # print(' frame elo...')
                        # frames[i_sample, :, :] = frame
                        data[i_run,i_scanpt_cam,i_sample,:,:] = frame
                    else:
                        print(" No frame...")
                # Now the loop exits and clear the event, notify all waiting threads...
                else:
                    frame = hdcamcon.get_lastframedata()
                    # print(i_sample, end='')
                    if frame is not False:
                        # print(' frame elo...')
                        # frames[i_sample, :, :] = frame
                        data[i_run,i_scanpt_cam,i_sample,:,:] = frame
                    else:
                        print(" No frame...")
                    result = hdcamcon.wait_capevent_frameready(timeout_ms)
            
            # print('Single scan point end...')
            scan_start_time = time.perf_counter()
            # clear the event and notify all waiting threads...
            with condition:
                trigger_event.clear()   # clear the event
                condition.notify_all()  # notify the waiting threads, here the PB thread..
            scan_end_time = time.perf_counter()
            scan_time_list.append(scan_end_time - scan_start_time)

    else:
        # dcamcon.simulate_frame
        None

    return [data, scan_time_list, i_scanpt_cam]


def process_data(i_max, roi, n_frames, data):
    """Process data for plotting. Find the mean signal, reference and the contrast
    Note: This is done after the whole scan is performed."""
    hsize = roi[2];
    vsize = roi[3];
    mean_sig = np.zeros(i_max);
    mean_ref = np.zeros(i_max);
    contrast = np.zeros(i_max);
    for i_scanpt in range(0, i_max):
        # signal_frames = np.array([data[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize] for j in range(0,n_frames,2)])
        # reference_frames = np.array([data[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize] for j in range(1,n_frames,2)])
        # calculate the sum of pixels for each signal and ref frame and form an array
        # sum_of_signal_frames = np.array([np.sum(data[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize]) for j in range(0,n_frames,2)])
        # sum_of_reference_frames = np.array([np.sum(data[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize]) for j in range(1,n_frames,2)])

        signal_frames = data[i_scanpt, 0::2, :, :];
        reference_frames = data[i_scanpt, 1::2, :, :];

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
    xValues = param[1:i_max]
    if live:
        plt.plot([x / x_unit for x in xValues], contrast[1:i_max], 'b--')
        plt.show()
        plt.pause(0.0001)
    else:
        # plt.figure("Signal, Reference & Contrast Plots: %d points" %(i_max))
        # plt.figure()
        plt.subplot(121)
        plt.plot([x / x_unit for x in xValues], mean_sig[1:i_max], '.-', [x / x_unit for x in xValues], mean_ref[1:i_max], '.-')
        plt.legend(['Avg Sig', 'Avg Ref'])
        plt.xlabel(x_label)

        plt.subplot(122)
        plt.plot([x / x_unit for x in xValues], contrast[1:i_max], '.-b')
        plt.xlabel(x_label);
        plt.ylabel('Contrast')
        plt.show()


# ----------------------------------------------------------------------------
class PBThread(threading.Thread):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.trigger_event = kwargs.get('trigger_event')
        self.condition = kwargs.get('condition')
        self.t_seq_total = kwargs.get('t_seq_total')
        self.parameters = kwargs.get('parameters')
        self.sequence = kwargs.get('sequence')
        self.seqArgList = kwargs.get('seqArgList')
        self.t_exposure = args[0]
        self.N_total = args[1]
        self.args = args
        self.kwargs = kwargs
        self.inst_set_time_list = []
        self.time_taken = []
        # self.sg_time = 
    
    def run(self):
        """ Control the scanpoint loop and pulse sequence for experiment as well as camera trigger
    
        """
        global run_sequence, trial_run#, sequence
        Nscanpts = len(self.parameters)
        # inst_set_time_list = []
        # time.sleep(0.02)
        
        
        for i_scanpt in range(0, Nscanpts):
            t_seq_total_i = self.t_seq_total[i_scanpt]
            param = self.parameters[i_scanpt]
            print("\x1b[38;2;255;150;10m>", i_scanpt + 1, ' / ', Nscanpts, ': ', param, "\x1b[0m")
            
            # TODO: thread this part as well for parallel operation..
            start_time = time.perf_counter()
            instructionList = []
            if self.sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                if trial_run[0] == 'n':
                    sgctrl.set_SG_freq(param)
            else:
                self.seqArgList[0] = param
            
            the_list = pbctrl.PB_program(instr, self.sequence, self.seqArgList)
            # print('the_list: ', the_list)
            for i in range(0, len(the_list)):
                instructionList.append(the_list[i][0])
            # print('instructionList: ', instructionList)

            if trial_run[1] == 'n':
                # print('Starting PB sequence...')
                if instr == 'cam':
                    pbctrl.run_sequence_for_camera(instructionList, self.t_exposure, self.t_seq_total, self.N_total)
                elif instr == 'cam_level1':
                    pbctrl.run_sequence_for_camera_level_trigger_1(instructionList)
                elif 'levelm' in instr: # changed to 'levelm' in instr; prev instr == 'cam_levelm': to accomodate rot_field_controlled ODMR into the code
                    t_align = self.args[2]      # incoming time in [ns]
                    if instr == 'cam_levelm_trigger_ao':
                        pbctrl.run_sequence_for_camera_level_trigger_many_ac(instructionList, self.t_exposure, t_align, t_seq_total_i, self.N_total, expCfg.Nsamples)
                    else:
                        pbctrl.run_sequence_for_camera_level_trigger_many_dc(instructionList, self.t_exposure, t_align, t_seq_total_i, self.N_total)
                elif instr == 'cam_timeseries':
                    t_align_dc = self.args[2]
                    pbctrl.custom_trigger(self.t_exposure, t_align_dc)
                elif instr == 'cam_timeseries_trigger_ao':
                    t_align_rot = self.args[2]
                    t_measurement = self.args[3]
                    t_align_rot_extended = self.args[4]
                    pbctrl.custom_trigger_rot_field(self.t_exposure, t_align_rot, t_measurement, t_align_rot_extended)
                
                self.trigger_event.set()  # Set the event to trigger the camera
                self.inst_set_time_list.append((time.perf_counter()-start_time)*1e3)
                with self.condition:
                    # check whether the event is 'set'... if yes, wait inside the loop - no further execution...
                    while self.trigger_event.is_set():   # this statement dictates the wait condition: wait until the condition evaluates to False = event cleared..
                    # for 1st execution, this is False
                        self.condition.wait()        # waits for some event notification via notify() or notify_all()..
                        # whenever notified, the wait ends and proceeds again to check the condition
                    # In this case: it is event clear notification: the condition evaluates to False (the event is cleared = not_set) and below statements are executed...
                self.time_taken.append((time.perf_counter()-start_time)*1e3)
        # return [i_scanpt, inst_set_time_list]

def start_dc_alignment_field(align_data, t_align, n_repetitions):
    """Perform the dc alignment before each run to initialize the propeller.
    Perform a [125 ms (depending on cutoff frequency) +/- Bz, 125 ms Bx] sequence 5-10 times before starting the rotating field sequence. For this, turn the x, z channel ON and let the pulseblaster control the timing.
    """
    print("> Alignment field started!")
    print(f"Set t_align = {t_align}")
    
    if t_align < 50:
        instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW ^ concfg.start_trig, Inst.CONTINUE, 0, t_align*ms/2],
                        [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align*ms/2]]
    else:
        # t_align is in ms
        number_of_repetitions = int(np.ceil(np.ceil((t_align/2/1e3)/(np.max(align_data.shape)/samp_rate))/2)*2)
        pattern_time = np.max(align_data.shape)/samp_rate     # in [s]
        print(f"pattern time = {pattern_time}")

        # instructionList = [[concfg.laser ^ concfg.bz ^ concfg.start_trig, Inst.LOOP, 10, pattern_time*1e9],
        #                     [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align*ms/2-pattern_time*1e9)],
        #                     [concfg.laser ^ concfg.bx ^ concfg.by, Inst.END_LOOP, 0, (t_align*ms/2)]]
        if direction == 'z':
            instructionList = [[concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.start_trig, Inst.LOOP, 10, pattern_time*1e9],
                                [concfg.laser ^ concfg.bx ^ concfg.by, Inst.CONTINUE, 0, (t_align*ms/2-pattern_time*1e9)],
                                [concfg.laser ^ concfg.bz, Inst.END_LOOP, 0, (t_align*ms/2)]]
        elif direction == 'x':
            instructionList = [[concfg.laser ^ concfg.bz ^ concfg.start_trig, Inst.LOOP, 10, pattern_time*1e9],
                                [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align*ms/2-pattern_time*1e9)],
                                [concfg.laser ^ concfg.bx ^ concfg.by, Inst.END_LOOP, 0, (t_align*ms/2)]]
    pbctrl.run_sequence_for_diode([instructionList])
    # time.sleep(10)
    print("writing dc alignment data")
    ao_task.write(align_data)
    

if __name__ == '__main__':
    # global data
    warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
    
    expCfg = import_module(expCfgFile)
    # expCfg.N_scanPts = len(expCfg.scannedParam)
    t_AOM = expCfg.t_AOM
    clk_cyc = 1e3 / concfg.PBclk  # One clk cycle of PB = inverse of the clk freq of PB
    print('\x10 \x1b[38;2;250;250;0mRunning ' + expCfg.saveFileName + ' sequence\x1b[0m')
    [SG, hdcamcon] = initialize_instr(expCfg.sequence)

    seqctrl.check_params(expCfgFile)
    print("\x10 Parameter checks completed...")

    [from_liveframes, roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr,expCfg)
    if trial_run[1]=='n':
        [exit_code, ao_task, t_exposure, align_aovoltage, align_field, _] = from_liveframes
        print(f"t_exposure = {t_exposure} s")
        # t_exposure in seconds
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
    # print(t_exposure)
    if trial_run[1]=='n':
        N_total = [] if expCfg.sequence == 'esr_seq' else [int(np.floor(t_exposure*1e9/t)) for t in t_seq_total[0]]
    else:
        N_total = [] if expCfg.sequence == 'esr_seq' else [int(np.floor(t_exposure*1e6/t)) for t in t_seq_total[0]]
    N_total = []
    print(f"N_total = {N_total}")
    print(f"Set t_exposure = {t_exposure} [s]")
    # -------------------------- 20062023-------------------------
    # prepare data for AO - either DC field or rotating field
    if 'trigger_ao' in instr:   # this takes care of triggered timeseries + ODMR type experiments
        # prepare rotational alignment data pattern
        [dt, t_array, b_rot] = daqctrl.triggered_ao_data(direction, rot_angle, align_field, rot_field_amp, rot_field_freq)
        t_align_rot = np.max(b_rot.shape)*dt       # t_align_rot in seconds: total time of the rotational alignment pattern
        calib_factor = daqctrl.coil_calibration([1,1,1])                            # coil calibration factor
        if instr == 'cam_timeseries_trigger_ao':
            # for timeseries acquisition, we acquire data for the whole duration of measurement, even when the rotating field is ON
            # hence need to adjust the field ON time (t_align_rot) so that frame acquistion (syncrhonous readout mode) can fit in even numbers (sig+ref)
            # make the number of frames that can fit in the alignment time an even number
            n_frames_rot_alignment = np.ceil(np.ceil(t_align_rot/t_exposure)/2)*2       # UNITS???? t_exposure is seconds, and t_align_rot is in seconds (calculated from dt in seconds)
            t_align_rot_extended = n_frames_rot_alignment*t_exposure       # extended/modified t_align_rot to fit even number of frames (signals and references)
            t_zero = t_align_rot_extended - t_align_rot
            ti = t_array[-1] + dt
            temp = ti + np.array(np.arange(0,t_zero,dt))
            samp_rate = 1/dt
            bzero_padding = np.zeros((np.max(temp.shape), b_rot.shape[-1]))
            t_array = np.hstack((t_array, temp))

            b_rot_padded = np.ascontiguousarray(np.vstack((b_rot, bzero_padding)))      # field array with zero padding at the end
            v_rot_padded = np.zeros_like(b_rot_padded)
            for i in range(0, b_rot_padded.shape[-1]):
                v_rot_padded[:,i] = b_rot_padded[:,i]*calib_factor[i]                   # voltage array from the field array, includes the calibration factor..
                
            v_rot_padded_prepared = daqctrl.prepare_data_for_write(v_rot_padded)
            plt.figure(); plt.plot(t_array, b_rot_padded)       # t_array is required just for plotting purposes to verify the actual output
            plt.figure(); plt.plot(t_array, v_rot_padded)
            print(f"Expected t_align_rot = {t_align_rot} seconds")
            print(f"Extended t_align_rot = {t_align_rot_extended} seconds")
            print(f"#frames in alignment time (sig+ref) = {n_frames_rot_alignment}")

            v_rot_pattern = v_rot_padded_prepared.copy()
            # also prepare data for DC alignment
            # align_data_dc = daqctrl.prepare_data_for_write(np.array([daqctrl.coil_calibration(align_field) for _ in range(0, np.max(v_rot_pattern.shape))]))

            # instructionList = [[concfg.laser ^ concfg.bz ^ concfg.start_trig, Inst.LOOP, 10, t_align_rot_extended*1e9],
            #             [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_dc*ms/2-t_align_rot_extended*1e9)],
            #             # [concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_rot_extended*1e9],
            #             # [concfg.laser ^ concfg.bx ^ concfg.by, Inst.END_LOOP, 0, (t_align_dc*ms/2-t_align_rot_extended*1e9)]
            #             [concfg.laser ^ concfg.bx ^ concfg.by, Inst.END_LOOP, 0, (t_align_dc*ms/2)]]
            #             # [concfg.laser ^ concfg.bz ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_rot_extended*1e9],
            #             # [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_dc*ms/2-t_align_rot_extended*1e9)],
            #             # [concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_rot_extended*1e9],
            #             # [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 4, (t_align_dc*ms/2-t_align_rot_extended*1e9)]]
            # pbctrl.run_sequence_for_diode([instructionList])
            # time.sleep(10)
        
        elif instr == 'cam_levelm_trigger_ao':
            # for measurement case, the frame acquisition happens after the field is OFF, this is more like the other cases with slight modifications in the timescales for the rotating field
            v_rot = np.zeros_like(b_rot)
            for i in range(0, b_rot.shape[-1]):
                v_rot[:,i] = b_rot[:,i]*calib_factor[i]                   # voltage array from the field array, includes the calibration factor..
            
            v_rot_prepared = daqctrl.prepare_data_for_write(v_rot)
            v_rot_pattern = v_rot_prepared.copy()
        
        # else:       # handle cases when it is not triggered AO, ie only DC field measurement
    if trial_run[1] == 'n' and hdcamcon is not None:
        if instr == 'cam_timeseries_trigger_ao':
            expCfg.Nsamples = int(n_frames_rot_alignment/2 + np.ceil(t_meas/1e3 *fps/2))      # there is a division by 2 since there will be a multiplication to automatically consider the signal and references in read_save_details(); 200 = FPS
        [frames_per_cyc, Nsamples] = read_save_details(roi, Nscanpts, expCfg.Nsamples)

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

        # t_exposure = hdcamcon.get_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME)  # mind the value... in seconds !!!

    # status = pbctrl.pb_stop(); pbctrl.errorCatcher(status) already stopped...
    # pbctrl.custom_trigger(t_exposure *1e9, t_align_dc *1e6)
    closed = False
    display_parameters = dialog.yesno_box(expCfg.saveFileName + ' Params', 'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' % (str(concfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples, param[0] / expCfg.plotXaxisUnits, param[-1] / expCfg.plotXaxisUnits))
    
    if display_parameters == 'yes':
        if trial_run[1] == 'n':
            try:
                # continue_run = 'yes'
                # save_flag = False
                run_start_time = time.perf_counter()    # entire measurement time
                camctrl.configure_camera(instr)
                # camctrl.query_cam_values()
                # t_exposure = hdcamcon.setget_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME, val=t_exposure) # in [s]
                print(f"Actual exposure time = {hdcamcon.get_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME)*1e3:0.2f} ms")
                PB_wait_time = []
                # set buffer and start capture sequence..
                while not camctrl.startingcapture(size_buffer=Nsamples, sequence=True):
                    print("Retrying..")
                
                print(f"Set t_align_DC = {t_align_dc} seconds")
                
                # if 'trigger_ao' not in instr:
                # t_align_dc is in ms
                if t_align_dc < 50:
                    instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_dc*ms/2],
                                    [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align_dc*ms/2]]
                else:
                    duty = 0.07
                    x = np.ceil(duty*t_align_dc/(1-duty)/10)*10       # in ms
                    instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, x],
                                        [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_dc*ms/2-x)],
                                        [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, (t_align_dc*ms/2)]]
                    
                pbctrl.run_sequence_for_diode([instructionList])

                print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)    # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
                print('\x10 %d frames in each cycles...' % frames_per_cyc[0])       # hence, total # frames = frames_per_cyc[0]*expCfg.Nsamples = Nsamples
                hsize = roi[2]; vsize = roi[3]
                print("------Acquiring %dx%d------" % (hsize, vsize))

                data = np.zeros((expCfg.Nruns, Nscanpts, Nsamples, vsize, hsize), dtype='uint16')
                # data_raw_time = np.zeros((expCfg.Nruns, Nscanpts, Nsamples), dtype='float64')

                for i_run in range(0, expCfg.Nruns):     # TODO: add Nruns loop as well
                    print("\x1b[38;2;255;150;10mRun: ", i_run + 1, ' / ', expCfg.Nruns, "\x1b[0m")
                    print('Starting threads...')
                    trigger_event = threading.Event()       # Event object to signal between threads
                    condition = threading.Condition()
                    
                    kw_args = {'t_seq_total': t_seq_total, 'parameters': param, 'sequence': expCfg.sequence, 'seqArgList': seqArgList, 'trigger_event': trigger_event, 'condition': condition}
                    # Create and start the scan thread: t_seq_total, parameters, sequence, seqArgList, trigger_event, condition
                    acq_start_time = time.perf_counter()
                    # print(f"t_align_rot = {t_align_rot}")
                    # print(f"t_align_rot_extended = {t_align_rot_extended}")
                    # print(f"t_meas = {t_meas}")
                    if instr == 'cam_timeseries_trigger_ao':
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, t_meas *1e6, t_align_rot_extended *1e9, **kw_args)
                    elif instr == 'cam_levelm_trigger_ao':
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, **kw_args)
                    else:
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_dc *1e6, **kw_args)

                    # scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, t_meas *1e6, t_align_rot_extended *1e9, **kw_args) if instr == 'cam_timeseries_trigger_ao' else PBThread(t_exposure *1e9, N_total, t_align_dc *1e6, **kw_args)     # there are >2 cases now, hence if-else block is used..
                    # camera_thread = threading.Thread(target=acquire_data, args=(trigger_event, condition, roi, Nsamples, Nscanpts))
                    # camera_thread.start()

                    # # DC alignment at the beginning of each run if it is a triggered AO based manipulation
                    # # start DC field
                    # output_aovoltage = daqctrl.coil_calibration(align_field)
                    # daqctrl.set_outputs_to_constant(ao_task, output_aovoltage)
                    # if t_align_dc < 50:
                    #         instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_dc*ms/2],
                    #                         [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align_dc*ms/2]]
                    # else:
                    #     duty = 0.07
                    #     # x = np.ceil(duty*t_align_dc/(1-duty)/10)*10       # in ms
                    #     # instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 15, x],
                    #     #                 [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_dc*ms/2-x)],
                    #     #                 [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, (t_align_dc*ms/2)]]
                    #     instructionList = [[concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.MW, Inst.CONTINUE, 15, x],
                    #                     [concfg.laser ^ concfg.bx ^ concfg.by, Inst.CONTINUE, 0, (t_align_dc*ms/2-x)],
                    #                     [concfg.laser ^ concfg.bz, Inst.BRANCH, 0, (t_align_dc*ms/2)]]
                    # pbctrl.run_sequence_for_diode([instructionList])
                    # time.sleep(15*t_align_dc/1e3)

                    if 'trigger_ao' in instr:
                        daqctrl.create_retriggerable_ao_task(ao_task, v_rot_pattern)
                        # if instr == 'cam_timeseries_trigger_ao':  # this may not be required now since the variable 'v_rot_pattern' is generalized from 'v_rot_padded_prepared'
                        print("Starting alignment pattern!")
                        daqctrl.start_retriggerable_ao_task(ao_task, v_rot_pattern, dt)
                        # time.sleep(0.04)
                    scan_thread.start()
                    
                    # Main thread handling camera acquisition
                    # acq_start_time = time.perf_counter()
                    [data, scan_time_list, i_scanpt_cam] = acquire_data(trigger_event, condition, roi, Nsamples, Nscanpts, i_run)
                    
                    # acq_time = time.perf_counter() - acq_start_time

                    # Wait for the pulse generation thread to complete
                    scan_thread.join()
                    # camera_thread.join()
                    # [i_scanpt_inst, inst_set_time_list] = scan_thread.result
                    acq_time = time.perf_counter() - acq_start_time

                    # check if i_scanpt from cam and inst are same..
                    # if i_scanpt_inst == i_scanpt_cam:
                    # i_scanpt = i_scanpt_cam
                    i_scanpt = Nscanpts
                    # else:
                    #     print("Error in acquisition: i_scanpt_inst != i_scanpt_cam!!")
                    # if expCfg.Nruns > 1:
                    #     camctrl.live_frames(ao_task, roi=[], exposure=t_exposure, t_align_dc=t_align_dc, field=align_field, running=True)
                    
                # all runs complete.. process for display...
                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print(f"Run time = {run_time:0.2f} s")
                plt.figure(num=time.strftime(" [%H:%M:%S]", time.localtime()))
                for i_run in range(0, expCfg.Nruns):
                    processed_data = process_data(i_scanpt, roi, Nsamples, data[i_run,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                    plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples, live=False)

                # camctrl.query_cam_values()

                # Close all after full acquisition
                closed = close_all(SG, hdcamcon, ao_task)

                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Data file saving
                    if expCfg.Nruns == 1:
                        if save_flag == False:  # Ask for file number iff no save was performed
                            [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
                        save_flag = save_data(datafilename, data[0,:,:,:,:], f_number)
                    else:
                        for i_run in range(0, expCfg.Nruns):
                            datafilename = savePath + expCfg.saveFileName + "_camera_" + str(i_run) + ".tiff"
                            save_flag = save_data(datafilename, data[i_run,:,:,:,:], i_run)

                else:
                    print("\x10 \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")
                # if (i_run + 1) < expCfg.Nruns:
                #     continue_run = dialog.yesno_box('Continue', "Continue Run?")
                #     if continue_run == 'yes':
                #         continue
                #     else:
                #         sys.exit("Run Interrupted...")
                # time_list = 

            except KeyboardInterrupt:
                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
                # Plot the last run data upto the point where it was interrupted...
                processed_data = process_data(i_scanpt, roi, Nsamples, data[0,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples, live=False)
                # Then ask whether to save it...
                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Ask for filename only if there was no save operation
                    if save_flag == False:
                        [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
                    save_flag = save_data(datafilename, data[0,:,:,:,:], f_number)
            # nicher ei 'else' ta kno dewa hoyechhilo??? 18-06-2023
            # else:
            #     if (i_run+1)==expCfg.Nruns: # Save data of the final run here
            #         save_data(len(concfg.input_terminals), datafilename, data, expCfg.Nsamples, reads_per_cyc)
            finally:
                if not closed:
                    closed = close_all(SG, hdcamcon, ao_task)

                print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))

                print("data = \x1b[38;2;250;150;50m%d x %d x %d x %d x%d\x1b[0m" % np.shape(data[:,:,:,:,:]))
                # Save parameters (only if there was one save operation)...
                if save_flag:
                    expParamList[1] = i_scanpt + 1  # expParamList[1] -> value of N_scanPts
                    expParamList[3] = i_run + 1  # expParamList[3] -> value of Nruns

                    if expCfg.Nruns>1:
                        paramfilename = savePath+"params_"+str(f_number)+".txt"

                    save_parameters(paramfilename, param_save_format, expParamList)

                    [next_params_format, expParamList] = extra_param_save_details(scan_time_list, run_time, roi, t_exposure)
                    save_parameters(paramfilename, next_params_format, expParamList)

                # sys.exit(exit_code)       # the exit code is the one returned from the pyqt app
                # scan_time_list in seconds; exec_time in seconds
                # ----------------------------------------------------------------------

                # matlab_engine = eng.start_matlab()

                # if expCfg.sequence == 'esr':
                #     analysis = eng.ESR_raw_analysis()
        elif trial_run[1] == 'y':
            # if ao_task is not None and ao_task.is_task_done() == False:
            #     print("Task is already initialized and running. No need to reinitialize.")
            # else:
            ao_task = daqctrl.config_ao("P6363")
            print("> Task defined...")
            
            # DC alignment at the beginning of each run if it is a triggered AO based manipulation
            # start DC field
            output_aovoltage = daqctrl.coil_calibration(align_field)
            daqctrl.set_outputs_to_constant(ao_task, output_aovoltage)
            if t_align_dc < 50:
                instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_dc*ms/2],
                                [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align_dc*ms/2]]
            else:
                duty = 0.07
                x = np.ceil(duty*t_align_dc/(1-duty)/10)*10       # in ms
                # instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 10, x],
                #                 [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_dc*ms/2-x)],
                #                 [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, (t_align_dc*ms/2)]]
                instructionList = [[concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.MW, Inst.CONTINUE, 15, x],
                                        [concfg.laser ^ concfg.bx ^ concfg.by, Inst.CONTINUE, 0, (t_align_dc*ms/2-x)],
                                        [concfg.laser ^ concfg.bz, Inst.BRANCH, 0, (t_align_dc*ms/2)]]
            pbctrl.run_sequence_for_diode([instructionList])
            time.sleep(10*t_align_dc/1e3)
            if 'trigger_ao' in instr:
                daqctrl.create_retriggerable_ao_task(ao_task, v_rot_pattern)
                if instr == 'cam_timeseries_trigger_ao':
                    print("Starting alignment pattern!")
                    daqctrl.start_retriggerable_ao_task(ao_task, v_rot_pattern, dt)
                print("> Alignment field started!")
            
            instructionList = []
            if expCfg.sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                None
            else:
                seqArgList[0] = param[seq_no_plot[-1]]
            # print(seqArgList)
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
                # print(t_exposure)
                pbctrl.run_sequence_for_camera_level_trigger_many(instructionList, t_exposure * 1e9, t_align_dc * 1e6, t_seq_total[seq_no_plot[-1]], N_total)      # t_exposure [s], t_align_dc [ms] at the beginning...
            # pbctrl.pb_close(); print("\x10 Pulse Blaster closed...\x1b[0m")
            elif instr == 'cam_timeseries':
                pbctrl.custom_trigger(t_exposure *1e9, t_align_dc *1e6)
            elif instr == 'cam_timeseries_trigger_ao':
                print("Run PB continuously...")
                pbctrl.custom_trigger_rot_field_continuous_uncorrected(t_exposure *1e9, t_align_rot *1e9, t_meas *1e6, t_align_rot_extended *1e9)
                # pbctrl.custom_trigger_rot_field_continuous(t_exposure *1e9, t_align_rot *1e9, t_meas *1e6, t_align_rot_extended *1e9)
                # pbctrl.custom_trigger_rot_field(t_exposure *1e9, t_align_rot_extended *1e9, t_meas*1e6, t_align_rot_extended *1e9)
            
            # time.sleep(4)
            # instructionList = [[concfg.start_trig, Inst.CONTINUE, 0, 1*ms],
            #                    [0, Inst.BRANCH, 0, 1*ms]]
            # pbctrl.run_sequence_for_diode([instructionList])
            # daqctrl.set_outputs_to_zero(ao_task)
            # ao_task.close()
            pbctrl.pb_close()
        
    else:
        print("Measurement aborted...")
        if not closed:
            closed = close_all(SG, hdcamcon, ao_task)
        
    # return instructionList

    # %%
