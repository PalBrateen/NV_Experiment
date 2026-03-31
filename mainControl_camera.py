# mainControl_camera -v2.0
# Initialization and Definition
# %reset
# stylesheets @ https://claude.ai/chat/4727170d-279a-4468-be67-3b6a9ac3fdc1
import  connectionConfig as concfg, matplotlib.pyplot as plt, numpy as np, sys, time, dialog, cv2, os, \
    tifffile as tifff, threading, logging, dcamcon, psutil
from PyQt5.QtWidgets import QApplication
from sequencecontrol import sequencecontrol
from SGcontrol import SignalGenerator, SignalGenerator_sim
from Camcontrol import CameraThread, CameraWorker, InputApp
from PBcontrol import PulseBlaster, ns, ms, us, s, Inst
from DAQcontrol import DAQ_write_pattern, AnalogOutputTask
# from spinapi import ns, us, ms, Inst
from os.path import isdir, isfile;
from importlib import import_module
import dcam, warnings, matplotlib
from queue import Queue, Empty
# from matlab import engine as eng
# %matplotlib qt5
plt.rcParams.update({'figure.max_open_warning': 0})

global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, sg, expCfg, f_number, hdcamcon, instr, data, pb, ao_task, camera_worker#, data_raw_time

instr = 'cam_levelm'
# options: cam, cam_level1, cam_levelm, cam_syncm, cam_timeseries, cam_timeseries_trigger_ao, cam_syncm_trigger_ao_ac, cam_syncm_trigger_ao_dc,  cam_levelm_trigger_ao_ac, cam_levelm_trigger_ao_dc,  cam_levelm_no_trigger_ao = cam_levelm
# WARNING: for sync-trigger, check the acquire_data() for discarding 1st frame
expCfgFile = 'rabi' + '_config'
# N_total = [2, 2]  # a 2-element list (sig and ref) for total number of repetitions.. if empty then N_total is allowed to change for each scanpt..
# N_total = [] #if expCfgFile == 'esr_config' else [12354,12354]

trial_run = ['n', 'n']  # 1st = sg, 2nd = camera
# all times here are in seconds
fps = 996.3
rot_field_amp = 50      # field amplitude in [gauss]
rot_field_freq = 20     # field frequency [Hz]

direction = 'z'
rot_angle = 60      # rotation angle in degrees
t_meas = 50        # measurement time [ms] for cam_timeseries_trigger_ao only

t_exposure = 20     # exposure [ms]; set to 5 ms for 'timeseries' measurements
align_field = [ 50 , 50, 50 ]          # [bx, by, bz] G
align_field = [ 0 , 0 , 0 ]

t_align_dc = 200         # [ms]; set to 250 ms for 'timeseries' measurements (time: Bz + Bx)
# propeller orientation
prop_theta = 0      # [degrees]
prop_phi = 0        # [degrees]

# TODO: Put the test field in the GUI for exposure control and return the values
# test_field = 10; test_theta = 0; test_phi = 0
# test_field = [-14.1*4, -20, -10*1]        # [bx, by, bz] G
test_field = [0, 0, 0]
min_focus_time = 5      # [s]

# cam_levelm: level trigger measurements
# cam_timeseries_trigger_ao: take the timeseries measurement at 5ms exposure time for some time interval after performing the 3x(Bz, Bx) + fractional rotating field alignment
# cam_syncm_trigger_ao_ac: perform ODMR with triggered rotating (ac) field for diffusion control with sync-readout triggered camnera acquisition
# cam_syncm_trigger_ao_dc

t_exposure /= 1e3       # [s]
if prop_phi == 0:
    align_field[1] = 0
elif prop_phi == 90:
    align_field[0] = 0
    
# ------------- Plotting details----------------------------
seq_no_plot = [0]
voltage_unit = 1  # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100  # The dpi of the displayed pulse sequence plot
plotPulseSequence = False
livePlotUpdate = False
ao_task = None


# TODO: set up the logger
fsplit = lambda b: (2870 - 2.8*b, 2870 + 2.8*b)
# test_field = test_field*np.array([np.sin(test_theta*np.pi/180)*np.cos(test_phi*np.pi/180), np.sin(test_theta*np.pi/180)*np.sin(test_phi*np.pi/180), np.cos(test_theta*np.pi/180)])

seqctrl = sequencecontrol()
def initialize_instr(sequence):
    global hdcamcon, pb, ao_task, camera_worker, sg
    # TODO: how to handle 'sequence' to PulseBlaster()
    pb = PulseBlaster()              # pb not configured in __init__()
    try:
        pb.configure()
    except Exception as e:
        print(f"Error PB Init: {e}")
    try:
        ao_task = AnalogOutputTask()    # AO alreay configured here in __init__()
    except Exception as e:
        print(f"Error DAQ Init: {e}")
    try:
        if sequence not in ['aom_timing', 'rodelay']: #trial_run[0] == 'n' and 
            # Do not initialize sg if it is a trial run or the sequence is present in the list ['aom_timing', 'rodelay']
            sg = SignalGenerator() if trial_run[0] == 'n' else SignalGenerator_sim()
            if sg != '':
                sg.enable_sg_output();
                print("SG Output Enabled...")
                sg.set_sg_amp(expCfg.MW_power)
                sg.set_sg_freq(2.87e9)
                sg.setup_sg_pulse_mod();
                print("SG Ext Pulse Mod Enabled...")
        # elif trial_run[0] == 'y':
        # sg = SignalGenerator_sim(simulate=True)
    except Exception as e:
        print(f"Error SG Init: {e}")
    # if trial_run[1] == 'n':
    try:
        camera_worker = CameraWorker(simulate=(trial_run[1]=='y'))
        # print(camera_worker.roi)
        camera_worker.init_cam()
        hdcamcon = camera_worker.hdcamcon
        # this is incorrect and hard to debug.. Camcontrol already has a configure_camera() method..
        # Keeping for now.. Change later
        if 'level' in instr:
            camera_worker.triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.LEVEL
        elif 'sync' in instr:
            camera_worker.triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.SYNCREADOUT
        else:
            sys.exit("Invalid instrument name...")
        hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE, val=camera_worker.triggeractive)
        # else:
        #     hdcamcon = None
    except Exception as e:
        print(f"Error Camera Init: {e}")
    return [sg, ao_task, hdcamcon]

def set_core_affinity():
    # Assuming you want cores 0,1 for instrument control
    process = psutil.Process(os.getpid())
    process.cpu_affinity([0, 1, 2])
    print(f"Process pinned to cores: {process.cpu_affinity()}")

def close_all(sg, hdcamcon, ao_task):
    """End the measurement
    Closes all the instruments.
    """
    try:
        # if (trial_run[1] == 'n') and (hdcamcon is not None):
            # stops capture, releases buffer, closes camera and uninitializes DCAM-API
        camera_worker.uninit_cam()
            # print("\x1b[38;2;10;250;50mCamera Closed...")

        if (ao_task is not None):
            print(f"AO closing:: status: {ao_task.task_state}")
            ao_task.set_outputs_to_constant([0,0,0])
            time.sleep(0.5)
            ao_task.__del__()
            print('AO stopped...')
        
        # if (trial_run[0] == 'n') and (sg is not None):
        sg.set_sg_freq(2.87e9)
        # sg.disable_sg_output()
        sg.uninit_sg()

        # pb.pb_init();
        # if t_align_dc < 50:
        #     pb.run_only_daq(250 *ms)        # let t_align_dc=250
        # else:
        #     pb.run_only_daq(t_align_dc *ms)
        closed = True
        pb.stop_sequence();
        pb.closePB();
        print("Pulse Blaster closed...\x1b[0m")
        return True
    except Exception as e:
        print(f"Error: {e}")

def measurement_to_focus(t_exposure, t_rot_alignment):
    """Necessary operations to transition from measurement state to focusing state
    """

    # PB settings
    pb.focus_adjustment_sequence(t_exposure=t_exposure, t_align=t_rot_alignment)
    # while camera_worker.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
    #     pb.start_sequence()
    #     time.sleep(pb.min_focus_time/s)
    #     pb.stop_sequence()

    # camera settings
    camera_worker.trigger_mode = dcamcon.DCAMPROP.TRIGGERSOURCE.INTERNAL
    camera_worker.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.TRIGGERSOURCE, camera_worker.trigger_mode) if not camera_worker.simulate else None

    # # DAQ AO settings
    # ao_task.set_outputs_to_constant(align_field)

def focus_to_measurement():
    """Necessary operations to transition from focusing state to measurement state
    """
    # camera settings
    camera_worker.configure_camera()

    # # DAQ AO settings
    # ao_task.create_retriggerable_ao_task(data_shape)
    # ao_task.start_retirggerable_ao_task(patern_data)

    # PB settings
    None

def focus_adjustment(t_exposure, t_rot_alignment):
    # transit from measurement to focus
    measurement_to_focus(t_exposure, t_rot_alignment)

    # start focusing
    print(f"Before adjustment.. {camera_worker.query_camera_status()}")
    camera_worker.start_capture()
    print(f"After adjustment.. {camera_worker.query_camera_status()}")

    # while camera_worker.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
    #     pass

    # camera_worker.stop_capture()
    focus_to_measurement()


def initialize_exp(instr, expCfg):
    """Initialize the experimental parameters and check whether the sequence can be sent to the PB...
    Also, start the initial PB sequence.

    Returns:
    -------
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
            # param = np.array([0, 200] + [202]*100).astype(np.float64)
            
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
    # if trial_run[0] == 'n':# or trial_run[1] == 'n':
    # instructionList = [start_initial_PB_seq()]
    # pb.run_sequence_for_diode(instructionList)
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
    pb.run_sequence_for_diode([instructionList])
            
    print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
    if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq', 'aom_timing', 'rodelay', 'timeseries_seq']:# and trial_run[0] == 'n':
        sg.set_sg_freq(expCfg.MW_freq)
    else:
        # print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")
        pass

    # Create the data save folder...
    print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')
    cwd = os.getcwd()  # current working directory - generally 'NV_Experiment' folder. Data folder is created a level up
    # rfind() gives the last occurence of '\\'
    savePath = cwd[0:cwd.rfind('\\')] + "\\Saved_Data\\" + time.strftime("%Y-%m-%d", time.localtime()) + '\\'
    if not (isdir(savePath)):
        os.makedirs(savePath)
        
    # if expCfg.Nruns > 1:
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
    
    roi = [968, 752, 644, 808]      # region of laser spot
    roi = [980, 1040, 48, 44]
    # roi = [1016, 1204, 88, 84]          # a part of diamond
    roi = [1020, 1124, 236, 240]          # part of diamond in confocal
    roi = [1116, 908, 52, 52]
    roi = []
    
    # viewing live image and adjust the exposure time
    # if trial_run[1] == 'n' and hdcamcon is not None:
    pbchannels = concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.bz
    # manage the duty cycle of MW
    instructionList = [[pbchannels ^ concfg.MW, Inst.CONTINUE, 0, 40 * ms],
                    [pbchannels, Inst.BRANCH, 0, 600 * ms]]
    pb.run_sequence_for_diode([instructionList])
    # ao_task = daqctrl.config_ao(dev='P6363')
    
    app = QApplication(sys.argv)
    camera_app = InputApp(ao_task, camera_worker, roi, t_exposure, t_align_dc, align_field)
    camera_app.show()
    exit_code = app.exec_()

    # from_liveframes = camera_worker.live_frames(ao_task, roi=roi, exposure=t_exposure, t_align=t_align_dc, field=align_field, running=False)
    from_liveframes = [camera_worker.exposure, camera_app.align_field, camera_worker.last_frame]
    # [camera_worker_obj.exposure, camera_worker_obj.align_voltage, camera_worker_obj.align_field, camera_worker_obj.last_frame]
    pb.run_only_daq(t_align_dc *ms)
    # else:
    #     from_liveframes = [None, None, None]
    # now select the ROI -- -- -- --> now inside camera_worker
    # if roi == [] and trial_run[1] == 'n':
    #     roi = camera_worker.select_roi(data=from_liveframes[-1], roi=roi) if (trial_run[1] == 'n' and hdcamcon is not None) else None
    # if trial_run[1] == 'n':
    roi = camera_worker.roi

    # stop the initial PB sequence after initializing the parameters
    # status = pb.pb_stop()
    # pb.errorCatcher(status)
    # instead of stopping PB, run the manipulation fields
    pb.run_only_daq(t_align_dc *ms)
    
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
        the_list = pb.PB_program(instr, sequence, seqArgList)
        for j in range(0, len(the_list)):
            instructionList = the_list[j][0]
            inst_times = the_list[j][4]

            # if instr in ['cam', 'cam_levelm', 'cam_syncm_trigger_ao', 'cam_levelm_trigger_ao']:   # this is generalized below...
            if 'level' in instr or 'sync' in instr:
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
    def wait_for_frame(timeout_ms, timeout_happened = 0):
        frame_ready = hdcamcon.wait_capevent_frameready(timeout_ms)
        timeout_happened = 0
        if frame_ready is not True:
            # frame does not come
            if frame_ready != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                print('\x1b[38;2;255;20;10m-NG: Dcam.wait_event() failed with error {}\x1b[0m'.format(frame_ready))
                # break
            # else TIMEOUT error happens
            timeout_happened += 1
            if timeout_happened == 1:
                print('\x1b[38;2;255;20;10mWaiting for a frame to arrive.', end='')
                if hdcamcon.get_propertyvalue(dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
                    print(' Check your trigger source.', end='')
                else:
                    print(' Check <timeout_ms>.', end='')
                print(' Press Ctrl+C to abort.\x1b[0m')
            else:
                print('.\x1b[0m')
                if timeout_happened > 5:
                    timeout_happened = 0
            # recursive_wait_for_frame(timeout_ms, timeout_happened)

        return frame_ready

    scan_time_list = []  # time (in seconds) for each scannedParam
    global trial_run, data#, data_raw_time

    timeout_ms = 1000#int(t_exposure*1e3+20)  # revisit...

    if trial_run[1] == 'n':

        for i_scanpt_cam in range(0, Nscanpts):
            # print("Camera: waiting for trigger...")
            trigger_event.wait()        # waiting for trigger to be set...
            # print("Camera: receives event trigger...")
            # print("PB started...")
            pb.start_sequence()
            t1 = time.perf_counter()
            frame_ready = wait_for_frame(timeout_ms=timeout_ms)
            # print(f"frame_ready={frame_ready}")
            # frame_ready = hdcamcon.wait_capevent_frameready(timeout_ms)
            # if frame_ready is not True:
                # frame_ready = wait_for_frame(timeout_ms=timeout_ms)

            # timeout_happened = 0
            # while frame_ready is not True:
            #     print("Timeout...")
            # # if (frame_ready) is not True:
            #     # frame does not come
            #     if frame_ready != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
            #         print('-NG: Dcam.wait_event() failed with error {}'.format(frame_ready))
            #         break

            #     # TIMEOUT error happens
            #     timeout_happened += 1
            #     if timeout_happened == 1:
            #         print('Waiting for a frame to arrive.', end='')
            #         if hdcamcon.get_propertyvalue(dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
            #             print(' Check your trigger source.', end='')
            #         else:
            #             print(' Check <timeout_ms>.', end='')
            #         print(' Press Ctrl+C to abort.')
            #     else:
            #         print('.')
            #         if timeout_happened > 5:
            #             timeout_happened = 0
            #     print("Waiting for frame...")
            #     frame_ready = hdcamcon.wait_capevent_frameready(timeout_ms)
            #     # data_raw_time[i_run,i_scanpt_cam,0] = time.perf_counter()
            #     # continue
            # -----------------------------------------------
            # # activate the below lines for sync-readout trigger ONLY
            # # frame discard is not required for level trigger since exposusre is direclty controlled by PB
            # if i_run>0 or (i_run==0 and i_scanpt_cam>0):
            #     discard_frame = hdcamcon.get_lastframedata()
            #     if discard_frame is False:
            #         print("No frame...")
            #         sys.exit()
            #     frame_ready = hdcamcon.wait_capevent_frameready(timeout_ms)
            # ------------------------------------------------
            # print(f"Cam Nsamples = {Nsamples}")
            for i_sample in range(0, Nsamples):
                # print(i_sample, end=' ')
                if frame_ready is True:     # wait_capevent_frameready() succeeded
                    if i_sample+1 == Nsamples:
                        # Last frame for a scanpt:: stop PB and take out the last frame
                        # pb.run_only_daq(t_align_dc *ms)
                        pb.stop_sequence()      # may / may not stop the sequence, PBcontrol already takes care..
                        # print("PB stopped EXT..")
                        frame = hdcamcon.get_lastframedata()
                        # pb.stop_sequence()
                        if frame is not False:
                            # print('frame elo...')
                            # frames[i_sample, :, :] = frame
                            data[i_run,i_scanpt_cam,i_sample,:,:] = frame
                            # pb.stop_sequence()
                        else:
                            print("\x1b[38;2;255;20;10mNo frame...\x1b[0m")
                    # Now the loop exits and clear the event, notify all waiting threads...
                    else:
                        frame = hdcamcon.get_lastframedata()
                        if frame is not False:
                            # print('frame elo...')
                            # frames[i_sample, :, :] = frame
                            data[i_run,i_scanpt_cam,i_sample,:,:] = frame
                        else:
                            print("\x1b[38;2;255;20;10mNo frame...\x1b[0m")
                        # frame_ready = hdcamcon.wait_capevent_frameready(timeout_ms)
                        frame_ready = wait_for_frame(timeout_ms=timeout_ms)
                        # print(f"frame_ready={frame_ready}")
                else:
                    print("\x1b[38;2;255;20;10mTimeout: No frame...\x1b[0m")
                    frame_ready = wait_for_frame(timeout_ms=timeout_ms)
                    # print(f"frame_ready={frame_ready}")
            
            # print('Single scan point end...')
            scan_start_time = time.perf_counter()
            # clear the event and notify all waiting threads...
            with condition:
                # print("Camera: clears trigger event...")
                trigger_event.clear()   # clear the event
                condition.notify_all()  # notify the waiting threads, here the PB thread..
                # print("Camera: notified PB for event clear...")
            scan_end_time = time.perf_counter()
            scan_time_list.append(scan_end_time - scan_start_time)

    else:
        # dcamcon.simulate_frame
        None

    return [data, scan_time_list, i_scanpt_cam]

def acquire_data_sim(trigger_event, condition, roi, Nsamples, Nscanpts, i_run):
    scan_time_list = []  # time (in seconds) for each scannedParam
    global data#, data_raw_time

    timeout_ms = 1000  # revisit...

    for i_scanpt_cam in range(0, Nscanpts):
        # print("Camera: waiting for trigger...")
        trigger_event.wait()        # waiting for trigger to be set...
        # print("Camera: receives event trigger...")
        # print("PB started...")
        pb.start_sequence()
        t1 = time.perf_counter()
        frame_ready = True
        # print(f"frame_ready={frame_ready}")

        if i_scanpt_cam>0:
            time.sleep(t_exposure)
            discard_frame = (np.random.rand(*camera_worker.roi[-2:])*(2**16)).astype(np.uint16).transpose()
            if discard_frame is False:
                # print("No frame...")
                sys.exit()
            frame_ready = True #hdcamcon.wait_capevent_frameready(timeout_ms)
        # ------------------------------------------------
        # print(f"Cam Nsamples = {Nsamples}")
        for i_sample in range(0, Nsamples):
            # print(i_sample, end=' ')
            if frame_ready is True:     # wait_capevent_frameready() succeeded
                if i_sample+1 == Nsamples:
                    # Last frame for a scanpt:: stop PB and take out the last frame
                    # pb.run_only_daq(t_align_dc *ms)
                    time.sleep(t_exposure)
                    frame = (np.random.rand(*camera_worker.roi[-2:])*(2**16)).astype(np.uint16).transpose()
                    # pb.stop_sequence()
                    if frame is not False:
                        # print('frame elo...')
                        # frames[i_sample, :, :] = frame
                        data[i_run,i_scanpt_cam,i_sample,:,:] = frame
                        # pb.stop_sequence()
                    else:
                        print("\x1b[38;2;255;20;10mNo frame...\x1b[0m")
                    # print("PB stopped EXT..")
                    pb.stop_sequence()      # may / may not stop the sequence, PBcontrol already takes care..
                # Now the loop exits and clear the event, notify all waiting threads...
                else:
                    time.sleep(t_exposure)
                    frame = (np.random.rand(*camera_worker.roi[-2:])*(2**16)).astype(np.uint16).transpose()
                    if frame is not False:
                        # print('frame elo...')
                        # frames[i_sample, :, :] = frame
                        data[i_run,i_scanpt_cam,i_sample,:,:] = frame
                    else:
                        print("\x1b[38;2;255;20;10mNo frame...\x1b[0m")
                    frame_ready = True #hdcamcon.wait_capevent_frameready(timeout_ms)
                    # print(f"frame_ready={frame_ready}")
            else:
                print("\x1b[38;2;255;20;10mTimeout: No frame...\x1b[0m")
                frame_ready = True #hdcamcon.wait_capevent_frameready(timeout_ms)
                # print(f"frame_ready={frame_ready}")
        
        # print('Single scan point end...')
        scan_start_time = time.perf_counter()
        # clear the event and notify all waiting threads...
        with condition:
            # print("Camera: clears trigger event...")
            trigger_event.clear()   # clear the event
            condition.notify_all()  # notify the waiting threads, here the PB thread..
            # print("Camera: notified PB for event clear...")
        scan_end_time = time.perf_counter()
        scan_time_list.append(scan_end_time - scan_start_time)

    return [data, scan_time_list, i_scanpt_cam]

# TODO: Something takes a lot of time: either processing or plotting.. Inspect and Optimize
def process_data(i_max, roi, n_frames, data):
    """Process data for plotting. Find the mean signal, reference and the contrast
    Note: This is done after the whole scan is performed."""
    hsize = roi[2];
    vsize = roi[3];
    mean_sig = np.zeros(i_max);
    mean_ref = np.zeros(i_max);
    contrast = np.zeros(i_max);
    for i_scanpt in range(0, i_max):
        signal_frames = data[i_scanpt, 0::2, :, :];
        reference_frames = data[i_scanpt, 1::2, :, :];
        mean_sig[i_scanpt] = np.sum(np.mean(signal_frames, 0))
        mean_ref[i_scanpt] = np.sum(np.mean(reference_frames, 0))

        contrast[i_scanpt] = mean_sig[i_scanpt] / mean_ref[i_scanpt]
    return [mean_sig, mean_ref, contrast]

def plot_raw_data(data):
    data_reshape = np.reshape(data, (np.prod(data.shape[0:-2]), *data.shape[-2:]))
    data_reshape1 = np.reshape(data, (data.shape[0], np.prod(data.shape[1:-2]), *data.shape[-2:]))
    fig, (ax1, ax2) = plt.subplots(2,1, num=time.strftime("Raw: [%H:%M:%S]", time.localtime()))
    ax1.plot(np.sum(data_reshape, axis=(1,2)))
    ax1.set_ylabel('Sum of px')
    ax2.plot(np.sum(data_reshape1, axis=(2,3)).transpose());
    ax2.set_xlabel('Frame #'); ax2.set_ylabel('Sum of px')
    plt.tight_layout()


def plot_data(i_max, param, processed_data, x_label, x_unit, roi, n_frames, i_run, live=False):
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
        # plt.figure()
        plt.subplot(121)
        plt.plot([x / x_unit for x in xValues], mean_sig[0:i_max], '.-', label=f'S{i_run}')
        plt.plot([x / x_unit for x in xValues], mean_ref[0:i_max], '.--',  label=f'R{i_run}')
        # plt.legend()#['Avg Sig', 'Avg Ref'])
        plt.xlabel(x_label)

        plt.subplot(122)
        plt.plot([x / x_unit for x in xValues], contrast[0:i_max], '.-', label=f'{i_run}')
        plt.xlabel(x_label);
        plt.ylabel('Contrast')
        # plt.legend()
        plt.show()
        plt.tight_layout()
# --------------------------------------------------------------------------

class focus_adjustment_thread(threading.Thread):
    
    def __init__(self, t_exposure, t_align, camera_worker, trigger_event):
        super().__init__()
        self.t_exposure = t_exposure
        self.t_align = t_align
        self.camera_worker = camera_worker
        self.trigger_event = trigger_event
        # self.trigger_condition = trigger_condition
        # self.trigger_focus_condition = trigger_focus_condition

    def run(self):
        # self.trigger_event.wait()        # waiting for trigger to be set...
        # # PB settings
        # pb.focus_adjustment_sequence(t_exposure=self.t_exposure, t_align=self.t_align)
        # camera_worker.query_cam_settings()
        # print("starting PB sequence loop")
        # while self.camera_worker.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
        #     pb.start_sequence()
        #     time.sleep(pb.min_focus_time/s)
        #     pb.stop_sequence()
        # print("PB sequence loop stopped...")
        # with self.trigger_condition:
        #     self.trigger_event.clear()   # clear the event
        #     print("Trigger event cleared...")
        #     self.trigger_condition.notify_all()  # notify the waiting threads, here the PB thread..


        pb.focus_adjustment_sequence(t_exposure=self.t_exposure, t_align=self.t_align)
        # camera_worker.query_cam_settings()
        # print("starting PB sequence loop")
        self.trigger_event.set()
        self.trigger_event._flag = True
        # print(f"Flag status: {self.trigger_event._flag}")
        while (self.camera_worker.camera_status == dcamcon.DCAMCAP_STATUS.BUSY) \
        and self.trigger_event.is_set():
        # while self.trigger_event.is_set():
            print(f"{pb.min_focus_time/s} [s], {pb.n_repeat}")
            pb.start_sequence()
            time.sleep(pb.min_focus_time/s)
            # pb.stop_sequence()
        # print("PB sequence loop stopped...")
        self.trigger_event.clear()   # clear the event
        self.trigger_event._flag = False
        # print("Trigger event cleared...")

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
            print(f"\x1b[38;2;255;150;10m> {i_scanpt + 1} / {Nscanpts}: {param:0.0f}\x1b[0m")
            
            # TODO: thread this part as well for parallel operation..
            start_time = time.perf_counter()
            instructionList = []
            if self.sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                # if trial_run[0] == 'n':
                sg.set_sg_freq(param)
            else:
                self.seqArgList[0] = param
            
            the_list = pb.PB_program(instr, self.sequence, self.seqArgList)
            # print('the_list: ', the_list)
            for i in range(0, len(the_list)):
                instructionList.append(the_list[i][0])
            # print('instructionList: ', instructionList)

            # if trial_run[1] == 'n':
            if instr == 'cam':
                pb.run_sequence_for_camera(instructionList, self.t_exposure, self.t_seq_total, self.N_total)
            elif instr == 'cam_level1':
                pb.run_sequence_for_camera_level_trigger_1(instructionList)
            elif instr == 'cam_levelm':
                # run this for level trigger and without any alignment field (ac/dc) control
                pb.run_sequence_for_camera_level_trigger_many(instructionList, self.t_exposure,
                                                                  t_seq_total_i, self.N_total, expCfg.Nsamples)
                
            # elif 'levelm' in instr: # changed to 'levelm' in instr; prev instr == 'cam_levelm': to accomodate rot_field_controlled ODMR into the code
            # elif instr == 'cam_syncm_trigger_ao':     # changed to below for generalization
            elif 'syncm_trigger_ao' in instr:
                t_align = self.args[2]      # incoming time in [ns]
                pb.run_sequence_for_camera_sync_trigger_many_bac(instructionList, self.t_exposure, t_align,
                                                                 t_seq_total_i, self.N_total, expCfg.Nsamples)
                # print(f"PB Nsamples = {expCfg.Nsamples}")
            
            # elif instr == 'cam_levelm_trigger_ao':    # changed to below for generalization
            elif 'levelm_trigger_ao' in instr:      # valid for both ac and dc triggers
                t_align = self.args[2]      # incoming time in [ns]
                pb.run_sequence_for_camera_level_trigger_many_bac(instructionList, self.t_exposure, t_align,
                                                                  t_seq_total_i, self.N_total, expCfg.Nsamples)
                # print(f"PB Nsamples = {expCfg.Nsamples}")
            
            # elif instr == '??':
            #     t_align = self.args[2]      # incoming time in [ns]
            #     pb.run_sequence_for_camera_level_trigger_many_bdc(instructionList, self.t_exposure,
            #                                                       t_align, t_seq_total_i, self.N_total)

            elif instr == 'cam_timeseries':
                t_align_dc = self.args[2]
                pb.custom_trigger(self.t_exposure, t_align_dc)

            elif instr == 'cam_timeseries_trigger_ao':
                t_align_rot = self.args[2]
                t_measurement = self.args[3]
                t_align_rot_extended = self.args[4]
                pb.custom_trigger_rot_field(self.t_exposure, t_align_rot, t_measurement, t_align_rot_extended)
            
            # print("PB loaded: sets trigger event...")
            self.trigger_event.set()  # Set the event to trigger the camera
            self.inst_set_time_list.append((time.perf_counter()-start_time)*1e3)
            with self.condition:
                # check whether the event is 'set'... if yes, wait inside the loop - no further execution...
                while self.trigger_event.is_set():   # this statement dictates the wait condition: wait until the condition evaluates to False = event cleared..
                # for 1st execution, this is False
                    # print("PB waiting: camera to complete...")
                    self.condition.wait()        # waits for some event notification via notify() or notify_all()..
                    # print("PB wait ends: next PB scanpt...")
                    # whenever notified, the wait ends and proceeds again to check the condition
                # In this case: it is event clear notification: the condition evaluates to False (the event is cleared = not_set) and below statements are executed...
            self.time_taken.append((time.perf_counter()-start_time)*1e3)
        # return [i_scanpt, inst_set_time_list]

def start_dc_alignment_field(align_data, t_align, n_repetitions):
    """Perform the dc alignment before each run to initialize the propeller.
    Perform a [125 ms (depending on cutoff frequency) +/- Bz, 125 ms Bx] sequence 5-10 times before starting the rotating field sequence. For this, turn the x, z channel ON and let the pulseblaster control the timing.
    """
    # print("> Alignment field started!")
    print(f"Set t_align = {t_align}")
    
    if t_align < 50:
        instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW ^ concfg.start_trig, Inst.CONTINUE, 0, t_align*ms/2],
                        [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align*ms/2]]
    else:
        # t_align is in ms
        number_of_repetitions = int(np.ceil(np.ceil((t_align/2/1e3)/(np.max(align_data.shape)/samp_rate))/2)*2)
        pattern_time = np.max(align_data.shape)/samp_rate     # in [s]
        # print(f"pattern time = {pattern_time}")

        if direction == 'z':
            instructionList = [[concfg.laser ^ concfg.bx ^ concfg.by ^ concfg.start_trig, Inst.LOOP, 10, pattern_time*1e9],
                                [concfg.laser ^ concfg.bx ^ concfg.by, Inst.CONTINUE, 0, (t_align*ms/2-pattern_time*1e9)],
                                [concfg.laser ^ concfg.bz, Inst.END_LOOP, 0, (t_align*ms/2)]]
        elif direction == 'x':
            instructionList = [[concfg.laser ^ concfg.bz ^ concfg.start_trig, Inst.LOOP, 10, pattern_time*1e9],
                                [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align*ms/2-pattern_time*1e9)],
                                [concfg.laser ^ concfg.bx ^ concfg.by, Inst.END_LOOP, 0, (t_align*ms/2)]]
    pb.run_sequence_for_diode([instructionList])
    # time.sleep(10)
    print("writing dc alignment data")
    ao_task.write(align_data)
    

if __name__ == '__main__':
    set_core_affinity()
    warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
    expCfg = import_module(expCfgFile)

    # expCfg.N_scanPts = len(expCfg.scannedParam)
    t_AOM = expCfg.t_AOM
    clk_cyc = 1e3 / concfg.PBclk  # One clk cycle of PB = inverse of the clk freq of PB
    print('\x10 \x1b[38;2;250;250;0mRunning ' + expCfg.saveFileName + ' sequence\x1b[0m')
    [sg, ao_task, hdcamcon] = initialize_instr(expCfg.sequence)
    # ao_task = AnalogOutputTask()

    seqctrl.check_params(expCfgFile)
    print("\x10 Parameter checks completed...")

    [from_liveframes, roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr,expCfg)
    # if trial_run[1]=='n':
    [t_exposure, align_field, _] = from_liveframes
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
    # if trial_run[1]=='n':
    N_total = [] if expCfg.sequence == 'esr_seq' else [int(np.floor(t_exposure*1e9/t)) for t in t_seq_total[0]]
    # else:
    #     N_total = [] if expCfg.sequence == 'esr_seq' else [int(np.floor(t_exposure*1e6/t)) for t in t_seq_total[0]]
    N_total = []
    print(f"N_total = {N_total}")
    print(f"Set t_exposure = {t_exposure} [s]")
    
    # -------------------------- 20062023-------------------------
    # prepare data for AO - either DC field or rotating field
    if 'trigger_ao' in instr:   # this takes care of triggered timeseries + ODMR type experiments
        # prepare data for AO - either DC field or rotating field
        if 'ac' in instr:
            # prepare rotational alignment data pattern
            [dt, t_array, b_rot] = DAQ_write_pattern.triggered_ao_data_ac(direction, rot_angle, align_field, rot_field_amp, rot_field_freq)
            t_align_rot = np.max(b_rot.shape)*dt       # t_align_rot in [s]: total time of the rotational alignment pattern
            # calib_factor = ao_task.vi_calibration           # coil calibration factor
            # daqctrl.coil_calibration([1,1,1])

        elif 'dc' in instr:
            print('DC pattern data loading...')
            # prepare for DC alignment data pattern
            [dt, t_array, b_rot] = DAQ_write_pattern.triggered_ao_data_dc(t_align=t_align_dc/1e3, rotation_theta=prop_theta, rotation_phi=prop_phi, amp=rot_field_amp)
            t_align_rot = np.max(b_rot.shape)*dt       # t_align_rot in [s]: total time of the DC alignment pattern

        # if instr == 'cam_timeseries_trigger_ao':    # this is changed to below.. more modular
        if 'timeseries' in instr:
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
            v_rot_padded = b_rot_padded/ao_task.vi_calibration
                
            v_rot_padded_prepared = AnalogOutputTask.prepare_data_for_write(v_rot_padded)
            plt.figure(); plt.plot(t_array, b_rot_padded)       # t_array is required just for plotting purposes to verify the actual output
            plt.figure(); plt.plot(t_array, v_rot_padded)
            print(f"Expected t_align_rot = {t_align_rot} seconds")
            print(f"Extended t_align_rot = {t_align_rot_extended} seconds")
            print(f"#frames in alignment time (sig+ref) = {n_frames_rot_alignment}")

            v_rot_pattern = v_rot_padded_prepared.copy()
        
        # elif instr == 'cam_syncm_trigger_ao' or instr == 'cam_levelm_trigger_ao':
        elif 'sync' in instr or 'level' in instr:      # generalized 'else'
            v_rot = b_rot/ao_task.vi_calibration.T
            
            test_voltage = (np.array(test_field)/ao_task.vi_calibration).T
            v_rot = np.hstack((v_rot, test_voltage))
            
            v_rot_prepared = AnalogOutputTask.prepare_data_for_write(v_rot)
            v_rot_pattern = v_rot_prepared.copy()
            # below is not needed.. anyways this is not a part of alignment but test field
            # t_align_rot = np.max(v_rot_pattern.shape)*dt       # t_align_rot in seconds: total time of the rotational alignment pattern
    
    elif instr in ['cam_levelm', 'cam_syncm']:
        ao_task.set_outputs_to_constant(test_field)
        # else:       # handle cases when it is not triggered AO, ie only DC field measurement
    
    # if trial_run[1] == 'n' and hdcamcon is not None:
    if instr == 'cam_timeseries_trigger_ao':
        expCfg.Nsamples = int(n_frames_rot_alignment/2 + np.ceil(t_meas/1e3 *fps/2))      # there is a division by 2 since there will be a multiplication to automatically consider the signal and references in read_save_details(); 200 = FPS
    [frames_per_cyc, Nsamples] = read_save_details(roi, Nscanpts, expCfg.Nsamples)

    # -------------------------- 20062023-------------------------

    # expCfg.Nsamples = 1 mane holo at each scan pt 1 frame (like N APD readouts for averaging). So multiply by frames_per_cyc[0]
    # thle Nsamples = 2 hbe... at each scan pt 2 frames.. then go to next scan pt..
    # Nsamples = ekta scanpt e total kotogulo signal (ba reference) frames...
    # Nsamples = frames_per_cyc[0]*expCfg.Nsamples        # feed this as the 'n_frames' in capture() of Camcontrol.py... eta read_save_details() theke ashbe ebar...
    save_flag = False
    # if trial_run[1] == 'n':  # same as if hdcamcon is not None...
        # check this part and simplify..
        # just coverting everything to new library for now..----------------------

    camera_worker.configure_camera(instr)  # configure the camera for triggered acquisition

    # status = pb.pb_stop(); pb.errorCatcher(status) already stopped...
    # pb.custom_trigger(t_exposure *1e9, t_align_dc *1e6)
    
    closed = False
    display_parameters = dialog.yesno_box(expCfg.saveFileName + ' Params', 'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' % (str(concfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples, param[0] / expCfg.plotXaxisUnits, param[-1] / expCfg.plotXaxisUnits))
    paramfilename = savePath+"params_"+str(f_number)+".txt"

    if display_parameters == 'yes':
        if trial_run[1] == 'n':
            try:
                run_start_time = time.perf_counter()    # entire measurement time
                camera_worker.configure_camera(instr)
                camera_worker.query_cam_settings()
                print(f"Actual exposure time = {hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME)*1e3:0.2f} ms")
                PB_wait_time = []
                # set buffer and start capture sequence..
                while not camera_worker.start_capture(buffer_size=Nsamples, sequence=True, prep_only=True):
                    print("Retrying start_capture()..")
                
                if 'trigger_ao' in instr:
                    print(f"t_align_rot = {t_align_rot} ms")
                    # t_align_dc is in ms
                    if t_align_rot < 50:
                        instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW ^ concfg.start_trig, Inst.CONTINUE, 0, t_align_rot*ms/2],
                                        [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, t_align_rot*ms/2]]
                    else:
                        # duty = 0.07
                        duty = (t_exposure+8e-3)/(t_align_rot+2*(t_exposure+8e-3))
                        x = np.ceil(duty*t_align_rot/(1-duty)/10)*10       # in ms
                        instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, x],
                                            [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (t_align_rot*ms/2-x)],
                                            [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, (t_align_rot*ms/2)]]
                        
                    pb.run_sequence_for_diode([instructionList])

                print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)    # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
                print('\x10 %d frames in each cycles...' % frames_per_cyc[0])       # hence, total # frames = frames_per_cyc[0]*expCfg.Nsamples = Nsamples
                hsize = roi[2]; vsize = roi[3]
                print("------Acquiring %dx%d------" % (hsize, vsize))

                data = np.zeros((expCfg.Nruns, Nscanpts, Nsamples, vsize, hsize), dtype='uint16')
                # data_raw_time = np.zeros((expCfg.Nruns, Nscanpts, Nsamples), dtype='float64')

                # trigger_focus_event = threading.Event()
                # # trigger_focus_condition = threading.Condition()
                # focus_thread = focus_adjustment_thread(t_exposure*1e9, t_align_rot*1e9, camera_worker, trigger_focus_event)
                # focus_thread.start()
                # # focus_thread.join()

                # main loop
                for i_run in range(0, expCfg.Nruns):     # TODO: add Nruns loop as well
                    print("\x1b[38;2;255;150;10mRun: ", i_run + 1, ' / ', expCfg.Nruns, "\x1b[0m")
                    if (i_run+1)%1 == 0:
                        print("Adjusting focus...")
                        # camera settings
                        # camera_worker.query_cam_settings()
                        camera_worker.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.TRIGGERSOURCE, dcamcon.DCAMPROP.TRIGGERSOURCE.INTERNAL)
                        # camera_worker.query_cam_settings()
                        if instr == 'cam_levelm' or instr == 'cam_syncm':
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9)
                            
                        elif 'levelm_trigger' in instr:
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9, t_align=t_align_rot*1e9)
                            
                        else:
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9, t_align=t_align_rot*1e9,)
                            
                        print(f"{pb.min_focus_time/s} [s], {pb.n_repeat}")
                        pb.start_sequence()

                        time_start = time.perf_counter()
                        camera_worker.cv_window_status = 0
                        while (time.perf_counter()-time_start) < pb.min_focus_time/s:
                            timeout_happened = 0
                            
                            # print("Eta run hochhe??")
                            res = camera_worker.hdcamcon.wait_capevent_frameready(timeout_millisec=int(t_exposure*1e3+2))
                            if res is not True:
                                # frame does not come
                                if res != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                                    print('-NG: Dcam.wait_event() failed with error {}\x1b[0m'.format(res))
                                    break

                                # TIMEOUT error happens
                                timeout_happened += 1
                                if timeout_happened == 1:
                                    print('Waiting for a frame to arrive.', end='')
                                    if camera_worker.hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE) \
                                        == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
                                        print(' Check your trigger source.', end ='')
                                    else:
                                        print(' Check your <timeout_millisec> calculation in the code.', end='')
                                    print(' Press Ctrl+C to abort.\x1b[0m')
                                else:
                                    print('.')
                                    if timeout_happened > 5:
                                        timeout_happened = 0
                                # continue

                            # wait_capevent_frameready() succeeded
                            camera_worker.last_frame = camera_worker.hdcamcon.get_lastframedata()
                            # print('frame elo...')
                            if camera_worker.last_frame is not False:
                                if not camera_worker.display_frame(str(i_run+1), camera_worker.last_frame):
                                    # if q | Q is pressed on the cv2 window
                                    # self.stop()
                                    cv2.destroyWindow(str(i_run+1))
                                    # self.query_camera_status()       # expecting BUSY
                                    camera_worker.camera_status = dcamcon.DCAMCAP_STATUS.READY       # this is done to indicate the class that the capture has stopped
                                    
                                    print("Live View stopped...")
                                    break
                        cv2.destroyWindow(str(i_run+1))

                        pb.stop_sequence()
                        
                        camera_worker.configure_camera(instr)
                        # camera_worker.query_cam_settings()
                        print("Focus adjusted...")
                        # for 'cam_levelm', focus adjustment leaves a frame in the buffer, which needs to be cleared...
                        # if instr == 'cam_levelm':
                        # this might be for all 'level' triggered applications, hence modified below...
                        if 'level' in instr:
                            res = camera_worker.hdcamcon.wait_capevent_frameready(timeout_millisec=int(t_exposure*1e3+2))
                            camera_worker.last_frame = camera_worker.hdcamcon.get_lastframedata()
                        # time.sleep(2)
                    
                    print('Starting threads...')
                    trigger_event = threading.Event()       # Event object to signal between threads
                    condition = threading.Condition()
                    
                    kw_args = {'t_seq_total': t_seq_total, 'parameters': param, 'sequence': expCfg.sequence,\
                               'seqArgList': seqArgList, 'trigger_event': trigger_event, 'condition': condition}
                    # Create and start the scan thread: t_seq_total, parameters, sequence, seqArgList, trigger_event, condition
                    acq_start_time = time.perf_counter()
                    # print(f"t_align_rot = {t_align_rot}")
                    # print(f"t_align_rot_extended = {t_align_rot_extended}")
                    # print(f"t_meas = {t_meas}")

                    if instr == 'cam_timeseries_trigger_ao':
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, t_meas *1e6,\
                                               t_align_rot_extended *1e9, **kw_args)
                        
                    # elif instr == 'cam_syncm_trigger_ao' or instr == 'cam_levelm_trigger_ao':   # changed to below for generalization
                    elif 'syncm_trigger' in instr or 'levelm_trigger' in instr:
                        print('PB thread...')
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, **kw_args)

                    else:
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_dc *1e6, **kw_args)

                    if 'trigger_ao' in instr:
                        ao_task.create_retriggerable_ao_task(v_rot_pattern.shape)
                        # if instr == 'cam_timeseries_trigger_ao':  # this may not be required now since the variable 'v_rot_pattern' is generalized from 'v_rot_padded_prepared'
                        print("Starting alignment pattern!")
                        ao_task.start_retriggerable_ao_task(v_rot_pattern)
                        # time.sleep(0.04)
                        
                    scan_thread.start()
                    
                    # Main thread handling camera acquisition
                    # acq_start_time = time.perf_counter()
                    [data, scan_time_list, i_scanpt_cam] = acquire_data(trigger_event, condition, roi,\
                                                                        Nsamples, Nscanpts, i_run)
                    
                    # acq_time = time.perf_counter() - acq_start_time

                    # Wait for the pulse generation thread to complete
                    scan_thread.join()
                    # camera_thread.join()
                    # [i_scanpt_inst, inst_set_time_list] = scan_thread.result
                    acq_time = time.perf_counter() - acq_start_time

                    i_scanpt = Nscanpts
                    
                # all runs complete.. process for display...
                # focus_thread.join()
                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print(f"Run time = {run_time:0.2f} s")
                plot_raw_data(data=data)
                plt.figure(num=time.strftime(" [%H:%M:%S]", time.localtime()))
                for i_run in range(0, expCfg.Nruns):
                    processed_data = process_data(i_scanpt, roi, Nsamples, data[i_run,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                    plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel,
                              expCfg.plotXaxisUnits, roi, Nsamples, i_run, live=False)

                # camera_worker.query_cam_settings()

                # Close all after full acquisition
                closed = close_all(sg, hdcamcon, ao_task)

                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Data file saving
                    # if expCfg.Nruns == 1:
                    # if save_flag == False:  # Ask for file number iff no save was performed
                    #     [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
                    datafilename = savePath + expCfg.saveFileName + "_" + str(f_number)
                    # save_flag = save_data(datafilename, data[0,:,:,:,:], f_number)
                    for i_run in range(0, expCfg.Nruns):
                        datafilename = savePath + expCfg.saveFileName + "_camera_" + str(i_run) + ".tiff"
                        save_flag = save_data(datafilename, data[i_run,:,:,:,:], i_run)

                else:
                    print("\x10 \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")

            except KeyboardInterrupt:
                # focus_thread.join()
                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
                # Plot the last run data upto the point where it was interrupted...
                processed_data = process_data(i_scanpt, roi, Nsamples, data[0,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel,\
                          expCfg.plotXaxisUnits, roi, Nsamples, i_run, live=False)
                # Then ask whether to save it...
                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Ask for filename only if there was no save operation
                    if save_flag == False:
                        [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
                    save_flag = save_data(datafilename, data[0,:,:,:,:], f_number)
            finally:
                # focus_thread.join()
                if not closed:
                    closed = close_all(sg, hdcamcon, ao_task)

                print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))

                print("data = \x1b[38;2;250;150;50m%d x %d x %d x %d x%d\x1b[0m" % data.shape)
                # Save parameters (only if there was one save operation)...
                if save_flag:
                    expParamList[1] = i_scanpt + 1  # expParamList[1] -> value of N_scanPts
                    expParamList[3] = i_run + 1  # expParamList[3] -> value of Nruns

                    # if expCfg.Nruns>1:
                    #     paramfilename = savePath+"params_"+str(f_number)+".txt"

                    save_parameters(paramfilename, param_save_format, expParamList)

                    [next_params_format, expParamList] = extra_param_save_details(scan_time_list,\
                                                                                  run_time, roi, t_exposure)
                    save_parameters(paramfilename, next_params_format, expParamList)
                # plt.close('all')

                # sys.exit(exit_code)       # the exit code is the one returned from the pyqt app
                # scan_time_list in seconds; exec_time in seconds
        elif trial_run[1] == 'y':
            # TODO: #6 the laser is not ON during trial run --  why is that??
            try:
                run_start_time = time.perf_counter()    # entire measurement time
                print(f"Actual exposure time = {camera_worker.exposure*1e3:0.2f} ms")
                PB_wait_time = []
                
                if 'trigger_ao' in instr:
                    print(f"Set t_align_DC = {t_align_dc} ms")
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
                        
                    pb.run_sequence_for_diode([instructionList])

                print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)    # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
                print('\x10 %d frames in each cycles...' % frames_per_cyc[0])       # hence, total # frames = frames_per_cyc[0]*expCfg.Nsamples = Nsamples
                hsize = roi[2]; vsize = roi[3]
                print("------Acquiring %dx%d------" % (hsize, vsize))

                data = np.zeros((expCfg.Nruns, Nscanpts, Nsamples, vsize, hsize), dtype='uint16')

                # main loop
                for i_run in range(0, expCfg.Nruns):     # TODO: add Nruns loop as well
                    print("\x1b[38;2;255;150;10mRun: ", i_run + 1, ' / ', expCfg.Nruns, "\x1b[0m")
                    if (i_run+1)%5 == 0:
                        print("Adjusting focus...")
                        camera_worker.query_cam_settings()
                        if instr == 'cam_levelm' or instr == 'cam_syncm':
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9)
                            
                        elif 'levelm_trigger' in instr:
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9, t_align=t_align_rot*1e9)
                            
                        else:
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9, t_align=t_align_rot*1e9,)
                        
                        print(f"{pb.min_focus_time/s} [s], {pb.n_repeat}")
                        pb.start_sequence()

                        # camera_worker.start_capture(run_only=True)
                        time_start = time.perf_counter()
                        # while camera_worker.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
                        camera_worker.cv_window_status = 0
                        while (time.perf_counter()-time_start) < pb.min_focus_time/s:
                            timeout_happened = 0
                            
                            res = True #camera_worker.hdcamcon.wait_capevent_frameready(2000)
                            
                            # wait_capevent_frameready() succeeded
                            camera_worker.last_frame = (np.random.rand(*camera_worker.roi[-2:]) \
                                                        *(2**16)).astype(np.uint16).transpose()
                            # print('frame elo...')
                            if camera_worker.last_frame is not False:
                                if not camera_worker.display_frame(camera_worker.device_title, camera_worker.last_frame):
                                    # if q | Q is pressed on the cv2 window
                                    # self.stop()
                                    cv2.destroyWindow(camera_worker.device_title)
                                    # self.query_camera_status()       # expecting BUSY
                                    # next is done to indicate the class that the capture has stopped
                                    camera_worker.camera_status = dcamcon.DCAMCAP_STATUS.READY
                                    
                                    print(f"Live View stopped...")
                                    break
                        cv2.destroyWindow(camera_worker.device_title)
                        
                        pb.stop_sequence()

                        # camera settings
                        # camera_worker.configure_camera(instr)
                        # camera_worker.query_cam_settings()
                        print("Focus adjusted...")
                    
                    print('Starting threads...')
                    trigger_event = threading.Event()       # Event object to signal between threads
                    condition = threading.Condition()
                    
                    kw_args = {'t_seq_total': t_seq_total, 'parameters': param, 'sequence': expCfg.sequence,\
                               'seqArgList': seqArgList, 'trigger_event': trigger_event, 'condition': condition}
                    # Create and start the scan thread: t_seq_total, parameters, sequence, seqArgList, trigger_event, condition
                    acq_start_time = time.perf_counter()

                    if instr == 'cam_timeseries_trigger_ao':
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, t_meas *1e6,\
                                               t_align_rot_extended *1e9, **kw_args)
                    
                    # elif instr == 'cam_syncm_trigger_ao' or instr == 'cam_levelm_trigger_ao':   # changed to below for generalization
                    elif 'syncm_trigger' in instr or 'levelm_trigger' in instr:
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, **kw_args)
                    
                    else:
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_dc *1e6, **kw_args)
                    
                    if 'trigger_ao' in instr:
                        ao_task.create_retriggerable_ao_task(v_rot_pattern.shape)
                        print("Starting alignment pattern!")
                        ao_task.start_retriggerable_ao_task(v_rot_pattern)
                        # time.sleep(0.04)

                    scan_thread.start()
                    
                    # Main thread handling camera acquisition
                    # acq_start_time = time.perf_counter()
                    [data, scan_time_list, i_scanpt_cam] = acquire_data_sim(trigger_event, condition, roi,\
                                                                            Nsamples, Nscanpts, i_run)
                    
                    # acq_time = time.perf_counter() - acq_start_time

                    # Wait for the pulse generation thread to complete
                    scan_thread.join()

                    acq_time = time.perf_counter() - acq_start_time

                    i_scanpt = Nscanpts
                    
                # all runs complete.. process for display...

                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print(f"Run time = {run_time:0.2f} s")
                plot_raw_data(data=data)
                plt.figure(num=time.strftime(" [%H:%M:%S]", time.localtime()))
                for i_run in range(0, expCfg.Nruns):
                    processed_data = process_data(i_scanpt, roi, Nsamples, data[i_run,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                    plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel,
                              expCfg.plotXaxisUnits, roi, Nsamples, live=False)
                # camera_worker.query_cam_settings()

                # Close all after full acquisition
                closed = close_all(sg, hdcamcon, ao_task)

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

            except KeyboardInterrupt:
                # focus_thread.join()
                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
                # Plot the last run data upto the point where it was interrupted...
                processed_data = process_data(i_scanpt, roi, Nsamples, data[0,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel,\
                          expCfg.plotXaxisUnits, roi, Nsamples, live=False)
                # Then ask whether to save it...
                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Ask for filename only if there was no save operation
                    if save_flag == False:
                        [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
                    save_flag = save_data(datafilename, data[0,:,:,:,:], f_number)
            finally:
                # focus_thread.join()
                if not closed:
                    closed = close_all(sg, hdcamcon, ao_task)

                print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))

                print("data = \x1b[38;2;250;150;50m%d x %d x %d x %d x%d\x1b[0m" % data.shape)
                # Save parameters (only if there was one save operation)...
                if save_flag:
                    expParamList[1] = i_scanpt + 1  # expParamList[1] -> value of N_scanPts
                    expParamList[3] = i_run + 1  # expParamList[3] -> value of Nruns

                    if expCfg.Nruns>1:
                        paramfilename = savePath+"params_"+str(f_number)+".txt"

                    save_parameters(paramfilename, param_save_format, expParamList)

                    [next_params_format, expParamList] = extra_param_save_details(scan_time_list,\
                                                                                  run_time, roi, t_exposure)
                    save_parameters(paramfilename, next_params_format, expParamList)
                # plt.close('all')

                # sys.exit(exit_code)       # the exit code is the one returned from the pyqt app
                # scan_time_list in seconds; exec_time in seconds
        elif trial_run[1] == 'y':
            # TODO: #6 the laser is not ON during trial run --  why is that??
            try:
                run_start_time = time.perf_counter()    # entire measurement time
                print(f"Actual exposure time = {camera_worker.exposure*1e3:0.2f} ms")
                PB_wait_time = []
                
                if 'trigger_ao' in instr:
                    print(f"Set t_align_DC = {t_align_dc} ms")
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
                        
                    pb.run_sequence_for_diode([instructionList])

                print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)    # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
                print('\x10 %d frames in each cycles...' % frames_per_cyc[0])       # hence, total # frames = frames_per_cyc[0]*expCfg.Nsamples = Nsamples
                hsize = roi[2]; vsize = roi[3]
                print("------Acquiring %dx%d------" % (hsize, vsize))

                data = np.zeros((expCfg.Nruns, Nscanpts, Nsamples, vsize, hsize), dtype='uint16')

                # main loop
                for i_run in range(0, expCfg.Nruns):     # TODO: add Nruns loop as well
                    print("\x1b[38;2;255;150;10mRun: ", i_run + 1, ' / ', expCfg.Nruns, "\x1b[0m")
                    if (i_run+1)%5 == 0:
                        print("Adjusting focus...")
                        camera_worker.query_cam_settings()
                        if instr == 'cam_levelm' or instr == 'cam_syncm':
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9)
                            
                        elif 'levelm_trigger' in instr:
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9, t_align=t_align_rot*1e9)
                            
                        else:
                            pb.focus_adjustment_sequence(instr, t_exposure=t_exposure*1e9, Nsamples=Nsamples, \
                                                         focus_time=min_focus_time*1e9, t_align=t_align_rot*1e9,)
                        
                        print(f"{pb.min_focus_time/s} [s], {pb.n_repeat}")
                        pb.start_sequence()

                        # camera_worker.start_capture(run_only=True)
                        time_start = time.perf_counter()
                        # while camera_worker.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
                        camera_worker.cv_window_status = 0
                        while (time.perf_counter()-time_start) < pb.min_focus_time/s:
                            timeout_happened = 0
                            
                            res = True #camera_worker.hdcamcon.wait_capevent_frameready(2000)
                            
                            # wait_capevent_frameready() succeeded
                            camera_worker.last_frame = (np.random.rand(*camera_worker.roi[-2:]) \
                                                        *(2**16)).astype(np.uint16).transpose()
                            # print('frame elo...')
                            if camera_worker.last_frame is not False:
                                if not camera_worker.display_frame(camera_worker.device_title, camera_worker.last_frame):
                                    # if q | Q is pressed on the cv2 window
                                    # self.stop()
                                    cv2.destroyWindow(camera_worker.device_title)
                                    # self.query_camera_status()       # expecting BUSY
                                    # next is done to indicate the class that the capture has stopped
                                    camera_worker.camera_status = dcamcon.DCAMCAP_STATUS.READY
                                    
                                    print(f"Live View stopped...")
                                    break
                        cv2.destroyWindow(camera_worker.device_title)
                        
                        pb.stop_sequence()

                        # camera settings
                        # camera_worker.configure_camera(instr)
                        # camera_worker.query_cam_settings()
                        print("Focus adjusted...")
                    
                    print('Starting threads...')
                    trigger_event = threading.Event()       # Event object to signal between threads
                    condition = threading.Condition()
                    
                    kw_args = {'t_seq_total': t_seq_total, 'parameters': param, 'sequence': expCfg.sequence,\
                               'seqArgList': seqArgList, 'trigger_event': trigger_event, 'condition': condition}
                    # Create and start the scan thread: t_seq_total, parameters, sequence, seqArgList, trigger_event, condition
                    acq_start_time = time.perf_counter()

                    if instr == 'cam_timeseries_trigger_ao':
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, t_meas *1e6,\
                                               t_align_rot_extended *1e9, **kw_args)
                    
                    # elif instr == 'cam_syncm_trigger_ao' or instr == 'cam_levelm_trigger_ao':   # changed to below for generalization
                    elif 'syncm_trigger' in instr or 'levelm_trigger' in instr:
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_rot *1e9, **kw_args)
                    
                    else:
                        scan_thread = PBThread(t_exposure *1e9, N_total, t_align_dc *1e6, **kw_args)
                    
                    if 'trigger_ao' in instr:
                        ao_task.create_retriggerable_ao_task(v_rot_pattern.shape)
                        print("Starting alignment pattern!")
                        ao_task.start_retriggerable_ao_task(v_rot_pattern)
                        # time.sleep(0.04)

                    scan_thread.start()
                    
                    # Main thread handling camera acquisition
                    # acq_start_time = time.perf_counter()
                    [data, scan_time_list, i_scanpt_cam] = acquire_data_sim(trigger_event, condition, roi,\
                                                                            Nsamples, Nscanpts, i_run)
                    
                    # acq_time = time.perf_counter() - acq_start_time

                    # Wait for the pulse generation thread to complete
                    scan_thread.join()

                    acq_time = time.perf_counter() - acq_start_time

                    i_scanpt = Nscanpts
                    
                # all runs complete.. process for display...

                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print(f"Run time = {run_time:0.2f} s")
                plot_raw_data(data=data)
                plt.figure(num=time.strftime(" [%H:%M:%S]", time.localtime()))
                for i_run in range(0, expCfg.Nruns):
                    processed_data = process_data(i_scanpt, roi, Nsamples, data[i_run,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                    plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel,
                              expCfg.plotXaxisUnits, roi, Nsamples, live=False)
                # camera_worker.query_cam_settings()

                # Close all after full acquisition
                closed = close_all(sg, hdcamcon, ao_task)

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

            except KeyboardInterrupt:
                # focus_thread.join()
                run_end_time = time.perf_counter()
                run_time = run_end_time - run_start_time  # entire measurement time
                print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
                # Plot the last run data upto the point where it was interrupted...
                processed_data = process_data(i_scanpt, roi, Nsamples, data[0,:,:,:,:])  # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
                plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel,\
                          expCfg.plotXaxisUnits, roi, Nsamples, live=False)
                # Then ask whether to save it...
                savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
                if savefile_yn == 'yes':
                    # Ask for filename only if there was no save operation
                    if save_flag == False:
                        [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
                    save_flag = save_data(datafilename, data[0,:,:,:,:], f_number)
            finally:
                # focus_thread.join()
                if not closed:
                    closed = close_all(sg, hdcamcon, ao_task)

                print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))

                print("data = \x1b[38;2;250;150;50m%d x %d x %d x %d x%d\x1b[0m" % data.shape)
                # Save parameters (only if there was one save operation)...
                if save_flag:
                    expParamList[1] = i_scanpt + 1  # expParamList[1] -> value of N_scanPts
                    expParamList[3] = i_run + 1  # expParamList[3] -> value of Nruns

                    if expCfg.Nruns>1:
                        paramfilename = savePath+"params_"+str(f_number)+".txt"

                    save_parameters(paramfilename, param_save_format, expParamList)

                    [next_params_format, expParamList] = extra_param_save_details(scan_time_list,\
                                                                                  run_time, roi, t_exposure)
                    save_parameters(paramfilename, next_params_format, expParamList)
                # plt.close('all')
    else:
        print("Measurement aborted...")
        if not closed:
            closed = close_all(sg, hdcamcon, ao_task)
