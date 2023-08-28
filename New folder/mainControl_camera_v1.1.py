#%% Initialization and Definition
# reset

import Camcontrol as camctrl; import PBcontrol as PBctrl; import sequencecontrol as seqctrl; import SGcontrol as SGctrl; import connectionConfig as conCfg
from spinapi import ns, us, ms, Inst
import matplotlib.pyplot as plt; import numpy as np; from os.path import isdir, isfile; from os import makedirs; from importlib import import_module
import sys, time, dialog, cv2
# %matplotlib qt5

# from matlab import engine as eng

# def main():
print("\x10 \x1b[0mImports Successful...")
plt.rcParams.update({'figure.max_open_warning': 0})     # No warnings on opening mult fig windows
global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, SG, expCfg
expCfgFile = 'esr'+'_config'

trial_run = 'n'
seq_no_plot = [0]
voltage_unit = 1      # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100                      # The dpi of the displayed pulse sequence plot
plotPulseSequence = True
# t_exposure = 0.1     # in seconds.. eta dorkar o nei.. user will define it..
livePlotUpdate = False

def initialize_instr(sequence):
    try:
        PBctrl.configurePB()
        print('\x10 PB: \x1b[38;2;250;250;0mv'+PBctrl.pb_get_version()+'\x1b[0m')   # Display the PB board version using pb_get_version()
    except:
        print("Error Initializing PB !!")
    # finally:
    #     sys.exit()
    if trial_run in ['n','N'] and sequence not in ['aom_timing',  'rodelay']:       # Do not initialize SG if it is a trial run or the sequence is present in the list ['aom_timing', 'T1ms0', 'rodelay']
        SG = SGctrl.initSG(conCfg.serialaddr, conCfg.model_name)      # Initialize SG using RS-232 address and the model name
        if not(SG.query('*idn?') == ''):
            print("\x10 \x1b[38;2;250;250;0mSG384 Initialized...\x1b[0m")
            SG.write('remt')
            SG.write('disp2')
            SGctrl.enable_SG_op(SG)
            print("SG Output Enabled...")
            SGctrl.set_SG_amp(SG,expCfg.MW_power)
            SGctrl.setup_SG_pulse_mod(SG)
            print("SG Ext Pulse Mod Enabled...")

        cam = camctrl.open_camera()
        cam["trigger_source"] = 1       # Internal(1), external(2), software(3), master_pulse(4)        
        return [SG, cam]
    else:
        print("\x10 Trial run... \x1b[38;2;250;250;0mNo SG")
        return [None, None]
    
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
    
    param = expCfg.scannedParam;    n_error=0
    param_save_format = expCfg.formattingSaveString
    sequenceArgs = expCfg.updateSequenceArgs()      # Variables used in the pulse sequence
    expParamList = expCfg.updateExpParamList()      # List of experimental parameters
    
    # ------------------------------------------------------------------
        # est_time = expCfg.t_tot * expCfg.Nsamples * expCfg.N_scanPts ??
    # check the scanned parameters for errors, if errors are found, remove those params
    if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:   # for sequences except ESR, set MW frequency
        seqArgList = [param[seq_no_plot[-1]]]
        seqArgList.extend(sequenceArgs)     # Make a dummy seqArgList just to create it
        [n_error, param] = seqctrl.param_err_check(instr, expCfg.sequence, expCfg.PBchannels, seqArgList,  expCfg.scannedParam, expCfg.N_scanPts)
        if n_error>0:
            print('\x1b[1;37;41m'+'Err: Check Sequences...\x1b[0m')
            print('\x1b[38;2;250;0;0m'+str(n_error)+'\x1b[0m parameters removed...')

            # Close Error plots??
            close_plots = dialog.yesno_box('Close Plots', 'Close the Error Plots?')
            if close_plots == 'yes':
                plt.close('all')

            print("\x10 Sequences checked... \x1b[38;2;100;250;0mErrors removed...\x1b[0m")
            # est_time = (2*expCfg.t_AOM*Nscanpts + sum(param))*expCfg.Nsamples
        else:
            print('\x1b[38;2;100;250;0m----No Errors----\x1b[0m')
    else:       # for ESR exp, set MW freq to start pt of the scan -> there is no scannedParam
        seqArgList = sequenceArgs
        
    # define 'Nscanpts' as the length of the 'param' variable...
    Nscanpts = len(param)
    if Nscanpts>0:
        print("\x10 \x1b[38;2;250;100;10m%d\x1b[0m scan pts" % Nscanpts)
    else:
        print("\x1b[38;2;200;200;10mCheck sequences for subtle problems...\x1b[0m")
    
    # Plotting the pulse sequences -------------------------------------------
    if plotPulseSequence:
        view_sequence(expCfg.sequence, expCfg.PBchannels, seqArgList, False, param, seq_no_plot, seq_plot_dpi)
        
    # Start the initial sequence now ------------
    print("Starting Initial Sequence...")
    if trial_run in ['n','N']:
        instructionList = [start_initial_PB_seq()]
        PBctrl.run_sequence_for_diode(instructionList)
        print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
        if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq','aom_timing','rodelay'] and trial_run in ['n','N']:
            SGctrl.set_SG_freq(SG, expCfg.MW_freq)
    else:
        print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    # Create the data save folder...
    savePath = expCfg.savePath + time.strftime("%Y-%m-%d", time.localtime()) + '\\'
    if not (isdir(savePath)):
        makedirs(savePath)
    print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')

    # adjust the exposure time and select ROI for experiment
    camctrl.live_view(cam)      # see image live and adjust the exposure time
    # now select the ROI
    roi = camctrl.select_roi(cam)
    # roi = [1048, 1136, 464, 460]
    # cam["exposure_time"] = 0.05

    return [roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList]
    
def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[conCfg.laser, Inst.CONTINUE, 0, 500*ms],
                       [conCfg.laser, Inst.BRANCH, 0, 500*ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...
    
    return instructionList

def close_all(*args):       # arguments are 'cam' and 'SG'
    cam = args[0]
    if trial_run in ['n','N'] and len(args)==2:
        SG = args[1]
        if 'SG' in vars():
            SG.write('freq2.87ghz')
            SG.write('lcal')
            SGctrl.disable_SG_op(SG)
            # SG.close()
        camctrl.close_camera(cam); print("\x1b[38;2;10;250;50m\x10 Camera Closed...")
    # PBctrl.pb_init();
    PBctrl.pb_stop(); PBctrl.pb_close(); print("\x10 Pulse Blaster closed...\x1b[0m")
    # return DAQclosed

def read_save_details(roi, N_scanpts, Nsamples_expCfg):
    """defines the 'frame_per_cyc' to be captured, and the write formats of the 'scannedParam' and the 'data'
    Returns: 
    frames_per_cyc            - define the number of frames captured in one cycle of the sequence, 1 signal, 1 ref, so 2. 
    scannedparam_write_format - write format of the scannedparams in the data file
    data_write_format         - write format of the data, eta line by line er format ta.. not the whole at a time..
    """
    hsize = roi[2]; vsize = roi[3];
    frames_per_cyc = [2]
    Nsamples = frames_per_cyc[0]*Nsamples_expCfg

    scannedparam_write_format = "%g\t" * N_scanpts
    scannedparam_write_format = scannedparam_write_format[0:-1] + "\n"

    data_write_format = ""
    # niche n_frames holo number of frames acquired in each cycle.. (dekhe tai mone hochhe, bhule gechhi)
    # replacing the command range(0, n_frames*hsize) with range(0, frames_per_cyc*hsize)
    data_write_format = "%d\t" * Nsamples*hsize    # n_frames*hsize... (n_frames*vsize na)
    # uporer line edit korlam... frames_per_cyc[0]*hsize...20230701...
    # # puro file ta ekbar e likhe debo... niche... memory error dekhachhe at full resolution for 1000 pts...
    # data_write_format = (data_write_format[0:-1] + "\n") * (N_scanpts*vsize) # remove \t at the end and add \n
    data_write_format = (data_write_format[0:-1] + "\n")    # puro file ekbar e save kora jachhe na... tai single line format
    
    return [frames_per_cyc, Nsamples, scannedparam_write_format, data_write_format]


def extra_param_save_details(scan_time_list, exec_time):
    extra_params_format = ' %s\t%0.2f\n %s\t%0.2f\n %s\t%0.4g\n %s\t%0.4g\n %s\t%g\n %s\t%d\n %s\t%d\n %s\t%d\n %s\t%d\n'  # %s\t%f\n
    # scan_time_list in ms; exec_time in seconds
    param_list = ['Max_scan_time(us):',max(scan_time_list)*1e6, 'Min_scan_time(us):',min(scan_time_list)*1e6, 'Total_scan_time(s):',np.sum(scan_time_list), 'Total_run_time(s):',exec_time, 'Step:', (expCfg.scannedParam[1]-expCfg.scannedParam[0]), 'X0:',roi[0], 'Y0:',roi[1], 'W:',roi[2], 'H:',roi[3]]
    # param_list[1] = i_scanpt+1..... Ekhane ki hbe?? Kon parameter save korbo??
    # ekhane ekta if kore, jodi 'i'=1 hoe, thle first param_list ta return korbe, nhle porer param_list ta return korbe.
    # Eta jodi kora hoe, thle runs>1 hole ba multi-param scan hole, notun param er sathe tar details save kora jabe...
    # Config file e ekta variable lagbe jeta dekhabe je kon variable ta scan hoechhe, other than scannedParam
    param_list = tuple(param_list)
    return [extra_params_format, param_list]
    
def prepare_for_saving(savePath):
    global file_number
    if 'file_number' in vars():       # Returns a dict of all local variables
        print("Previous file number: \x1b[38;2;250;150;0m"+file_number+'\x1b[0m')
    file_number = input("File name:"+expCfg.saveFileName+"_camera_#: ")
    
    paramfilename = savePath + expCfg.saveFileName + "_camera_" + "params_" + file_number +".txt"
    datafilename = savePath + expCfg.saveFileName + "_camera_" + file_number +".txt"
    if (isfile(datafilename)):
        print('\x1b[38;2;250;250;0mAppending file: '+expCfg.saveFileName+'_'+file_number+'.txt\x1b[0m')
        
    return [paramfilename, datafilename, file_number]

def save_data(datafilename, data_raw, data_write_format, scannedparam_write_format, param):
    print("\x10 Saving to file.....")
    with open(datafilename, 'a') as datafile:
        # first e file tar first line e sob kota parameters save kora hok...
        datafile.write(scannedparam_write_format % tuple(param)) # param ekta list
        # ebar line-by-line data (image frames) gulo lekha hok...
        for line in data_raw:
            datafile.write(data_write_format % tuple(line))
        # file lekha completed...
    
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
            if  sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                # expCfg.t_AOM = t_AOM*parameter[seq_no_plot[i]]/parameter[0]
                # sequenceArgs = expCfg.updateSequenceArgs()
                # seqArgList = [parameter[seq_no_plot[i]]]
                # seqArgList.extend(sequenceArgs)
                # expCfg.t_AOM = t_AOM
                seqArgList[0] = parameter[seq_no_plot[i]]
        plt.figure(dpi = plot_dpi)
        the_list = PBctrl.PB_program(instr,sequence,seqArgList)
        for j in range(0, len(the_list)):
            instructionList = the_list[j][0]
            inst_times = the_list[j][4]
            
            plt.subplot(1,2,(j+1))
            [t_us,channelPulses,yTicks] = seqctrl.plot_sequence(instructionList, PBchannels)
            for channel in channelPulses:
                plt.plot(t_us, list(channel))
            plt.yticks(yTicks, PBchannels.keys())          # Include the names of the PB channels
            plt.xlabel('Time (us)')
            plt.ylabel('Channel')
            plt.title(sequence + ' Plot. Param @ '+str(seq_no_plot[i])+': ' + str(parameter[seq_no_plot[i]]) + 'ns\nTransitions [us]: '+str([vals for vals in inst_times.values()]), fontsize=10)
    # Way to display the total no of instructions in the plot, since it will vary with different scanPts:
        # print("Total %d inst" % len(inst_times.keys()))

def acquire_data(t_exposure, t_seq_total, roi, Nsamples, parameter, sequence, seqArgList, trial):
    # setup next scan iteration (e.g. for ESR experiment, change microwave frequency; for T2 experiment, reprogram pulseblaster with new delay)
    instructionList=[]
    if sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
        if trial in ['n','N']:
            SGctrl.set_SG_freq(SG, parameter)
    else:
        seqArgList[0] = parameter
    the_list = PBctrl.PB_program(instr,sequence,seqArgList)
    # print(the_list)
    for i in range(0, len(the_list)):
        instructionList.append(the_list[i][0])
    # print(instructionList)
    #--------------------eta aage chhilo... get_frames() diye acquisition--------------------
    # scan_start_time = time.perf_counter()   # time in seconds
    # PBctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total)
    # frames = camctrl.get_frames(cam,roi,Nsamples)    # capture frames from camera
    # PBctrl.pb_stop()
    #--------------------niche notun ta... stream() ke bhenge--------------------
    hsize = roi[2]; vsize = roi[3];
    frames=np.zeros((vsize,hsize,Nsamples),dtype=np.uint16);
    stream = camctrl.ham.Stream(cam,Nsamples)   # the statement assigns the buffer according to the no of frames to acquire..
    stream.__enter__()
    # PBctrl.pb_stop()        # eta dewa ta thik noi... already ekta pb_stop() achhe sesh e... eta dewa hoyechhe initial sequence ta ke stop korar jonne.. that is NOT correct...
    camctrl.cam.start()
    time.sleep(0.28)
    scan_start_time = time.perf_counter()   # time in seconds
    PBctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total)
    # what shows that a sample/frame has been acquired???
    for i, frame_buffer in enumerate(stream):
        if i==Nsamples-1:
            scan_end_time = time.perf_counter()
        # ekhane 'run_sequence' dile hbe na.. trigger er jonno wait korchhe.. tar mane frame gulo start() korar porei 'stream' hoye jaye..
        # print(frame_buffer)
        frame = camctrl.ham.copy_frame(frame_buffer)
        frames[:,:,i] = frame
    
    stream.__exit__()
    PBctrl.pb_stop()
    
    scan_time = (scan_end_time - scan_start_time)     # in seconds
    # if ith_scan_pt==0:
    #     time.sleep(100)
    return [frames, scan_time, instructionList]

def process_data(i_max, roi, n_frames, data_raw):
    hsize = roi[2]; vsize = roi[3];
    mean_sig = np.zeros(i_max);     mean_ref = np.zeros(i_max);     contrast = np.zeros(i_max);
    for i in range(0,i_max):
        # signal_frames = np.array([data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize] for j in range(0,n_frames,2)])
        # reference_frames = np.array([data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize] for j in range(1,n_frames,2)])
        sum_of_signal_frames = np.array([np.sum(data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize]) for j in range(0,n_frames,2)])
        sum_of_reference_frames = np.array([np.sum(data_raw[i*vsize:(i+1)*vsize,j*hsize:(j+1)*hsize]) for j in range(1,n_frames,2)])

        mean_sig[i] = np.mean(sum_of_signal_frames)
        # print(mean_sig[i])
        mean_ref[i] = np.mean(sum_of_reference_frames)
        # print(mean_ref[i])
        # define the contrast
        contrast[i] = mean_sig[i]/mean_ref[i]
    return [mean_sig, mean_ref, contrast]

def plot_data(i_max, param, processed_data, x_label, x_unit, roi, n_frames, live=False):
    """Plot the signal, reference and the contrast"""
    # n_frames = total frames in one cycle of sequence..
    [mean_sig, mean_ref, contrast] = processed_data
    xValues = param[0:i_max]
    if live:
        plt.plot([x/x_unit for x in xValues], contrast[0:i_max], 'b--')
        plt.pause(0.0001)
    else:
        # plt.figure("Signal, Reference & Contrast Plots: %d points" %(i_max))
        plt.figure()
        plt.subplot(121)
        plt.plot([x/x_unit for x in xValues], mean_sig[0:i_max], '.-',[x/x_unit for x in xValues], mean_ref[0:i_max], '.-')
        plt.legend(['Avg Sig','Avg Ref'])
        plt.xlabel(x_label)

        plt.subplot(122)
        plt.plot([x/x_unit for x in xValues],  contrast[0:i_max], '.-b')
        plt.xlabel(x_label); plt.ylabel('Contrast')
    
# ----------------------------------------------------------------------------

# if __name__ == '__main__':
    # instructionList = main()

expCfg = import_module(expCfgFile)
# expCfg.N_scanPts = len(expCfg.scannedParam)
t_AOM = expCfg.t_AOM
clk_cyc = 1e3/conCfg.PBclk      # One clk cycle of PB = inverse of the clk freq of PB
print('\x10 \x1b[38;2;250;250;0mRunning '+expCfg.saveFileName+' sequence\x1b[0m')
[SG, cam] = initialize_instr(expCfg.sequence)

seqctrl.check_params(expCfgFile)
print("\x10 Parameter checks completed...")

instr='cam'
[roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr, expCfg)
roi = [int(element/4)*4 for element in roi];
# roi[2] = 8; roi[3] = 12;
hsize = roi[2]; vsize = roi[3];
# hsize = 8; vsize =12;


if not (expCfg.sequence == 'esr_seq'):
    # t_manip = [param for ]
    t_manip = np.array(param)
else:
    t_manip = np.zeros(len(param))

# below 2 lines are for T1 with varying 'init' duration such that the duty cycle remains same...
# t_seq_total_sig = [t_AOM*param[i]/param[0] for i in range(0,len(param))] + t_manip
# # t_seq_total_ref = [t_AOM*param[i]/param[0] for i in range(0,len(param))] + t_manip
# t_seq_total_ref = np.zeros(len(t_manip),dtype='float64')        # for T1 seq with laser ALWAYS ON

t_seq_total_sig = expCfg.t_AOM + t_manip
t_seq_total_ref = expCfg.t_AOM + t_manip                      # for other sequences..
# t_seq_total_ref = np.zeros(len(t_manip),dtype='float64')        # for T1 seq with laser ALWAYS ON
t_seq_total = np.transpose(np.array([t_seq_total_sig, t_seq_total_ref]))

#-------------------------- 20062023-------------------------

[frames_per_cyc, Nsamples, scannedparam_write_format, data_write_format] = read_save_details(roi, Nscanpts, expCfg.Nsamples)

#-------------------------- 20062023-------------------------

# expCfg.Nsamples = 1 mane holo at each scan pt 1 frame (like N APD readouts for averaging). So multiply by frames_per_cyc[0]
# thle Nsamples = 2 hbe... at each scan pt 2 frames.. then go to next scan pt..
# Nsamples = ekta scanpt e total kotogulo signal (ba reference) frames...
# Nsamples = frames_per_cyc[0]*expCfg.Nsamples        # feed this as the 'n_frames' in capture() of Camcontrol.py... eta read_save_details() theke ashbe ebar...
save_flag=False
# configure the camera for triggered acquisition
camctrl.configure_camera(cam)
# set the ROI for the acquisition
camctrl.set_roi(cam,roi,True)
t_exposure = cam["exposure_time"].value     # mind the value... in seconds !!!
# t_exposure = 0.1
display_parameters = dialog.yesno_box(expCfg.saveFileName+' Params', 'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' %(str(conCfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples, param[0]/expCfg.plotXaxisUnits, param[-1]/expCfg.plotXaxisUnits))
PBctrl.pb_stop()
if display_parameters == 'yes':
    try:
        # continue_run = 'yes'
        # DAQclosed = False;
        # save_flag = False
        # DAQtask = DAQctrl.configure_daq(Nsamples)        # configure_daq() accepts no. of samples to read from DAQ-AI & returns the task created; variable 'readTask'
        print("Exposure time: %0.1f ms" %(cam["exposure_time"].value * 1e3))
        print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)  # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
        print('\x10 %d frames in each cycles...' % frames_per_cyc[0])
        # from above, the total number of frames to acquire is thus frames_per_cyc[0]*expCfg.Nsamples = Nsamples
        
        print("------Acquiring %dx%d------" % tuple([roi[2], roi[3]]))
        for i_run in range(0, expCfg.Nruns):
            data_raw = np.zeros((vsize*Nscanpts,hsize*Nsamples));   # check the dimension... hsize, vsize, Nscanpts, Nsamples... OK?
            scan_time_list = []   # time (in seconds) for each scannedParam
            # print("\x10",i_run+1,'/',expCfg.Nruns)
            start_time = time.perf_counter()        # in seconds
            for i_scanpt in range (0, Nscanpts):
                print("\x10 ",i_scanpt+1,' / ',Nscanpts,': ',param[i_scanpt])
                # print("Exp [ms]: ",(cam["exposure_time"].value)*1e3)
                # expCfg.t_AOM = t_AOM*param[i_scanpt]/param[0]
                # print("t_AOM [ns] ",expCfg.t_AOM)
                # sequenceArgs = expCfg.updateSequenceArgs()
                # seqArgList = [param[i_scanpt]]
                # seqArgList.extend(sequenceArgs)
                # print(seqArgList)
                [frames, scan_time, instructionList] = acquire_data(t_exposure*1e9, t_seq_total[i_scanpt], roi, Nsamples, param[i_scanpt], expCfg.sequence, seqArgList, trial_run)
                data_raw[i_scanpt*vsize:(i_scanpt+1)*vsize,:] = np.concatenate(tuple(frames[:,:,j] for j in range(0,Nsamples)), axis=1)       # 2d-array... ekhane for j in range(0,frames_per_cyc[0]) chhilo... otake range(0,Nsamples) kore dilam.. think.. concat korte gele j should run over all the frames.. jodi total 4 frames thake (2 sig, 2 ref), then j=0,1,2,3...
                scan_time_list.append(scan_time)    # in seconds

            end_time = time.perf_counter()          # in seconds
            exec_time = end_time - start_time       # in seconds
            processed_data = process_data(i_scanpt, roi, Nsamples, data_raw)    # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
            plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples, live=False)
            print("\x10 Execution time = %.2fs" % exec_time)        # in seconds
            print("\x10 Scan time = %.2fs" % np.sum(scan_time_list))   # in seconds; scan_time_list in seconds
            # Close all if this is the last run
            if i_run+1 == expCfg.Nruns:
                close_all(cam, SG) if trial_run in ['n','N'] and expCfg.sequence not in ['T1ms0','aom_timing','rodelay'] else close_all(cam)
                # a=1
            
            savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
            if savefile_yn == 'yes':
                # Data file saving
                if save_flag == False:          # Ask for file number iff no save was performed
                    [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
                save_flag = save_data(datafilename, data_raw, data_write_format, scannedparam_write_format, param)
            else:
                print("\x10 \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")
            if (i_run+1) < expCfg.Nruns:
                continue_run = dialog.yesno_box('Continue',"Continue Run?")
                if continue_run == 'yes':
                    continue
                else:
                    sys.exit("Run Interrupted...")
                
    except KeyboardInterrupt:
        end_time = time.perf_counter()
        exec_time = end_time - start_time       # in seconds
        print('\x1b[38;2;250;100;0mUser Interrupted. Quitting...\x1b[0m')
        # Plot the last run data upto the point where it was interrupted...
        processed_data = process_data(i_scanpt, roi, Nsamples, data_raw)    # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
        plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples, live=False)
        # Then ask whether to save it...
        savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
        if savefile_yn == 'yes':
            # Ask for filename only if there was no save operation
            if save_flag == False:
                [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
            save_flag = save_data(datafilename, data_raw, data_write_format, scannedparam_write_format, param)
    # nicher ei 'else' ta kno dewa hoyechhilo??? 18-06-2023
    # else:
    #     if (i_run+1)==expCfg.Nruns: # Save data of the final run here
    #         save_data(len(conCfg.input_terminals), datafilename, data_raw, data_write_format, expCfg.Nsamples, reads_per_cyc)
    finally:
        close_all(cam, SG) if trial_run in ['n','N'] and expCfg.sequence not in ['T1ms0','aom_timing','rodelay'] else close_all(cam)
        
        print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))
        
        print("\x10 len(frames) = \x1b[38;2;250;150;50m%d\x1b[0m" % len(frames))
        # Save parameters (only if there was one save operation)...
        if save_flag:
            expParamList[1] = i_scanpt+1        # expParamList[1] -> value of N_scanPts
            expParamList[3] = i_run+1           # expParamList[3] -> value of Nruns
            save_parameters(paramfilename, param_save_format, expParamList)
            
            [next_params_format, expParamList] = extra_param_save_details(scan_time_list, exec_time)
            save_parameters(paramfilename, next_params_format, expParamList)
        
        # scan_time_list in seconds; exec_time in seconds
        #----------------------------------------------------------------------
        
        # matlab_engine = eng.start_matlab()
        
        # if expCfg.sequence == 'esr':
        #     analysis = eng.ESR_raw_analysis()
        
else:
    print("Dialog closed...")
    close_all(cam, SG) if trial_run in ['n','N'] and expCfg.sequence not in ['T1ms0','aom_timing','rodelay'] else close_all(cam)
    # sys.exit()
# return instructionList

