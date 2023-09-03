#%% Initialization and Definition
# reset

import Camcontrol as camctrl; import PBcontrol_v2 as PBctrl; import sequencecontrol as seqctrl; import SGcontrol as SGctrl; import connectionConfig as conCfg
from spinapi import ns, us, ms, Inst
import matplotlib.pyplot as plt; import numpy as np; from os.path import isdir, isfile; from importlib import import_module
import sys, time, dialog, cv2, os
# %matplotlib qt

# from matlab import engine as eng

# def main():
print("\x10 \x1b[0mImports Successful...")
plt.rcParams.update({'figure.max_open_warning': 0})     # No warnings on opening mult fig windows
global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, SG, expCfg, scan_time
expCfgFile = 'rabi'+'_config'
trial_run = 'y'
seq_no_plot = [-1]
voltage_unit = 1      # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100                      # The dpi of the displayed pulse sequence plot
plotPulseSequence = True
exposure_time = 0.1     # in seconds
livePlotUpdate = False

def initialize_instr(sequence):
    # if trial_run not in ['n']:
    try:
        PBctrl.configurePB()
        print('\x10 PB: \x1b[38;2;250;250;0mv'+PBctrl.pb_get_version()+'\x1b[0m')   # Display the PB board version using pb_get_version()
    except:
        print("Error Initializing PB !!")
    # finally:
    #     sys.exit()
    if trial_run in ['n','N'] and sequence not in ['aom_timing', 'T1ms0', 'rodelay']:       # Do not initialize SG if it is a trial run or the sequence is present in the list ['aom_timing', 'T1ms0', 'rodelay']
        SG = SGctrl.initSG(conCfg.serialaddr, conCfg.model_name)      # Initialize SG using RS-232 address and the model name
        print('\x10 Signal Generator: \x1b[38;2;250;250;0m' + SG.query('*idn?') + '\x1b[0m')    # Query instrument
        print("SG Initialized...")
        SG.write('remt')
        SG.write('disp2')
        SGctrl.enable_SG_op(SG)
        print("SG Output Enabled...")
        SGctrl.set_SG_amp(SG,expCfg.MW_power)
        SGctrl.setup_SG_pulse_mod(SG)
        print("SG Ext Pulse Mod Enabled...")

        # cam = camctrl.open_camera()
        # cam["trigger_source"] = 1       # Internal(1), external(2), software(3), master_pulse(4)
        # camctrl.adjust_roi(cam)         
        # return [SG, cam]
        return SG
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
    # if trial_run in ['n','N']:
    instructionList = [start_initial_PB_seq()]
    PBctrl.run_sequence_for_diode(instructionList)
    print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
    print("Initial Sequence Started...")
    if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq','aom_timing','rodelay'] and trial_run in ['n','N']:
        SGctrl.set_SG_freq(SG, expCfg.MW_freq)
    # else:
    #     print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    # Create the data save folder...
    savePath = expCfg.savePath + time.strftime("%Y-%m-%d", time.localtime()) + '\\'
    if not (isdir(savePath)):
        os.makedirs(savePath)
    print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')

    # adjust the exposure time and select ROI for experiment
    # camctrl.live_view(cam)      # see image live and adjust the exposure time
    # now select the ROI
    # [frame, roi] = camctrl.select_roi(cam)
    roi = [100, 100, 500, 600]

#-------------------------- 20062023-------------------------
    return [roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList]
#-------------------------- 20062023-------------------------
    
def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[conCfg.laser, Inst.CONTINUE, 0, 500*ms],
                       [conCfg.laser, Inst.BRANCH, 0, 500*ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...
    
    return instructionList

def close_all(*args):       # arguments are 'cam' and 'SG'
    # cam = args[0]
    if trial_run in ['n','N'] and len(args)==3:
        SG = args[2]
        if 'SG' in vars():
            SG.write('freq2.87ghz')
            SG.write('lcal')
            SGctrl.disable_SG_op(SG)
        # camctrl.close_camera(cam); print("\x1b[38;2;10;250;50m\x10 Camera Closed...")
    # PBctrl.pb_init(); PBctrl.pb_stop(); PBctrl.pb_close(); print("\x10 Pulse Blaster closed...\x1b[0m")
    # return DAQclosed

def read_save_details(roi, N_scanpts):
    """ defines the 'frame_per_cyc' to be captured, and the write formats of the 'scannedParam' and the 'data'
    Returns: 
    frames_per_cyc            - define the number of frames captured in one cycle of the sequence, 1 signal, 1 ref, so 2. 
    scannedparam_write_format - write format of the scannedparams in the data file
    data_write_format         - write format of the data, eta line by line er format ta.. not the whole at a time..
    """
    width = roi[2]; height = roi[3];
    frames_per_cyc = [2]
    
    scannedparam_write_format = ""
    for i in range(0, N_scanpts):
        scannedparam_write_format += "%g\t"
    scannedparam_write_format = scannedparam_write_format[0:-1] + "\n"

    data_write_format = ""
    # niche n_frames holo number of frames acquired in each cycle.. (dekhe tai mone hochhe, bhule gechhi)
    # replacing the command range(0, n_frames*width) with range(0, frames_per_cyc*width)
    for i in range(0, frames_per_cyc[0]*height):   # n_frames*width... na n_frames*height...
        data_write_format += "%d\t"
    data_write_format = data_write_format[0:-1] + "\n" # remove \t at the end and add \n
    # nicher 4 line e sample ta roilo.. 
    # write_format = ""
    # for i in range(0, n_frames*height):
    #     write_format += "%d\t"
    # write_format = write_format[0:-1]+"\n"
    return [frames_per_cyc, scannedparam_write_format, data_write_format]


def param_save_details(scan_time, exec_time):
    next_params_format = ' %s\t%0.2f\n %s\t%0.2f\n %s\t%0.4g\n %s\t%0.4g\n %s\t%g\n'  # %s\t%f\n
    # scan_time in ms; exec_time in seconds
    param_list = ['Max_scan_time(us):',max(scan_time)*1e3, 'Min_scan_time(us):',min(scan_time)*1e3, 'Total_scan_time(s):',(np.sum(scan_time)/1e3), 'Total_run_time(s):',exec_time, 'Step:', (expCfg.scannedParam[1]-expCfg.scannedParam[0])]
    # param_list[1] = i_scanpt+1..... Ekhane ki hbe?? Kon parameter save korbo??
    # ekhane ekta if kore, jodi 'i'=1 hoe, thle first param_list ta return korbe, nhle porer param_list ta return korbe.
    # Eta jodi kora hoe, thle runs>1 hole ba multi-param scan hole, notun param er sathe tar details save kora jabe...
    # Config file e ekta variable lagbe jeta dekhabe je kon variable ta scan hoechhe, other than scannedParam
    param_list = tuple(param_list)
    return [next_params_format, param_list]
    
def prepare_for_saving(savePath):
    global file_number
    if 'file_number' in vars():       # Returns a dict of all local variables
        print("Previous file number: \x1b[38;2;250;150;0m"+file_number+'\x1b[0m')
    file_number = input("File name:"+expCfg.saveFileName+"_#: ")
    
    paramfilename = savePath + expCfg.saveFileName + "_" + "params_" + file_number +".txt"
    datafilename = savePath + expCfg.saveFileName + "_" + file_number +".txt"
    if (isfile(datafilename)):
        print('\x1b[38;2;250;250;0mAppending file: '+expCfg.saveFileName+'_'+file_number+'.txt\x1b[0m')
        
    return [paramfilename, datafilename, file_number]

def save_data(datafilename, data_raw, data_write_format, scannedparam_write_format, param):
    print("\x10 Saving to file.....")
    # try:
    datafile = open(datafilename, 'a')
    # first e file tar first line e sob kota parameters save kora hok...
    datafile.write(scannedparam_write_format % tuple(i for i in param)) # param ekta list, protita element ke alada-alada kore ekta tuple
    # ebar line-by-line data (image frames) gulo lekha hok...
    for line in data_raw:
        datafile.write(data_write_format % tuple(line))
    # file lekha completed...
    datafile.close()
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
                seqArgList[0] = parameter[seq_no_plot[i]]
        plt.figure(dpi = plot_dpi)
        the_list = PBctrl.PB_program(instr,sequence,seqArgList)
        for j in range(0, len(the_list)):
            instructionList = the_list[j][0]
            inst_times = the_list[j][4]
            # instructionList=PBctrl.PB_program_camera(sequence,seqArgList)[0][0]
            # inst_times = PBctrl.PB_program_camera(sequence,seqArgList)[0][4]
            plt.subplot(1,2,(j+1))
            [t_us,channelPulses,yTicks] = seqctrl.plot_sequence(instructionList, PBchannels)
            for channel in channelPulses:
                plt.plot(t_us, list(channel))
            plt.yticks(yTicks, PBchannels.keys())          # Include the names of the PB channels
            plt.xlabel('Time (us)')
            plt.ylabel('Channel')
            plt.title(sequence + ' Plot. Param @ '+str(seq_no_plot[i])+': ' + str(parameter[seq_no_plot[i]]) + 'ns\nTransitions at: '+str([vals for vals in inst_times.values()]), fontsize=10)
    # Way to display the total no of instructions in the plot, since it will vary with different scanPts:
        # print("Total %d inst" % len(inst_times.keys()))

def acquire_data(t_exposure, t_seq_total, roi, Nsamples, parameter, ith_scan_pt, sequence, seqArgList, trial):
    # setup next scan iteration (e.g. for ESR experiment, change microwave frequency; for T2 experiment, reprogram pulseblaster with new delay)
    width=roi[2]; height=roi[3];
    instructionList=[]
    if sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
        if trial in ['n','N']:
            SGctrl.set_SG_freq(SG, parameter[ith_scan_pt])
    else:
        seqArgList[0] = parameter[ith_scan_pt]
    the_list = PBctrl.PB_program(sequence,seqArgList)
    # print(the_list)
    for i in range(0, len(the_list)):
        instructionList.append(the_list[i][0])
    # print(instructionList)
    PBctrl.run_sequence_for_camera(instructionList, t_exposure, t_seq_total)
    
    scan_start_time = time.perf_counter()   # time in seconds
    # frames = camctrl.capture(cam, Nsamples,width,height)    #read DAQ
    frames = np.zeros((width,height,Nsamples));
    for i in range(0,Nsamples):
        frame = np.random.rand(width, height)*2**16
        frame = frame.astype(int)
        frames[:,:,i] = frame       # numpy nd array
    
    scan_end_time = time.perf_counter()
    
    scan_time.append((scan_end_time - scan_start_time)*1e3)     # in ms
    # if ith_scan_pt==0:
    #     time.sleep(100)
    return [frames, scan_time, instructionList]

def plot_data(i_max, param, data_raw, x_label, x_unit, roi, n_frames, live=False):
    """Plot the signal, reference and the contrast"""
    # n_frames = total frames in one cycle of sequence..
    width = roi[2]; height = roi[3];
    mean_sig = np.zeros(i_max); mean_ref = np.zeros(i_max);
    contrast = np.zeros(i_max);
    for i in range(0,i_max):
        # signal_frames = np.array([data_raw[i*width:(i+1)*width,j*height:(j+1)*height] for j in range(0,n_frames,2)])
        # reference_frames = np.array([data_raw[i*width:(i+1)*width,j*height:(j+1)*height] for j in range(1,n_frames,2)])
        sum_of_signal_frames = np.array([np.sum(data_raw[i*width:(i+1)*width,j*height:(j+1)*height]) for j in range(0,n_frames,2)])
        sum_of_reference_frames = np.array([np.sum(data_raw[i*width:(i+1)*width,j*height:(j+1)*height]) for j in range(1,n_frames,2)])

        mean_sig[i] = np.mean(sum_of_signal_frames)
        mean_ref[i] = np.mean(sum_of_reference_frames)
        # define the contrast
        contrast[i] = mean_sig[i]/mean_ref[i];
    xValues = param[0:i_max]
    if live:
        plt.plot([x/x_unit for x in xValues], contrast[0:i_max], 'b--')
        plt.pause(0.0001)
    else:
        plt.figure()

        plt.subplot(121)
        plt.plot([x/x_unit for x in xValues], mean_sig[0:i_max], '.-',[x/x_unit for x in xValues], mean_ref[0:i_max], '.-')
        plt.legend(['Avg Sig','Avg Ref'])
        plt.xlabel(x_label)

        plt.subplot(122)
        plt.plot([x/x_unit for x in xValues],  contrast[0:i_max], '.-b')
        plt.xlabel(x_label); plt.ylabel('Contrast')
        plt.title('Plotting %d points' %(i_max))
    
# ----------------------------------------------------------------------------

# if __name__ == '__main__':
    # instructionList = main()

expCfg = import_module(expCfgFile)
# expCfg.N_scanPts = len(expCfg.scannedParam)
clk_cyc = 1e3/conCfg.PBclk      # One clk cycle of PB = inverse of the clk freq of PB
print('\x10 \x1b[38;2;250;250;0mRunning '+expCfg.saveFileName+' sequence\x1b[0m')
[SG, cam] = initialize_instr(expCfg.sequence)
# cam["exposure_time"] = exposure_time
seqctrl.check_params(expCfgFile)
print("\x10 Parameter checks completed...")
instr = 'cam'
#-------------------------- 20062023-------------------------
[roi, savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr,expCfg)
roi = [int(element/4)*4 for element in roi];
# width = roi[2]; height = roi[3];
width = 8; height =12;
roi[2] = 8; roi[3] = 12;
#-------------------------- 20062023-------------------------
#%%
t_manip=0
if not (expCfg.sequence == 'esr_seq'):
    t_manip = np.array(param)
t_seq_total0 = expCfg.t_AOM + t_manip
t_seq_total1 = expCfg.t_AOM + t_manip
t_seq_total = np.transpose(np.array([t_seq_total0, t_seq_total1]))

# i am putting this in the initialize_exp() function --- may delete this later .. 18-06-2023
# savePath = expCfg.savePath + time.strftime("%Y-%m-%d", time.localtime()) + '\\'
# if not (isdir(savePath)): # "Saved Data" folder is created in the pwd, if it does not exist
#     makedirs(savePath)
#     print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')

#-------------------------- 20062023-------------------------

[frames_per_cyc, scannedparam_write_format, data_write_format] = read_save_details(roi, Nscanpts)

#-------------------------- 20062023-------------------------

# expCfg.Nsamples = 1 mane holo at each scan pt 1 frame (like N APD readouts for averaging). So multiply by frames_per_cyc[0]
# thle Nsamples = 2 hbe... at each scan pt 2 frames.. then go to next scan pt..
# Nsamples = ekta scanpt e total kotogulo signal (ba reference) frames...
Nsamples = frames_per_cyc[0]*expCfg.Nsamples        # feed this as the 'n_frames' in capture() of Camcontrol.py
save_flag=False
# configure the camera for triggered acquisition
# camctrl.configure_camera(cam)
# set the ROI for the acquisition
# camctrl.adjust_roi(cam,roi,True)
# t_exposure = cam["exposure_time"].value     # min the value... in seconds !!!
t_exposure = 0.1
display_parameters = dialog.yesno_box(expCfg.saveFileName+' Params', 'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' %(str(conCfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples, param[0]/expCfg.plotXaxisUnits, param[-1]/expCfg.plotXaxisUnits))

if display_parameters == 'yes':
    try:
        # continue_run = 'yes'
        # DAQclosed = False;
        # save_flag = False
        # DAQtask = DAQctrl.configure_daq(Nsamples)        # configure_daq() accepts no. of samples to read from DAQ-AI & returns the task created; variable 'readTask'
        print("\x10 Camera configured for %d frames..." % expCfg.Nsamples)  # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
        print('\x10 %d frames in each cycles...' % frames_per_cyc[0])
        
        # <<<<<<<<---------------------Run experiment--------------------->>>>>>>>
        print("------Acquisition Started------")
        for i_run in range(0, expCfg.Nruns):
            data_raw = np.zeros((width*Nscanpts,height*Nsamples));      # check the dimension... width, height, Nscanpts, Nsamples...?
            scan_time = []   # time (in ms) for each scannedParam
            # print("\x10",i_run+1,'/',expCfg.Nruns)
            start_time = time.perf_counter()        # in seconds
            for i_scanpt in range (0, Nscanpts):
                print(i_scanpt+1,' / ',Nscanpts,': ',param[i_scanpt])
                 
                [frames, scan_time, instructionList] = acquire_data(t_exposure*1e9, t_seq_total[i_scanpt], roi, Nsamples, param, i_scanpt, expCfg.sequence, seqArgList, trial_run)
                data_raw[i_scanpt*width:(i_scanpt+1)*width,:] = np.concatenate(tuple(frames[:,:,j] for j in range(0,frames_per_cyc[0])), axis=1)       # 2d-array

                # nicher line gulo thaklo sudhu reference er jonne.. diode ar camera r code unify korar somoy lagbe...------------
                # [mean_sig, mean_bg, contrast, data_raw] = process_data(frame, reads_per_cyc, param, i_scanpt, mean_sig, mean_bg, contrast, data_raw, len(conCfg.input_terminals))

                # def process_data(cts, reads_per_cyc, parameter, ith_scan_pt, mean_sig, mean_bg, contrast, data_raw, n_channels):
                #     # Process data for saving    
                #     data_raw = np.concatenate(tuple(np.concatenate(tuple(all_frames[j][:,:,i] for i in range(0,n_frames)), axis=1) for j in range(0,N_scanpts)), axis=0)
                #     # Process data for plotting
                #     mean_sig[ith_scan_pt] = np.mean(cts[0][0::reads_per_cyc[0]]) if n_channels==2 else np.mean(cts[0::reads_per_cyc[0]])
                #     mean_bg[ith_scan_pt] = np.mean(cts[0][1::reads_per_cyc[0]]) if n_channels==2 else np.mean(cts[1::reads_per_cyc[0]])
                #     contrast.append(calc_contrast(mean_sig[ith_scan_pt], mean_bg[ith_scan_pt], 's/r'))
                #     return [mean_sig, mean_bg, contrast, data_raw]
                # above lines required ony for diode measurements----------------------------------------------------------

                if livePlotUpdate == 'yes':
                    plot_data(i_scanpt+1, param, data_raw, expCfg.xAxisLabel, expCfg.plotXaxisUnits, live=True)                
            end_time = time.perf_counter()          # in seconds
            exec_time = end_time - start_time       # in seconds
            plot_data(i_scanpt, param, data_raw, expCfg.xAxisLabel, expCfg.plotXaxisUnits, roi, Nsamples, live=False)
            print("\x10 Execution time = %.2fs" % exec_time)        # in seconds
            print("\x10 Scan time = %.2fs" % (np.sum(scan_time)/1e3))   # in seconds; scan_time in ms
            # Close all if this is the last run
            if i_run+1 == expCfg.Nruns:
                # close_all(cam, SG) if trial_run in ['n','N'] and expCfg.sequence not in ['T1ms0','aom_timing','rodelay'] else close_all(cam)
                a=1
            
            savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
            if savefile_yn == 'yes':
                # Data file saving
                if save_flag == False:          # Ask for file number iff no save was performed
                    [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
                save_flag = save_data(datafilename, data_raw, data_write_format, scannedparam_write_format, param)
            else:
                print("\x10 \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")
            if i_run+1 < expCfg.Nruns:    
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
        plot_data(i_scanpt, param, data_raw, expCfg.xAxisLabel, expCfg.plotXaxisUnits, live=False)
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
        # close_all(cam, SG) if trial_run in ['n','N'] and expCfg.sequence not in ['T1ms0','aom_timing','rodelay'] else close_all(cam)
        
        print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m frames at each pt." % (frames_per_cyc[0], expCfg.Nsamples))
        
        print("\x10 len(frames) = \x1b[38;2;250;150;50m%d\x1b[0m" % len(frames))
        # Save parameters (only if there was one save operation)...
        if save_flag:
            expParamList[1] = i_scanpt+1        # expParamList[1] -> value of N_scanPts
            expParamList[3] = i_run+1           # expParamList[3] -> value of Nruns
            save_parameters(paramfilename, param_save_format, expParamList)
            
            [next_params_format, expParamList] = param_save_details(scan_time, exec_time)
            save_parameters(paramfilename, next_params_format, expParamList)
        
        # scan_time in ms; exec_time in seconds
        #----------------------------------------------------------------------
        
        # matlab_engine = eng.start_matlab()
        
        # if expCfg.sequence == 'esr':
        #     analysis = eng.ESR_raw_analysis()
        
else:
    print("Dialog closed...")
    # close_all(cam, SG) if trial_run in ['n','N'] and expCfg.sequence not in ['T1ms0','aom_timing','rodelay'] else close_all(cam)
    # sys.exit()
# return instructionList

# %%
