#%% Initialization and Definition
# reset
import DAQcontrol as daqctrl, PBcontrol_v2 as pbctrl, sequencecontrol as seqctrl, SGcontrol as sgctrl, connectionConfig as concfg, matplotlib.pyplot as plt, numpy as np, sys, time, dialog, os, ametekcontrol as amctrl
from spinapi import ns, us, ms, Inst
from os.path import isdir, isfile; from os import makedirs; from importlib import import_module

# from matlab import engine as eng
# %matplotlib qt5

# def main():
print("\x10 \x1b[0mImports Successful...")
plt.rcParams.update({'figure.max_open_warning': 0})     # No warnings on opening mult fig windows

global expCfgFile, trial_run, seq_no_plot, voltage_unit, seq_plot_dpi, plotPulseSequence, clk_cyc, SG, expCfg
expCfgFile = 'rabi'+'_config'
test_field = [0,0,10]


trial_run = ['n','n']       # 1st=SG, 2nd=PB, 3rd=ametek
seq_no_plot = [-1]
voltage_unit = 1      # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100                      # The dpi of the displayed pulse sequence plot
plotPulseSequence = True

# sweep_value = [1*us]
# sweep_step = 100*ns
# scannedParam2_max
# scannedParam2_step
# scannedParam2 = np.linspace(scannedParam2_min, scannedParam2_max, scannedParam2_step)
# save_raw_data = True            # Try shot-by-shot normalization

def initialize_instr(sequence):
    if trial_run[1] == 'n':
        try:
            pbctrl.configurePB()
            print('\x10 PB: \x1b[38;2;250;250;0mv'+pbctrl.pb_get_version()+'\x1b[0m')   # Display the PB board version using pb_get_version()
        except:
            print("Error Initializing PB !!")
    
    if trial_run[0] == 'n' and sequence not in ['aom_timing', 'rodelay']:
        # Do not initialize SG if it is a trial run or the sequence is present in the list ['aom_timing', 'rodelay']
        SG = sgctrl.init_sg(concfg.serialaddr,concfg.model_name)
        if SG != '':
            sgctrl.enable_SG_op();
            print("SG Output Enabled...")
            sgctrl.set_SG_amp(expCfg.MW_power)
            sgctrl.setup_SG_pulse_mod();
            print("SG Ext Pulse Mod Enabled...")
    else:
        SG = None
    # if trial_run[1] == 'n':
    #     ametek = amctrl.init_ametek();
    #     if not(ametek.query("ID") == ''):
    #         print("\x10 \x1b[38;2;250;250;0m7270 Initialized...\x1b[0m")
    return [SG]


def close_all(ai_task, SG, ao_task):   # arguments are closed, ai_task and SG
    """End the measurement
    Closes all the instruments.
    """
    if (ai_task is not None):    # Close DAQ task
        daqctrl.close_daq_task(ai_task)
        print("AI stopped...")
    
    if (ao_task is not None):
        daqctrl.start_ao(ao_task,[0,0,0])
        daqctrl.close_daq_task(ao_task)
        print('AO stopped...')
    
    if (trial_run[0] == 'n') and (SG is not None):
        SG.write('freq2.87ghz')
        sgctrl.disable_SG_op()
        sgctrl.uninit_sg()

    # pbctrl.pb_init();
    closed = True
    pbctrl.pb_stop();
    pbctrl.pb_close();
    print("Pulse Blaster closed...\x1b[0m")
    return True


def initialize_exp(instr, expCfg):
    ###### Include trial run check so that the error plots are only displayed when it is not a trial run
    # expCfg.N_scanPts = len(expCfg.scannedParam) # protection against non-integer user inputs for Nscanpts.
    # --------------------------------------------------------------
    # sequenceArgs = []
    # expParamList = []
    
    # Nscanpts = expCfg.N_scanPts
    param = expCfg.scannedParam;    n_error=0
    param_save_format = expCfg.formattingSaveString
    sequenceArgs = expCfg.updateSequenceArgs()      # Variables used in the pulse sequence
    expParamList = expCfg.updateExpParamList()      # List of experimental parameters
    
    # ------------------------------------------------------------------
        # est_time = expCfg.t_tot * expCfg.Nsamples * expCfg.N_scanPts ??
    if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:   # for sequences except ESR, set MW frequency
        seqArgList = [param[seq_no_plot[-1]]]
        seqArgList.extend(sequenceArgs)     # Make a seqArgList with dummy 1st element... just to create it.. pore change hoye jabe..
        [n_error, param] = seqctrl.param_err_check(instr, expCfg.sequence, expCfg.PBchannels, seqArgList,  expCfg.scannedParam, expCfg.N_scanPts)
        if n_error>0:
            print('\x1b[1;37;41m'+'Err: Check Sequences...\x1b[0m')
            print('\x1b[38;2;250;0;0m'+str(n_error)+' parameters removed...\x1b[0m')

            # # Close Error plots??
            # close_plots = dialog.yesno_box('Close Plots', 'Close the Error Plots?')
            # if close_plots == 'yes':
            #     plt.close('all')
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
    
    # instructionList = pbctrl.PB_program(expCfg.sequence,seqArgList)[0]        # ekhane change chhilo.. (24/02/23)

    # Start the initial sequence now ------------
    print("Starting Initial Sequence...")
    if trial_run[1] == 'n':
        instructionList = [start_initial_PB_seq()]
        pbctrl.run_sequence_for_diode(instructionList)
        print("\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
        if expCfg.sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq','rodelay'] and trial_run[0] == 'n':
            sgctrl.set_SG_freq(expCfg.MW_freq)
    else:
        # instructionList = []
        print("\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    # Create the data save folder...
    cwd = os.getcwd()
    savePath = cwd[0:cwd.rfind('\\')]+ "\\Saved_Data\\" + time.strftime("%Y-%m-%d", time.localtime()) + '\\'         # Path to folder where data will be saved
    if not (isdir(savePath)):
        os.makedirs(savePath)
    print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')
    return [savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList]
    
def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[concfg.laser, Inst.CONTINUE, 0, 500*ms],
                       [concfg.laser, Inst.BRANCH, 0, 500*ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...
    
    return instructionList

def read_save_details(sequence, n_channels, Nsamples_expCfg):
    # Enter no of [<apd>,<pd>] signals acquired in each cycle
    if sequence in ['modesr']:#, 'esr_seq']:
        reads_per_cyc = [4,1] if n_channels==2 else [4]
    elif sequence in ['drift_seq']:
        reads_per_cyc = [1,1] if n_channels==2 else [1]
    elif 'dig_mod' in sequence:
        reads_per_cyc = [expCfg.N_laser*2*2]      # not acquiring from 2 channels (APD, PD) at the moment
    else:
        reads_per_cyc = [2,2] if n_channels==2 else [2]
    Nsamples = sum(reads_per_cyc)*Nsamples_expCfg

    data_write_format = "%g\t"
    if 'dig_mod' in sequence:
        data_write_format += "%0.8f\t%0.8f\n"
    else:
        data_write_format += "%0.3f\t" * Nsamples
        # for i in range(0, sum(reads_per_cyc)):
        #     data_write_format += "%0.3f\t"
        data_write_format = data_write_format[0:-1] + "\n" # remove \t at the end and add \n
    return [reads_per_cyc, Nsamples, data_write_format]

def extra_param_save_details(scan_time_list, exec_time):
    extra_params_format = ' %s\t%0.2f\n %s\t%0.2f\n %s\t%0.4g\n %s\t%0.4g\n %s\t%g\n'  # %s\t%f\n
    # scan_time_list in seconds; exec_time in seconds
    param_list = ['Max_scan_time(us):',max(scan_time_list)*1e6, 'Min_scan_time(us):',min(scan_time_list)*1e6, 'Total_scan_time(s):',np.sum(scan_time_list), 'Total_run_time(s):',exec_time, 'Step:', (expCfg.scannedParam[1]-expCfg.scannedParam[0])]
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
    file_number = input("File name:"+expCfg.saveFileName+"_#: ")
    
    paramfilename = savePath + expCfg.saveFileName + "_" + "params_" + file_number +".txt"
    datafilename = savePath + expCfg.saveFileName + "_" + file_number +".txt"
    if (isfile(datafilename)):
        print('\x1b[38;2;250;250;0mAppending file: '+expCfg.saveFileName+'_'+file_number+'.txt\x1b[0m')
        
    return [paramfilename, datafilename, file_number]

def save_data_txt_txt(datafilename, data_raw, data_write_format):
    print("\x10 Saving.....")
    with open(datafilename, 'a') as datafile:
        datafile.write("%s\n" % (expCfg.saveFileName + file_number))
        for line in data_raw:
            datafile.write(data_write_format % tuple(line))
    
    print("\x10 Data saved to\x1b[38;2;100;250;50m %s_%s\x1b[0m !!!" % (expCfg.saveFileName, file_number))
    return True
       
# def save_data_txt_txt_for_mod(n_channels, datafilename, data_raw, data_write_format, Nsamples, reads_per_cyc):
#     print("\x10 Saving.....")
#     # try:
#     datafile = open(datafilename, 'a')
#     datafile.write("--%s--\n" % datafilename[-5:])
#     for lines in data_raw:      # Each 'lines' in data_raw stands for a specific scanpt...
#         for i in range(0,expCfg.N_laser*expCfg.Nsamples,1):
#             line = [lines[0], lines[1][i], lines[2][i]]
#             datafile.write(data_write_format % tuple(line))
#         for i in range(0, expCfg.N_laser*expCfg.Nsamples,1):
#             line = [lines[0], lines[3][i], lines[4][i]]
#             datafile.write(data_write_format % tuple(line))
#     datafile.close()
#     print("\x10 Data saved to\x1b[38;2;100;250;50m %s_%s\x1b[0m !!!" % (expCfg.saveFileName, file_number))
#     return True
    
def save_parameters(paramfilename, param_save_format, param_list):
    paramFile = open(paramfilename, 'a')
    paramFile.write(param_save_format % tuple(param_list))
    paramFile.close()

def calc_contrast(signal, reference, op):
    if op == '+-':
        contrast = signal - reference
    elif op == '-+':
        contrast = -signal + reference
    elif op == 's/r':
        contrast = signal/reference
    elif op == 'r/s':
        contrast = reference/signal
    return contrast
    
def view_sequence(sequence, PBchannels, seqArgList, only_plot=False, parameter=[0], seq_no_plot=[0], plot_dpi=100):
    for i in range(0, len(seq_no_plot)):
        if only_plot == False:
            if  sequence not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
                seqArgList[0] = parameter[seq_no_plot[i]]
        plt.figure(num=f"{sequence} sequence plot", dpi = plot_dpi)
        the_list = pbctrl.PB_program(instr,sequence,seqArgList)
        for j in range(0, len(the_list)):
            instructionList = the_list[j][0]
            inst_times = the_list[j][4]
            # instructionList=pbctrl.PB_program_camera(sequence,seqArgList)[0][0]
            # inst_times = pbctrl.PB_program_camera(sequence,seqArgList)[0][4]
            # plt.subplot(1,2,(j+1))
            [t_us,channelPulses,yTicks] = seqctrl.plot_sequence(instructionList, PBchannels)
            for channel in channelPulses:
                plt.plot(t_us, list(channel))
            plt.yticks(yTicks, PBchannels.keys())          # Include the names of the PB channels
            plt.xlabel('Time (us)')
            plt.ylabel('Channel')
            plt.title(sequence + ' Pulse Seq plot. Param val @ '+str(seq_no_plot[i])+': ' + str(parameter[seq_no_plot[i]]) + 'ns\nTransitions at: '+str([vals for vals in inst_times.values()]), fontsize=12)
    # Way to display the total no of instructions in the plot, since it will vary with different scanPts:
        # print("Total %d inst" % len(inst_times.keys()))

def acquire_data(Nsamples, parameter, ai_task, sequence, seqArgList, trial):
    # setup next scan iteration (e.g. for ESR experiment, change microwave frequency; for T2 experiment, reprogram pulseblaster with new delay)
    if sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
        if trial[0] == 'n':
            sgctrl.set_SG_freq(parameter)
    else:
        seqArgList[0] = parameter
    # instructionList = pbctrl.PB_program(instr,sequence,seqArgList)[0]

    # #notun ... eta lagbe
    instructionList=[]
    the_list = pbctrl.PB_program(instr,sequence,seqArgList)
    # print(the_list)
    for i in range(0, len(the_list)):
        instructionList.append(the_list[i][0])
    # print(instructionList)

    # start = time.perf_counter_ns()
    # print(i_scanpt+1,' / ',Nscanpts)
    # if i_scanpt == 1:
    #     time.sleep(20)
    # stop = time.perf_counter_ns()

    pbctrl.run_sequence_for_diode(instructionList)
    # cts = []
    scan_start_time = time.perf_counter()   # time in seconds
    cts = daqctrl.read_daq(ai_task, Nsamples,61*60)    #read DAQ
    scan_end_time = time.perf_counter()
    
    scan_time = (scan_end_time - scan_start_time)     # in seconds
    return [cts, scan_time]

def process_data(i_max, reads_per_cyc, data_raw):
    # note that this function can be made simpler for multi-channel acquisition
    # Process data for plotting
    mean_sig = np.zeros((i_max, len(concfg.input_terminals)));
    mean_ref = np.zeros((i_max, len(concfg.input_terminals)));
    contrast = np.zeros((i_max, len(concfg.input_terminals)));
    
    signals_x = data_raw[0:i_max,[i for i in range(0, expCfg.Nsamples*reads_per_cyc[0], reads_per_cyc[0])]]
    references_x = data_raw[0:i_max,[i for i in range(1, expCfg.Nsamples*reads_per_cyc[0], reads_per_cyc[0])]]
    mean_sig[:,0] = np.mean(signals_x,axis=1)
    mean_ref[:,0] = np.mean(references_x,axis=1)
    
    contrast[:,0] = calc_contrast(mean_sig[:,0], mean_ref[:,0], 's/r')
    
    if len(concfg.input_terminals)>1:
        signals_y = data_raw[0:i_max,[i for i in range(expCfg.Nsamples*reads_per_cyc[1], expCfg.Nsamples*reads_per_cyc[1]*2, reads_per_cyc[1])]]
        references_y = data_raw[0:i_max,[i for i in range(expCfg.Nsamples*reads_per_cyc[1]+1, expCfg.Nsamples*reads_per_cyc[1]*2, reads_per_cyc[1])]]
        mean_sig[:,1] = np.mean(signals_y,axis=1)
        mean_ref[:,1] = np.mean(references_y,axis=1)
        
        contrast[:,1] = calc_contrast(mean_sig[:,1], mean_ref[:,1], 's/r')
    
    return [mean_sig, mean_ref, contrast]

# def process_data_for_mod(cts, reads_per_cyc, parameter, ith_scan_pt, mean_sig, mean_bg, contrast, data_raw):
#     MW_off_data = []; MW_on_data = [];
#     for i in range(0,expCfg.Nsamples,1):        # expCfg.Nsamples = N_MW
#         MW_off_data += cts[4*i*expCfg.N_laser:(4*i+2)*expCfg.N_laser:1]   # laser mod data at 0 to 2*N_laser, this includes the reference data also...
#         MW_on_data += cts[(4*i+2)*expCfg.N_laser:(4*i+4)*expCfg.N_laser:1]      # MW mod data at 2*N_laser to 4*N_laser, including the reference data also...
    
#     sig_wo_mw = MW_off_data[0::2];  ref_wo_mw = MW_off_data[1::2];
#     sig_w_mw = MW_on_data[0::2];    ref_w_mw = MW_on_data[1::2];
    
#     data_raw.append([parameter[ith_scan_pt], sig_wo_mw, ref_wo_mw, sig_w_mw, ref_w_mw])
    
#     # Process for plotting
#     mean_sig[ith_scan_pt] = np.mean(sig_w_mw)
#     mean_bg[ith_scan_pt] = np.mean(sig_wo_mw)
#     contrast.append(calc_contrast(mean_sig[ith_scan_pt], mean_bg[ith_scan_pt], 's/r'))
    
#     return [mean_sig, mean_bg, contrast, data_raw]

def plot_data(i_max, param, processed_data, x_label, x_unit, live=False):
    [mean_sig, mean_ref, contrast] = processed_data
    xValues = param[0:i_max]
    if live:
        plt.plot([x/x_unit for x in xValues], contrast[0:i_max], 'b--')
        plt.pause(0.0001)
    else:
        plt.figure()
        plt.subplot(121)
        plt.plot([x/x_unit for x in xValues], mean_sig[0:i_max], '.-',[x/x_unit for x in xValues], mean_ref[0:i_max], '.-')
        plt.legend(['Sig X', 'Sig Y', 'Ref X', 'Ref Y']) if len(concfg.input_terminals)>1 else plt.legend(['Sig','Ref']) 
        plt.xlabel(x_label)
        plt.subplot(122)
        plt.plot([x/x_unit for x in xValues],  contrast[0:i_max])
        plt.xlabel(x_label);    plt.ylabel('Contrast')
        plt.legend(['Contrast X', 'Contrast Y']) if len(concfg.input_terminals)>1 else plt.legend() 
        plt.title('Plotting %d points' %(i_max))

# ----------------------------------------------------------------------------

# if __name__ == '__main__':
    # instructionList = main()

expCfg = import_module(expCfgFile)
# expCfg.N_scanPts = len(expCfg.scannedParam)
clk_cyc = 1e3/concfg.PBclk      # One clk cycle of PB = inverse of the clk freq of PB
print('\x10 \x1b[38;2;250;250;0mRunning '+expCfg.saveFileName+' sequence\x1b[0m')
[SG] = initialize_instr(expCfg.sequence)

seqctrl.check_params(expCfgFile)
print("\x10 Parameter checks completed...")
instr = 'diode'

# expCfg.Nruns = Nruns #..........eta notun file e acche... kno??
[savePath, param_save_format, seqArgList, expParamList, Nscanpts, param, instructionList] = initialize_exp(instr, expCfg)

# Experiment operations (with all hardwares working)... ki korbo eta???
# if trial_run == 'y':
#     print("Cannot proceed further")

[reads_per_cyc, Nsamples, data_write_format] = read_save_details(expCfg.sequence, len(concfg.input_terminals), expCfg.Nsamples)

Nsamples = reads_per_cyc[0]*expCfg.Nsamples     # Nsamples is controlled by the no of apd signals i.e. reads_per_cyc[0] (even though 1 pd signal is saved, we have to acquire all, at the time when apd signal is acquired); this is sent to configure_DAQ

display_parameters = dialog.yesno_box(expCfg.saveFileName+' Params', 'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' %(str(concfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples, param[0]/expCfg.plotXaxisUnits, param[-1]/expCfg.plotXaxisUnits))
closed = False
if display_parameters == 'yes':    
    try:
        # if trial_run == 'n':
        ao_task = daqctrl.config_ao(dev="U9263")
        daqctrl.start_ao(ao_task, daqctrl.coil_calibration(test_field))
        save_flag = False
        pbctrl.run_sequence_for_diode(instructionList)        # Run the PB hardware
        print("Initial Sequence Started...")
        # continue_run = 'yes'
        closed = False;  ai_task = daqctrl.config_ai(Nsamples)        # configure_daq() accepts no. of samples to read from DAQ-AI & returns the task created; variable 'readTask'
        print("\x10 DAQ configured for %d samples..." % expCfg.Nsamples)
        print('\x10 %d samples in each cycles...' % reads_per_cyc[0])
        
        # <<<<<<<<---------------------Run experiment--------------------->>>>>>>>
        print("------Acquisition Started------")
        for i_run in range(0, expCfg.Nruns):

            data_raw = np.zeros((Nscanpts,Nsamples*len(concfg.input_terminals)))
            scan_time_list = []   # time (in seconds) for each scannedParam
            # print("\x10",i_run+1,'/',expCfg.Nruns)
            start_time = time.perf_counter()        # in seconds
            for i_scanpt in range (0, Nscanpts):
                # if i_scanpt>0:
                #     break
                print(i_scanpt+1,' / ',Nscanpts,': ',param[i_scanpt])
                [cts, scan_time] = acquire_data(Nsamples, param[i_scanpt], ai_task, expCfg.sequence, seqArgList, trial_run)
                # in general for multi-channel acquisition, cts is a list of lists... this needs to be taken care
                # the sig1,ref1,sig2,ref2,... order has not been changed cts and data_raw
                data_raw[i_scanpt,:] = np.ravel(np.array(cts))
                scan_time_list.append(scan_time)
                # time.sleep(100e-3)
                # savefile_yn = dialog.yesno_box('Proceed',"Proceed with scan?")
                # if savefile_yn == 'yes':
                #     continue
                # else:
                #     break      
            end_time = time.perf_counter()          # in seconds
            exec_time = end_time - start_time       # in seconds
            processed_data = process_data(i_scanpt+1, reads_per_cyc, data_raw)    # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
            plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, live=False)
            print("\x10 Execution time = %.2fs" % exec_time)        # in seconds
            print("\x10 Scan time = %.2fs" % (np.sum(scan_time_list)))   # in seconds

            data_raw = np.insert(data_raw,0,param[0:i_scanpt+1],axis=1)     # insert the parameter value at the head of the array

            if i_run+1 == expCfg.Nruns: # Close all if this is the last run
                closed = close_all(ai_task, SG, ao_task)
            savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
            if savefile_yn == 'yes':
                # Data file saving
                if save_flag == False:          # Ask for file number iff no save was performed
                    [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
                save_flag = save_data_txt_txt(datafilename, data_raw, data_write_format)
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
        processed_data = process_data(i_scanpt, reads_per_cyc, data_raw)    # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
        plot_data(i_scanpt, param, processed_data, expCfg.xAxisLabel, expCfg.plotXaxisUnits, live=False)
        # Then ask whether to save it...
        savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
        if savefile_yn == 'yes':
            # Ask for filename only if there was no save operation
            if save_flag == False:
                [paramfilename, datafilename, file_number] = prepare_for_saving(savePath)
            # A conditional save statement.. Different for dig_mod sequences...
            save_flag = save_data_txt_txt(datafilename, data_raw, data_write_format)
        # save_data_txt(len(concfg.input_terminals), datafilename, data_raw, data_write_format, expCfg.Nsamples, reads_per_cyc)
    # else:        
    #     if (i_run+1)==expCfg.Nruns: # Save data of the final run here
    #         save_data_txt(len(concfg.input_terminals), datafilename, data_raw, data_write_format, expCfg.Nsamples, reads_per_cyc)
    finally:
        if not closed:
            closed = close_all(ai_task, SG, ao_task)
        # close AMETEK
         # amctrl.stop_oscillator(ametek)
        
        print("\x10 Read \x1b[38;2;250;150;50m%d*%d\x1b[0m samples at each pt." % (reads_per_cyc[0], expCfg.Nsamples))
        if save_flag:   # below line should "not" run if there was "no" save operation / data acquisition...
            print("\x10 len(cts[0]) = \x1b[38;2;250;150;50m%d\x1b[0m" % len(cts[0])) if len(concfg.input_terminals) == 2 else print("\x10 len(cts) = \x1b[38;2;250;150;50m%d\x1b[0m" % len(cts))
        # Save parameters only if there was one save operation
        if save_flag:
            expParamList[1] = i_scanpt+1        # expParamList[1] -> value of N_scanPts
            expParamList[3] = i_run+1           # expParamList[3] -> value of Nruns
            save_parameters(paramfilename, param_save_format, expParamList)
            
            [next_params_format, expParamList] = extra_param_save_details(scan_time_list, exec_time)
            save_parameters(paramfilename, next_params_format, expParamList)
        
        # scan_time in seconds; exec_time in seconds
        #----------------------------------------------------------------------
        
        # matlab_engine = eng.start_matlab()
        
        # if expCfg.sequence == 'esr':
        #     analysis = eng.ESR_raw_analysis()
        
else:
    print("Measurement aborted...")
    # if trial_run[0] == 'n' and expCfg.sequence not in ['rodelay']:
    print(closed)
    if not closed:
        closed = close_all(None, SG, None)
    # sys.exit()
# return instructionList
