#%% Initialization and Definition
# reset
from fileinput import filename
import connectionConfig as concfg, matplotlib.pyplot as plt, numpy as np, time
import dialog, psutil, json, yaml, os, logging, matplotlib as mpl, h5py
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogInputTask, AnalogOutputTask#, DAQ_write_pattern
from sequencecontrol import sequencecontrol
from SGcontrol import SignalGenerator, SignalGenerator_sim
from spinapi import ns, us, ms, Inst
from pathlib import Path
from importlib import import_module
# %matplotlib qt5
plt.rcParams.update({'figure.max_open_warning': 0})     # No warnings on opening mult fig windows
# def main():
plt.style.use('dark_background')
plt.rcParams['axes.prop_cycle'] = mpl.rcParamsOrig['axes.prop_cycle']

global trial_run, f_number, instr, data, pb, ao_task, ai_task#, data_array_time
expCfgFile = 'esr'+'_config'
expCfg = import_module(expCfgFile)
params = expCfg.params_dict

params['test_field'] = [25, 75, 0]        # 0.001 G = 100 nT
# params['test_field'] = [0, 0, 0]

trial_run = ['y','n']       # 1st=SG, 2nd=PB, 3rd=ametek
seq_no_plot = [-1]
voltage_unit = 1      # mV voltage... Convert the voltages in cts to mV unit
seq_plot_dpi = 100                      # The dpi of the displayed pulse sequence plot
plotPulseSequence = True

load_pb_all_params = True if 'train' in expCfgFile else False
reload_pb = False if 't1ms0_train' in expCfgFile else True

# TODO: set up the logger ???
fsplit = lambda b: (2870 - 2.8*b, 2870 + 2.8*b)
# test_field = test_field*np.array([np.sin(test_theta*np.pi/180)*np.cos(test_phi*np.pi/180), np.sin(test_theta*np.pi/180)*np.sin(test_phi*np.pi/180), np.cos(test_theta*np.pi/180)])

# TODO: save a png of the plot in a common folder with file name
# TODO: save a h5 file with the parameter - far better than putting the parameters in the yaml file!!
seqctrl = sequencecontrol(params)
ai_task: AnalogInputTask
def initialize_instr(sequence):
    global ao_task, pb #,ai_task, 
    # TODO: how to handle 'sequence' to PulseBlaster()
    sg = SignalGenerator_sim()
    
    pb = PulseBlaster(params)               # pb not configured in __init__()
    # pb.set_sequencecontrol(seqctrl)         # Set the sequence control object to PulseBlaster
    # seqctrl.set_pbcontrol(pb)            # Set the pulseblaster object to sequence control

    try:
        pb.configure()
    except Exception as e:
        logging.exception(f"❌ Error in PB Init: {e}")

    try:
        ao_task = AnalogOutputTask(dev="P6363", channels=[0,1,2], coil='small_confocal')    # AO alreay configured here in __init__()
        # ao_task = None
        # TODO: need to do this properly - how to init the ai task
        # ai_task = AnalogInputTask(sampling_rate=concfg.daq_max_samp_rate, voltage_range=(-10, 10), channels=concfg.input_terminals, trigger_source=concfg.samp_clk_terminal)
        
    except Exception as e:
        logging.exception(f"❌ Error in DAQ Init: {e}")
    try:
        if sequence not in ['aom_timing', 'T1ms0_train']: #trial_run[0] == 'n' and 
            # Do not initialize sg if it is a trial run or the sequence is present in the list ['aom_timing',]
            sg = SignalGenerator() if trial_run[0] == 'n' else sg
            if sg != '':
                sg.enable_ntype(1)
                print("✔ SG Output Enabled...")
                sg.set_amp(expCfg.MW_power)
                sg.set_freq(params['mw']['freq'][0])
                # sg.set_freq(2.7327e9)
                # sg.setup_ext_pulse_mod()
                sg.setup_sg_fm('external', dev=5e5)
                # sg.setup_sg_am('external', dev=100)
                sg.enable_modulation(1)
                print("✔ SG Mod Enabled...")
        # else:
        #     sg = None
        # elif trial_run[0] == 'y':
        # sg = SignalGenerator_sim(simulate=True)
    except Exception as e:
        # print(f"❌ Error in SG Init: {e}")
        logging.exception(f"❌ Error in SG Init: {e}")
        sg = SignalGenerator_sim()
    # control_sequences = 
    devs = sg, ao_task#, ai_task
    return devs

def set_core_affinity():
    # Assuming you want cores 0,1 for instrument control
    process = psutil.Process(os.getpid())
    process.cpu_affinity([0, 1])
    print(f"✔ Process pinned to cores: {process.cpu_affinity()}")

def close_all(sg=None, ao_task:AnalogOutputTask=None, ai_task:AnalogInputTask=None):
    """
    End the measurement. Closes all the instruments.
    """
    try:
        if (ao_task is not None):
            print(f"🔃 AO closing:: status: {ao_task._task_state}")
            ao_task.set_outputs_to_constant([0,0,0])
            time.sleep(0.5)
            ao_task.__del__()
            print('✔ AO stopped...')
        if ai_task is not None:
            ai_task.__del__()
            print('✔ AI stopped...')
        
        # if (trial_run[0] == 'n') and (sg is not None):
        if sg is not None:
            # sg.set_freq(2.87e9)
            # sg.enable_ntype(0)
            sg.uninit()

        # pb.pb_init();
        # if t_align_dc < 50:
        #     pb.run_only_daq(250 *ms)        # let t_align_dc=250
        # else:
        #     pb.run_only_daq(t_align_dc *ms)
        closed = True
        # pb.stop_sequence();
        pb.closePB();
        print("✔ Pulse Blaster closed...\x1b[0m")
        return True
    except Exception as e:
        logging.exception(f"❌ Error: {e}")

def initialize_exp(instr):
    ###### Include trial run check so that the error plots are only displayed when it is not a trial run
    
    # Nscanpts = expCfg.N_scanPts
    param = params['scan'][params['scan']['names'][0]]['values'];    n_error=0
    sequenceArgs = params['seq']['args']
    instructionList = []
    
    # ------------------------------------------------------------------
        # est_time = expCfg.t_tot * expCfg.Nsamples * expCfg.N_scanPts ??
    if params['seq']['sequence'] not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:   # for sequences except ESR, set MW frequency
        seqArgList = [param[seq_no_plot[-1]]]
        seqArgList.extend(sequenceArgs)     # Make a seqArgList with dummy 1st element... just to create it.. pore change hoye jabe..

        [n_error, param] = seqctrl.param_err_check(instr, params['seq']['sequence'], seqArgList, param, len(param))
        if n_error>0:
            print(f"❌ \x1b[1;37;41mErr: Check Sequences...\x1b[0m")
            print(f"✔ \x1b[38;2;250;0;0m'{str(n_error)} parameters removed...\x1b[0m")

            # Close Error plots??
            # close_plots = dialog.yesno_box('Close Plots', 'Close the Error Plots?')
            # if close_plots == 'yes':
            #     plt.close('all')
            print("✔ Sequences checked... \x1b[38;2;100;250;0mErrors removed...\x1b[0m")
            # est_time = (2*expCfg.t_AOM*params['scan']['Nscanpts'] + sum(param))*expCfg.Nsamples
        else:
            print('✔ \x1b[38;2;100;250;0m----No Errors----\x1b[0m')
    else:       # for ESR exp, set MW freq to start pt of the scan -> there is no scannedParam
        seqArgList = sequenceArgs
    
    # define 'Nscanpts' as the length of the 'param' variable...
    params['scan'][params['scan']['names'][0]]['Nscanpts'] = len(param)
    params['scan'][params['scan']['names'][0]]['values'] = param

    if len(param)>0:
        print(f"▶ \x1b[38;2;250;100;10m{len(param)}\x1b[0m scan pts")
    else:
        print("❌ \x1b[38;2;200;200;10mCheck sequences for subtle problems...\x1b[0m")
    
    # Plotting the pulse sequences -------------------------------------------
    if plotPulseSequence:
        # seqArgList.append(params['pb']['channels'])
        # print(seqArgList)
        sequencecontrol.view_sequence(instr, params['seq']['sequence'], seqArgList+[params['pb']['channels']],
                                      False, param, seq_no_plot, seq_plot_dpi)
    
    # instructionList = pbctrl.PB_program(params['seq']['sequence'],seqArgList)[0]        # ekhane change chhilo.. (24/02/23)

    # Start the initial sequence now ------------
    print("▶ Starting Initial Sequence...")
    if trial_run[1] == 'n':
        instructionList = [start_initial_PB_seq()]
        print(instructionList)
        pb.run_sequence_for_diode(instructionList)
        print("▶\x1b[38;2;50;250;50m----------PB Running----------\x1b[0m")
        if params['seq']['sequence'] not in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq', 'aom_timing', 'T1ms0_train'] and trial_run[0] == 'n':
            sg.set_freq(params['mw']['freq'])
    else:
        # instructionList = []
        print("▶\x1b[38;2;250;250;0m----------PB NOT Running----------\x1b[0m")

    ## Create the data save folder...
    
    return [seqArgList, instructionList]
    
def start_initial_PB_seq():
    instructionList = []
    # the default initial sequence... LASER ta sob somoy ON thakbe ekhane...
    instructionList = [[concfg.laser, Inst.CONTINUE, 0, 500*ms],
                       [concfg.laser, Inst.BRANCH, 0, 500*ms]]
    # jodi onnyo kono initial sequence lage, eg some sequene involving a parameter, thle 'sequenceArgs' ke edit korte hbe...
    
    return instructionList

def read_details(sequence, n_channels, Nsamples_expCfg):
    # Enter no of [<apd>,<pd>] signals acquired in each cycle
    if sequence in ['modesr']:#, 'esr_seq']:
        reads_per_cyc = [4,1] if n_channels==2 else [4]
    elif sequence in ['drift_seq', 'T1ms0_train']:
        reads_per_cyc = [1,1] if n_channels==2 else [1]
    # elif 'dig_mod' in sequence:
    #     reads_per_cyc = [expCfg.N_laser*2*2]      # not acquiring from 2 channels (APD, PD) at the moment
    elif 'lia' in sequence:
        reads_per_cyc = [1,1] if n_channels==2 else [1]
    else:
        reads_per_cyc = [2,2] if n_channels==2 else [1]
    daq_Nsamples = sum(reads_per_cyc)*Nsamples_expCfg

    return [reads_per_cyc, daq_Nsamples]

# def extra_param_save_details(scan_time_list, exec_time):
#     extra_params_format = ' %s\t%0.2f\n %s\t%0.2f\n %s\t%0.4g\n %s\t%0.4g\n %s\t%g\n'  # %s\t%f\n
#     # scan_time_list in seconds; exec_time in seconds
#     param_list = ['Max_scan_time(us):',max(scan_time_list)*1e6, 'Min_scan_time(us):',min(scan_time_list)*1e6, 'Total_scan_time(s):',np.sum(scan_time_list), 'Total_run_time(s):',exec_time, 'Step:', (expCfg.scannedParam[1]-expCfg.scannedParam[0])]
#     # param_list[1] = i_scanpt+1..... Ekhane ki hbe?? Kon parameter save korbo??
#     # ekhane ekta if kore, jodi 'i'=1 hoe, thle first param_list ta return korbe, nhle porer param_list ta return korbe.
#     # Eta jodi kora hoe, thle runs>1 hole ba multi-param scan hole, notun param er sathe tar details save kora jabe...
#     # Config file e ekta variable lagbe jeta dekhabe je kon variable ta scan hoechhe, other than scannedParam
#     param_list = tuple(param_list)
#     return [extra_params_format, param_list]
    
def prepare_for_saving():
    
    # if 'f_number' in vars():       # Returns a dict of all local variables
    #     print("Previous file number: \x1b[38;2;250;150;0m"+f_number+'\x1b[0m')
    
    # expt code dir
    code_dir = Path.cwd()  # current working directory - generally 'NV_Experiment' folder. Data folder is created a level up
    data_dir = code_dir.parent / "Saved_Data" / time.strftime("%Y-%m-%d", time.localtime())
    if not (Path.is_dir(data_dir)):
        Path.mkdir(data_dir, parents=True, exist_ok=True)
    
    print(f"▶ Save folder: \x1b[38;2;100;250;30m{time.strftime('%Y-%m-%d', time.localtime())}\x1b[0m")
    
    global folder_number
    folder_number = input(f"▶ Enter folder number:: {params['seq']['sequence']}_#: ")
    
    # file saving directory
    while True:
        file_dir = data_dir / f"{params['seq']['sequence']}_{folder_number}"
        
        if not Path.is_dir(file_dir):
            Path.mkdir(file_dir)
            break
        elif not any(file_dir.iterdir()):  # Check if the directory is empty
            break
        else:
            folder_number = input(f"❌ Folder exists and is not empty. Re-enter folder number:: {params['seq']['sequence']}_#: ")
    
    if not (Path.is_dir(file_dir)):
        Path.mkdir(file_dir)
    
    datafilename = file_dir / f"{params['save']['savefileprefix']}_{folder_number}.h5"

    params['save']['savepath'] = file_dir
    params['save']['datafilename'] = datafilename
    params['save']['paramfilename'] = datafilename.with_suffix('.yaml')
        
    return folder_number

def savefile(filename: Path, data, data_type:str='data') -> bool:
    """
    Saves data to a file in a format determined by the file extension.
    Supported formats: .npy, .json, .yaml, .txt, .csv
    Args:
        filename (Path): The path to the file where data will be saved.
        data: The data to be saved.
    Returns:
        bool: True if the file was saved successfully, False otherwise.
    """

    def convert_numpy_to_python(obj):
        """Recursively convert numpy types to native Python types"""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.generic):
            return obj.item()
        elif isinstance(obj, dict):
            return {key: convert_numpy_to_python(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy_to_python(item) for item in obj]
        elif isinstance(obj, Path):
            return str(obj)
        else:
            return obj
    
    try:
        if 'npy' in filename.suffix:
            np.save(filename, data, allow_pickle=False)
            
        elif 'json' in filename.suffix:
            with open(filename, 'w') as f:
                json.dump(convert_numpy_to_python(data), f, indent=4)

        elif 'yaml' in filename.suffix:
            with open(filename, 'w') as f:
                yaml.dump(convert_numpy_to_python(data), f, indent=4, default_flow_style=False)
            
        elif 'txt' in filename.suffix:
            with open(filename, 'w') as f:
                f.write(data)

        elif 'csv' in filename.suffix:
            np.savetxt(filename, data, delimiter=',')

        elif 'h5' in filename.suffix or 'hdf5' in filename.suffix:
            with h5py.File(filename, 'a') as f:
                grp = f.require_group("data")
                for run in range(data.shape[0]):
                    grp.create_dataset(f"run_{run}", data=data[run])
        
        print(f"✔  Saved to\x1b[38;2;100;250;50m {filename.name}\x1b[0m !!!")
        return True
    except Exception as e:
        logging.exception(f"❌ Error saving file {filename.name}: {e}")
        return False


def calc_contrast(signal, reference, op):
    if op == '+-':
        contrast = signal - reference
    elif op == '-+':
        contrast = -signal + reference
    elif op == 's/r':
        contrast = signal/reference
    elif op == 'r/s':
        contrast = reference/signal
    else:
        contrast = (signal - reference)/(signal + reference)
    return contrast


# TODO: Check
def acquire_data(Nsamples, parameter, sequence, seqArgList, trial):
    the_list: list
    # setup next scan iteration (e.g. for ESR experiment, change microwave frequency; for T2 experiment, reprogram pulseblaster with new delay)
    if sequence in ['esr_dig_mod_seq', 'esr_seq', 'pesr_seq', 'modesr', 'drift_seq']:
        if trial[0] == 'n':
            sg.set_freq(parameter)
    else:
        seqArgList[0] = parameter
    # _, instructionList = pbctrl.PB_program(instr,sequence,seqArgList)[0]
    if sequence == 'esr_seq' and i_scanpt==0:
        # print(seqArgList)
        # #notun ... eta lagbe
        instructionList=[]
        name, the_list = PulseBlaster.PB_program(instr, sequence, seqArgList)
        # print(name)
        # print(the_list)
        for i in range(0, len(the_list)):
            instructionList.append(the_list[i][0])
        # print(instructionList)

        # start = time.perf_counter_ns()
        # print(i_scanpt+1,' / ',params['scan']['Nscanpts'])
        # if i_scanpt == 1:
        #     time.sleep(20)
        # stop = time.perf_counter_ns()
        # print('In acquire_data().. Starting sequence...')
        pb.run_sequence_for_diode(instructionList)
        # cts = []
    scan_start_time = time.perf_counter()   # time in seconds
    # cts = daqctrl.read_daq(ai_task, Nsamples,61*60)    #read DAQ
    # print(f'Nsamples = {Nsamples}')
    # print('Starting capture...')

    cts = ai_task.read_daq(Nsamples, timeout=params['seq']['t_total(s)']*params['seq']['Nsamples']+5)
    # cts = [np.mean(cts[int(2e6*1e-3)*i:int(2e6*1e-3)*(i+1)]) for i in range(params['daq']['ai']['daq_Nsamples'])]
    
    scan_end_time = time.perf_counter()
    
    scan_time = (scan_end_time - scan_start_time)     # in seconds
    return [cts, scan_time]

# TODO: Check
def acquire_data_all_params(Nsamples, parameters, sequence, seqArgList, trial):
    print("Acquiring..")
    instructionList=[]
    seqArgList[0] = parameters
    _, the_list = PulseBlaster.PB_program(instr, sequence, seqArgList)
    # print(the_list)
    for i in range(0, len(the_list)):
        instructionList.append(the_list[i][0])
    # print(instructionList)

    Nsamples = Nsamples if reload_pb else Nsamples*params['scan'][params['scan']['names'][0]]['Nscanpts']*params['scan']['Nruns']
    # print(f"Nsamples = {Nsamples}")
    pb.run_sequence_for_diode(instructionList)
    
    scan_start_time = time.perf_counter()   # time in seconds
    cts = ai_task.read_daq(Nsamples, timeout=60*10)
    scan_end_time = time.perf_counter()
    scan_time = (scan_end_time - scan_start_time)     # in seconds
    return [cts, scan_time]

# TODO: Check
def process_data(i_max, reads_per_cyc, data):
    # TODO: this function can be made simpler for multi-channel acquisition
    # Process data for plotting
    mean_sig = np.zeros((i_max, len(concfg.input_terminals)));
    mean_ref = np.zeros((i_max, len(concfg.input_terminals)));
    contrast = np.zeros((i_max, len(concfg.input_terminals)));
    
    signals_x = data[0:i_max,[i for i in range(0, expCfg.Nsamples*reads_per_cyc[0], reads_per_cyc[0])]].T
    mean_sig[:,0] = np.mean(signals_x,axis=0).T
    
    if params['seq']['sequence'].lower() not in ['t1ms0_train']:
        references_x = data[0:i_max,[i for i in range(1, expCfg.Nsamples*reads_per_cyc[0], reads_per_cyc[0])]].T
        mean_ref[:,0] = np.mean(references_x,axis=0).T
    
        contrast[:,0] = calc_contrast(mean_sig[:,0], mean_ref[:,0], 's/r')
    
    if len(concfg.input_terminals)>1:
        signals_y = data[0:i_max,[i for i in range(expCfg.Nsamples*reads_per_cyc[1], expCfg.Nsamples*reads_per_cyc[1]*2, reads_per_cyc[1])]].T
        references_y = data[0:i_max,[i for i in range(expCfg.Nsamples*reads_per_cyc[1]+1, expCfg.Nsamples*reads_per_cyc[1]*2, reads_per_cyc[1])]].T
        mean_sig[:,1] = np.mean(signals_y,axis=0).T
        mean_ref[:,1] = np.mean(references_y,axis=0).T
        
        contrast[:,1] = calc_contrast(mean_sig[:,1], mean_ref[:,1], 's/r')
    
    if params['seq']['sequence'].lower() not in ['t1ms0_train']:
        return [mean_sig, mean_ref, contrast]
    else:
        return [mean_sig]


def plot_raw_data(data):
    plt.figure(figsize=(10,5))
    plt.plot(np.ravel(data))
    plt.title('Raw Data')
    plt.xlabel('Samples')
    plt.ylabel('Voltage')
    plt.tight_layout()

def plot_data(i_max, param, processed_data, live=False, ax=None):
    mean_sig = processed_data[0]
    xValues = param[0:i_max]
    x_unit = params['plot']['plotXaxisUnits']
    x_label = params['plot']['plotXaxisLabel']

    # if live:
    #     plt.plot([x/x_unit for x in xValues], contrast[0:i_max], 'b--')
    #     plt.pause(0.0001)
    # else:
    if ax is None:
        fig, ax = plt.subplots(1,2, num=time.strftime(" [%H:%M:%S]", time.localtime()), figsize=(10,5))
        
    ax[0].plot([x/x_unit for x in xValues], mean_sig[0:i_max], '.-', label='Sig')
    if params['seq']['sequence'].lower() not in ['t1ms0_train']:
        mean_ref, contrast = processed_data[1:]
        ax[0].plot([x/x_unit for x in xValues], mean_ref[0:i_max], '.-', label='Ref')
        # plt.legend(['Sig X', 'Sig Y', 'Ref X', 'Ref Y']) if len(concfg.input_terminals)>1 else plt.legend(['Sig','Ref']) 
        # ax[0].legend()
        ax[0].set_xlabel(x_label)
        
        ax[1].plot([x/x_unit for x in xValues],  contrast[0:i_max], label=[i for i in range(len(concfg.input_terminals))])
        ax[1].set_xlabel(x_label);    ax[1].set_ylabel('Contrast')
        # ax[1].legend(['Contrast X', 'Contrast Y']) if len(concfg.input_terminals)>1 else plt.legend()
        ax[1].set_title('Plotting %d points' %(i_max))
    plt.tight_layout()
    # plt.savefig(args, kwargs)

# ----------------------------------------------------------------------------


# clk_cyc = 1e3/concfg.PBclk      # One clk cycle of PB = inverse of the clk freq of PB
if trial_run[0] == 'n':
    sg: SignalGenerator
else:
    sg: SignalGenerator_sim

print(f'🔄\x1b[38;2;250;250;0m Running {params['save']['savefileprefix']} sequence @ {time.strftime("%H:%M:%S", time.localtime())}\x1b[0m')
[sg, ao_task] = initialize_instr(params['seq']['sequence']) # type: ignore

seqctrl.check_params()
print("✔ Parameter checks completed...")
instr = 'diode'

# params['scan']['Nruns'] = Nruns #..........eta notun file e acche... kno??
[seqArgList, instructionList] = initialize_exp(instr)
# pb.stop_sequence()

# Experiment operations (with all hardwares working)... ki korbo eta???
# if trial_run == 'y':
#     print("Cannot proceed further")

[reads_per_cyc, daq_Nsamples] = read_details(params['seq']['sequence'], len(concfg.input_terminals), params['seq']['Nsamples'])

# Nsamples = int(2* 2e6* expCfg.t_AOM/1e9 * expCfg.Nsamples)
print(f"▶ DAQ-Nsamples = {daq_Nsamples}")

display_parameters = dialog.yesno_box(f"{params['save']['savefileprefix']} Parameters",
                                      (f"Channels\t: {str(concfg.input_terminals)} \n"
                                        f"Runs\t: {params['scan']['Nruns']}\n"
                                        f"ScanPts\t: {params['scan'][params['scan']['names'][0]]['Nscanpts']}\n"
                                        f"Samples\t: {params['seq']['Nsamples']}\n"
                                        f"Start\t: {min(params['scan'][params['scan']['names'][0]]['values'])/params['plot']['plotXaxisUnits']}\n"
                                        f"End\t: {max(params['scan'][params['scan']['names'][0]]['values'])/params['plot']['plotXaxisUnits']}\n"
                                        "Proceed ?")
                                    )
closed = False
# paramfilename = f"{params['save']['savepath']}params_{str(folder_number)}.txt"
# params['save']['paramfilename'] = paramfilename

params['pb']['reload_pb'] = reload_pb
params['pb']['load_pb_all_params'] = load_pb_all_params

params['daq']['ai']['daq_Nsamples'] = daq_Nsamples
params['daq']['reads_per_cyc'] = reads_per_cyc

Nruns = params['scan']['Nruns']
Nscanpts = params['scan'][params['scan']['names'][0]]['Nscanpts']
param = params['scan'][params['scan']['names'][0]]['values']

scan_time: list = []
run_start_time: float = 0.; save_flag: bool = False
i_scanpt: int = 0; i_run: int = 0
data_array: np.ndarray = np.array([])
if display_parameters == 'yes':    
    try:
        # if trial_run == 'n':
        ao_task.set_outputs_to_constant(output_field_in_gauss=params['test_field'])
        # save_flag = False
        pb.run_sequence_for_diode(instructionList)        # Run the PB hardware
        print("✔ Initial Sequence Started...")
        # continue_run = 'yes'
        closed = False
        # ai_task = daqctrl.config_ai(Nsamples)        # configure_daq() accepts no. of samples to read from DAQ-AI & returns the task created; variable 'readTask'
        
        #TODO Clear up the mess in AnalogInputTask initialization
        if reload_pb:
            # ai_task = AnalogInputTask(sampling_rate=concfg.daq_max_samp_rate, voltage_range=(-0.1, 0.1),
            #                       channels=concfg.input_terminals, trigger_source=concfg.start_trig_terminal, samps_per_chan=int(Nsamples))
            # without LIA
            # ai_task = AnalogInputTask(voltage_range=(-0.2, 0.2), channels=concfg.input_terminals, sampling_source=concfg.samp_clk_terminal,
            #                           start_trigger_source=concfg.start_trig_terminal, samps_per_chan=int(daq_Nsamples))
            
            # ai_task = AnalogInputTask(channels=concfg.input_terminals, voltage_range=(-0.2, 0.2), 
            #                           sampling_source='internal', sampling_rate=2e6, samps_per_chan=int(daq_Nsamples*2e6*1e-3),
            #                           start_trigger_source=concfg.start_trig_terminal,
            #                           pause_trigger_source=concfg.samp_clk_terminal,
            #                           )
            
            # with LIA one sample per trigger
            # ai_task = AnalogInputTask(voltage_range=(-10, 10), channels=concfg.input_terminals,
            #                           sampling_rate=10e3, sampling_source=concfg.samp_clk_terminal,
            #                           start_trigger_source=concfg.start_trig_terminal, samps_per_chan=int(daq_Nsamples))
            # with LIA multiple samples
            ai_task = AnalogInputTask(voltage_range=(-10, 10), channels=concfg.input_terminals,
                                      sampling_rate=10e3, sampling_source='',
                                      start_trigger_source=concfg.start_trig_terminal, samps_per_chan=int(daq_Nsamples))
        else:
            if load_pb_all_params:
                ai_task = AnalogInputTask(sampling_rate=concfg.daq_max_samp_rate, voltage_range=(-0.2, 0.2),
                                  channels=concfg.input_terminals, start_trigger_source=concfg.start_trig_terminal,
                                  samps_per_chan=int(daq_Nsamples*Nscanpts*Nruns))
        # print(f"Samples per channel = {ai_task._task.timing.samp_quant_samp_per_chan}")

        print(f"▶ DAQ configured for {params['seq']['Nsamples']} samples...")
        print(f'▶ {reads_per_cyc[0]} samples in each cycles...')
        data_array = np.zeros((params['scan']['Nruns'], params['scan'][params['scan']['names'][0]]['Nscanpts'],
                               params['daq']['ai']['daq_Nsamples']*len(concfg.input_terminals)))
        # data_list_runs = []

        # <<<<<<<<---------------------Run experiment--------------------->>>>>>>>
        print(f"🔄 Acquisition Started @ {time.strftime("%H:%M:%S", time.localtime())}")
        run_start_time = time.perf_counter()        # in seconds
        if reload_pb:
            for i_run in range(0, params['scan']['Nruns']):
                print("▶ \x1b[38;2;255;150;10mRun: ", i_run + 1, ' / ', params['scan']['Nruns'], "\x1b[0m")

                # data_array = np.zeros((params['scan']['Nscanpts'],Nsamples*len(concfg.input_terminals)))
                scan_time_list = []   # time (in seconds) for each scannedParam
                # print("",i_run+1,'/',params['scan']['Nruns'])
                # data_list = []
                for i_scanpt in range (0, params['scan'][params['scan']['names'][0]]['Nscanpts']):
                    # if i_scanpt>0:
                    #     break
                    if i_scanpt % 10 == 0:
                        print(i_scanpt+1,' / ',params['scan'][params['scan']['names'][0]]['Nscanpts'],': ',param[i_scanpt])
                        print(f"Waiting: {params['seq']['t_total(s)']*params['seq']['Nsamples']}s")
                    [cts, scan_time] = acquire_data(daq_Nsamples, param[i_scanpt], params['seq']['sequence'],
                                                    seqArgList+[params['pb']['channels']], trial_run)
                    # in general for multi-channel acquisition, cts is a list of lists... this needs to be taken care
                    # the sig1,ref1,sig2,ref2,... order has not been changed cts and data_array
                    data_array[i_run, i_scanpt,:] = np.ravel(np.array(cts))         # type: ignore
                    # data_list.append(np.array(cts))
                    scan_time_list.append(scan_time)
                # data_list_runs.append(data_list)
        else:
            scan_time_list = []
            if load_pb_all_params:
                [cts, scan_time] = acquire_data_all_params(daq_Nsamples, param, params['seq']['sequence'], seqArgList, trial_run)
                data_array = np.array(cts).reshape(Nruns, Nscanpts, daq_Nsamples)
                i_scanpt = params['scan'][params['scan']['names'][0]]['Nscanpts']
            scan_time_list.append(scan_time) # type: ignore

        run_end_time = time.perf_counter()          # in seconds
        exec_time = run_end_time - run_start_time       # in seconds
        if i_run+1 == params['scan']['Nruns']: # type: ignore # Close all if this is the last run
            closed = close_all(sg, ao_task, ai_task) # type: ignore
        # data_list_runs = np.array(data_list_runs)
        print(" Run time = %.2fs" % exec_time)        # in seconds
        print(" Scan time = %.2fs" % (np.sum(scan_time_list)))   # type: ignore # in seconds

        data_array += 0.0035
        # plot_raw_data(data=data_array)
        fig, axs = plt.subplots(1,2, num=f"{time.strftime(' [%H:%M:%S]', time.localtime())}", figsize=(10,5))
        for i_run in range(0, Nruns):
            processed_data = process_data(i_scanpt, reads_per_cyc, data_array[i_run])    # type: ignore # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
            plot_data(i_scanpt, param, processed_data, ax=axs) if 'train' not in expCfgFile else None # type: ignore
        
        # data_array = np.insert(data_array,0,param[0:i_scanpt+1],axis=2)     # type: ignore # insert the parameter value at the head of the array

        savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
        if savefile_yn == 'yes':
            if save_flag == False:          # Ask for file number iff no save was performed
                f_number = prepare_for_saving()
            save_flag = savefile(params['save']['datafilename'], data_array)
        else:
            print("▶ \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")
                # if i_run+1 < params['scan']['Nruns']:    
                #     continue_run = dialog.yesno_box('Continue',"Continue Run?") 
                #     if continue_run == 'yes':
                #         continue
                #     else:
                #         sys.exit("Run Interrupted...")

        # savefile_yn = dialog.yesno_box('Data Saving', "Save data to file?")
        # if savefile_yn == 'yes':
        #     # Data file saving
        #     if params['scan']['Nruns'] == 1:
        #         if save_flag == False:  # Ask for file number iff no save was performed
        #             [paramfilename, datafilename, f_number] = prepare_for_saving(savePath)
        #         save_flag = save_data(datafilename, data_array[0,:,:,:,:], f_number)
        #     else:
        #         for i_run in range(0, params['scan']['Nruns']):
        #             # for loop not required... save as numpy array
        #             datafilename = savePath + expCfg.savefileprefix + "_" + str(f_number)
        #             save_flag = save_data(datafilename, data_array)
        #     # np.save(datafilename+".npy", data_list_runs, allow_pickle=False)
        #     save_flag = True
        # else:
        #     print(" \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")
    

    except KeyboardInterrupt:
        run_end_time = time.perf_counter()
        exec_time = run_end_time - run_start_time       # in seconds
        print('🛑\x1b[38;2;250;100;0m User Interrupted. Quitting...\x1b[0m')
        # Plot the last run data upto the point where it was interrupted...
        processed_data = process_data(i_scanpt, reads_per_cyc, data_array[i_run])    # ekhane i_scanpt newa hochhe.. (i_scanpt+1) noi.. expt majhe stop korle last scan pt ta baad dewa hochhe...
        fig, axs = plt.subplots(1,2, num=f"{time.strftime(' [%H:%M:%S]', time.localtime())}", figsize=(10,5))
        plot_data(i_scanpt, param, processed_data, ax=axs)
        
        # Then ask whether to save it...
        savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
        if savefile_yn == 'yes':
            # Ask for filename only if there was no save operation
            if save_flag == False:
                f_number = prepare_for_saving()
            # A conditional save statement.. Different for dig_mod sequences...
            save_flag = savefile(params['save']['datafilename'], data_array)
        
        # save_data(len(concfg.input_terminals), datafilename, data_array, data_write_format, expCfg.Nsamples, reads_per_cyc)
    # else:        
    #     if (i_run+1)==params['scan']['Nruns']: # Save data of the final run here
    #         save_data(len(concfg.input_terminals), datafilename, data_array, data_write_format, expCfg.Nsamples, reads_per_cyc)
    finally:
        if not closed:
            closed = close_all(sg, ao_task, ai_task)
        
        print(f"▶ Read \x1b[38;2;250;150;50m{reads_per_cyc[0]}*{expCfg.Nsamples}\x1b[0m samples at each pt.")
        if save_flag:   # below line should "not" run if there was "no" save operation / data acquisition...
            print(f"▶ len(cts[0]) = \x1b[38;2;250;150;50m{len(cts[0])}\x1b[0m") if len(concfg.input_terminals) == 2\
                else print(f"▶ len(cts) = \x1b[38;2;250;150;50m{len(cts)}\x1b[0m")
        # Save parameters only if there was one save operation
        if save_flag:
            params['scan'][params['scan']['names'][0]]["Nscanpts"] = i_scanpt+1        # expParamList[1] -> value of N_scanPts
            params['scan']["Nruns"] = i_run+1           # expParamList[3] -> value of Nruns
            params['daq']['ai'].update({key: value for key, value in vars(ai_task).items()
                                        if not key.startswith('_')})
            
            seqctrl_name, the_list = PulseBlaster.PB_program(instr, params['seq']['sequence'],
                                                             seqArgList+[params['pb']['channels']])
            params['pb']['sequence_name'] = seqctrl_name
            params['pb']['sequence_list'] = [the_list[0][0], the_list[0][-1]]
            params['scan']['exec_time(s)'] = exec_time
            params['scan']['scan_time(ms)'] = [f"{num*1e3:0.2f}" for num in [np.mean(scan_time_list), np.std(scan_time_list)]]
            params['pb'].update(concfg.params)
            for name in params['scan']['names']:
                del params['scan'][name]['values']
            savefile(params['save']['paramfilename'], params)

        # scan_time in seconds; exec_time in seconds      
else:
    print("❌ Measurement aborted...")
    # if trial_run[0] == 'n' and params['seq']['sequence'] not in ['rodelay']:
    print(closed)
    if not closed:
        # closed = close_all(sg, ao_task, ai_task)
        closed = close_all(sg, ao_task)
    # sys.exit()
# return instructionList

# TODO: Put scan times in the parameter save file

# def process_data_for_mod(cts, reads_per_cyc, parameter, ith_scan_pt, mean_sig, mean_bg, contrast, data_array):
#     MW_off_data = []; MW_on_data = [];
#     for i in range(0,expCfg.Nsamples,1):        # expCfg.Nsamples = N_MW
#         MW_off_data += cts[4*i*expCfg.N_laser:(4*i+2)*expCfg.N_laser:1]   # laser mod data at 0 to 2*N_laser, this includes the reference data also...
#         MW_on_data += cts[(4*i+2)*expCfg.N_laser:(4*i+4)*expCfg.N_laser:1]      # MW mod data at 2*N_laser to 4*N_laser, including the reference data also...
    
#     sig_wo_mw = MW_off_data[0::2];  ref_wo_mw = MW_off_data[1::2];
#     sig_w_mw = MW_on_data[0::2];    ref_w_mw = MW_on_data[1::2];
    
#     data_array.append([parameter[ith_scan_pt], sig_wo_mw, ref_wo_mw, sig_w_mw, ref_w_mw])
    
#     # Process for plotting
#     mean_sig[ith_scan_pt] = np.mean(sig_w_mw)
#     mean_bg[ith_scan_pt] = np.mean(sig_wo_mw)
#     contrast.append(calc_contrast(mean_sig[ith_scan_pt], mean_bg[ith_scan_pt], 's/r'))
    
#     return [mean_sig, mean_bg, contrast, data_array]

# def save_data_txt_txt_for_mod(n_channels, datafilename, data_array, data_write_format, Nsamples, reads_per_cyc):
#     print(" Saving.....")
#     # try:
#     datafile = open(datafilename, 'a')
#     datafile.write("--%s--\n" % datafilename[-5:])
#     for lines in data_array:      # Each 'lines' in data_array stands for a specific scanpt...
#         for i in range(0,expCfg.N_laser*expCfg.Nsamples,1):
#             line = [lines[0], lines[1][i], lines[2][i]]
#             datafile.write(data_write_format % tuple(line))
#         for i in range(0, expCfg.N_laser*expCfg.Nsamples,1):
#             line = [lines[0], lines[3][i], lines[4][i]]
#             datafile.write(data_write_format % tuple(line))
#     datafile.close()
#     print(" Data saved to\x1b[38;2;100;250;50m %s_%s\x1b[0m !!!" % (expCfg.savefileprefix, f_number))
#     return True
# %%
