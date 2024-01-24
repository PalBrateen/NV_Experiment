#%% Initialize instruments
from os.path import isdir, isfile; from os import makedirs; from importlib import import_module
import sys, time, dialog, os
import pyvisa as visa, numpy as np, matplotlib.pyplot as plt
kHz=1e3; MHz=1e6; GHz=1e9

rm = visa.ResourceManager()
rm.list_resources()

SG = rm.open_resource('ASRL5::INSTR')
SA = rm.open_resource('USB0::0x0699::0x0408::C049247::INSTR')
try:
    if (not SA.query("*IDN?") == '') and (not SG.query("*IDN?") == ''):
        print("\x1b[38;2;250;250;0mINIT Success!\x1b[0m")
except visa.errors.VisaIOError:
    print('\x1b[38;2;250;10;10mINIT Error!')

#%% User inputs
# specan
n_averages = 4
scan_start_freq = 200 *kHz
scan_stop_freq = 2 *MHz

# ODMR
odmr_start_freq = int( 2.5 *GHz)
odmr_stop_freq = int( 3.3 *GHz)
odmr_step_freq = int( 2 *MHz)
mw_ampl = 0

#%% Set the SG and SA with the parameters for scan
odmr_n_freq = round((odmr_stop_freq - odmr_start_freq)/odmr_step_freq + 1)

odmr_freq = np.linspace(odmr_start_freq,odmr_stop_freq,odmr_n_freq,endpoint=True)
# odmr_n_freq = len(odmr_freq)
# print(odmr_n_freq)

SG.write('ampr '+str(mw_ampl)+'dbm')
SG.write('enbr1')

SA.write("RF:DETECT:MOD AUTO")
SA.write("RF:RBW 200")
SA.write("SEL:RF_AVE 1")
SA.write("RF:STAR "+str(scan_start_freq))   
SA.write("RF:STOP "+str(scan_stop_freq))    
SA.write("RF:RF_AVE:NUMAV "+str(n_averages))         # set the number of averaging runs
# SA.write("SENS1:SWE:TIME?")                         # query the sweep time and set the timeout time accordingly
# sweep_time = eval(SA.read())
# SA.timeout = (sweep_time*n_averages+10)*1e3       # extra 10 seconds

# specan parameters for saving
SA.write(":DAT:SOU RF_NORM")
SA.write(":DAT:STAR 1")
SA.write(":DAT:STOP 2e7")
SA.write(":WFMO:ENC ASC")
SA.write(":WFMO:BYT_N 4")
SA.write(":HEAD 1")
SA.write(":VERB 1")

#%% Function definitions
def trace_array(trace_string):
    trace = []      # individual freq er trace; a list
    pos = trace_string.find(',',0)       # find syntax: find(char, start, end)
    wordlength = pos+1
    trace.append(eval(trace_string[0:pos]))
    while trace_string.find(',',pos+1)>0:
        trace.append(eval(trace_string[pos+1:(trace_string.find(',',pos+1))]))
        pos = trace_string.find(',',pos+1)
    trace.append(eval(trace_string[pos+1:(trace_string.find(',',pos+1))]))
    pos = trace_string.find(',',pos+1)
    return trace

def read_val():
    value_string = SA.read()
    value_string_wo_info = value_string[value_string.find(" ")+1:]
    value = eval(value_string_wo_info)
    return [value,value_string_wo_info]

def acquire_specan_trace():
    error = True
    SA.write("RF:SPECTRUMT RESET")
    
    SA.write("RF:RF_AVE:COUN?"); avg_runs = read_val()[0]
    i=0;
    while n_averages==avg_runs:
        SA.write("RF:RF_AVE:COUN?"); avg_runs = read_val()[0]
        i+=1
        continue
    # print(i)
    
    start_time = time.perf_counter()
    # time.sleep(2)         eta toh ar lagbe na.. checking toh aagei hoye jachhe
    # SA.write("INIT1:CONT OFF")
    # SA.write("INIT1; *WAI")      # begin trace measurement and wait for completion
    SA.write("RF:RF_AVE:COUN?"); avg_runs = read_val()[0]
    i=0
    while n_averages>avg_runs:
        # time.sleep(0.100)
        SA.write("RF:RF_AVE:COUN?"); avg_runs = read_val()[0]
        # print(avg_runs)
        i+=1
        continue
    stop_time = time.perf_counter()
    # print(i)
    i=0
    while True:
        try:
            i+=1
            SA.write(":CURV?");     trace = read_val()
            # error = False
            break
        except visa.errors.VisaIOError:
            print('Timeout occured!')
            # error = True
            # continue
    # print(i)
    
    trace_array = np.array(trace[0])     # trace[0] = tuple of the float values
    trace_string = trace[1]         # trace[1] = trace string wo info
    
    # trace unit conversion: values are in watt
    # trace_array = 10*np.log10(1e3*trace_array)        # convert values to dBm
    return [trace_array, trace_string, stop_time-start_time]

def prepare_for_saving(savePath):
    global file_number
    if 'file_number' in vars():       # Returns a dict of all local variables
        print("Previous file number: \x1b[38;2;250;150;0m"+file_number+'\x1b[0m')
    file_number = input("File name: specan_ESR"+"_#: ")
    
    # paramfilename = savePath + expCfg.saveFileName + "_" + "params_" + file_number +".txt"
    datafilename = savePath + "specan_ESR_" + file_number +".txt"
    if (isfile(datafilename)):
        print('\x1b[38;2;250;250;0mAppending file: '+'specan_ESR_'+file_number+'.txt\x1b[0m')
        
    return [datafilename, file_number]

def save_data_txt(datafilename, data_string):
    print("\x10 Saving.....")
    with open(datafilename, 'a') as datafile:
        datafile.write("%s\n" % ("specan_ESR_" + file_number))
        datafile.write(data_string)        
    
    print("\x10 Data saved to\x1b[38;2;100;250;50m %s_%s\x1b[0m !!!" % ("specan_ESR", file_number))
    return True

#%% measurement run
# number of averages thik set hoyechhe?
SA.write("RF:RF_AVE:NUMAV?"); val = read_val()[0]
while not(val==n_averages):
    SA.write("RF:RF_AVE:NUMAV "+str(n_averages))
    SA.write("RF:RF_AVE:NUMAV?"); val = read_val()[0]

# odmr_n_freq = 6
set_freq_simulation = np.linspace(int(200e3),int(2e6),odmr_n_freq,endpoint=True)
# SA.write("SENS1:SWE:POIN?")
sweep_points = 751 #-------------------------------
traces = np.zeros([odmr_n_freq,2*sweep_points])     # sob freq er trace store korar jonno
contrast_traces = np.zeros([odmr_n_freq,sweep_points])
contrast_traces1 = np.zeros([odmr_n_freq,sweep_points])
acq_time = np.zeros([2,odmr_n_freq])
save_data_string=""
save_flag = False

# data acquisition starts
start_time = time.perf_counter()
for i in range(0,odmr_n_freq):
    print(i+1,' / ',odmr_n_freq,': ',odmr_freq[i])
    SG.write("FREQ "+str(odmr_freq[i]))
    SG.write("enbr1")           # MW ON
    
    # SG.write("FREQ "+str(set_freq_simulation[i]))
    # SG.write("ampl-30dbm")
    # SG.write("enbl1")           # MW ON
    
    # time.sleep(0.1)
    [sig_trace, trace_str, sig_time] = acquire_specan_trace()
    acq_time[0,i] = sig_time
    traces[i,0:sweep_points] = sig_trace
    save_data_string += str(odmr_freq[i]) + "," + trace_str

    SG.write("enbr0")           # MW OFF
    
    # SG.write("enbl0")           # MW OFF
    
    # time.sleep(0.1)
    [ref_trace, trace_str, ref_time] = acquire_specan_trace()
    acq_time[1,i] = ref_time
    traces[i,sweep_points:2*sweep_points] = ref_trace
    save_data_string = save_data_string[0:-1] + "," + trace_str   # There is a \n at the end of the SA-returned string; make use of that

    contrast_traces[i] = sig_trace/ref_trace
    # contrast_traces[i] = sig_trace
    # contrast_traces1[i] = ref_trace

end_time = time.perf_counter()
print("Runtime = %.2fs" % (end_time-start_time))

# release SG
SG.write("FREQ 2.87GHz")
SG.write('lcal')

#%% Plotting
plt.figure(); plt.plot(np.transpose(acq_time))

SA.write('RF:STAR?'); specan_start_freq = read_val()[0]
SA.write('RF:STOP?'); specan_stop_freq = read_val()[0]
specan_step_freq = (specan_stop_freq - specan_start_freq)/(sweep_points-1)
specan_freq = np.linspace(specan_start_freq,specan_stop_freq,sweep_points,endpoint=True)

# # Plot the spectrum
# plt.figure()
# for i in range(0,odmr_n_freq):
#     plt.plot(specan_freq, np.transpose(contrast_traces[i]))
# plt.axis([min(specan_freq), max(specan_freq),-120,-20])

# plt.figure()
# for i in range(0,odmr_n_freq):
#     plt.plot(specan_freq, np.transpose(contrast_traces1[i]))

# plt.axis([min(specan_freq), max(specan_freq),-120,-20])
# # plt.legend('')

# Plot the spectrogram
# generate 2 2d grids for the x & y bounds
spectrogramy, spectrogramx = np.meshgrid(odmr_freq, specan_freq, indexing='ij')
z = contrast_traces

fig, ax = plt.subplots(dpi=200)

c = ax.pcolor(spectrogramx, spectrogramy, z, vmin=np.min(z), vmax=np.max(z))
# ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax.axis([spectrogramx.min(), spectrogramx.max(), spectrogramy.min(), spectrogramy.max()])
fig.colorbar(c, ax=ax)
plt.xlabel("Fourier Frequency [Hz]")
plt.ylabel("MW Frequency [GHz]")

#%% file saving - saving all the traces---------------------
cwd = os.getcwd()
savePath = cwd[0:cwd.rfind('\\')]+ "\\Saved_Data\\" + time.strftime("%Y-%m-%d", time.localtime()) + '\\'         # Path to folder where data will be saved
if not (isdir(savePath)): 
    os.makedirs(savePath)
print("\x10 Save folder: \x1b[38;2;100;250;30m" + time.strftime("%Y-%m-%d", time.localtime()) + '\x1b[0m')

savefile_yn = dialog.yesno_box('Data Saving',"Save data to file?")
if savefile_yn == 'yes':
    # Data file saving
    if save_flag == False:          # Ask for file number iff no save was performed
        [datafilename, file_number] = prepare_for_saving(savePath)
    save_flag = save_data_txt(datafilename, save_data_string)
else:
    print("\x10 \x1b[38;2;250;50;10mData NOT saved !!!\x1b[0m")


# %%
