#%% Initialize instruments
from os.path import isdir, isfile; from os import makedirs; from importlib import import_module
import sys, time, dialog, os
import pyvisa as visa, numpy as np, matplotlib.pyplot as plt
import time
kHz=1e3; MHz=1e6; GHz=1e9

rm = visa.ResourceManager()
rm.list_resources()

SG = rm.open_resource('ASRL5::INSTR')
SA = rm.open_resource('TCPIP0::169.254.179.26::inst0::INSTR')
if (not SA.query("*IDN?") == '') and (not SG.query ("*IDN?") == ''):
    print("Success!")

#%% User inputs
# specan
n_averages = 2
scan_start_freq = 200 *kHz
scan_stop_freq = 2 *MHz

# ODMR
odmr_start_freq = int( 2.8 *GHz)
odmr_stop_freq = int( 2.85 *GHz)
odmr_step_freq = int( 10 *MHz)
mw_ampl = 0

#%% Set the SG and SA with the parameters for scan
odmr_n_freq = round((odmr_stop_freq - odmr_start_freq)/odmr_step_freq + 1)

freq = np.linspace(odmr_start_freq,odmr_stop_freq,odmr_n_freq,endpoint=True)
# odmr_n_freq = len(freq)
print(odmr_n_freq)

SG.write('ampr '+str(mw_ampl)+'dbm')
SG.write('enbr1')
SA.write("SENS1:FREQ:STAR "+str(scan_start_freq))   
SA.write("SENS1:FREQ:STOP "+str(scan_stop_freq))    
SA.write("SENS1:SWE:COUN "+str(n_averages))         # set the number of averaging runs
SA.write("SENS1:SWE:TIME?")                         # query the sweep time and set the timeout time accordingly
sweep_time = eval(SA.read())
SA.timeout = (sweep_time*n_averages+10)*1e3       # extra 10 seconds

#%% measurement run
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

def acquire_specan_trace(set_mw_freq):
    SG.write("FREQ "+str(set_mw_freq))

    SA.write("INIT1:CONT OFF")
    SA.write("INIT1; *WAI")      # begin trace measurement and wait for completion
    SA.write('TRAC? TRACE1')    # acquire the trace, TRACE1 here
    trace_string = SA.read()               # the trace values are returned as a string

    return trace_string

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

# odmr_n_freq = 6

SA.write("SENS1:SWE:POIN?")
sweep_points = eval(SA.read())
traces = np.zeros([odmr_n_freq,2*sweep_points])     # sob freq er trace store korar jonno
contrast_traces = np.zeros([odmr_n_freq,sweep_points])
save_data_string=""
save_flag = False

# data acquisition starts
for i in range(0,odmr_n_freq):
    SG.write("enbr1")           # MW ON
    trace_str = acquire_specan_trace(freq[i])
    sig_trace = np.array(trace_array(trace_str))
    traces[i,0:sweep_points] = sig_trace
    save_data_string += str(freq[i]) + "," + trace_str

    SG.write("enbr0")           # MW OFF
    trace_str = acquire_specan_trace(freq[i])
    ref_trace = np.array(trace_array(trace_str))
    traces[i,sweep_points:2*sweep_points] = trace_array(trace_str)
    save_data_string = save_data_string[0:-1] + "," + trace_str   # There is a \n at the end of the SA-returned string; make use of that

    # contrast_traces[i] = sig_trace/ref_trace
    contrast_traces[i] = sig_trace
SG.write("FREQ 2.87GHz")
SG.write('lcal')

# file saving - saving all the traces------------------------------------------------
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
