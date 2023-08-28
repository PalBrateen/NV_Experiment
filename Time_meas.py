# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 12:41:52 2022
@author: PC
"""

import connectionConfig as concfg
import numpy as np
import DAQcontrol as daq
import time, os.path
from os import makedirs
import nidaqmx, dialog
from matplotlib import pyplot as plt
plt.close('all')

# Change the source of Sample clock and conv clock to default ONBOARD CLOCK
concfg.samp_clk_terminal = ''
concfg.conv_clk_terminal = ''
concfg.start_trig_terminal = ''
concfg.input_terminals = ['P6363/ai15']

rate = 2e6      #concfg.cal_samp_rate() (enter value in Hz)
print('Sampling Rate:', rate)
# Nsamples = int( 2e6 )
# est_time = Nsamples/rate *1e3    # in ms
est_time = 10 *1e3   # in ms (enter value in seconds)
Nsamples = int( est_time/1e3 * rate )

# cts=[]
read_task = daq.configure_daq(Nsamples, rate)
# while True:
print("Start...")
read_task.start()
start = time.perf_counter_ns()
cts = daq.read_daq(read_task, Nsamples, est_time/1e3+10)    # in volts
# cts = read_task.read(Nsamples, 10)
stop = time.perf_counter_ns()
read_task.stop()

# read_task.timing.samp_clk_src?
act_time = (stop - start) /1e6  # in ms
# plt.figure()
# plt.plot(cts)
print("Mean voltage = %f\nStnd Dev = %f\n" % (np.mean(cts), np.std(cts)))

daq.close_daq_task(read_task)

saving = dialog.yesno_box('File Saving', 'Save File?')
if saving == 'yes':
    file_number = input('File name: Time_#: ')
    filename = "D:\\Brateen\\Python_codes\\qdSpectro-active\\Saved_Data\\2023-08-18\\Time_"+file_number+".txt"
    datafile = open(filename, 'a')
    datafile.write("%g\t%g\n%g\t%g\n" %(rate, Nsamples, est_time, act_time))
    data_write_format = "%0.3f\n"
    for i in range(len(cts)):
        datafile.write(data_write_format %(cts[i]))
    datafile.close()
