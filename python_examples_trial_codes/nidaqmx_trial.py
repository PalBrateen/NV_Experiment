# -*- coding: utf-8 -*-
#%% Imports
import nidaqmx
from  nidaqmx.constants import VoltageUnits, TerminalConfiguration, AcquisitionType, Edge, DigitalWidthUnits, MIOAIConvertTimebaseSource
import sequenceControl as seqctrl
import PBcontrol as PBctrl
from connectionConfig import laser, samp_clk, conv_clk, start_trig
from matplotlib import pyplot as plt
import time, os
#%% Create a task and assign channels
channels = ["Dev1/ai7", "Dev1/ai0"]#, "Dev1/ai20", "Dev1/ai16"]
if 'j' not in vars():
    j=1
else:
    j+=1

# folder = 
filename = os.getcwd()+"\\Saved_Data\\"+time.strftime("%Y-%m-%d", time.localtime()) + '\\'+'sampling_check'+str(j)+'.txt'
task = nidaqmx.Task()
for channel in channels:
    task.ai_channels.add_ai_voltage_chan(channel,"", TerminalConfiguration.DIFFERENTIAL, -0.1, 0.1)
    
# #%% Dislay the channels
# print(task.channels)
# print("Max sampling rate for %d task(s) = %g MSa/s" % (len(channels),
#            (task.timing.samp_clk_max_rate/1e6)))
# print("Max ADC conv rate for %d task(s) = %g /s" % (len(channels),(task.timing.ai_conv_max_rate)))

#%% Configure the sample clock and convert clock
src_samp_clk = 'PFI14'
src_conv_clk = 'PFI15'
Nsamples = 10
# samp_rate_per_chan_per_s = task.timing.samp_clk_max_rate
task.timing.cfg_samp_clk_timing(0.5e6, src_samp_clk, Edge.RISING, AcquisitionType.FINITE, Nsamples)
# task.timing.samp_clk_src = src_samp_clk
# task.timing.samp_clk_active_edge = Edge.RISING

# print('Samp Rate for the task: %g MSa/s' %(samp_rate_per_chan_per_s/1e6))

task.timing.ai_conv_src = src_conv_clk
# conv_clk_rate = task.timing.samp_clk_max_rate
# task.timing.ai_conv_rate = 1e6
task.timing.ai_conv_active_edge = Edge.RISING
# print('ADC Conv Rate for the task: %g MSa/s' %(conv_clk_rate/1e6))

#%% Configure pulse sequence
sequence = 'simult_samp'
instructionList = PBctrl.PB_program(sequence, [])[0]
PBchannels = {'Conv\nCLK':conv_clk,'Samp\nCLK':samp_clk,'Signal':laser}
# [t_us,channelPulses,yTicks] = seqctrl.plot_sequence(instructionList, PBchannels)
# for channel in channelPulses:
    # plt.plot(t_us, list(channel))
# plt.yticks(yTicks, PBchannels.keys())          # Include the names of the PB channels
# plt.xlabel('Time (us)')
# plt.ylabel('Channel')
#%% Acquire samples
PBctrl.run_sequence(instructionList)
cts = task.read(Nsamples,60)
task.close()
# print('1st channel:\n'+str(cts[0]))
# print('2nd channel:\n'+str(cts[1]))

# File saving
data_file = open(filename,'a')
data_file.write("%s\t%s\n" % tuple(channels))
for i in range(0, Nsamples):
    line = [cts[0][i], cts[1][i]]
    data_file.write("%0.3f\t%0.3f\n" % tuple(line))
PBctrl.pb_init(); PBctrl.pb_stop(); PBctrl.pb_close()
data_file.close()
#%% Close task
task.close()

#%% Table for sampling rate vs # channels

print("# tasks \t Max Samp Rate(MSa/s) \t Max ADC Conv Rate(/s)")
chan_num = [i for i in range(0,8)]
chan_num.extend([i for i in range(16,24)])
for j in range(0,16):
    channels = []
    for i in chan_num[0:j+1]:
        channels.append("Dev1/ai"+str(i))
    task = nidaqmx.Task()
    for channel in channels:
        task.ai_channels.add_ai_voltage_chan(channel,"", TerminalConfiguration.DIFFERENTIAL, -0.1,0.1, VoltageUnits.VOLTS)
    print("%d\t%g\t%g" %(j+1, (task.timing.samp_clk_max_rate/1e6),(task.timing.ai_conv_max_rate)))
    task.close()
