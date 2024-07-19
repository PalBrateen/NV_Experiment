# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:00:53 2024

@author: PC
"""
# import DAQcontrol as daqctrl
ao_task = daqctrl.config_ao("U9263")
daqctrl.start_ao(ao_task, daqctrl.coil_calibration([50,0,0]))
# sgctrl.set_sg_freq(freq2.)

#%%
daqctrl.start_ao(ao_task,[0,0,0])
daqctrl.close_daq_task(ao_task)

#%%
pbctrl.pb_stop()
pbctrl.pb_close()
