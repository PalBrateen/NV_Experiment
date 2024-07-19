# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:06:44 2024

@author: PC
"""

import DAQcontrol as daqctrl, time
ao_task = daqctrl.config_ao("U9263")

#%%

daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,50]))
time.sleep(1)
daqctrl.start_ao(ao_task, daqctrl.coil_calibration([50,0,0]))
time.sleep(1)
daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,0]))

#%%
daqctrl.close_daq_task(ao_task)

#%%
for i in range(0,5):
    daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,-50]))
    time.sleep(1)
    daqctrl.start_ao(ao_task, daqctrl.coil_calibration([+50,0,0]))
    time.sleep(1)
    daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,0]))
    time.sleep(5)
    
#%%
daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,50]))
time.sleep(1)
daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,0]))

#%%
daqctrl.start_ao(ao_task, daqctrl.coil_calibration([50,0,0]))
time.sleep(1)
# daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,0]))

#%%
for i in range(0,5):
    daqctrl.start_ao(ao_task, daqctrl.coil_calibration([50,0,0]))
    time.sleep(1)
    daqctrl.start_ao(ao_task, daqctrl.coil_calibration([0,0,0]))
    time.sleep(5)