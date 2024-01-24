# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:26:03 2021

@author: PC
"""
#%%
from spinapi import *
import connectionConfig as concfg

pb_init()
pb_get_version()
pb_core_clock(500)
ONE_PERIOD = 0x200000
# shortpulseFLAG = ONE_PERIOD

pb_start_programming(PULSE_PROGRAM)
on_time  = 500 *ms
off_time = 500 *ms
# start = pb_inst_pbonly(shortpulseFLAG^11, Inst.CONTINUE, 0, on_time)
start = pb_inst_pbonly(concfg.laser, Inst.CONTINUE,0,on_time)
s = pb_inst_pbonly(concfg.laser, Inst.BRANCH, start, off_time)
pb_stop_programming()

pb_reset()
pb_start()

# print("Enter any key to stop")
# input("Enter any key to stop: ") 

# pb_stop()
pb_close()

#%% Reset the PB
pb_init(); pb_stop(); pb_close()