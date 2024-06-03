# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 13:20:48 2024

@author: PC
"""

import time, PBcontrol_v2 as PBctrl, matplotlib.pyplot as plt
time1=[]
time2=[]
for i in range(0,100000):
    t1 = time.perf_counter()
    PBctrl.pb_start()
    t2 = time.perf_counter()
    PBctrl.pb_stop()
    t3 = time.perf_counter()
    time1.append(t2-t1)
    time2.append(t3-t2)

#%%
plt.figure(1)
plt.plot(time1)
plt.plot(time2)
plt.figure(2)
plt.hist(time1)
plt.hist(time2)