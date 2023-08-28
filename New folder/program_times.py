# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 12:06:22 2021

@author: PC
"""
#%%
import time
a = input('a = ')
start = time.perf_counter_ns()
if int(a)>10:
    step = time.perf_counter_ns()
    a=10
    stop = time.perf_counter_ns()
else:
    step = time.perf_counter_ns()
    print('a<10')
    stop = time.perf_counter_ns()

print("%f\t%f\t%f" % (step-start, stop-start, stop-step))

#%%
def func(a,b):
    print(a+b)
a = 5; b = 10;
start = time.perf_counter()
start = time.perf_counter()
func(a,b)
step1 = time.perf_counter()
print(a^b)
step2 = time.perf_counter()
stop = time.perf_counter()

print((step1-start)*1e6)
print((step2-step1)*1e6)
print((stop-step2)*1e6)

#%%
import connectionConfig as conCfg
import time

start = time.perf_counter()
print(len(conCfg.input_terminals))
step1 = time.perf_counter()
n_channels = len(conCfg.input_terminals)
print(n_channels)
step2 = time.perf_counter()

print((step1-start)*1e6)
print((step2-step1)*1e6)