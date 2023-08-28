# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 19:25:39 2022

@author: PC
"""

import time

def func(a,b):
    print(a+b)
a = 5; b = 10;
start = time.perf_counter()
print(a+b)
step1 = time.perf_counter()
func(a,b)
step2 = time.perf_counter()

print(step1-start)
print(step2-step1)