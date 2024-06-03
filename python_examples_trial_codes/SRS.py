# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 12:21:58 2022

@author: PC
"""

import pyvisa as visa
import time
rm = visa.ResourceManager()
SG = rm.open_resource('asrl3::instr')
#%%
SG.write('enbr0')
SG.write('freq2.87ghz')
# time.sleep(5)
for i in range(0,10):
    SG.write('enbr1')
    time.sleep(2)
    SG.write('enbr0')
    time.sleep(2)