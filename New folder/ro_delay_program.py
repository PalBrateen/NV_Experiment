# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 12:34:51 2022

@author: PC
"""

intervals = [0.8, 1, 1.2, 2, 5, 10, 50, 100, 500, 1000, 2000, 5000]

import rodelay_config as expCfg
import mainControl as main

for interval in intervals:
    # print(interval)
    expCfg.interval = interval
    