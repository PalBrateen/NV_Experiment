# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 18:58:24 2021

@author: brate
"""
from matplotlib import pyplot as plt
import numpy as np
from random import shuffle

plt.figure()

for i in range(2,10):
    x = list(range(0,i))
    shuffle(x)
    y = list(range(0,i))
    shuffle(y)
    plt.plot(x,y,'b.')
    plt.pause(1)
    print(x,y)