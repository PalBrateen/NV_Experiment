# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 13:23:21 2024

@author: PC
"""

import pyvisa as visa

def init_ametek():
    rm = visa.ResourceManager()
    # rm.list_resources()
    ametek = rm.open_resource("USB0::0x0A2D::0x001B::10236409::RAW")
    return ametek

def stop_oscillator(ametek):
    ametek.write("OA 0")
    
# stop_oscillator(ametek)