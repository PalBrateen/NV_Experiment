# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 16:50:32 2022

@author: PC
"""

import pyvisa as visa

rm = visa.ResourceManager()
serialaddr = '3'
SRSaddr = 'asrl'+str(serialaddr)+'::instr'
SRS = rm.open_resource(SRSaddr)
deviceID = SRS.query('*idn?')
print(deviceID)

SRS.write('freq10khz')
SRS.write('enbl1')
SRS.write('ampl2vpp')

SRS.write('modl1')
print('\x10 Modluation enabled...')
SRS.write('type1')
print('\x10 Freq Modulation enabled...')
SRS.write('mfnc1')
print('\x10 FM function set to sine...')
SRS.write('rate10khz')
print('\x10 Modulation rate set to 100 kHz...')
SRS.write('fdev0.5kHz')
print('x10 Freq deviation set to 10 kHz')