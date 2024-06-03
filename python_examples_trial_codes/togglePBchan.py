#togglePBchan.py

import sequenceControl as seqCtl
from connectionConfig import *
from PBcontrol import *
import numpy as np

PBchanDict = {'A': AOM, 'M': uW, 'D': DAQ, 'I': I, 'Q':Q, 'S':STARTtrig}
deviceNameDict = {'A': 'AOM', 'M': 'Microwaves', 'D': 'DAQ gate', 'I': 'I', 'Q':'Q','S':'Start trigger'}

#Function definitions ------------------------------------------
def PBtoggle(device,currentState):
    configurePB()
    PBchanAddressOfToggledDevice = PBchanDict[device]
    newState = PBchanAddressOfToggledDevice^currentState
    pb_start_programming(PULSE_PROGRAM);
    start = pb_inst_pbonly(newState, Inst.CONTINUE, 0, 200.0*ms);
    pb_inst_pbonly(newState, Inst.BRANCH, start, 100.0*ms);
    pb_stop_programming();
    pb_start()
    pb_close()
    if PBchanDict[device]&newState:
        print(deviceNameDict[device], 'ON')
    else:
        print(deviceNameDict[device], 'OFF')
    return newState
    
def printStatus(state):
    devicesThatAreON =[]
    for device in deviceNameDict.keys():
        if PBchanDict[device]&state:
            devicesThatAreON.append(deviceNameDict[device])
    if devicesThatAreON ==[]:
        print('\n No devices are on.\n')
    else:
        print('\n Devices currently on:')
        print(*devicesThatAreON, sep = ", ")
        
# Main loop ---------------------------------------------------
#Print instructions
print('Press A to toggle AOM')
print('Press M to toggle microwaves')
print('Press D to toggle DAQ gate/sample clock')
print('Press S to toggle DAQ start trigger')
print('Press I to toggle I phase control')
print('Press Q to toggle Q phase control \n')

#Wait for user command:
state= 0 # bit flag of devices which are currently on
while True:
        print('\n To quit, press E. For a list on devices that are on, press W')
        device = input("Please enter key to toggle a device...\n")
        if device not in ['A','M','D', 'I', 'Q', 'S', 'E', 'W']:
            print('Device not recognised. Please choose from A, M, D, I, Q or S or press E to quit.')
            continue
        if device == 'E':
            break
        if device == 'W':
            printStatus(state)
        else:
            state = PBtoggle(device,state)