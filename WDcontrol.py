# WDcontrol

import pyvisa as visa
import sys
# Frequency unit multiplier definitions
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

# unicode = lambda s: str(s)

##-------------------- Function definitions--------------------
def initSG(serialaddr='3'):
    SGaddr = 'asrl'+str(serialaddr)+'::instr'      # Construct instrument identifier from serial address
    rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager
    SG = rm.open_resource(SGaddr)         # Start 
    
    # try:
    #     deviceID = SG.query('*IDN?')       # Try quering SG identity:
    # except Exception as excpt:
    #     print('Error: could not query SG...') # ' Please check Serial address is correct and SG Serial communication is enabled. Exception details:', type(excpt).__name__,'.',excpt)
    #     sys.exit()

    return SG

def SG_err_check(SG):
    # err = SG.query('LERR?')
    # if int(err) != 0:
    #     print('SG error: error code', int(err),'. Please refer to SG manual for a description of error codes.')
    #     sys.exit()
    return None
            
def enable_SG_op(SG):
    SG.write('o1')      # o) set RF On(1) or Off(0) 1
    # SG_err_check(SG)

def disable_SG_op(SG):
    SG.write('o0')      # o) set RF On(1) or Off(0) 1
    # print("SG Output disabled...")
    
def set_SG_amp(SG,RFamplitude=3):
    SG.write('h1')      # h) set RF High(1) or Low(0) Power 1
    SG.write('a'+str(RFamplitude))      # a) set RF Power (0=mimimum, 3=maximum) 2
    
def set_SG_freq(SG, freq):    
    SG.write('f'+str(freq))     # f) RF Frequency Now (MHz) 1000.000

def set_SG_disp(SG, disp):
    return None

def setup_SG_mod(SG,sequence):
    #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
    #and disables modulation for ESR, Rabi and T1 sequences.
    # if sequence in ['ESGeq', 'RabiSeq', 'T1seq']:
    #     disable_mod(SG)
    # elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
    #     enable_iq_mod(SG)
    # else:
    #     print('Error in SGcontrol.py: unrecognised sequence name passed to setupSGmodulation.')
    #     sys.exit()
    return None
    
def setup_SG_pulse_mod(SG):
    # SG.write('modl1')
    # SG_err_check(SG)
    # SG.write('type4')
    # SG_err_check(SG)
    # SG.write('pfnc5')
    # SG_err_check(SG)
    return None
    

def enable_iq_mod(SG):
    # SG_err_check(SG)
    # SG.write('MODL 1')     #Enable modulation
    # SG_err_check(SG)
    # SG.write('TYPE 6')     #Set modulation type to IQ
    # SG_err_check(SG)
    # SG.write('QFNC 5')     #Set IQ modulation function to external
    return None

def disable_mod(SG):
    # SG.write('MODL 0')
    # SG_err_check(SG)
    return None

def query_mod_status(SG):
    # mod_status = SG.query('MODL?')
    # SG_err_check(SG)
    # if mod_status=='1\r\n':
    #     print('SG modulation is on...')
    #     mod_type = SG.query('TYPE?')
    #     SG_err_check(SG)
    #     if mod_type =='6\r\n':
    #         print('...and is set to IQ')
    #     elif mod_type == '4\r\n':
    #         print('... and is set to Pulse modulation')
    #         mod_type = SG.query('PFNC?')
    #         if mod_type == '5\r\n':
    #             print('... External')
    #         else:
    #             print('... Square' if mod_type == '3\r\n' else '... Noise (PRBS)')
    #     else:
    #         print('Modulation is set to '+mod_status+'. Set either 4 (for IQ) or 6 (for Pulse)')
    # else:
    #     print('SG modulation is off.')
    # return mod_status
    return None