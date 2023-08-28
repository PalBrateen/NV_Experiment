# SGcontrol

import pyvisa as visa
import sys
# Frequency unit multiplier definitions
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

# unicode = lambda s: str(s)

##-------------------- Function definitions--------------------
def initSG(serialaddr='3',modelName='SG384'):
    """
    Opens a RS-232 communication channel with the SG.
    It also clears the ESR and INSR registers as well as the LERR error buffer.
    ESR = Standard Event Status Register
    INSR = Instrument Status Register
    LERR = Last Error
    
    Parameters
    ----------
    serialaddr : int
        serial addess of SG connected via RS-232 interface.
    modelName : string
        Model name of SG, e.g.'SG384'

    Returns
    -------
    SG : object of SerialInstrument class
        DESCRIPTION.

    """
    SGaddr = 'asrl'+str(serialaddr)+'::instr'      # Construct instrument identifier from serial address
    rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager
    SG = rm.open_resource(SGaddr)         # Start 
    
    try:
        deviceID = SG.query('*IDN?')       # Try quering SG identity:
    except Exception as excpt:
        print('Error: could not query SG...') # ' Please check Serial address is correct and SG Serial communication is enabled. Exception details:', type(excpt).__name__,'.',excpt)
        sys.exit()
        
    # if 'Stanford Research Systems,'+modelName not in deviceID:
    #     print('Error: Instrument at this Serial address, (',serialaddr,') is not an SG '+modelName+'. When sent an identity query, \'*IDN?\', it returned ',deviceID,'. Please check the your SG signal generator\'s Serial address and/or model name.\n') 
    #     sys.exit()
    # Clear the ESR (Standard Event Status Register), INSR (Instrument Status Register) and LERR (Last Error Buffer):
    
    SG.write('*CLS')
    return SG

def SG_err_check(SG):
    err = SG.query('LERR?')
    if int(err) != 0:
        print('SG error: error code', int(err),'. Please refer to SG manual for a description of error codes.')
        sys.exit()
            
def enable_SG_op(SG):
    SG.write('ENBR1')
    SG_err_check(SG)

def disable_SG_op(SG):
    SG.write('ENBR0')
    SG_err_check(SG)
    print("SG Output disabled...")
    
def set_SG_amp(SG,RFamplitude, units='dBm'):
    SG.write('AMPR'+str(RFamplitude)+''+units)
    SG_err_check(SG)
    
def set_SG_freq(SG, freq, units='Hz'):
    """
    Sets frequency of the SG output.

    Parameters
    ----------
    SG : TYPE
        DESCRIPTION.
    freq : TYPE
        DESCRIPTION.
    units : TYPE, optional
        DESCRIPTION. The default is 'Hz'.

    Returns
    -------
    None.

    """
    
    SG.write('FREQ'+str(freq)+''+units)
    SG_err_check(SG)

def set_SG_disp(SG, disp):
    SG.write('disp'+str(disp))
    SG_err_check(SG)

def setup_SG_mod(SG,sequence):
    #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
    #and disables modulation for ESR, Rabi and T1 sequences.
    if sequence in ['ESGeq', 'RabiSeq', 'T1seq']:
        disable_mod(SG)
    elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
        enable_iq_mod(SG)
    else:
        print('Error in SGcontrol.py: unrecognised sequence name passed to setupSGmodulation.')
        sys.exit()
    
def setup_SG_pulse_mod(SG):
    """
    Setup external pulse modulation for producing MW pulses without MW switches.

    Parameters
    ----------
    SG : TYPE
        Object of SerialInstrument class.

    Returns
    -------
    None.

    """
    SG.write('modl1')
    SG_err_check(SG)
    SG.write('type4')
    SG_err_check(SG)
    SG.write('pfnc5')
    SG_err_check(SG)
    

def enable_iq_mod(SG):
    SG_err_check(SG)
    SG.write('MODL 1')     #Enable modulation
    SG_err_check(SG)
    SG.write('TYPE 6')     #Set modulation type to IQ
    SG_err_check(SG)
    SG.write('QFNC 5')     #Set IQ modulation function to external

def disable_mod(SG):
    SG.write('MODL 0')
    SG_err_check(SG)
    
def query_mod_status(SG):
    mod_status = SG.query('MODL?')
    SG_err_check(SG)
    if mod_status=='1\r\n':
        print('SG modulation is on...')
        mod_type = SG.query('TYPE?')
        SG_err_check(SG)
        if mod_type =='6\r\n':
            print('...and is set to IQ')
        elif mod_type == '4\r\n':
            print('... and is set to Pulse modulation')
            mod_type = SG.query('PFNC?')
            if mod_type == '5\r\n':
                print('... External')
            else:
                print('... Square' if mod_type == '3\r\n' else '... Noise (PRBS)')
        else:
            print('Modulation is set to '+mod_status+'. Set either 4 (for IQ) or 6 (for Pulse)')
    else:
        print('SG modulation is off.')
    return mod_status