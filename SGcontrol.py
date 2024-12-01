# SGcontrol
#%%
import pyvisa as visa
import sys, time, matplotlib.pyplot as plt
# Frequency unit multiplier definitions
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

# unicode = lambda s: str(s)
global SG

##-------------------- Function definitions--------------------
def init_sg(serialaddr='5',modelName='SG384'):
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
    global SG
    # SGaddr = 'asrl'+str(serialaddr)+'::instr'      # Construct instrument identifier from serial address
    SGaddr = "TCPIP0::169.254.59.63::inst0::INSTR"
    rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager
    SG = rm.open_resource(SGaddr)         # Start 
    
    try:
        deviceID = SG.query('*IDN?')       # Try quering SG identity
        print("\x10 \x1b[38;2;250;250;0mSG384 Init'd...\x1b[0m")
        # SG.write('*CLS')
        # SG.write('remt')
        SG.write('disp2')
        return SG
        
    except Exception as excpt:
        print('Error: could not query SG...') # ' Please check Serial address is correct and SG Serial communication is enabled. Exception details:', type(excpt).__name__,'.',excpt)
        sys.exit()
        
def uninit_sg():
    SG.close()

def SG_err_check():
    global SG
    err = SG.query('LERR?')
    if int(err) != 0:
        print('SG error: error code', int(err),'. Please refer to SG manual for a description of error codes.')
        sys.exit()
            
def enable_SG_op():
    global SG
    SG.write('ENBR1')
    SG_err_check()

def disable_SG_op():
    global SG
    SG.write('ENBR0')
    SG_err_check()
    print("SG Output disabled...")
    
def set_SG_amp(RFamplitude, units='dBm'):
    global SG
    SG.write('AMPR'+str(RFamplitude)+''+units)
    SG_err_check()
    
def set_SG_freq(freq, units='Hz'):
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
    global SG
    SG.write('FREQ'+str(freq)+''+units)
    # SG_err_check()

def set_SG_disp(disp):
    global SG
    SG.write('disp'+str(disp))
    SG_err_check()

def setup_SG_mod(sequence):
    global SG
    #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
    #and disables modulation for ESR, Rabi and T1 sequences.
    if sequence in ['ESGeq', 'RabiSeq', 'T1seq']:
        disable_mod(SG)
    elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
        enable_iq_mod(SG)
    else:
        print('Error in SGcontrol.py: unrecognised sequence name passed to setupSGmodulation.')
        sys.exit()
    
def setup_SG_pulse_mod():
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
    global SG
    SG.write('modl1')
    SG_err_check()
    SG.write('type4')
    SG_err_check()
    SG.write('pfnc5')
    SG_err_check()
    

def enable_iq_mod():
    global SG
    
    SG.write('MODL 1')     #Enable modulation
    SG_err_check()
    SG.write('TYPE 6')     #Set modulation type to IQ
    SG_err_check()
    SG.write('QFNC 5')     #Set IQ modulation function to external
    SG_err_check()

def disable_mod():
    global SG
    SG.write('MODL 0')
    SG_err_check()
    
def query_mod_status():
    global SG
    mod_status = SG.query('MODL?')
    SG_err_check()
    if mod_status=='1\r\n':
        print('SG modulation is on...')
        mod_type = SG.query('TYPE?')
        SG_err_check()
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
#%%
if __name__ == '__main__':
    # sg = init_sg()
    t=[]
    for i in range(0,int(1e1)):
        freq = 2.87e9 + 10*i
        t1 = time.perf_counter()
        # print(freq)
        set_SG_freq(freq)
        t2 = time.perf_counter()
        t.append((t2-t1)*1e3)

    plt.figure(); plt.hist(t)
    plt.figure(); plt.plot(t)
    #%%
    sg = init_sg()
# %%
    uninit_sg()
# %%
