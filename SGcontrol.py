# self._instrcontrol
#%%
import pyvisa as visa, sys, time, logging
from typing import Union, List, Tuple
from enum import IntEnum
from experiment_config import SG_ADDR
# Frequency unit multiplier definitions
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

def printt(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")

class ModulationType(IntEnum):
    NONE = -1
    AMPLITUDE = 0
    FREQUENCY = 1
    PHASE = 2
    SWEEP = 3
    PULSE = 4
    BLANK = 5
    IQ = 6

class ModulationFunction(IntEnum):
    SINE = 0
    RAMP = 1
    TRIANGLE = 2
    SQUARE = 3
    NOISE = 4
    EXTERNAL = 5

class SGDISPLAY(IntEnum):
    MODULATION_TYPE = 0
    MODULATION_FUNCTION = 1
    FREQUENCY = 2
    PHASE = 3
    MODULATION_RATE = 4
    MODULATION_PERIOD = 4
    MODULATION_DEVIATION = 5
    DUTY_CYCLE = 5
    AMPLITUDE_NTYPE = 6
    AMPLITUDE_BNC= 7
    AMPLITUDE_RF_DOUBLER = 8
    AMPLITUDE_CLOCK = 9
    OFFSET_BNC = 10
    OFFSET_REAR_DC = 11
    OFFSET_CLOCK = 12

class ErrorCodes(IntEnum):
    NO_ERROR = 0
    ILLEGAL_VALUE = 10
    ILLEGAL_MODE = 11
    NOT_ALLOWED = 12
    RECALL_FAILED = 13
    NO_CLOCK_OPTION = 14
    NO_RF_DOUBLER_OPTION = 15
    NO_IQ_OPTION = 16
    FAILED_SELF_TEST = 17

    LOST_DATA = 30
    NO_LISTENER = 32

    FAILED_ROM_CHECK = 40
    FAILED_EEPROM_CHECK = 42
    FAILED_FPGA_CHECK = 43
    FAILED_SRAM_CHECK = 44
    FAILED_GPIB_CHECK = 45
    FAILED_LF_DDS_CHECK = 46
    FAILED_RF_DDS_CHECK = 47
    FAILED_20MHZ_PLL = 48
    FAILED_100MHZ_PLL = 49
    FAILED_19MHZ_PLL = 50
    FAILED_1GHZ_PLL = 51
    FAILED_4GHZ_PLL = 52
    FAILED_DAC = 53

    ILLEGAL_COMMAND = 110
    UNDEFINED_COMMAND = 111
    ILLEGAL_QUERY = 112
    ILLEGAL_SET = 113
    NULL_PARAMETER = 114
    EXTRA_PARAMETERS = 115
    MISSING_PARAMETERS = 116
    PARAMETER_OVERFLOW = 117
    INVALID_FLOAT = 118
    INVALID_INTEGER = 120
    INTEGER_OVERFLOW = 121
    INVALID_HEXADECIMAL = 122
    SYNTAX_ERROR = 126
    ILLEGAL_UNITS = 127
    MISSING_UNITS = 128

    COMMUNICATION_ERROR = 170
    OVERRUN = 171

    TOO_MANY_ERRORS = 254

    @classmethod
    def get_description(cls, code: int) -> str:
        """Get human-readable description for error code."""
        descriptions = {
            0: "No error",
            10: "Illegal value - parameter out of range",
            11: "Illegal mode - command not valid in current mode",
            12: "Not allowed - operation not permitted",
            13: "Recall failed - settings recall error",
            14: "No clock option installed",
            15: "No RF doubler option installed",
            16: "No IQ option installed",
            17: "Failed self test",
            30: "Lost data - buffer overflow",
            32: "No listener - GPIB communication error",
            40: "Failed ROM check",
            42: "Failed EEPROM check",
            43: "Failed FPGA check",
            44: "Failed SRAM check",
            45: "Failed GPIB check",
            46: "Failed LF DDS check",
            47: "Failed RF DDS check",
            48: "Failed 20MHz PLL",
            49: "Failed 100MHz PLL",
            50: "Failed 19MHz PLL",
            51: "Failed 1GHz PLL",
            52: "Failed 4GHz PLL",
            53: "Failed DAC",
            110: "Illegal command",
            111: "Undefined command",
            112: "Illegal query",
            113: "Illegal set - cannot set this parameter",
            114: "Null parameter - missing value",
            115: "Extra parameters provided",
            116: "Missing parameters",
            117: "Parameter overflow",
            118: "Invalid floating point number",
            120: "Invalid integer",
            121: "Integer overflow",
            122: "Invalid hexadecimal",
            126: "Syntax error",
            127: "Illegal units",
            128: "Missing units",
            170: "Communication error",
            171: "Overrun - too much data",
            254: "Too many errors - error buffer full",
        }
        return descriptions.get(code, f"Unknown error (code {code})")

#%
class SignalGenerator():
    Hz = 1
    kHz = 1e3
    MHz = 1e6
    GHz = 1e9

    def __init__(self, auto_init_hardware=True) -> None:
        self._instr: visa.Resource
        self._instr_status: str
        self.addr: str# = ''
        self.modelname: str = 'SG384'
        # self._instr: visa.Resource
        self.amp_rf: float = -4
        self.amp_bnc: float = -4
        
        self.status_ntype: bool|int# = 1
        self.status_bnc: bool|int# = 0
        self.freq: float# = 2.87e9
        self.display: str = 'FREQUENCY'#'AMPLITUDE_NTYPE'
        self.mod_status: bool|int# = 0
        self.mod_type: str# = 'NONE'
        self.mod_func: str# = 'external'
        self.mod_rate: float# = 1e3
        self.mod_dev: float# = 0.0

        self.init(SG_ADDR) if auto_init_hardware else None
    
    def init(self, addr=SG_ADDR):
        """
        Opens a RS-232 communication channel with the self._instr.
        It also clears the Standard Event Status Register (ESR) and Instrument Status Register (INSR) registers as well as the Last Error (LERR) error buffer.

        Returns
        -------
        """
        found = False
        rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager
        printt("Searching SG384...")
        while not found:
            for ad in addr:
                try:
                    res = rm.open_resource(ad)
                    if self.modelname not in res.query('*IDN?'):
                        printt('❌ Error: could not query SG... Retrying')
                    else:
                        self.addr = ad
                        found = True
                        self._instr = rm.open_resource(self.addr)
                        break
                except:
                    pass
            time.sleep(1)
        self._instr_status = 'Connected' if found else 'Not connected'
        printt(f"✔\x1b[38;2;250;250;0m SG384 Init via {self.addr}!\x1b[0m")
        self._instr.write('*CLS')
        self.set_display(self.display)
        # self.enable_ntype(1)
        
    def uninit(self):
        self._instr.close()

    def write(self, cmd: str):
        self._instr.write(cmd)
        # self._check_err()
    
    def query(self, cmd: str) -> str:
        res = self._instr.query(cmd)
        # self._check_err()
        return res
    
    def _check_err(self):
        res = self._instr.query('LERR?')
        if int(res) != 0:
            printt(f'❌\x1b[38;2;0;100;250m SG error:: Code {int(res)}: {ErrorCodes.get_description(int(res))}\x1b[0m')
            # sys.exit()
            self._instr_status = 'Error'
    
    def query_all(self):
        try:
            self.amp_rf = float(self.query('AMPR?'))
            self.amp_bnc = float(self.query('AMPL?'))
            self.status_ntype = int(self.query('ENBR?'))
            self.status_bnc = int(self.query('ENBL?'))
            self.freq = float(self.query('FREQ?'))
            self.display = SGDISPLAY(int(self.query('DISP?'))).name
            self.mod_status = int(self.query('MODL?'))
            if self.mod_status:
                self.mod_type = ModulationType(int(self.query('TYPE?'))).name
                if self.mod_type == 'PULSE':
                    self.mod_func = ModulationFunction(int(self.query('PFNC?'))).name
                elif self.mod_type == 'SWEEP':
                    self.mod_func = ModulationFunction(int(self.query('SFNC?'))).name
                elif self.mod_type == 'IQ':
                    self.mod_func = ModulationFunction(int(self.query('QFNC?'))).name
                else:
                    self.mod_func = ModulationFunction(int(self.query('MFNC?'))).name
                self.mod_rate = float(self.query('RATE?'))
                if self.mod_type == 'FREQUENCY':
                    self.mod_dev = float(self.query('FDEV?'))
                elif self.mod_type == 'AMPLITUDE':
                    self.mod_dev = float(self.query('ADEP?'))
                elif self.mod_type == 'PULSE':
                    self.mod_dev = float(self.query('ADEP?'))   # TODO: correct this

        except Exception as e:
            logging.exception(f"❌ Error querying SG parameters: {e}")
            pass
    
    def set_display(self, display: str|int='FREQUENCY'):
        cmd = 'DISP'
        if isinstance(display, str):
            self.write(f'{cmd} {SGDISPLAY[str(display).upper()].value}')
        elif isinstance(display, int):
            self.write(f'DISP {SGDISPLAY(display)}')
        
        self.display = SGDISPLAY(int(self.query(f'{cmd}?'))).name

    def enable_ntype(self, enable:Union[bool, int]=True):
        cmd = 'ENBR'
        self.write(f'{cmd} {str(int(enable))}')

        self.status_ntype = int(self.query(f'{cmd}?'))
        if self.status_ntype:
            print("N-type ON...")
        else:
            print("N-type OFF...")

    def enable_bnc(self, enable:Union[bool, int]=True):
        cmd = 'ENBL'
        self.write(f'{cmd} {str(int(enable))}')

        self.status_bnc = int(self.query(f'{cmd}?'))
        if self.status_bnc:
            print("BNC ON...")
        else:
            print("BNC OFF...")
        
    def set_amp_rf(self, amp:float, units:str='dBm'):
        cmd = 'AMPR'
        self.write(f'{cmd} {str(amp)}{units}')
        
        self.amp_rf = float(self.query(f'{cmd}?'))

    set_amp = set_amp_rf
    
    def set_amp_bnc(self, amp:float, units='dBm'):
        cmd = 'AMPL'
        self.write(f'{cmd} {str(amp)}{units}')
        
        self.amp_bnc = float(self.query(f'{cmd}?'))
        
    def set_freq(self, freq:float, unit='Hz'):
        cmd = 'FREQ'
        self.write(f'{cmd} {str(freq)}{unit}')
        # 
        self.freq = freq
        # printt(self.freq)
    
    def enable_modulation(self, enable:Union[bool, int]=True):
        cmd = 'MODL'
        self.write(f'{cmd} {str(int(enable))}')
        
        self.mod_status = int(self.query(f'{cmd}?'))

    def set_mod_type(self, mod_type:str):
        cmd = 'TYPE'
        # if isinstance(mod_type, ModulationType) or isinstance(mod_type, str):
        self.write(f'{cmd} {ModulationType[str(mod_type).upper()].value}')   
        # elif isinstance(mod_type, int):
        #     self.write(f'TYPE {str(mod_type)}')
        
        self.mod_type = ModulationType(int(self.query(f'{cmd}?'))).name

    def set_mod_func(self, fm_func:str):
        cmd = None
        if self.mod_type in ['AMPLITUDE', 'FREQUENCY', 'PHASE']:
            cmd = 'MFNC'
            # if isinstance(fm_func, ModulationType) or isinstance(fm_func, str):
            self.write(f'{cmd} {ModulationFunction[str(fm_func).upper()].value}')
            # elif isinstance(fm_func, int):
            #     self.write(f'MFNC {str(fm_func)}')
        elif self.mod_type == 'PULSE':
            cmd = 'PFNC'
            # if isinstance(fm_func, ModulationType) or isinstance(fm_func, str):
            self.write(f'{cmd} {ModulationFunction[str(fm_func).upper()].value}')
            # elif isinstance(fm_func, int):
            #     self.write(f'PFNC {str(fm_func)}')
        elif self.mod_type == 'SWEEP':
            cmd = 'SFNC'
            # if isinstance(fm_func, ModulationType) or isinstance(fm_func, str):
            self.write(f'{cmd} {ModulationFunction[str(fm_func).upper()].value}')
            # elif isinstance(fm_func, int):
            #     self.write(f'SFNC {str(fm_func)}')
        elif self.mod_type == 'IQ':
            cmd = 'QFNC'
            # if isinstance(fm_func, ModulationType) or isinstance(fm_func, str):
            self.write(f'{cmd} {ModulationFunction[str(fm_func).upper()].value}')
            # elif isinstance(fm_func, int):
            #     self.write(f'QFNC {str(fm_func)}')
        assert cmd
        self.mod_func = ModulationFunction(int(self.query(f'{cmd}?'))).name

    def set_mod_rate(self, mod_rate:float=1e3):
        cmd = 'RATE'
        self.write(f'{cmd} {str(mod_rate)}')       #Set modulation rate
        
        self.mod_rate = float(self.query(f'{cmd}?'))
    
    set_mod_freq = set_mod_rate
    
    def set_mod_dev(self, mod_dev:float):
        cmd = None
        if self.mod_type.upper() == 'FREQUENCY':
            cmd = 'FDEV'
            self.write(f'{cmd} {str(mod_dev)}')       #Set frequency deviation
        elif self.mod_type.upper() == 'AMPLITUDE':
            cmd = 'ADEP'
            self.write(f'{cmd} {str(mod_dev)}')       #Set amplitude depth
        # elif
        assert cmd
        self.mod_dev = float(self.query(f'{cmd}?'))

    def setup_sg_mod(self, sequence:str):
        #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
        #and disables modulation for ESR, Rabi and T1 sequences.
        if sequence in ['Eself._instreq', 'RabiSeq', 'T1seq']:
            self.enable_modulation(1)
        elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
            self.enable_iq_mod()
        else:
            printt('Error in self._instrcontrol.py: unrecognised sequence name passed to setupself._instrmodulation.')
            sys.exit()
        
    def setup_ext_pulse_mod(self):
        """
        Setup external pulse modulation for producing MW pulses without MW switches.

        Parameters
        ----------
        self._instr : TYPE
            Object of SerialInstrument class.

        Returns
        -------
        None.

        """
        self.set_mod_type('PULSE')
        self.set_mod_func('external')
        self.enable_modulation()
    
    def setup_sg_fm(self, func:str="EXTERNAL", dev:float=1e3, modfreq:float=1e3,):
        self.set_mod_type('FREQUENCY')
        self.set_mod_func(func)
        self.set_mod_dev(dev)
        if func!="EXTERNAL":
            self.set_mod_rate(modfreq)

    def setup_sg_am(self, func:str="EXTERNAL", dev:float=1e3, modfreq:float=1e3,):
        self.set_mod_type('amplitude')
        self.set_mod_func(func)
        self.set_mod_dev(dev)
        if func!="EXTERNAL":
            self.set_mod_rate(modfreq)

    def enable_iq_mod(self, func='EXTERNAL'):
        self.set_mod_type('iq')
        self.set_mod_func(func)
        self.enable_modulation(1)
        
    def query_mod_status(self):
        self.mod_status = self._instr.query('MODL?')
        
        if self.mod_status=='1\r\n':
            printt('self._instr modulation is on...')
            self.mod_type = self._instr.query('TYPE?')
            
            if self.mod_type =='6\r\n':
                printt('...and is set to IQ')
            elif self.mod_type == '4\r\n':
                printt('... and is set to Pulse modulation')
                self.mod_type = self._instr.query('PFNC?')
                if self.mod_type == '5\r\n':
                    printt('... External')
                else:
                    printt('... Square' if self.mod_type == '3\r\n' else '... Noise (PRBS)')
            else:
                printt(f'Modulation is set to {self.mod_status}. Set either 4 (for IQ) or 6 (for Pulse)')
        else:
            printt('self._instr modulation is off.')


class SignalGenerator_sim():

    def __init__(self) -> None:
        self._instr_status = None
        self.addr: str = ''
        self.modelname: str = 'SG384'
        # self._instr: visa.Resource
        self.amp_rf: float = -4
        self.amp_bnc: float = -4
        
        self.status_ntype: bool|int = 1
        self.status_bnc: bool|int = 0
        self.freq: float = 2.87e9
        self.display: SGDISPLAY|str = SGDISPLAY.FREQUENCY
        self.mod_status: bool|int = 0
        self.mod_type: ModulationType|str = ModulationType.NONE
        self.mod_func: ModulationFunction|str = ModulationFunction.EXTERNAL
        self.mod_rate: float = 1e3
        self.mod_dev: float = 0.0

        self.init(SG_ADDR)
    
    def init(self, addr):
        printt("Enabled")
        
    def uninit(self):
        printt("Closed")

    def _check_err(self):
        pass
    
    def set_display(self, display:Union[SGDISPLAY, str]=SGDISPLAY.FREQUENCY):
        print("Enabled")
        self.display = display

    def enable_ntype(self, enable:Union[bool, int]=True):
        if enable:
            print("N-type ON...")
        else:
            print("N-type OFF...")
        self.status_ntype = int(enable)

    def enable_bnc(self, enable:Union[bool, int]=True):
        if enable:
            print("BNC ON...")
        else:
            print("BNC OFF...")
        self.status_bnc = int(enable)
    
    def set_amp_rf(self, amp_rf:float, units='dBm'):
        print("Enabled")
        self.amp_rf = amp_rf

    set_amp = set_amp_rf
    
    def set_amp_bnc(self, amp:float, units='dBm'):
        print("Enabled")
        self.amp_bnc = amp
        
    def set_freq(self, freq:float, unit='Hz'):
        print("Enabled")
        self.freq = freq
    
    def enable_modulation(self, enable:Union[bool, int]=True):
        print("Enabled")
        self.mod_status = int(enable)

    def set_mod_type(self, mod_type:Union[ModulationType, str]):
        print("Enabled")
        self.mod_type = mod_type

    def set_mod_func(self, fm_func:Union[ModulationFunction, str]):
        print("Enabled")
        self.mod_func = fm_func

    def set_mod_rate(self, mod_rate:float=1e3):
        print("Enabled")
        self.mod_rate = mod_rate
    
    set_mod_freq = set_mod_rate
    
    def set_mod_dev(self, mod_dev:float):
        print("Enabled")
        self.mod_dev = mod_dev

    def setup_sg_mod(self, sequence:str):
        print("Enabled")
        
    def setup_ext_pulse_mod(self):
        print("Enabled")
    
    def setup_sg_fm(self, func:str="EXTERNAL", dev:float=1e3, modfreq:float=1e3,):
        print("Enabled")

    def setup_sg_am(self, func:str="EXTERNAL", dev:float=1e3, modfreq:float=1e3,):
        print("Enabled")

    def enable_iq_mod(self, func='EXTERNAL'):
        print("Enabled")
        
    def query_mod_status(self):
        pass
#%%
# if __name__ == '__main__':
#     # sg = init()
#     t=[]
#     for i in range(0,int(1e1)):
#         freq = 2.87e9 + 10*i
#         t1 = time.perf_counter()
#         # print(freq)
#         set_self._instr_freq(freq)
#         t2 = time.perf_counter()
#         t.append((t2-t1)*1e3)

#     plt.figure(); plt.hist(t)
#     plt.figure(); plt.plot(t)
#     #%%
#     sg = init()
# # %%
#     uninit()
#%%
if __name__ == '__main__':
    sg = SignalGenerator()
    sg.set_amp(5)
    sg.set_freq(2.87e9)
    sg.set_display('amplitude_ntype')
    sg.enable_ntype(0)
    sg.setup_ext_pulse_mod()
    sg.uninit()
