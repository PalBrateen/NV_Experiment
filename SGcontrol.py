# self.sgcontrol
#%%
import pyvisa as visa, sys, time
from connectionConfig import sg_addr, sg_model_name
from parameter_system import Instrument

# Frequency unit multiplier definitions
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

class SignalGenerator(Instrument):
    rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager

    def __init__(self, name: str = "sg1") -> None:
        # Initialize instrument attributes BEFORE calling super().__init__()
        # because super() will call _register_parameters() which needs these
        self.sg_status = None
        self.addr = ''
        self.modelname = sg_model_name
        # self.sg: visa.Resource
        self.rfamp = 3
        self.mod_status = ''
        self.mod_type = ''
        self.freq = 2.87e9

        # Initialize Instrument base class (calls _register_parameters())
        super().__init__(name)

        # Now initialize hardware connection
        self.init_sg(sg_addr)
        

    def init_sg(self, addr):
        """
        Opens a RS-232 communication channel with the self.sg.
        It also clears the Standard Event Status Register (ESR) and Instrument Status Register (INSR) registers as well as the Last Error (LERR) error buffer.

        Returns
        -------
        """

        searched = False
        print("Searching SG384...")
        while not searched:
            for ad in addr:
                try:
                    res = self.rm.open_resource(ad)
                    if self.modelname not in res.query('*IDN?'):
                        print('❌ Error: could not query SG... Retrying')
                    else:
                        self.addr = ad
                        searched = True
                        break
                except:
                    pass
            time.sleep(1)
        print(f"SG connected via {self.addr}")        
        self.sg = self.rm.open_resource(self.addr)

        print(">>\x1b[38;2;250;250;0mSG384 Init'd...\x1b[0m")
        self.sg.write('*CLS')
        self.sg.write('disp2')
        self.sg.write("enbr1")
            
    def uninit_sg(self):
        self.sg.close()

    # ========================================================================
    # PARAMETER SYSTEM INTEGRATION (New functionality)
    # ========================================================================

    def _register_parameters(self):
        """Register all signal generator parameters with the parameter system."""
        self.register_parameter("frequency", 2.87e9, (0, 20e9), "Hz")
        self.register_parameter("power", 3, (-130, 25), "dBm")
        self.register_parameter("phase", 0, (0, 360), "degrees")

    def _update_instrument(self, parameter_name: str, new_value):
        """
        Update hardware with new parameter value (safe method with validation).

        This is called by change_parameter() after validation.
        Use for setup and outer loop parameters.
        """
        if parameter_name == "frequency":
            self.freq = new_value
            self.sg.write(f'FREQ{new_value}Hz')
        elif parameter_name == "power":
            self.rfamp = new_value
            self.sg.write(f'AMPR{new_value}dBm')
        elif parameter_name == "phase":
            # Phase control if needed
            pass

    def _set_frequency_direct(self, value: float):
        """
        Fast frequency setting for ESR inner loop optimization.

        This bypasses validation and is called directly via function pointer
        in ParameterSweep for maximum performance (<10μs overhead).

        Args:
            value: Frequency in Hz
        """
        print("inside _set_frequency_direct() ??")
        self.freq = value
        self.sg.write(f'FREQ{value}Hz')

    def _set_power_direct(self, value: float):
        """
        Fast power setting for power sweep optimization.

        Args:
            value: Power in dBm
        """
        self.rfamp = value
        self.sg.write(f'AMPR{value}dBm')

    # ========================================================================
    # EXISTING METHODS (Backward compatibility - kept unchanged)
    # ========================================================================

    def sg_err_check(self):
        err = self.sg.query('LERR?')
        if int(err) != 0:
            print('SG error: error code', int(err),'. Please refer to self.sg manual for a description of error codes.')
            sys.exit()
                
    def enable_sg_output(self):
        self.sg.write('ENBR1')
        self.sg_err_check()

    def disable_sg_output(self):
        self.sg.write('ENBR0')
        self.sg_err_check()
        print("SG output disabled...")
        
    def set_sg_amp(self, rfamp, units='dBm'):
        self.rfamp = rfamp
        self.sg.write('AMPR'+str(self.rfamp)+''+units)
        self.sg_err_check()
        
    def set_sg_freq(self, freq, unit='Hz'):
        """
        Sets frequency of the self.sg output.

        Parameters
        ----------
        freq : TYPE
            DESCRIPTION.
        units : TYPE, optional
            DESCRIPTION. The default is 'Hz'.

        Returns
        -------
        None.

        """
        self.freq = freq
        self.sg.write('FREQ'+str(self.freq)+''+unit)
        # self.sg_err_check()

    def set_sg_disp(self, disp):
        self.disp = disp
        self.sg.write('disp'+str(self.disp))
        self.sg_err_check()

    def setup_sg_mod(self, sequence):
        
        #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
        #and disables modulation for ESR, Rabi and T1 sequences.
        if sequence in ['Eself.sgeq', 'RabiSeq', 'T1seq']:
            self.disable_sg_mod()
        elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
            self.enable_iq_mod()
        else:
            print('Error in self.sgcontrol.py: unrecognised sequence name passed to setupself.sgmodulation.')
            sys.exit()
        
    def setup_sg_pulse_mod(self):
        """
        Setup external pulse modulation for producing MW pulses without MW switches.

        Parameters
        ----------
        self.sg : TYPE
            Object of SerialInstrument class.

        Returns
        -------
        None.

        """
        self.sg.write('modl1')
        self.sg_err_check()
        self.sg.write('type4')
        self.sg_err_check()
        self.sg.write('pfnc5')
        self.sg_err_check()
        
    def enable_iq_mod(self):
        self.sg.write('MODL 1')     #Enable modulation
        self.sg_err_check()
        self.sg.write('TYPE 6')     #Set modulation type to IQ
        self.sg_err_check()
        self.sg.write('QFNC 5')     #Set IQ modulation function to external
        self.sg_err_check()

    def disable_sg_mod(self):
        self.sg.write('MODL 0')
        self.sg_err_check()
        
    def query_mod_status(self):
        self.mod_status = self.sg.query('MODL?')
        self.sg_err_check()
        if self.mod_status=='1\r\n':
            print('self.sg modulation is on...')
            self.mod_type = self.sg.query('TYPE?')
            self.sg_err_check()
            if self.mod_type =='6\r\n':
                print('...and is set to IQ')
            elif self.mod_type == '4\r\n':
                print('... and is set to Pulse modulation')
                self.mod_type = self.sg.query('PFNC?')
                if self.mod_type == '5\r\n':
                    print('... External')
                else:
                    print('... Square' if self.mod_type == '3\r\n' else '... Noise (PRBS)')
            else:
                print('Modulation is set to '+self.mod_status+'. Set either 4 (for IQ) or 6 (for Pulse)')
        else:
            print('self.sg modulation is off.')


class SignalGenerator_sim(Instrument):
    rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager

    def __init__(self, name: str = "sg1") -> None:
        # Initialize instrument attributes BEFORE calling super().__init__()
        self.sg_status = None
        self.addr = ''
        self.modelname = 'SG384'
        self.sg = None
        self.rfamp = 8
        self.mod_status = ''
        self.mod_type = ''
        self.freq = 2.87e9

        # Initialize Instrument base class (calls _register_parameters())
        super().__init__(name)

        # Now initialize simulated hardware
        self.init_sg()
        

    def init_sg(self):
        """
        Opens a RS-232 communication channel with the self.sg.
        It also clears the Standard Event Status Register (ESR) and Instrument Status Register (INSR) registers as well as the Last Error (LERR) error buffer.

        Returns
        -------
        """
        print(">>\x1b[38;2;250;250;0mSimSG: Init'd...\x1b[0m")
        self.sg = None

    def uninit_sg(self):
        print(">>\x1b[38;2;250;250;0mSimSG: Uninit'd...\x1b[0m")

    # ========================================================================
    # PARAMETER SYSTEM INTEGRATION (New functionality)
    # ========================================================================

    def _register_parameters(self):
        """Register all signal generator parameters with the parameter system."""
        self.register_parameter("frequency", 2.87e9, (0, 20e9), "Hz")
        self.register_parameter("power", 8, (-130, 25), "dBm")
        self.register_parameter("phase", 0, (0, 360), "degrees")

    def _update_instrument(self, parameter_name: str, new_value):
        """
        Update simulated hardware with new parameter value.

        This is called by change_parameter() after validation.
        """
        if parameter_name == "frequency":
            self.freq = new_value
            print(f"SimSG: Frequency set to {new_value} Hz")
        elif parameter_name == "power":
            self.rfamp = new_value
            print(f"SimSG: Power set to {new_value} dBm")
        elif parameter_name == "phase":
            print(f"SimSG: Phase set to {new_value} degrees")

    def _set_frequency_direct(self, value: float):
        """
        Fast frequency setting for ESR inner loop optimization (simulated).

        Args:
            value: Frequency in Hz
        """
        self.freq = value
        # Silently update (no print for performance)

    def _set_power_direct(self, value: float):
        """
        Fast power setting for power sweep optimization (simulated).

        Args:
            value: Power in dBm
        """
        self.rfamp = value
        # Silently update (no print for performance)

    # ========================================================================
    # EXISTING METHODS (Backward compatibility - kept unchanged)
    # ========================================================================
                
    def enable_sg_output(self):
        print("SimSG: output enabled...")
    
    def disable_sg_output(self):
        print("SimSG: output disabled...")
        
    def set_sg_amp(self, rfamp, units='dBm'):
        print(f"SimSG: RF Amplitude set to {rfamp} dBm...")
        
    def set_sg_freq(self, freq, units='Hz'):
        print(f"SimSG: Frequency set to {freq} Hz...")

    def set_sg_disp(self, disp):
        print(f"SimSG: Display set to {disp}...")

    def setup_sg_mod(self, sequence):
        pass
        # #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
        # #and disables modulation for ESR, Rabi and T1 sequences.
        # if sequence in ['Eself.sgeq', 'RabiSeq', 'T1seq']:
        #     self.disable_sg_mod(self.sg)
        # elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
        #     self.enable_iq_mod(self.sg)
        # else:
        #     print('Error in self.sgcontrol.py: unrecognised sequence name passed to setupself.sgmodulation.')
        #     sys.exit()
        
    def setup_sg_pulse_mod(self):
        print("SimSG: Pulse modulation setup...")
        
    def enable_iq_mod(self):
        # self.sg.write('MODL 1')     #Enable modulation
        # self.sg_err_check()
        # self.sg.write('TYPE 6')     #Set modulation type to IQ
        # self.sg_err_check()
        # self.sg.write('QFNC 5')     #Set IQ modulation function to external
        # self.sg_err_check()
        pass

    def disable_sg_mod(self):
        print("SimSG: Modulation disabled...")
        
    def query_mod_status(self):
        # self.mod_status = self.sg.query('MODL?')
        # self.sg_err_check()
        # if self.mod_status=='1\r\n':
        #     print('self.sg modulation is on...')
        #     self.mod_type = self.sg.query('TYPE?')
        #     self.sg_err_check()
        #     if self.mod_type =='6\r\n':
        #         print('...and is set to IQ')
        #     elif self.mod_type == '4\r\n':
        #         print('... and is set to Pulse modulation')
        #         self.mod_type = self.sg.query('PFNC?')
        #         if self.mod_type == '5\r\n':
        #             print('... External')
        #         else:
        #             print('... Square' if self.mod_type == '3\r\n' else '... Noise (PRBS)')
        #     else:
        #         print('Modulation is set to '+self.mod_status+'. Set either 4 (for IQ) or 6 (for Pulse)')
        # else:
        #     print('self.sg modulation is off.')
        pass
#%%
# if __name__ == '__main__':
#     # sg = init_sg()
#     t=[]
#     for i in range(0,int(1e1)):
#         freq = 2.87e9 + 10*i
#         t1 = time.perf_counter()
#         # print(freq)
#         set_self.sg_freq(freq)
#         t2 = time.perf_counter()
#         t.append((t2-t1)*1e3)

#     plt.figure(); plt.hist(t)
#     plt.figure(); plt.plot(t)
#     #%%
#     sg = init_sg()
# # %%
#     uninit_sg()
# %%
