# self.sgcontrol
#%%
import pyvisa as visa, sys, time
# Frequency unit multiplier definitions
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

class SignalGenerator():
    rm = visa.ResourceManager()             # Instantiate a resource manager; rm = object of type ResourceManager

    def __init__(self) -> None:
        self.sg_status = None
        self.addr = ''
        self.modelname = 'SG384'
        self.sg = None
        self.rfamp = 8
        self.mod_status = ''
        self.mod_type = ''
        self.freq = 2.87e9
        
        self.init_sg()
        

    def init_sg(self):
        """
        Opens a RS-232 communication channel with the self.sg.
        It also clears the Standard Event Status Register (ESR) and Instrument Status Register (INSR) registers as well as the Last Error (LERR) error buffer.

        Returns
        -------
        """

        addr = ["TCPIP0::169.254.59.63::inst0::INSTR", "ASRL5::INSTR"]
        searching = True
        while searching:
            for ad in addr:
                try:
                    res = self.rm.open_resource(ad)
                    if self.modelname not in res.query('*IDN?'):
                        print('Error: could not query SG... Retrying')
                    else:
                        print(f"SG connected via {ad}")
                        searching = False
                        break
                except:
                    pass
            time.sleep(1)
        self.addr = ad
        self.sg = self.rm.open_resource(ad)

        print(">>\x1b[38;2;250;250;0mSG384 Init'd...\x1b[0m")
        self.sg.write('*CLS')
        self.sg.write('disp2')
            
    def uninit_sg(self):
        self.sg.close()

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
        
    def set_sg_freq(self, freq, units='Hz'):
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
        self.sg.write('FREQ'+str(self.freq)+''+units)
        # self.sg_err_check()

    def set_sg_disp(self, disp):
        self.disp = disp
        self.sg.write('disp'+str(self.disp))
        self.sg_err_check()

    def setup_sg_mod(self, sequence):
        
        #Enables IQ modulation with an external source for T2, XY8 and correlation spectroscopy sequences
        #and disables modulation for ESR, Rabi and T1 sequences.
        if sequence in ['Eself.sgeq', 'RabiSeq', 'T1seq']:
            self.disable_sg_mod(self.sg)
        elif sequence in ['T2seq','XY8seq','correlSpecSeq']:
            self.enable_iq_mod(self.sg)
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
