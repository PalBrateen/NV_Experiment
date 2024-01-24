#%%-------------------------  USER INPUT  ---------------------------------------------------#
#PB clock frequency (in MHz):
PBclk = 500

#PB Connections ----------------------------------------------
#Enter below the bit numbers of the PB channels to which you connect your instruments, according to the definitions below.
#  PB_STARTtrig is the bit number of the PB channel used to generate the pulses fed to the Data Acquisition Card (DAQ) to trigger the start of data aquisition at each experiment scan point.
#  PB_DAQ is the bit number of the PB channel used to generate the pulses fed to the DAQ to gate/act as a sample clock to time the data aquisition.

# Start Trig -- PFI 15
# Conv CLK -- PFI 9
# Samp CLK -- PFI 14
PB_conv_clk = 0
PB_samp_clk = 1
PB_MW = 2
PB_AOM = 3
PB_start_trig = 4
PB_Q = 5
PB_I = 6
PB_camera = 7

# DAQ Connections-------------------------------------------------------
input_terminals = ["P6363/ai15"]#, "P6363/ai7"]       # Differential connection
# Detector connected to AI7 and PD to AI6
conv_clk_terminal = "PFI9"     # ADC conversion pulses
samp_clk_terminal = "PFI14"     # Start sampling from the channels in the scan list
start_trig_terminal = "PFI15"   # Trigger the start of data acquisiton
def cal_samp_rate():
    if len(input_terminals)==1:
        daq_max_samp_rate = 2e6    # Max samp rate in samp/CH/sec
    elif len(input_terminals)>1:
        daq_max_samp_rate = 1e6/len(input_terminals)    # Max samp rate in samp/CH/sec
    return daq_max_samp_rate

daq_max_samp_rate = cal_samp_rate()
min_voltage=-10  # Max/min voltage range of the photodiode signal
max_voltage=10

#SRS Connections-------------------------------------------------------
serialaddr = 5
model_name='SG384'

#------------------------- END OF USER INPUT ----------------------------------#

#%% Convert PB bit number to PB register address:
# To switch ON the samp_clk PB channel, a flag of 2^2=4 must be raised... which equals 0b100/0x4
# Similarly, 0x10(=16=2^4) passed as flag would turn ON the 4th bit (PB_I, here)
laser = 2**PB_AOM
start_trig = 2**PB_start_trig
samp_clk = 2**PB_samp_clk
MW = 2**PB_MW
conv_clk = 2**PB_conv_clk
I = 2**PB_I
Q = 2**PB_Q
camera = 2**PB_camera

# def update_connections():
    