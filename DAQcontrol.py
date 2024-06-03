#DAQ control
#%%
import connectionConfig as concfg
import nidaqmx
from  nidaqmx.constants import VoltageUnits, TerminalConfiguration, AcquisitionType, Edge
import sys

# def init_daq():

def close_daq_task(task):
    task.close()
    
def config_ai(Nsamples, rate=concfg.daq_max_samp_rate):
    """
    Create DAQmx task and assign channels with the configurations
    
    
    Parameters
    ----------
    Nsamples : int
        No. of samples to acquire from the detector.

    Returns
    -------
    read_task : nidaqmx.task.Task
        DESCRIPTION.

    """
    try:
        #Create and configure an analog input voltage task
        # NsampsPerDAQread=2*Nsamples
        read_task = nidaqmx.Task('NV_FL')
        for channels in concfg.input_terminals:
            read_task.ai_channels.add_ai_voltage_chan(channels,"", TerminalConfiguration.RSE, concfg.min_voltage, concfg.max_voltage, VoltageUnits.VOLTS)      # Create analog input channel for apd input
        print("\x10 Using channels: "+str(read_task.channels))
        
        read_task.timing.cfg_samp_clk_timing(rate, concfg.samp_clk_terminal, Edge.RISING, AcquisitionType.FINITE, Nsamples)                      # Configure sample clock: cfg_samp_clk_terminal_timing(rate, source=u'', active_edge=<Edge.RISING: 10280>, sample_mode=<AcquisitionType.FINITE: 10178>, samps_per_chan=1000)
        
        #Configure convert clock
        if len(concfg.input_terminals)>1:
            read_task.timing.ai_conv_src = concfg.conv_clk_terminal
        else:
            read_task.timing.ai_conv_src = concfg.samp_clk_terminal
        read_task.timing.ai_conv_active_edge = Edge.RISING
        #Configure start trigger
        readStartTrig = read_task.triggers.start_trigger
        if len(concfg.start_trig_terminal) != 0:
            readStartTrig.cfg_dig_edge_start_trig(concfg.start_trig_terminal,Edge.RISING)
    except Exception as excpt:
        print('Error configuring DAQ.\n Exception details:'+'\x1b[38;2;250;37;41m'+str(type(excpt).__name__)+'.'+str(excpt)+'\x1b[0m')
        close_daq_task(read_task)
        sys.exit()
    return read_task

def read_daq(task,Nsamples,timeout=120):
    try:
        counts = task.read(Nsamples, timeout)
    except Exception as excpt:
        print('Error: could not Read DAQ. Please check your DAQ\'s connections.\nException details:'+'\x1b[38;2;250;37;41m'+str(type(excpt).__name__)+'.'+str(excpt)+'\x1b[0m')
        sys.exit()
    return counts
    
# def read_daq_counter(task):
    
def coil_calibration(output_field: list):
    calibration = [26.5, 26.8, 43.5]
    output_aovoltage = [field/calib for field, calib in zip(output_field, calibration)]
    return output_aovoltage

def config_ao(dev='P6363'):
    """AO voltage sw operation."""
    try:
        ao_task = nidaqmx.Task()

        ao_task.ao_channels.add_ao_voltage_chan(dev+"/ao0")
        ao_task.ao_channels.add_ao_voltage_chan(dev+"/ao1")
        ao_task.ao_channels.add_ao_voltage_chan(dev+"/ao2")
    except Exception as excpt:
        print('Error configuring DAQ.\n Exception details:'+'\x1b[38;2;250;37;41m'+str(type(excpt).__name__)+'.'+str(excpt)+'\x1b[0m')
        close_daq_task(ao_task)
        sys.exit()
    return ao_task

def start_ao(ao_task, data):
    """Start AO"""
    print(ao_task.write(data, auto_start=True))

    return True

# examples...
# print(coil_calibration([60,60,60]))