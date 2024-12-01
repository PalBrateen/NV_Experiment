#DAQ control
#%%
import connectionConfig as concfg, nidaqmx, sys, os, numpy as np, matplotlib.pyplot as plt
from  nidaqmx.constants import VoltageUnits, TerminalConfiguration, AcquisitionType, Edge, Level, TriggerType, RegenerationMode
# def init_daq():
global samp_rate

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
    # calibration = [19,19,37]
    output_aovoltage = [field/calib for field, calib in zip(output_field, calibration)]
    return output_aovoltage

def config_ao(dev: str):
    """AO voltage sw operation."""
    try:
        ao_task = nidaqmx.Task("Coil_control")

        ao_task.ao_channels.add_ao_voltage_chan(dev+"/ao0")
        ao_task.ao_channels.add_ao_voltage_chan(dev+"/ao1")
        ao_task.ao_channels.add_ao_voltage_chan(dev+"/ao2")
    except Exception as excpt:
        print('Error configuring DAQ.\n Exception details:'+'\x1b[38;2;250;37;41m'+str(type(excpt).__name__)+'.'+str(excpt)+'\x1b[0m')
        close_daq_task(ao_task)
        # sys.exit()
    return ao_task

def start_ao(ao_task, data):
    """Start AO"""
    print(f"AO data = {data}")
    ao_task.write(data, auto_start=True)

    return True

def triggered_ao_data(direction, rot_angle, align_field, amp, freq):
    if direction == 'x':
        var = triggered_ao_data_x(rot_angle, align_field, amp, freq)
    elif direction == 'z':
        var = triggered_ao_data_z(rot_angle, align_field, amp, freq)
    return var
    
def triggered_ao_data_z(rot_angle, align_field, amp, freq=10):
    theta = 0       
    phi = 0
    phase = 90
    # amp = 2.63
    ti = 0
    tf = 1/freq     # tf is one time period

    fig = plt.figure()
    title = f'{rot_angle} deg of {freq} Hz theta_{theta} phi_{phi}'
    fig.suptitle(title)
    fig.canvas.manager.set_window_title(title)
    #
    global samp_rate
    samp_rate = 100e3     # daq sampling rate
    dt = 1/samp_rate    # time diff between two samples, sampling time
    # n_points = (tf*freq - ti)/dt     # total points for the waveform

    linespec = '-'
    fraction = rot_angle/360
    theta = theta *np.pi/180
    phi = phi *np.pi/180
    phase = phase *np.pi/180
    
    # t = np.linspace(ti,tf,int(n_points)+1, endpoint=True)      # time array
    t = np.arange(ti, tf, dt)
    t_fraction = t[0:int(fraction*len(t))+1]    # fT
    t_rem_fraction = t[int(fraction*len(t))+1:]
    t_array = t_fraction.copy()      # time array for actual b operation

    bx = amp*(np.cos(phi)*np.cos(2*np.pi*freq*t_fraction + np.pi/2) + np.sin(theta)*np.sin(phi)*np.sin(2*np.pi*freq*t_fraction + np.pi/2))
    by = amp*(np.sin(phi)*np.cos(2*np.pi*freq*t_fraction + np.pi/2) - np.sin(theta)*np.cos(phi)*np.sin(2*np.pi*freq*t_fraction + np.pi/2))
    bz = (align_field[-1]/np.abs(align_field[-1])) * amp*np.cos(theta)*np.sin(2*np.pi*freq*t_fraction + np.pi/2)

    b1 = np.transpose(np.array([bx, by, bz]))

    # plt.figure(); plt.plot(t_fraction, b1)
    ax = fig.add_subplot(2,4,1)
    ax.plot(t_fraction*1e3, b1,linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b1,linespec)

    # 2nd cycle RCP 1
    # with the same dt and shifted starting time point (t_fraction[-1]), form the time array of length same as that of prev cycle

    ti = t_array[-1] + dt
    temp = ti + t_fraction[0:-1]
    # temp = temp[0:-1]
    t_array = np.hstack((t_array, temp))

    b2 = np.flip(b1, axis=0)    # b1[1:len(t_fraction)+1]
    # print(b2.shape)
    b2 = b2[1:,:]
    # b2 = np.transpose(np.array([np.nan_to_num(bx*bx[-1]/np.absolute(bx[-1]), nan=0.0), np.nan_to_num(by*by[-1]/np.absolute(by[-1]), nan=0.0), -bz]))
    # plt.figure(); plt.plot(b2)
    ax = fig.add_subplot(2,4,2)
    ax.plot(temp*1e3, b2,linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b2,linespec)

    # 3rd cycle RCP 2

    ti = t_array[-1] + dt
    temp = ti + t_fraction[0:-1]
    t_array = np.hstack((t_array, temp))

    b3 = np.transpose(np.array([-bx, -by, bz]))
    b3 = b3[1:,:]
    # plt.figure(); plt.plot(b3)
    ax = fig.add_subplot(2,4,3)
    ax.plot(temp*1e3, b3,linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b3,linespec)

    # 4th cycle LCP

    ti = t_array[-1] + dt
    temp = ti + t_fraction[0:-2]
    t_array = np.hstack((t_array, temp))

    # b4 = b2     # this is wrong... It creates a reference to the same array in the memory
    b4 = b2.copy()
    b4[:,0:1] = -b2[:,0:1]
    b4 = b4[:-1,:]
    # plt.figure(); plt.plot(b4)
    ax = fig.add_subplot(2,4,4)
    ax.plot(temp*1e3, b4, linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot( b4,linespec)

    #
    # tzero = 100 /1e3        # enter value in ms
    # ti = t_array[-1] + dt
    # temp = ti + np.array(np.arange(0,tzero,dt))

    # Bzero = np.zeros((int(samp_rate*tzero), b1.shape[-1]))
    # t_array = np.hstack((t_array, temp))

    # b = np.ascontiguousarray(np.vstack((b1, b2, b3, b4, Bzero)))
    # b_before = np.array([b1[0,:] for _ in range(0,100)])
    # b_after = np.array([b4[-1,:] for _ in range(0,100)])
    b = np.ascontiguousarray(np.vstack((b1, b2, b3, b4)))
    # b = np.ascontiguousarray(np.vstack((b_before, b1, b2, b3, b4, b_after)))
    # t_array = 
    # plt.figure(); plt.plot(b)
    ax = fig.add_subplot(2,4,(5,8))
    ax.plot(t_array*1e3, b, linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b,linespec)
    ax.legend(['X', 'Y', 'Z'])

    return [dt, t_array, b]

def triggered_ao_data_x(rot_angle, align_field, amp, freq=10):
    theta = 0       
    phi = 0
    phase = 90
    # amp = 2.63
    ti = 0
    tf = 1/freq     # tf is one time period

    fig = plt.figure()
    title = f'{rot_angle} deg of {freq} Hz theta_{theta} phi_{phi}'
    fig.suptitle(title)
    fig.canvas.manager.set_window_title(title)
    #
    global samp_rate
    samp_rate = 100e3     # daq sampling rate
    dt = 1/samp_rate    # time diff between two samples, sampling time
    # n_points = (tf*freq - ti)/dt     # total points for the waveform

    linespec = '-'
    fraction = rot_angle/360
    theta = theta *np.pi/180
    phi = phi *np.pi/180
    phase = phase *np.pi/180
    
    # t = np.linspace(ti,tf,int(n_points)+1, endpoint=True)      # time array
    t = np.arange(ti, tf, dt)
    t_fraction = t[0:int(fraction*len(t))+1]    # fT
    t_rem_fraction = t[int(fraction*len(t))+1:]
    t_array = t_fraction.copy()      # time array for actual b operation

    bx = amp*(np.cos(phi)*np.cos(2*np.pi*freq*t_fraction) + np.sin(theta)*np.sin(phi)*np.sin(2*np.pi*freq*t_fraction))
    by = amp*(np.sin(phi)*np.cos(2*np.pi*freq*t_fraction) - np.sin(theta)*np.cos(phi)*np.sin(2*np.pi*freq*t_fraction))
    bz = (align_field[-1]/np.abs(align_field[-1])) * amp*np.cos(theta)*np.sin(2*np.pi*freq*t_fraction)

    b1 = np.transpose(np.array([bx, by, bz]))

    # plt.figure(); plt.plot(t_fraction, b1)
    ax = fig.add_subplot(2,4,1)
    ax.plot(t_fraction*1e3, b1,linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b1,linespec)

    # 2nd cycle RCP 1
    # with the same dt and shifted starting time point (t_fraction[-1]), form the time array of length same as that of prev cycle

    ti = t_array[-1] + dt
    temp = ti + t_fraction[0:-1]
    # temp = temp[0:-1]
    t_array = np.hstack((t_array, temp))

    b2 = np.flip(b1, axis=0)    # b1[1:len(t_fraction)+1]
    # print(b2.shape)
    b2 = b2[1:,:]
    # b2 = np.transpose(np.array([np.nan_to_num(bx*bx[-1]/np.absolute(bx[-1]), nan=0.0), np.nan_to_num(by*by[-1]/np.absolute(by[-1]), nan=0.0), -bz]))
    # plt.figure(); plt.plot(b2)
    ax = fig.add_subplot(2,4,2)
    ax.plot(temp*1e3, b2,linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b2,linespec)

    # 3rd cycle RCP 2

    ti = t_array[-1] + dt
    temp = ti + t_fraction[0:-1]
    t_array = np.hstack((t_array, temp))

    b3 = np.transpose(np.array([bx, by, -bz]))
    b3 = b3[1:,:]
    # plt.figure(); plt.plot(b3)
    ax = fig.add_subplot(2,4,3)
    ax.plot(temp*1e3, b3,linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b3,linespec)

    # 4th cycle LCP

    ti = t_array[-1] + dt
    temp = ti + t_fraction[0:-2]
    t_array = np.hstack((t_array, temp))

    # b4 = b2     # this is wrong... It creates a reference to the same array in the memory
    b4 = b2.copy()
    b4[:,-1] = -b2[:,-1]
    b4 = b4[:-1,:]
    # plt.figure(); plt.plot(b4)
    ax = fig.add_subplot(2,4,4)
    ax.plot(temp*1e3, b4, linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot( b4,linespec)

    #
    # tzero = 100 /1e3        # enter value in ms
    # ti = t_array[-1] + dt
    # temp = ti + np.array(np.arange(0,tzero,dt))

    # Bzero = np.zeros((int(samp_rate*tzero), b1.shape[-1]))
    # t_array = np.hstack((t_array, temp))

    # b = np.ascontiguousarray(np.vstack((b1, b2, b3, b4, Bzero)))
    # b_before = np.array([b1[0,:] for _ in range(0,100)])
    # b_after = np.array([b4[-1,:] for _ in range(0,100)])
    b = np.ascontiguousarray(np.vstack((b1, b2, b3, b4)))
    # b = np.ascontiguousarray(np.vstack((b_before, b1, b2, b3, b4, b_after)))
    # t_array = 
    # plt.figure(); plt.plot(b)
    ax = fig.add_subplot(2,4,(5,8))
    ax.plot(t_array*1e3, b, linespec)
    ax.set_ylim(np.floor(-amp), np.ceil(amp))
    # ax.plot(b,linespec)
    ax.legend(['X', 'Y', 'Z'])

    return [dt, t_array, b]

def prepare_data_for_write(data):
    # Ensure data is in the correct shape (samples as rows, channels as columns)
    # print(f"data shape = {data.shape}")
    if data.shape[0] != 3:  # Assuming 3 channels
        data = data.T
    # print(f"data shape = {data.shape}")
    # Create a new array that is guaranteed to be C_CONTIGUOUS and WRITEABLE
    data_copy = np.array(data, dtype=np.float64, order='C', copy=True)
    # print(data_copy)
    # Verify flags
    assert data_copy.flags['C_CONTIGUOUS'], "Array is not C_CONTIGUOUS"
    assert data_copy.flags['WRITEABLE'], "Array is not WRITEABLE"
    # print(f"data_copy shape = {data_copy.shape}")
    return data_copy

def create_retriggerable_ao_task(ao_task, data):
    """Make data compatible to DAQ and create and configure the output task"""
    global samp_rate
    # data = prepare_data_for_write(data)
    try:
        ao_task.stop()
        ao_task.out_stream.regen_mode = RegenerationMode.ALLOW_REGENERATION
        # ao_task.out_stream.offset = 0
        # ao_task.out_stream.relative_to = WriteRelativeTo.FIRST_SAMPLE
        # Configure timing
        ao_task.timing.cfg_samp_clk_timing(
            rate = samp_rate,
            sample_mode = AcquisitionType.FINITE,
            samps_per_chan = 1*np.max(data.shape)
        )
        # print(f"Host Buff size = {ao_task.out_stream.output_buf_size}")
        ao_task.out_stream.output_buf_size = np.max(data.shape)
        # print(f"Host Buff size = {ao_task.out_stream.output_buf_size}")
        # print(f"Onbrd Buff size = {ao_task.out_stream.output_onbrd_buf_size}")
        
        # configure digital start trigger
        aoStartTrig = ao_task.triggers.start_trigger      # get the start trigger configuration for the task
        aoStartTrig.cfg_dig_edge_start_trig(concfg.start_trig_terminal, Edge.RISING)
        aoStartTrig.retriggerable = True
        # print(f"Retriggerable start trigger = {aoStartTrig.retriggerable}")
        
        # configure digital pause trigger
        # ao_pause_trig = ao_task.triggers.pause_trigger
        # ao_pause_trig.trig_type = TriggerType.DIGITAL_LEVEL
        # ao_pause_trig.dig_lvl_src = concfg.start_trig_terminal
        # ao_pause_trig.dig_lvl_when = Level.LOW
        print("> Retriggerable AO configured!!")
    except Exception as excpt:
        print('Error configuring AO.\n Exception details:'+'\x1b[38;2;250;37;41m\n'+str(type(excpt).__name__)+'.'+str(excpt)+'\x1b[0m')
        close_daq_task(ao_task)
        status = False
    else:
        status = True

    # try:
    #     # Write data
    #     ao_task.write(data_copy, auto_start=True)
    #     # ao_task.wait_until_done()
    #     # ao_task.stop()
    # except nidaqmx.errors.DaqError as e:
    #     print(f"DAQmx Error: {e}")
    #     return False
    return status

def set_outputs_to_constant(task, field_data):
    print(f"> Setting to constant: {field_data}")
    # Stop the task
    task.stop()
    
    # Reconfigure timing for on-demand output
    task.timing.cfg_samp_clk_timing(
        1000,  # This rate doesn't matter for on-demand
        sample_mode=AcquisitionType.CONTINUOUS,
        samps_per_chan=1000
    )
    
    # Remove the start trigger configuration
    task.triggers.start_trigger.retriggerable = False
    task.triggers.start_trigger.disable_start_trig()
    
    # Write zeros to all channels
    # num_channels = len(task.ao_channels.channel_names)
    # print(f"Data = {field_data}")
    task.write(field_data, auto_start=True)
    
    # Wait for the write operation to complete
    # task.wait_until_done()
    # task.stop()

def start_retriggerable_ao_task(task, b, dt):
    try:
        # write_start_time = time.perf_counter_ns()
        samples_written = task.write(b, auto_start=True)
        # task.wait_until_done()
        # write_end_time = time.perf_counter_ns()
        # start_time = time.perf_counter_ns()
        # task.start()
        # end_time = time.perf_counter_ns()
        actual_sampling_rate = task.timing.samp_clk_rate

        # print(f"Actual samples written: {samples_written}")
        # print(f" Set sampling rate = {1/dt:g} Sa/s")
        # print(f"Actual sampling rate: {actual_sampling_rate:g} Sa/s")
    except nidaqmx.errors.DaqError as e:
        print(f"\x1b[38;2;250;37;41mDAQmx Error: {e}")
    except KeyboardInterrupt:
        pass

# this goes to main()
def start_rotating_field_task(task, b):
    rot_angle = 20
    freq = 10
    amp = 1
    [dt, t_array, b] = triggered_ao_data(amp, rot_angle, freq)
    b_rot_prepared = prepare_data_for_write(b)
    ao_prepared = create_retriggerable_ao_task(task, b_rot_prepared) if concfg.start_trig_terminal not in task.triggers.start_trigger.dig_edge_src else True
    if ao_prepared:
        start_retriggerable_ao_task(task, b, dt)

# # examples...
# if __name__ == '__main__':
# #%%
#     import time, matplotlib.pyplot as plt
#     t=[]
#     dev = "P6363"
#     ao_task = config_ao(dev)
#     time.sleep(1)
#     for i in range(0,int(1e1)):
#         t1 = time.perf_counter()
#         start_ao(ao_task, [i,i,i])
#         t2 = time.perf_counter()
#         t.append((t2-t1)*1e3)
#     plt.figure(); plt.plot(t)
#     plt.figure(); plt.hist(t)
#     start_ao(ao_task,[0,0,0])
#     close_daq_task(ao_task)
#     # print(coil_calibration([60,60,60]))
#     # %%
#     ao_task = config_ao(dev)
#     #%%
#     start_ao(ao_task, [0,0,0])
#     # %%
