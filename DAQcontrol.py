#DAQ control
#%%
import connectionConfig as concfg, nidaqmx, sys, os, numpy as np, matplotlib.pyplot as plt
from  nidaqmx.constants import VoltageUnits, TerminalConfiguration, AcquisitionType, Edge, Level, TriggerType, RegenerationMode
from typing import Union, List, Tuple, Optional
from abc import ABC, abstractmethod
from typing import Optional, Union, List  # Import at the top of file

# Base abstract class for common functionality
class DAQTask(ABC):
    def __init__(self):
        # self.task_name = task_name
        self.task_state = ""
    
    @abstractmethod
    def configure(self):
        pass
    
    @abstractmethod
    def start(self):
        pass
    
    @abstractmethod
    def stop(self):
        pass

# Specialized classes for each type
class AnalogInputTask(DAQTask):
    def __init__(self,
                 sampling_rate: float = 1000,
                 channels: Optional[List[int]] = [7],
                 voltage_range: tuple = (-10, 10),
                 trigger_source: Optional[str] = None,  # Optional trigger
                 trigger_edge: str = "rising",       # Only used if trigger_source exists
                 trigger_level: Optional[float] = None,  # Optional level for analog trigger
                 sampling_mode: str = "finite",
                 samples_to_acquire: Optional[int] = None
                 ):
        super().__init__()
        self.dev = "P6363"
        self.sampling_rate = sampling_rate
        self.sampling_mode = sampling_mode
        self.samples_to_acquire = samples_to_acquire
        self.channels = channels or []
        self.min_voltage, self.max_voltage = voltage_range
        
        # Store trigger configuration if provided
        self.trigger_config = None
        if trigger_source:
            self.trigger_config = {
                "source": trigger_source,
                "edge": trigger_edge,
                "level": trigger_level
            }
        self._task = None

    def configure(self):
        try:
            print("configuring..")
            # Basic configuration
            self._task = nidaqmx.Task()
            for channel in self.channels:
                self._task.ai_channels.add_ai_voltage_chan(
                    f"{self.dev}/ai{channel}",
                    min_val=self.min_voltage,
                    max_val=self.max_voltage
                )

            self._task.timing.cfg_samp_clk_timing(
                rate=self.sampling_rate,
                sample_mode=DAQConfiguration.get_sampling_mode(mode_name=self.sampling_mode),
                samps_per_chan=self.samples_to_acquire
            )

            # Configure trigger only if specified
            if self.trigger_config:
                self._configure_trigger()

            self.task_state = "configured"

        except nidaqmx.errors.Error as e:
            self.task_state = "error"
            raise RuntimeError(f"Configuration failed: {str(e)}")

    def _configure_trigger(self):
        """Private method to handle trigger configuration"""
        try:
            source = self.trigger_config["source"]
            if source.startswith("ai"):  # Analog trigger
                self._task.triggers.start_trigger.cfg_anlg_edge_start_trig(
                    trigger_source=source,
                    trigger_slope=DAQConfiguration.get_trigger_edge(edge_name=self.trigger_config["edge"]),
                    trigger_level=self.trigger_config["level"] or 0.0
                )
            else:  # Digital trigger
                self._task.triggers.start_trigger.cfg_dig_edge_start_trig(
                    trigger_source=source,
                    trigger_edge=DAQConfiguration.get_trigger_edge(edge_name=self.trigger_config["edge"])
                )
        except nidaqmx.errors.Error as e:
            raise RuntimeError(f"Trigger configuration failed: {str(e)}")

    def start(self):
        try:
            if self.task_state != "configured":
                raise RuntimeError("Task must be configured before starting")
            self._task.start()
            self.task_state = "running"
        except Exception as e:
            self.task_state = "error"
            raise RuntimeError(f"Start failed: {str(e)}")

    def read_daq(self,Nsamples,timeout=120):
        try:
            counts = self._task.read(Nsamples, timeout)
        except Exception as excpt:
            print('Error: could not Read DAQ. Please check your DAQ\'s connections.\nException details:'+'\x1b[38;2;250;37;41m'+str(type(excpt).__name__)+'.'+str(excpt)+'\x1b[0m')
            sys.exit()
        return counts

    def stop(self):
        try:
            if self.task_state == "running":
                self._task.stop()
                self.task_state = "stopped"
        except Exception as e:
            self.task_state = "error"
            raise RuntimeError(f"Stop failed: {str(e)}")
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, '_task') and self._task is not None:
                if self.task_state == "running":
                    self.stop()
                self._task.close()
        except Exception as e:
            print(f"Warning: Error during cleanup: {str(e)}")

class AnalogOutputTask(DAQTask):
    def __init__(self,
                 sampling_rate: float = 1000,
                 channels: Optional[List[int]] = [0,1,2],
                 voltage_range: tuple = (-10, 10),
                 trigger_source: Optional[str] = '',  # Optional trigger
                 trigger_edge: str = "rising",       # Only used if trigger_source exists
                 trigger_level: Optional[float] = None,  # Optional level for analog trigger
                 sampling_mode: str = "continuous",
                 samples_to_generate: Optional[int] = 1000       # this will define the output buffer
                 ):
        super().__init__()
        self.dev = "P6363"
        self.sampling_rate = sampling_rate
        self.sampling_mode = sampling_mode
        self.samples_to_generate = samples_to_generate
        self.channels = channels or []
        self.min_voltage, self.max_voltage = voltage_range
        self.vi_calibration = np.array([[26.5, 26.8, 43.5]])
        # self.vi_calibration = np.array([[19,19,37]])

        # Store trigger configuration if provided
        self.trigger_config = None
        if trigger_source != '':
            self.trigger_config = {
                "source": trigger_source,
                "edge": trigger_edge,
                "level": trigger_level
            }
        self._task = None  # to hold NI-DAQ task object
        
        self.configure()
    
    def configure(self):
        # Example NI-DAQ configuration
        try:
            print("configuring...")
            self._task = nidaqmx.Task()
            for channel in self.channels:
                self._task.ao_channels.add_ao_voltage_chan(
                    f"{self.dev}/ao{channel}",
                    min_val=self.min_voltage,
                    max_val=self.max_voltage
                )

            self._task.timing.cfg_samp_clk_timing(
                rate=self.sampling_rate,
                source='',
                active_edge=DAQConfiguration.get_trigger_edge(edge_name="rising"),
                sample_mode=DAQConfiguration.get_sampling_mode(mode_name=self.sampling_mode),
                samps_per_chan=self.samples_to_generate
            )

            # Configure trigger only if specified
            if self.trigger_config:
                self._configure_trigger()

            self.task_state = "configured"
        
        except nidaqmx.errors.Error as e:
            self.task_state = "error"
            raise RuntimeError(f"Configuration failed: {str(e)}")

    def _configure_trigger(self):
        """Private method to handle trigger configuration"""
        try:
            source = self.trigger_config["source"]
            if source.startswith("ai"):  # Analog trigger
                self._task.triggers.start_trigger.cfg_anlg_edge_start_trig(
                    trigger_source=source,
                    trigger_slope=DAQConfiguration.get_trigger_edge(edge_name=self.trigger_config["edge"]),
                    trigger_level=self.trigger_config["level"] or 0.0
                )
            elif source == '':
                # no trigger source
                pass
            else:  # Digital trigger
                self._task.triggers.start_trigger.cfg_dig_edge_start_trig(
                    trigger_source=source,
                    trigger_edge=DAQConfiguration.get_trigger_edge(edge_name=self.trigger_config["edge"])
                )
        except nidaqmx.errors.Error as e:
            raise RuntimeError(f"Trigger configuration failed: {str(e)}")

    def start(self, data):
        try:
            if self.task_state != "configured":
                raise RuntimeError("Task must be configured before starting")
            # self._task.start()
            self._task.write(data, auto_start=True)
            self.task_state = "running"
        except Exception as e:
            self.task_state = "error"
            raise RuntimeError(f"Start failed: {str(e)}")

    def stop(self):
        try:
            if self.task_state == "running" or self.task_state == "error":
                self._task.stop()
                self.task_state = "stopped"
        except Exception as e:
            self.task_state = "error"
            raise RuntimeError(f"Stop failed: {str(e)}")
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, '_task') and self._task is not None:
                if self.task_state == "running":
                    self.stop()
                if self.task_state != "closed":
                    self._task.close()
        except Exception as e:
            print(f"Warning: Error during cleanup: {str(e)}")

    @staticmethod
    def prepare_data_for_write(data):
        """ Ensure data is in the correct shape (samples as rows, channels as columns) for DAQmx write
        """

        # print(f"data shape before transpose = {data.shape}")
        if data.shape[0] != 3:  # Assuming 3 channels
            data = data.T
        # print(f"data shape after transpose = {data.shape}")

        # Create a new array that is guaranteed to be C_CONTIGUOUS and WRITEABLE
        data_copy = np.array(data, dtype=np.float64, order='C', copy=True)
        
        # Verify flags
        assert data_copy.flags['C_CONTIGUOUS'], "Array is not C_CONTIGUOUS"
        assert data_copy.flags['WRITEABLE'], "Array is not WRITEABLE"

        # print(f"data_copy shape = {data_copy.shape}")
        return data_copy

    def set_outputs_to_constant(self, output_field_in_gauss: Union[Tuple[float, float], List[float]]):
        """Configure for and Set outputs to a constant value
        Configures the task without a trigger and CONTINUOUS sample clock as it is in the default configuration. Then writes the data.
        Supply the output field in [gauss].

        If using ao_task._task.write, use the 'vi_calibration' attribute of AnalogOutputtask() to convert into voltage.
        Use np.array(<field_list>)/vi_calibration
        """
        # Configures the task without a trigger and FINITE sample clock. Then writes the data. -- this was before.
        print(f"> Setting to constant: {output_field_in_gauss} [gauss]")
        output_voltage = AnalogOutputTask.prepare_data_for_write(np.array(output_field_in_gauss)/self.vi_calibration)
        # Stop the task
        if self.task_state == "running":
            self._task.stop()
            self.task_state == "stopped"
        
        # Reconfigure timing for on-demand output
        self._task.timing.cfg_samp_clk_timing(
            1000,  # This rate doesn't matter for on-demand
            sample_mode=AcquisitionType.CONTINUOUS,
            samps_per_chan=1000
        )
        
        # Remove the start trigger configuration
        self._task.triggers.start_trigger.retriggerable = False
        self._task.triggers.start_trigger.disable_start_trig()
        
        # Write zeros to all channels
        # num_channels = len(task.ao_channels.channel_names)
        # print(f"Data = {field_data}")
        self._task.write(output_voltage, auto_start=True)
        
        # Wait for the write operation to complete
        # task.wait_until_done()
        # task.stop()

    def create_retriggerable_ao_task(self, data_shape: tuple):
        """configure the retriggerable output task"""

        self.samples_to_generate = np.max(data_shape)
        self.sampling_rate = 100e3
        # configure the retriggerable task
        try:
            if self.task_state == "running":
                self.stop()

            # set to regeneration mode
            self._task.out_stream.regen_mode = RegenerationMode.ALLOW_REGENERATION
            # ao_task.out_stream.offset = 0
            # ao_task.out_stream.relative_to = WriteRelativeTo.FIRST_SAMPLE

            # Configure sample clock timing
            self._task.timing.cfg_samp_clk_timing(
                rate = self.sampling_rate,
                source = '',
                active_edge = Edge.RISING,
                sample_mode = AcquisitionType.FINITE,
                samps_per_chan = 1*self.samples_to_generate
            )

            # check output buffer size
            # print(f"Host Buff size = {self._task.out_stream.output_buf_size}")
            self._task.out_stream.output_buf_size = self.samples_to_generate
            # print(f"Host Buff size = {self._task.out_stream.output_buf_size}")
            # print(f"Onbrd Buff size = {self._task.out_stream.output_onbrd_buf_size}")
            
            # configure digital start trigger
            aoStartTrig = self._task.triggers.start_trigger      # get the start trigger configuration for the task
            aoStartTrig.cfg_dig_edge_start_trig(concfg.start_trig_terminal, Edge.RISING)
            aoStartTrig.retriggerable = True
            # print(f"Retriggerable start trigger = {aoStartTrig.retriggerable}")
            
            # configure digital pause trigger
            # ao_pause_trig = ao_task.triggers.pause_trigger
            # ao_pause_trig.trig_type = TriggerType.DIGITAL_LEVEL
            # ao_pause_trig.dig_lvl_src = concfg.start_trig_terminal
            # ao_pause_trig.dig_lvl_when = Level.LOW
            # print("> Retriggerable AO configured!!")
            self.task_state = "configured"
        except Exception as e:
            self.task_state = "error"
            raise RuntimeError(f"Retrigger Config failed: {str(e)}")

    def start_retriggerable_ao_task(self, pattern_data):
        try:
            # write_start_time = time.perf_counter_ns()
            samples_written = self._task.write(pattern_data, auto_start=True)
            self.task_state = "running"
            # task.wait_until_done()
            # write_end_time = time.perf_counter_ns()
            # start_time = time.perf_counter_ns()
            # task.start()
            # end_time = time.perf_counter_ns()
            actual_sampling_rate = self._task.timing.samp_clk_rate

            # print(f"Actual samples written: {samples_written}")
            # print(f" Set sampling rate = {1/dt:g} Sa/s")
            # print(f"Actual sampling rate: {actual_sampling_rate:g} Sa/s")
        except nidaqmx.errors.DaqError as e:
            print(f"\x1b[38;2;250;37;41mDAQmx Error: {e}")
        except KeyboardInterrupt:
            pass

    # def voltage_to_current(self, output_field: list):
    #     self.vi_calibration = np.array([[26.5, 26.8, 43.5]])
    #     # calibration = [19,19,37]
    #     output_aovoltage = output_field/self.calibration
    #     return output_aovoltage

class CounterInputTask(DAQTask):
    def __init__(self, edge_config: str):
        super().__init__()
        self.edge_config = edge_config
    
    def configure(self):
        # Specific counter input configuration
        pass

# Configuration mapping
class DAQConfiguration:
    SAMPLING_MODES = {
        "continuous": nidaqmx.constants.AcquisitionType.CONTINUOUS,
        "finite": nidaqmx.constants.AcquisitionType.FINITE,
        "hwtimed": nidaqmx.constants.AcquisitionType.HW_TIMED_SINGLE_POINT
    }

    TRIGGER_EDGE = {
        "rising": nidaqmx.constants.Edge.RISING,
        "falling": nidaqmx.constants.Edge.FALLING
    }

    @classmethod
    def get_sampling_mode(cls, mode_name: str):
        return cls.SAMPLING_MODES.get(mode_name)

    @classmethod
    def get_trigger_edge(cls, edge_name: str):
        return cls.TRIGGER_EDGE.get(edge_name)

# Command pattern using dict
class DAQCommander:
    def __init__(self, task):
        self.task = task
        self.commands = {
            "start": self.task.start,
            "stop": self.task.stop,
            "configure": self.task.configure,
        }

    def execute(self, command: str):
        if command in self.commands:
            self.commands[command]()
        else:
            raise ValueError(f"Unknown command: {command}")

# Factory class to create appropriate task types
class DAQTaskFactory:
    @staticmethod
    def create_task(task_type: str, **kwargs):
        if task_type == "analog_input":
            return AnalogInputTask(**kwargs)
        elif task_type == "analog_output":
            return AnalogOutputTask(**kwargs)
        elif task_type == "counter_input":
            return CounterInputTask(**kwargs)
        else:
            raise ValueError(f"Unknown task type: {task_type}")
    
class DAQ_write_pattern():

    # def __init__(self):
    #     self.pattern_data = np.array([])
    #     self.time_array = np.array([])

    @staticmethod
    def triggered_ao_data(direction: str, rot_angle: float, align_field: list, amp: float, freq: float=10):
        if direction == 'x':
            return DAQ_write_pattern.triggered_ao_data_x(rot_angle, align_field, amp, freq)
        elif direction == 'z':
            return DAQ_write_pattern.triggered_ao_data_z(rot_angle, align_field, amp, freq)

    @staticmethod
    def triggered_ao_data_z(rotation_angle:float, align_field:list, amp:float, freq:float=10, ax:plt.axes=None):
        """Generate the rotating field data pattern for control about z-axis
        
        Parameters
        ----------
            rotation angle : float
                The degree of rotation of the propeller about the coil-z-axis
            align_field : list
                The DC alignment field operating at the beginning of the acquisition
            amp : float
                Amplitude of the applied control field
            freq : float
                Frequency of the applied control field
        """
        # self.frequency = freq
        # self.rotation_angle = rot_angle
        theta = 0       
        phi = 0
        phase = 90
        # amp = 2.63
        ti = 0
        tf = 1/freq     # tf is one time period

        fig = plt.figure()
        title = f'{rotation_angle} deg of {freq} Hz theta_{theta} phi_{phi}'
        fig.suptitle(title)
        fig.canvas.manager.set_window_title(title)
        
        #
        sampling_rate = 100e3     # daq sampling rate
        dt = 1/sampling_rate    # time diff between two samples, sampling time
        # n_points = (tf*freq - ti)/dt     # total points for the waveform

        linespec = '-'
        fraction = rotation_angle/360
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

        time_array = t_array
        pattern_data = AnalogOutputTask.prepare_data_for_write(b)
        return [dt, time_array, pattern_data]

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
        sampling_rate = 100e3     # daq sampling rate
        dt = 1/sampling_rate    # time diff between two samples, sampling time

        linespec = '-'
        fraction = rot_angle/360
        theta = theta *np.pi/180
        phi = phi *np.pi/180
        phase = phase *np.pi/180
        
        t = np.arange(ti, tf, dt)
        t_fraction = t[0:int(fraction*len(t))+1]    # fT
        t_rem_fraction = t[int(fraction*len(t))+1:]
        t_array = t_fraction.copy()      # time array for actual b operation

        bx = amp*(np.cos(phi)*np.cos(2*np.pi*freq*t_fraction) + np.sin(theta)*np.sin(phi)*np.sin(2*np.pi*freq*t_fraction))
        by = amp*(np.sin(phi)*np.cos(2*np.pi*freq*t_fraction) - np.sin(theta)*np.cos(phi)*np.sin(2*np.pi*freq*t_fraction))
        bz = (align_field[-1]/np.abs(align_field[-1])) * amp*np.cos(theta)*np.sin(2*np.pi*freq*t_fraction)

        b1 = np.transpose(np.array([bx, by, bz]))

        ax = fig.add_subplot(2,4,1)
        ax.plot(t_fraction*1e3, b1,linespec)
        ax.set_ylim(np.floor(-amp), np.ceil(amp))
        # ax.plot(b1,linespec)

        # 2nd cycle RCP 1
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
        
        b = np.ascontiguousarray(np.vstack((b1, b2, b3, b4)))
        ax = fig.add_subplot(2,4,(5,8))
        ax.plot(t_array*1e3, b, linespec)
        ax.set_ylim(np.floor(-amp), np.ceil(amp))
        ax.legend(['X', 'Y', 'Z'])

        time_array = t_array
        pattern_data = AnalogOutputTask.prepare_data_for_write(b)
        return [dt, time_array, pattern_data]
    
# def read_daq_counter(task):

# # this goes to main()
# def start_rotating_field_task(task, b):
#     rot_angle = 20
#     freq = 10
#     amp = 1
#     [dt, t_array, b] = triggered_ao_data(amp, rot_angle, freq)
#     b_rot_prepared = prepare_data_for_write(b)
#     ao_prepared = create_retriggerable_ao_task(task, b_rot_prepared) if concfg.start_trig_terminal not in task.triggers.start_trigger.dig_edge_src else True
#     if ao_prepared:
#         start_retriggerable_ao_task(task, b, dt)

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

a, b, c = DAQ_write_pattern.triggered_ao_data_z(60, [50,0,50], 50, 20)