#DAQ control
import nidaqmx, sys, os, numpy as np, logging, time
import matplotlib.pyplot as plt, matplotlib as mpl
from  nidaqmx.constants import VoltageUnits, TerminalConfiguration, AcquisitionType, Edge, Level, TriggerType, RegenerationMode, ProductCategory, CountDirection
from typing import Union, List, Tuple, Optional
from abc import ABC, abstractmethod
from typing import Optional, Union, List
plt.style.use(['dark_background'])
plt.rcParams['axes.prop_cycle'] = mpl.rcParamsOrig['axes.prop_cycle']

def printt(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")

# Base abstract class for common functionality
class DAQTask(ABC):
    def __init__(self, task_name=''):
        self.task_name = task_name
        self._task_state = ""
    
    @abstractmethod
    def configure(self):
        pass
    
    @abstractmethod
    def start(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def stop(self):
        pass

# Specialized classes for each type
class AnalogInputTask(DAQTask):
    def __init__(self,
                 dev = "P6363",
                 channels: list[int]|None = [21],
                 voltage_ranges: list[tuple[float,float]] = [(-1., 1.)],
                 sample_source: str = 'internal',
                 sample_rate: float = 10e3,
                 sample_mode: str = "finite",
                 samps_per_chan: int = 100,
                 start_trigger_source: str = '',  # Optional trigger
                 start_trigger_edge: str = "rising",       # Only used if start_trigger_source exists
                 pause_trigger_source: str = '',
                 name: str = '',
        ):
        super().__init__(name)
        self.dev = dev
        self.channels = channels or []
        self.voltage_ranges = voltage_ranges
        self.sample_clock = {
            'source': '' if sample_source.lower() in ['internal', ''] else sample_source,
            'edge': Edge.RISING,
            'rate': sample_rate,
            'mode': sample_mode,
            'samps_per_chan': samps_per_chan,
        }

        # Store trigger configuration if provided
        self.start_trigger:dict = {}
        if start_trigger_source != '':
            self.start_trigger = {
                "source": start_trigger_source,
                "edge": start_trigger_edge,
                # "level": trigger_level,
            }
        
        self.pause_trigger:dict = {}
        if pause_trigger_source != '':
            self.pause_trigger = {
                "source": pause_trigger_source,
                "type": "digital",
                "level": Level.LOW,
            }

        self._task: nidaqmx.Task
        self.configure()

    def configure(self):
        try:
            printt("Configuring Analog Input..")
            self._task = nidaqmx.Task()
            for idx, channel in enumerate(self.channels):
                self._task.ai_channels.add_ai_voltage_chan(
                    physical_channel=f"{self.dev}/ai{channel}",                    
                    terminal_config=TerminalConfiguration.RSE,
                    min_val=self.voltage_ranges[idx][0],
                    max_val=self.voltage_ranges[idx][1],
                )
            self._task.timing.cfg_samp_clk_timing(
                rate = self.sample_clock['rate'],
                source = self.sample_clock['source'],
                active_edge = self.sample_clock['edge'],
                sample_mode = DAQProperty.get_sampling_mode(mode=self.sample_clock['mode']),
                samps_per_chan = int(self.sample_clock['samps_per_chan']),
            )
            self._configure_triggers()
            self._task_state = "configured"
            printt("✔ Analog Input configured!")
            printt(f"Input      : {self.channels}, range: {self.voltage_ranges}")
            printt(f"Sampling   : {self.sample_clock}")
            printt(f"Start      : {self.start_trigger}")
            printt(f"Pause      : {self.pause_trigger}")

        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Configuration failed: {str(e)}")
            raise #RuntimeError(f"❌ Configuration failed: {str(e)}")

    def _configure_triggers(self):
        """Private method to handle trigger configuration"""
        try:
            # Digital START trigger
            if self.start_trigger:
                self._task.triggers.start_trigger.cfg_dig_edge_start_trig(
                    trigger_source=self.start_trigger["source"],
                    trigger_edge=DAQProperty.get_trigger_edge(edge=self.start_trigger["edge"])
                )
            
            # Pause trigger
            if self.pause_trigger:
                self._task.triggers.pause_trigger.trig_type = TriggerType.DIGITAL_LEVEL
                self._task.triggers.pause_trigger.dig_lvl_src = self.pause_trigger["source"]
                self._task.triggers.pause_trigger.dig_lvl_when = self.pause_trigger["level"]

        except Exception as e:
            logging.exception(f"❌ Trigger configuration failed: {str(e)}")
            raise #RuntimeError(f"Trigger configuration failed: {str(e)}")

    def start(self):
        try:
            if self._task_state != "configured":
                raise RuntimeError("Task must be configured before starting")
            self._task.start()
            self._task_state = "running"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Start failed: {str(e)}")
            raise #RuntimeError(f"❌ Start failed: {str(e)}")

    def read_daq(self,Nsamples,timeout=120.):
        try:
            counts = self._task.read(Nsamples, timeout)
        except Exception as excpt:
            logging.exception(f'❌ Error: could not Read DAQ. Please check your DAQ\'s connections.\
                  \nException details: \x1b[38;2;250;37;41m{str(type(excpt).__name__)}. {str(excpt)}\x1b[0m')
            sys.exit()
        return counts

    def stop(self):
        try:
            if self._task_state == "running":
                self._task.stop()
                self._task_state = "stopped"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Stop failed: {str(e)}")
            raise #RuntimeError(f"❌ Stop failed: {str(e)}")
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, '_task') and self._task is not None:
                if self._task_state == "running":
                    self.stop()
                self._task.close()
        except Exception as e:
            logging.exception(f"⚠ Warning: Error during cleanup: {str(e)}")

class  AnalogOutputTask(DAQTask):
    def __init__(self,
                 dev = "U9263",
                 channels: List[int] = [0,1,2],
                 voltage_range: tuple = (-10, 10),
                 sample_rate: float = 1000,
                 sample_mode: str = "continuous",
                 samps_per_chan: int = 1000,       # this will define the output buffer
                 coil: str = 'default',
                 start_trigger_source: Optional[str] = '',  # Optional trigger
                 start_trigger_edge: str = "rising",       # Only used if start_trigger_source exists
                 trigger_level: Optional[float] = None,  # Optional level for analog trigger
        ):
        super().__init__()
        self.dev = dev
        self.channels = channels or []
        self.min_voltage, self.max_voltage = voltage_range
        self.sample_rate = sample_rate
        self.sample_mode = sample_mode
        self.samps_per_chan = samps_per_chan
        self.coil = coil
        self.params:dict = {}
        
        # Store trigger configuration if provided
        self.trigger_config:dict = {}
        if start_trigger_source != '':
            self.trigger_config = {
                "source": start_trigger_source,
                "edge": start_trigger_edge,
                "level": trigger_level
            }
        self._task: nidaqmx.Task  # to hold NI-DAQ task object
        
        self.configure()
    
    def configure(self):
        # Example NI-DAQ configuration
        try:
            printt("Configuring Analog Output...")
            self._task = nidaqmx.Task()
            for channel in self.channels:
                self._task.ao_channels.add_ao_voltage_chan(
                    f"{self.dev}/ao{channel}",
                    min_val=self.min_voltage,
                    max_val=self.max_voltage,
                )

            self._task.timing.cfg_samp_clk_timing(
                rate=self.sample_rate,
                source='',
                active_edge=DAQProperty.get_trigger_edge(edge="rising"),
                sample_mode=DAQProperty.get_sampling_mode(mode=self.sample_mode),
                samps_per_chan=self.samps_per_chan
            )

            # Configure trigger only if specified
            if self.trigger_config:
                self._configure_trigger()

            self._task_state = "configured"
            printt(f"✔ Analog Output configured: channels = {self.channels}!")

        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Configuration failed: {str(e)}")
            raise #RuntimeError(f"❌ Configuration failed: {str(e)}")

    def _configure_trigger(self):
        """Private method to handle trigger configuration"""
        try:
            source = self.trigger_config["source"]
            # if source.startswith("ai"):  # Analog trigger
            #     self._task.triggers.start_trigger.cfg_anlg_edge_start_trig(
            #         start_trigger_source=source,
            #         trigger_slope=DAQProperty.get_trigger_edge(edge=self.trigger_config["edge"]),
            #         trigger_level=self.trigger_config["level"] or 0.0
            #     )
            if source == '':
                # no trigger source
                pass
            else:  # Digital trigger
                self._task.triggers.start_trigger.cfg_dig_edge_start_trig(
                    trigger_source=source,
                    trigger_edge=DAQProperty.get_trigger_edge(edge=self.trigger_config["edge"])
                )
        except Exception as e:
            logging.exception(f"❌ Trigger configuration failed: {str(e)}")
            raise #RuntimeError(f"Trigger configuration failed: {str(e)}")

    def start(self, data: Union[np.ndarray, List[float], Tuple[float, ...]]):
        try:
            if self._task_state != "configured":
                raise RuntimeError("Task must be configured before starting")
            # self._task.start()
            self._task.write(data, auto_start=True)
            self._task_state = "running"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Start failed: {str(e)}")
            raise #RuntimeError(f"❌ Start failed: {str(e)}")

    def stop(self):
        try:
            if self._task_state == "running" or self._task_state == "error":
                self._task.stop()
                self._task_state = "stopped"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Stop failed: {str(e)}")
            raise #RuntimeError(f"❌ Stop failed: {str(e)}")
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, '_task') and self._task is not None:
                if self._task_state == "running":
                    self.stop()
                if self._task_state != "closed":
                    self._task.close()
        except Exception as e:
            logging.exception(f"⚠ Warning: Error during cleanup: {str(e)}")

    @staticmethod
    def coil_calibration(coil):
        """ Calibration values for converting DAQ voltage [V] to field [gauss] for the coils
        Used in set_outputs_to_constant() method
        """
        if coil == 'hylum':
            vi_calibration = [26.5, 26.8, 43.5]         # hylum coil
        elif coil == 'confocal':
            vi_calibration = [100., 100., 64.]          # confocal coil
        elif coil == 'small_confocal':
            vi_calibration = [100., 100., 0.7]
        elif coil == 'propulsion':
            vi_calibration = [19., 19., 37.]            # propulsion coil
        else:
            vi_calibration = [1.]*3                     # no coil
        
        return np.array(vi_calibration)[:, np.newaxis]

    @staticmethod
    def prepare_data_for_write(data):
        """ Ensure data is in the correct shape (samples as rows, channels as columns) for DAQmx write
        """

        # printt(f"data shape before transpose = {data.shape}")
        if data.shape[0] != 3:  # Assuming 3 channels
            data = data.T
        # printt(f"data shape after transpose = {data.shape}")

        # Create a new array that is guaranteed to be C_CONTIGUOUS and WRITEABLE
        data_copy = np.array(data, dtype=np.float64, order='C', copy=True)
        
        # Verify flags
        assert data_copy.flags['C_CONTIGUOUS'], "Array is not C_CONTIGUOUS"
        assert data_copy.flags['WRITEABLE'], "Array is not WRITEABLE"

        # printt(f"data_copy shape = {data_copy.shape}")
        return data_copy

    def set_outputs_to_constant(self, output_field_in_gauss: Union[Tuple[float], List[float]]):
        """Configure for and Set outputs to a constant value.

        Configures the task without a trigger and CONTINUOUS sample clock as it is in the default configuration. Then writes the data.
        Supply the output field in [gauss].

        If using ao_task._task.write, use the 'vi_calibration' attribute of AnalogOutputtask() to convert into voltage.
        Use np.array(<field_list>)/vi_calibration
        """
        # the convention for the input data structure (typically a 2D NumPy array) is [Channels, Samples]

        # Configures the task without a trigger and FINITE sample clock. Then writes the data. -- this was before.
        printt(f"> Setting to constant: {output_field_in_gauss} [gauss]")

        output_field_in_gauss = [[i]*2 for i in output_field_in_gauss]  # create 2 samples for each channel
        # printt(np.array(output_field_in_gauss))
        # printt(AnalogOutputTask.coil_calibration(self.coil).T)
        output_voltage = AnalogOutputTask.prepare_data_for_write(np.array(output_field_in_gauss))
        output_voltage = np.divide(np.array(output_field_in_gauss), AnalogOutputTask.coil_calibration(self.coil))
        # printt(output_voltage)

        # Stop the task
        printt(f"task state = {self._task_state}")
        if self._task_state == "running":
            self._task.stop()
            self._task_state = "stopped"
        
        # Reconfigure timing for on-demand output
        self._task.timing.cfg_samp_clk_timing(
            1,  # This rate doesn't matter for on-demand
            sample_mode=AcquisitionType.CONTINUOUS,
            samps_per_chan=2
        )
        
        # Remove the start trigger configuration
        if "6363" in self._task.devices[-1].product_type or self._task.devices[-1].product_category == ProductCategory.X_SERIES_DAQ:
            self._task.triggers.start_trigger.retriggerable = False
            self._task.triggers.start_trigger.disable_start_trig()
        
        # Write zeros to all channels
        # num_channels = len(task.ao_channels.channel_names)
        # printt(f"Data = {field_data}")
        try:
            self._task.write(output_voltage, auto_start=True)
        except Exception as e:
            logging.exception(f"❌ Error writing constant output: {str(e)}")
            # raise
        finally:
            self._task.stop()
            pass
        
        # Wait for the write operation to complete
        # task.wait_until_done()
        # task.stop()
        return output_voltage

    def create_retriggerable_ao_task(self, data_shape: tuple):
        """configure the retriggerable output task"""

        self.samps_per_chan = np.max(data_shape)
        self.sample_rate = 100e3
        # configure the retriggerable task
        try:
            if self._task_state == "running":
                self.stop()

            # set to regeneration mode
            self._task.out_stream.regen_mode = RegenerationMode.ALLOW_REGENERATION
            # ao_task.out_stream.offset = 0
            # ao_task.out_stream.relative_to = WriteRelativeTo.FIRST_SAMPLE

            # Configure sample clock timing
            self._task.timing.cfg_samp_clk_timing(
                rate = self.sample_rate,
                source = '',
                active_edge = Edge.RISING,
                sample_mode = AcquisitionType.FINITE,
                samps_per_chan = 1*self.samps_per_chan
            )

            # check output buffer size
            # printt(f"Host Buff size = {self._task.out_stream.output_buf_size}")
            self._task.out_stream.output_buf_size = self.samps_per_chan
            # printt(f"Host Buff size = {self._task.out_stream.output_buf_size}")
            # printt(f"Onbrd Buff size = {self._task.out_stream.output_onbrd_buf_size}")
            
            # configure digital start trigger
            aoStartTrig = self._task.triggers.start_trigger      # get the start trigger configuration for the task
            aoStartTrig.cfg_dig_edge_start_trig(concfg.start_trig_terminal, Edge.RISING)
            aoStartTrig.retriggerable = True
            # printt(f"Retriggerable start trigger = {aoStartTrig.retriggerable}")
            
            # configure digital pause trigger
            # ao_pause_trig = ao_task.triggers.pause_trigger
            # ao_pause_trig.trig_type = TriggerType.DIGITAL_LEVEL
            # ao_pause_trig.dig_lvl_src = concfg.start_trig_terminal
            # ao_pause_trig.dig_lvl_when = Level.LOW
            # printt("> Retriggerable AO configured!!")
            self._task_state = "configured"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Retrigger Config failed: {str(e)}")
            raise #RuntimeError(f"❌ Retrigger Config failed: {str(e)}")

    def start_retriggerable_ao_task(self, pattern_data):
        try:
            # write_start_time = time.perf_counter_ns()
            samples_written = self._task.write(pattern_data, auto_start=True)
            self._task_state = "running"
            # task.wait_until_done()
            # write_end_time = time.perf_counter_ns()
            # start_time = time.perf_counter_ns()
            # task.start()
            # end_time = time.perf_counter_ns()
            actual_sampling_rate = self._task.timing.samp_clk_rate

            # printt(f"Actual samples written: {samples_written}")
            # printt(f" Set sample_clock rate = {1/dt:g} Sa/s")
            # printt(f"Actual sample_clock rate: {actual_sampling_rate:g} Sa/s")
        except Exception as e:
            logging.exception(f"❌ \x1b[38;2;250;37;41mDAQmx Error: {e}")
        except KeyboardInterrupt:
            pass

class CounterInputTask(DAQTask):
    def __init__(self,
                 dev: str = 'P6363',
                 counter: str = "P6363/ctr0",
                 channel: str = "/P6363/PFI7",
                 sample_source: str = 'internal',
                 sample_rate: float = 10e3,
                 sample_mode: str = 'finite',
                 samps_per_chan: int = 100,
                 start_trigger_source: str = '',
                 start_trigger_edge: str = "rising",
                 pause_trigger_source: str = '',
                 name: str = '',
        ):
        if name == '': name = f"task_{time.strftime('%H%M%S')}"
        super().__init__(name)
        self.dev: str = dev
        self.counter: str = counter
        self.channel: str = channel
        self.sample_clock: dict = {
            'source': '' if sample_source.lower() in ['internal', ''] else sample_source,
            'edge': Edge.RISING.name,
            'rate': sample_rate,
            'mode': sample_mode,
            'samps_per_chan': samps_per_chan,
            }
        # Trigger configuration
        self.start_trigger: dict = {}
        if start_trigger_source != '':
            self.start_trigger = {
                "type": "digital",
                "source": start_trigger_source,
                "edge": start_trigger_edge,
                }
        self.pause_trigger: dict = {}
        if pause_trigger_source != '':
            self.pause_trigger = {
                "type": "digital",
                "source": pause_trigger_source,
                "level": Level.LOW,
                }

        self._task: nidaqmx.Task
        self.configure()

    def configure(self):
        try:
            printt("Configuring Counter Input..")
            self._task = nidaqmx.Task()
            # Assign channel
            channel = self._task.ci_channels.add_ci_count_edges_chan(
                counter=self.counter, edge=Edge.RISING,
                initial_count=0, count_direction=CountDirection.COUNT_UP
            )
            channel.ci_count_edges_term = self.channel
            # Configure sample clock
            self._task.timing.cfg_samp_clk_timing(
                rate = self.sample_clock['rate'],
                source = self.sample_clock['source'],
                active_edge = DAQProperty.get_trigger_edge(edge=self.sample_clock['edge']),
                sample_mode = DAQProperty.get_sampling_mode(mode=self.sample_clock['mode']),
                samps_per_chan = int(self.sample_clock['samps_per_chan'])
            )
            self._configure_triggers()
            
            self._task_state = "configured"
            self.get_attributes()
            printt(f"✔ Counter Input configured: counter = {self.counter}, channel = {self.channel}!")
            printt(f"Sampling config: {self.sample_clock}")
            printt(f"Start config: {self.start_trigger}")
            printt(f"Pause config: {self.pause_trigger}")
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Configuration failed: {str(e)}")
            raise #RuntimeError(f"❌ Configuration failed: {str(e)}")

    def _configure_triggers(self):
        """Private method to handle trigger configuration"""
        try:
            # Arm start trigger
            if self.start_trigger:
                self._task.triggers.arm_start_trigger.trig_type = TriggerType.DIGITAL_EDGE
                self._task.triggers.arm_start_trigger.dig_edge_src = self.start_trigger["source"]
                self._task.triggers.arm_start_trigger.dig_edge_edge = DAQProperty.get_trigger_edge(edge=self.start_trigger["edge"])
            
            # Pause trigger
            if self.pause_trigger:
                self._task.triggers.pause_trigger.trig_type = TriggerType.DIGITAL_LEVEL
                self._task.triggers.pause_trigger.dig_lvl_src = self.pause_trigger["source"]
                self._task.triggers.pause_trigger.dig_lvl_when = self.pause_trigger["level"]

        except Exception as e:
            logging.exception(f"❌ Trigger configuration failed: {str(e)}")
            raise #RuntimeError(f"Trigger configuration failed: {str(e)}")
    
    def get_attributes(self):
        self.counter = self._task.ci_channels[0].name                   # type: ignore
        self.channel = self._task.ci_channels[0].ci_count_edges_term    # type: ignore
        if self.sample_clock:
            self.sample_clock['source']         = self._task.timing.samp_clk_src
            self.sample_clock['edge']           = self._task.timing.samp_clk_active_edge.name
            self.sample_clock['rate']           = self._task.timing.samp_clk_rate
            self.sample_clock['mode']           = self._task.timing.samp_quant_samp_mode.name
            self.sample_clock['samps_per_chan'] = self._task.timing.samp_quant_samp_per_chan
        if self.start_trigger:
            self.start_trigger['type']      = self._task.triggers.arm_start_trigger.trig_type.name
            self.start_trigger['source']    = self._task.triggers.arm_start_trigger.dig_edge_src
            self.start_trigger['edge']      = self._task.triggers.arm_start_trigger.dig_edge_edge.name
        if self.pause_trigger:
            self.pause_trigger['type']      = self._task.triggers.pause_trigger.trig_type.name
            self.pause_trigger['source']    = self._task.triggers.pause_trigger.dig_lvl_src
            self.pause_trigger['level']     = self._task.triggers.pause_trigger.dig_lvl_when.name

    def start(self):
        try:
            if self._task_state != "configured":
                raise RuntimeError("Task must be configured before starting")
            self._task.start()
            self._task_state = "running"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Start failed: {str(e)}")
            raise #RuntimeError(f"❌ Start failed: {str(e)}")
    
    # TODO: CIChannel.ci_count: Indicates the current value of the count register: Need to reset the counter after 2**32
    def read_daq(self,Nsamples,timeout=120.):
        try:
            counts = self._task.read(Nsamples, timeout)
            counts = np.diff(counts, axis=0, prepend=0)
        except Exception as excpt:
            logging.exception(f'❌ Error: could not Read DAQ. Please check your DAQ\'s connections.\
                  \nException details: \x1b[38;2;250;37;41m{str(type(excpt).__name__)}. {str(excpt)}\x1b[0m')
            sys.exit()
        return counts

    def stop(self):
        try:
            if self._task_state == "running":
                self._task.stop()
                self._task_state = "stopped"
        except Exception as e:
            self._task_state = "error"
            logging.exception(f"❌ Stop failed: {str(e)}")
            raise #RuntimeError(f"❌ Stop failed: {str(e)}")
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, '_task') and self._task is not None:
                if self._task_state == "running":
                    self.stop()
                self._task.close()
        except Exception as e:
            logging.exception(f"⚠ Warning: Error during cleanup: {str(e)}")


# Configuration mapping
class DAQProperty:
    SAMPLING_MODES = {
        "continuous": AcquisitionType.CONTINUOUS,
        "finite": AcquisitionType.FINITE,
        "hwtimed": AcquisitionType.HW_TIMED_SINGLE_POINT,
    }
    TRIGGER_EDGE = {
        "rising": Edge.RISING,
        "falling": Edge.FALLING,
    }

    @classmethod
    def get_sampling_mode(cls, mode: str) -> AcquisitionType:
        return cls.SAMPLING_MODES.get(mode.lower(), AcquisitionType.FINITE)

    @classmethod
    def get_trigger_edge(cls, edge: str) -> Edge:
        return cls.TRIGGER_EDGE.get(edge.lower(), Edge.RISING)


class DAQ_write_pattern():

    # def __init__(self):
    #     self.pattern_data = np.array([])
    #     self.time_array = np.array([])

    @staticmethod
    def triggered_ao_data_ac(direction: str, rot_angle: float, align_field: list, amp: float, freq: float=10):
        if direction == 'x':
            return DAQ_write_pattern.triggered_ao_data_x(rot_angle, align_field, amp, freq)
        elif direction == 'z':
            return DAQ_write_pattern.triggered_ao_data_z(rot_angle, align_field, amp, freq)

    @staticmethod
    def triggered_ao_data_z(rotation_angle:float, align_field:list, amp:float, freq:float=10, ax=None):
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
        fig.canvas.manager.set_window_title(title)      # type: ignore
        
        #
        sample_rate = 100e3     # daq sample_clock rate
        dt = 1/sample_rate    # time diff between two samples, sample_clock time
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
        # printt(b2.shape)
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

        parameters = {
            'manip_field': {
                'dir': 'z',
                'rotation_angle': rotation_angle,
                'freq': freq,
                'amp': amp,
                'theta': theta,
                'phi': phi,
            }
        }
        return [dt, time_array, pattern_data, parameters]

    @staticmethod
    def triggered_ao_data_x(rot_angle:float, align_field:list, amp:float, freq:float=10.):
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
        sample_rate = 100e3     # daq sample_clock rate
        dt = 1/sample_rate    # time diff between two samples, sample_clock time

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
        # printt(b2.shape)
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

        parameters = {
            'manip_field': {
                'dir': 'z',
                'rotation_angle': rot_angle,
                'freq': freq,
                'amp': amp,
                'theta': theta,
                'phi': phi,
            }
        }
        
        return [dt, time_array, pattern_data]
    
    @staticmethod
    def triggered_ao_data_dc(t_align:float, rotation_theta:float, rotation_phi:float, amp:float, ax=None):
        """For demonstrating vector control of the propeller -> Bz - Bx - Bz - Bx -> rotation angle -> Test field
        The test field is written in mainControl.py

        Parameters
        ----------
            t_align : float
                Time duration for each Bz/Bx alignment field
            rotation_theta : float
                The degree of rotation of the propeller about the coil-z-axis
            align_field : list
                The DC alignment field operating at the beginning of the acquisition
            amp : float
                Amplitude of the applied control field
        """
        N_aligns = 3        # number of Bz - Bx rotations

        sample_rate = 100e3     # daq sample_clock rate [Hz]
        dt = 1/sample_rate    # time diff between two samples, sample_clock time [s]
        t = np.arange(0, t_align/2, dt)       # incoming t_align in [s]
        
        # TODO: #5 This is either (Bz-Bx) or (Bz-By) alignment field.. Need (Bz - (aBx + bBy)) alignment field
        if rotation_phi == 0:
            # (Bx - Bz) field alignment
            bx = amp*np.ones_like(t)
            by = np.zeros_like(t)
            bz = np.zeros_like(t)
            b1 = np.transpose(np.array([bx, by, bz]))

            bx = np.zeros_like(t)
            by = np.zeros_like(t)
            bz = amp*np.ones_like(t)
            b2 = np.transpose(np.array([bx, by, bz]))
        
        elif rotation_phi == 90:
            # (By - Bz) field alignment
            bx = np.zeros_like(t)
            by = amp*np.ones_like(t)
            bz = np.zeros_like(t)
            b1 = np.transpose(np.array([bx, by, bz]))

            bx = np.zeros_like(t)
            by = np.zeros_like(t)
            bz = amp*np.ones_like(t)
            b2 = np.transpose(np.array([bx, by, bz]))
        
        else:
            raise ValueError("Invalid phi angle. Only 0 and 90 degrees are allowed.")
        
        # t_array = np.hstack((t, t)*N_aligns)
        b = np.ascontiguousarray(np.vstack((b2, b1)*N_aligns))

        # rotation about the coil-z-axis
        # t = np.arange(0, 0.05, dt)     # 50 ms rotation time only
        theta = rotation_theta *np.pi/180
        phi = rotation_phi *np.pi/180
        bx = amp*np.cos(phi)*np.sin(theta)*np.ones_like(t)
        by = amp*np.sin(phi)*np.sin(theta)*np.ones_like(t)
        bz = amp*np.cos(theta)*np.ones_like(t)
        b3 = np.transpose(np.array([bx, by, bz]))

        # t_array = np.hstack((t_array, t))
        b = np.ascontiguousarray(np.vstack((b, b3)))
        t_array = np.arange(0, len(b)*dt, dt)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(t_array*1e3, b)
        # ax.set_ylim(np.floor(-amp), np.ceil(amp))
        ax.legend(['X', 'Y', 'Z'])

        pattern_data = AnalogOutputTask.prepare_data_for_write(b)
        return [dt, t_array, pattern_data]

    @staticmethod
    def check_pattern_output(pattern_data, coil: str='default'):
        ao_task = AnalogOutputTask()
        pattern_data = pattern_data/AnalogOutputTask.coil_calibration(coil).T
        AnalogOutputTask.prepare_data_for_write(pattern_data)
        ao_task.create_retriggerable_ao_task(pattern_data.shape)
        printt("Starting pattern!")
        ao_task.start_retriggerable_ao_task(pattern_data)

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
# if __name__ == '__main__':
    # a, b, c = DAQ_write_pattern.triggered_ao_data_z(60, [50,0,50], 50, 20)
# TODO: #4 Need an simple way to test any output pattern.. Write a function that takes the pattern and writes it to the DAQ