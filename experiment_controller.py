"""
Experiment Controller Classes for NV Experiment.

This module provides a unified class structure for both diode/APD-based
and camera-based experiments, enabling future unification while maintaining
clean separation of concerns.

Key Classes:
- ExperimentController (ABC): Base class for all experiment types
- DiodeExperimentController: For voltage/count-based experiments (mainControl_diode.py)
- CameraExperimentController: For image-based experiments (mainControl_camera.py)

Date: December 2024
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
import numpy as np
import time
import logging
from pathlib import Path

from parameter_system import ExperimentConfig, ParameterSweep
from save_manager import SaveManager
from SGcontrol import SignalGenerator, SignalGenerator_sim
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogInputTask, AnalogOutputTask
import connectionConfig as concfg


# ============================================================================
# BASE EXPERIMENT CONTROLLER
# ============================================================================

class ExperimentController(ABC):
    """
    Abstract base class for experiment control.

    Provides common interface for both diode and camera-based experiments.

    Benefits:
    - Common interface for both acquisition modes
    - Shared functionality in base class
    - Future unification straightforward
    - Clean separation of mode-specific logic
    """

    def __init__(self, params_dict: dict, acquisition_mode: str):
        """
        Initialize experiment controller.

        Args:
            params_dict: Parameter dictionary from config file
            acquisition_mode: 'diode' or 'camera'
        """
        self.params_dict = params_dict
        self.acquisition_mode = acquisition_mode
        self.instruments: Dict[str, Any] = {}
        self.save_manager = SaveManager()

        # Extract common parameters
        self.sequence = params_dict.get('seq', {}).get('sequence', 'unknown')
        self.num_runs = params_dict.get('scan', {}).get('Nruns', 1)

    @abstractmethod
    def initialize_instruments(self, trial_run: List[str]):
        """
        Initialize instruments.

        Must be implemented by subclasses to initialize mode-specific instruments.

        Args:
            trial_run: List indicating simulation mode [sg_trial, pb_trial, ...]
        """
        pass

    @abstractmethod
    def _acquire_single_point(self, param_values: dict) -> dict:
        """
        Acquire data at one parameter point.

        Must be implemented by subclasses for mode-specific acquisition.

        Args:
            param_values: Dict of current parameter values (e.g., {'frequency': 2.87e9})

        Returns:
            Dict with measurement results (must include all param_values)
        """
        pass

    def cleanup(self):
        """Close all instruments."""
        if 'sg' in self.instruments and self.instruments['sg'] is not None:
            try:
                self.instruments['sg'].uninit_sg()
                print("✔ SG closed")
            except Exception as e:
                logging.exception(f"⚠ Error closing SG: {e}")

        if 'pb' in self.instruments and self.instruments['pb'] is not None:
            try:
                self.instruments['pb'].closePB()
                print("✔ PB closed")
            except Exception as e:
                logging.exception(f"⚠ Error closing PB: {e}")

        if 'ai_task' in self.instruments and self.instruments['ai_task'] is not None:
            try:
                self.instruments['ai_task'].stop()
                print("✔ AI task stopped")
            except Exception as e:
                logging.exception(f"⚠ Error closing AI task: {e}")

        if 'ao_task' in self.instruments and self.instruments['ao_task'] is not None:
            try:
                self.instruments['ao_task'].stop()
                print("✔ AO task stopped")
            except Exception as e:
                logging.exception(f"⚠ Error closing AO task: {e}")

    def get_instruments_dict(self) -> Dict[str, Any]:
        """Get dictionary of initialized instruments."""
        return self.instruments


# ============================================================================
# DIODE EXPERIMENT CONTROLLER
# ============================================================================

class DiodeExperimentController(ExperimentController):
    """
    Controller for voltage/count-based experiments (APD/diode).

    Replaces procedural code in mainControl_diode.py with class-based approach.

    Instruments:
    - SignalGenerator (or simulation)
    - PulseBlaster
    - AnalogInputTask (DAQ AI)
    - AnalogOutputTask (DAQ AO) - optional for rotating fields
    """

    def __init__(self, params_dict: dict):
        """
        Initialize diode experiment controller.

        Args:
            params_dict: Parameter dictionary from config file
        """
        super().__init__(params_dict, 'diode')

    def initialize_instruments(self, trial_run: List[str] = ['n', 'n']):
        """
        Initialize instruments for diode-based experiments.

        Args:
            trial_run: [sg_mode, pb_mode] where 'n'=real, 'y'=simulation
                      e.g., ['n', 'n'] = both real, ['y', 'n'] = simulated SG
        """
        print("\n" + "="*70)
        print(" "*15 + "Initializing Diode Experiment Instruments")
        print("="*70)

        # Initialize Signal Generator
        try:
            if trial_run[0] == 'n' and self.sequence not in ['aom_timing', 'T1ms0_train']:
                self.instruments['sg'] = SignalGenerator(name="sg1")
                self.instruments['sg'].enable_sg_output()
                print("✔ SG Output Enabled")

                # Set MW power and frequency from params
                mw_params = self.params_dict.get('mw', {})
                if 'power' in mw_params:
                    self.instruments['sg'].set_sg_amp(mw_params['power'])
                if 'freq' in mw_params:
                    freq = mw_params['freq']
                    if isinstance(freq, np.ndarray):
                        freq = freq[0]  # Use first value if array
                    self.instruments['sg'].set_sg_freq(freq)

                self.instruments['sg'].setup_sg_pulse_mod()
                print("✔ SG Ext Pulse Mod Enabled")
            else:
                self.instruments['sg'] = SignalGenerator_sim(name="sg1")
                print("✔ Simulated SG initialized")
        except Exception as e:
            logging.exception(f"❌ Error initializing SG: {e}")
            self.instruments['sg'] = SignalGenerator_sim(name="sg1")

        # Initialize Pulse Blaster
        try:
            self.instruments['pb'] = PulseBlaster(self.params_dict, name="pb1")
            self.instruments['pb'].configure()
            print("✔ PulseBlaster configured")
        except Exception as e:
            logging.exception(f"❌ Error initializing PulseBlaster: {e}")
            raise

        # Initialize Analog Input (DAQ)
        try:
            daq_params = self.params_dict.get('daq', {})
            self.instruments['ai_task'] = AnalogInputTask(
                dev=concfg.daq_dev,
                channels=concfg.input_terminals,
                voltage_range=(-10, 10),
                sampling_rate=daq_params.get('sampling_rate', concfg.daq_max_samp_rate),
                trigger_source=concfg.samp_clk_terminal,
                name="ai_task"
            )
            print("✔ Analog Input configured")
        except Exception as e:
            logging.exception(f"❌ Error initializing AI task: {e}")
            raise

        # Initialize Analog Output (optional, for rotating fields)
        if self.params_dict.get('daq', {}).get('use_ao', False):
            try:
                self.instruments['ao_task'] = AnalogOutputTask(
                    dev=concfg.daq_dev,
                    channels=concfg.output_terminals,
                    sampling_rate=1000,
                    name="ao_task"
                )
                print("✔ Analog Output configured")
            except Exception as e:
                logging.exception(f"❌ Error initializing AO task: {e}")
                self.instruments['ao_task'] = None
        else:
            self.instruments['ao_task'] = None

        print("="*70)
        print("✔ All instruments initialized")
        print("="*70 + "\n")

    def _acquire_single_point(self, param_values: dict) -> dict:
        """
        Acquire data at one parameter point for diode-based experiments.

        This method would contain the acquisition logic from mainControl_diode.py.
        For now, it's a placeholder that needs to be filled with actual logic.

        Args:
            param_values: Dict of current parameter values

        Returns:
            Dict with measurement results
        """
        # TODO: Implement actual acquisition logic from mainControl_diode.py
        # This would include:
        # 1. Set parameters (already done by ParameterSweep if using param system)
        # 2. Program PulseBlaster sequence
        # 3. Start AI task
        # 4. Trigger sequence
        # 5. Read data
        # 6. Process data
        # 7. Return results

        # Placeholder implementation
        pb = self.instruments['pb']
        ai_task = self.instruments['ai_task']

        # Start sequence
        pb.start_sequence()

        # Acquire data (simplified - needs actual logic)
        daq_params = self.params_dict.get('daq', {})
        Nsamples = daq_params.get('Nsamples', 1000)
        counts = ai_task.read_daq(Nsamples)

        # Process and return
        return {
            **param_values,
            # 'signal': 0.0,  # Placeholder
            # 'contrast': 0.0,  # Placeholder
            'data': counts,
            'timestamp': time.time()
        }


# ============================================================================
# CAMERA EXPERIMENT CONTROLLER
# ============================================================================

class CameraExperimentController(ExperimentController):
    """
    Controller for image-based experiments (camera).

    Replaces procedural code in mainControl_camera.py with class-based approach.

    Instruments:
    - SignalGenerator (or simulation)
    - PulseBlaster
    - Camera (CameraWorker)
    - AnalogOutputTask (DAQ AO) for rotating fields
    """

    def __init__(self, params_dict: dict, camera_mode: str = 'level'):
        """
        Initialize camera experiment controller.

        Args:
            params_dict: Parameter dictionary from config file
            camera_mode: 'level', 'sync', or 'timeseries' trigger mode
        """
        super().__init__(params_dict, 'camera')
        self.camera_mode = camera_mode

    def initialize_instruments(self, trial_run: List[str] = ['n', 'y']):
        """
        Initialize instruments for camera-based experiments.

        Args:
            trial_run: [sg_mode, camera_mode] where 'n'=real, 'y'=simulation
        """
        print("\n" + "="*70)
        print(" "*15 + "Initializing Camera Experiment Instruments")
        print("="*70)

        # Initialize Signal Generator
        try:
            if trial_run[0] == 'n':
                self.instruments['sg'] = SignalGenerator(name="sg1")
                self.instruments['sg'].enable_sg_output()

                mw_params = self.params_dict.get('mw', {})
                if 'power' in mw_params:
                    self.instruments['sg'].set_sg_amp(mw_params['power'])
                if 'freq' in mw_params:
                    freq = mw_params['freq']
                    if isinstance(freq, np.ndarray):
                        freq = freq[0]
                    self.instruments['sg'].set_sg_freq(freq)

                print("✔ SG initialized")
            else:
                self.instruments['sg'] = SignalGenerator_sim(name="sg1")
                print("✔ Simulated SG initialized")
        except Exception as e:
            logging.exception(f"❌ Error initializing SG: {e}")
            self.instruments['sg'] = SignalGenerator_sim(name="sg1")

        # Initialize Pulse Blaster
        try:
            self.instruments['pb'] = PulseBlaster(self.params_dict, name="pb1")
            self.instruments['pb'].configure()
            print("✔ PulseBlaster configured")
        except Exception as e:
            logging.exception(f"❌ Error initializing PulseBlaster: {e}")
            raise

        # Initialize Camera
        # NOTE: Camera initialization requires Camcontrol module
        # This is a placeholder - actual implementation needs CameraWorker
        try:
            # from Camcontrol import CameraWorker
            # import dcamcon
            # self.instruments['camera'] = CameraWorker(simulate=trial_run[1]=='y')
            #
            # # Configure trigger mode
            # if 'level' in self.camera_mode:
            #     self.instruments['camera'].triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.LEVEL
            # elif 'sync' in self.camera_mode:
            #     self.instruments['camera'].triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.SYNCREADOUT
            #
            # print(f"✔ Camera initialized (mode: {self.camera_mode})")
            print("ℹ Camera initialization placeholder (needs Camcontrol module)")
            self.instruments['camera'] = None
        except Exception as e:
            logging.exception(f"❌ Error initializing Camera: {e}")
            self.instruments['camera'] = None

        # Initialize Analog Output for rotating fields
        try:
            self.instruments['ao_task'] = AnalogOutputTask(
                dev=concfg.daq_dev,
                channels=concfg.output_terminals,
                sampling_rate=1000,
                name="ao_task"
            )
            print("✔ Analog Output configured")
        except Exception as e:
            logging.exception(f"❌ Error initializing AO task: {e}")
            self.instruments['ao_task'] = None

        print("="*70)
        print("✔ All instruments initialized")
        print("="*70 + "\n")

    def _acquire_single_point(self, param_values: dict) -> dict:
        """
        Acquire data at one parameter point for camera-based experiments.

        This method would contain the acquisition logic from mainControl_camera.py.
        For now, it's a placeholder.

        Args:
            param_values: Dict of current parameter values

        Returns:
            Dict with measurement results
        """
        # TODO: Implement actual acquisition logic from mainControl_camera.py
        # This would include:
        # 1. Set parameters
        # 2. Program PulseBlaster sequence
        # 3. Prepare camera
        # 4. Trigger sequence
        # 5. Acquire images
        # 6. Process images
        # 7. Return results

        # Placeholder implementation
        return {
            **param_values,
            'mean_signal': 0.0,  # Placeholder
            'contrast': 0.0,  # Placeholder
            'timestamp': time.time()
        }


if __name__ == "__main__":
    print("\n" + "="*70)
    print(" "*15 + "Experiment Controller Module Loaded")
    print("="*70)
    print("\nAvailable classes:")
    print("  - ExperimentController (ABC): Base class")
    print("  - DiodeExperimentController: For APD/diode experiments")
    print("  - CameraExperimentController: For camera experiments")
    print("\nBenefits:")
    print("  - Common interface for both acquisition modes")
    print("  - Clean separation of mode-specific logic")
    print("  - Future unification straightforward")
    print("  - Reusable components across experiments")
    print("="*70 + "\n")
