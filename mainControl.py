"""
Unified Main Control for NV Experiment

This unified entry point handles both:
- Diode/APD-based experiments (voltage/count measurements)
- Camera-based experiments (image measurements)

The acquisition mode is automatically detected from the config file or can be
explicitly specified.

Usage:
    python mainControl.py --config esr_config --mode diode
    python mainControl.py --config rabi_config --mode camera
    python mainControl.py --config esr_config  # Auto-detect from config

Date: December 2024
"""

# ============================================================================
# IMPORTS
# ============================================================================

import connectionConfig as concfg
import matplotlib.pyplot as plt
import numpy as np
import time
import dialog
import psutil
import os
import logging
import argparse
from pathlib import Path
from importlib import import_module

# Unified controller imports
from experiment_controller import DiodeExperimentController, CameraExperimentController
from parameter_system import ParameterSweep
from save_manager import SaveManager
from config_loader import load_config as load_yaml_config

# Legacy imports for backward compatibility
from PBcontrol import PulseBlaster, ns, ms, us, s, Inst
from DAQcontrol import AnalogOutputTask, AnalogInputTask
from sequencecontrol import sequencecontrol
from SGcontrol import SignalGenerator, SignalGenerator_sim, GHz

plt.rcParams.update({'figure.max_open_warning': 0})

# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

# Acquisition mode detection
DIODE_SEQUENCES = ['esr_seq', 'rabi_seq', 't1_seq', 't2_seq', 'aom_timing']
CAMERA_SEQUENCES = ['cam_levelm', 'cam_syncm', 'cam_timeseries']

# Default settings
DEFAULT_TRIAL_RUN = ['n', 'n']  # [sg_mode, pb/camera_mode]
DEFAULT_PLOT_DPI = 100
DEFAULT_VOLTAGE_UNIT = 1  # mV

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def detect_acquisition_mode(params_dict: dict) -> str:
    """
    Auto-detect acquisition mode from config parameters.

    Args:
        params_dict: Parameter dictionary from config file

    Returns:
        'diode' or 'camera'
    """
    sequence = params_dict.get('seq', {}).get('sequence', '')

    # Check sequence name
    if any(seq in sequence for seq in CAMERA_SEQUENCES):
        return 'camera'
    elif any(seq in sequence for seq in DIODE_SEQUENCES):
        return 'diode'

    # Check for camera-specific parameters
    if 'camera' in params_dict or 'exposure_time' in params_dict.get('seq', {}):
        return 'camera'

    # Default to diode
    return 'diode'


def set_core_affinity(cores=[0, 1]):
    """Pin process to specific CPU cores for better performance."""
    try:
        process = psutil.Process(os.getpid())
        process.cpu_affinity(cores)
        print(f"✔ Process pinned to cores: {process.cpu_affinity()}")
    except Exception as e:
        print(f"⚠ Could not set core affinity: {e}")


def setup_plotting():
    """Configure matplotlib for experiment plotting."""
    plt.rcParams.update({
        'figure.max_open_warning': 0,
        'font.size': 10,
        'figure.dpi': DEFAULT_PLOT_DPI
    })


def load_experiment_config(config_path: str) -> dict:
    """
    Load experiment configuration from either Python (.py) or YAML (.yaml/.yml) file.

    Args:
        config_path: Path to config file (with or without extension)

    Returns:
        params_dict compatible with experiment system

    Examples:
        # Load Python config
        params = load_experiment_config('esr_config')
        params = load_experiment_config('esr_config.py')

        # Load YAML config
        params = load_experiment_config('configs/esr_config.yaml')
        params = load_experiment_config('configs/rabi_config.yml')
    """
    config_path:Path = Path(config_path)

    # Try different extensions if no extension provided
    if not config_path.suffix:
        # Try YAML first
        for ext in ['.yaml', '.yml']:
            yaml_path = config_path.with_suffix(ext)
            if yaml_path.exists():
                print(f"✓ Loading YAML config: {yaml_path}")
                return load_yaml_config(str(yaml_path))

        # Try Python
        py_path = config_path.with_suffix('.py')
        if py_path.exists():
            config_name = config_path.stem
            print(f"✓ Loading Python config: {py_path}")
            expCfg = import_module(config_name)
            return expCfg.params_dict

        # Try without extension (import by name)
        try:
            config_name = config_path.name
            print(f"✓ Loading Python config module: {config_name}")
            expCfg = import_module(config_name)
            return expCfg.params_dict
        except:
            pass

        raise FileNotFoundError(f"Config file not found: {config_path} (.yaml, .yml, or .py)")

    # Extension provided - load based on extension
    if config_path.suffix in ['.yaml', '.yml']:
        if not config_path.exists():
            raise FileNotFoundError(f"YAML config not found: {config_path}")
        print(f"✓ Loading YAML config: {config_path}")
        return load_yaml_config(str(config_path))

    elif config_path.suffix == '.py':
        if not config_path.exists():
            raise FileNotFoundError(f"Python config not found: {config_path}")
        config_name = config_path.stem
        print(f"✓ Loading Python config: {config_path}")
        expCfg = import_module(config_name)
        return expCfg.params_dict

    else:
        raise ValueError(f"Unknown config format: {config_path.suffix} (expected .py, .yaml, or .yml)")


# ============================================================================
# GENERIC MEASUREMENT FUNCTIONS
# ============================================================================

def diode_measurement(instruments, current_params):
    """
    Generic diode/APD measurement.
    Parameters already set by sweep (including PB reprogrammed if needed).

    Args:
        instruments: Dict of instrument instances (with '_config' key)
        current_params: Dict of current parameter values

    Returns:
        Dict with measurement results
    """
    pb:PulseBlaster = instruments['pb']
    ai_task:AnalogInputTask = instruments['ai_task']
    Nsamples:int = instruments['_config']['seq']['Nsamples']

    print(f"ai_task = {ai_task}")

    pb.start_sequence()
    data = ai_task.read_daq(Nsamples)
    print(f"data (len={len(data)})= {data[0:10]}")

    return {
        **current_params,
        'signal': np.mean(data[0]) if isinstance(data[0], (list, np.ndarray)) else data[0],
        'reference': np.mean(data[1]) if isinstance(data[1], (list, np.ndarray)) else data[1],
        'timestamp': time.time()
    }


def camera_measurement(instruments, current_params):
    """
    Generic camera measurement.
    Parameters already set by sweep (including PB reprogrammed if needed).

    Args:
        instruments: Dict of instrument instances (with '_config' key)
        current_params: Dict of current parameter values

    Returns:
        Dict with measurement results
    """
    pb = instruments['pb']
    camera = instruments['camera']
    Nsamples = instruments['_config']['seq']['Nsamples']

    pb.start_sequence()
    frames = camera.capture_frames(Nsamples)

    return {
        **current_params,
        'frames': frames,
        'mean_intensity': np.mean(frames),
        'timestamp': time.time()
    }


def run_sweep_experiment(params_dict, instruments, mode='diode'):
    """
    Universal sweep experiment runner.

    Args:
        params_dict: Configuration dict with 'sweep' key defining parameters
        instruments: Dict of instrument instances
        mode: 'diode' or 'camera'

    Returns:
        List of measurement results
    """
    from parameter_system import ParameterSweep

    # Store config reference for measurement functions
    instruments['_config'] = params_dict

    # Set PB instrument type
    if 'pb' in instruments:
        instruments['pb'].set_instr_type('diode' if mode == 'diode' else 'cam_levelm')

    # Initialize PB sequence BEFORE sweep (important for ESR where PB params aren't swept)
    if 'pb' in instruments:
        pb: PulseBlaster = instruments['pb']
        # Check if any PB parameters are being swept
        sweep_params = params_dict.get('sweep', {}).keys()
        pb_params = ['pulse_duration', 'tau', 't_AOM', 'ro_delay']
        pb_is_swept = any(param in pb_params for param in sweep_params)

        # If no PB parameters are swept, program it once now
        if not pb_is_swept:
            print("  ℹ Programming PulseBlaster (no PB params in sweep)")
            pb._reprogram_sequence(pb_channels=params_dict['pb'].get('chhanels', {}))

    # Build sweep
    sweep = ParameterSweep(instruments)

    sweep_config = params_dict.get('sweep', {})
    for param_name, values in sweep_config.items():
        sweep.add_sweep(param_name, values)

    sweep.print_sweep_info()

    # Select measurement function
    measure_func = diode_measurement if mode == 'diode' else camera_measurement

    # Run sweep
    results = sweep.run(measure_func)

    return results


# ============================================================================
# UNIFIED EXPERIMENT RUNNER
# ============================================================================

class UnifiedExperimentRunner:
    """
    Unified experiment runner for both diode and camera modes.

    This class:
    - Auto-detects acquisition mode
    - Creates appropriate controller (Diode or Camera)
    - Runs parameter sweeps
    - Saves data and parameters
    """

    def __init__(self, config_file: str, mode: str = 'auto', trial_run=None):
        """
        Initialize unified experiment runner.

        Args:
            config_file: Name of config file (without .py extension)
            mode: 'diode', 'camera', or 'auto' (auto-detect)
            trial_run: [sg_mode, pb/camera_mode] where 'n'=real, 'y'=simulation
        """
        self.config_file = config_file
        self.trial_run = trial_run or DEFAULT_TRIAL_RUN
        self.controller = None
        self.save_manager = SaveManager()

        # Load configuration
        print("\n" + "="*70)
        print(" "*20 + "Unified NV Experiment Control")
        print("="*70)
        print(f"\nLoading configuration: {config_file}")

        expCfg = import_module(config_file)
        self.params_dict = expCfg.params_dict

        # Detect or use specified mode
        if mode == 'auto':
            self.mode = detect_acquisition_mode(self.params_dict)
            print(f"✓ Auto-detected mode: {self.mode}")
        else:
            self.mode = mode
            print(f"✓ Using specified mode: {self.mode}")

        # Create appropriate controller
        if self.mode == 'diode':
            self.controller = DiodeExperimentController(self.params_dict)
        elif self.mode == 'camera':
            # Extract camera mode from params
            camera_mode = self._get_camera_mode()
            self.controller = CameraExperimentController(self.params_dict, camera_mode)
        else:
            raise ValueError(f"Unknown acquisition mode: {self.mode}")

        print(f"✓ Controller created: {type(self.controller).__name__}")

    def _get_camera_mode(self) -> str:
        """Extract camera trigger mode from params."""
        sequence = self.params_dict.get('seq', {}).get('sequence', '')

        if 'level' in sequence:
            return 'level'
        elif 'sync' in sequence:
            return 'sync'
        elif 'timeseries' in sequence:
            return 'timeseries'
        else:
            return 'level'  # Default

    def initialize(self):
        """Initialize instruments."""
        print("\n" + "="*70)
        print("Initializing instruments...")
        print("="*70)

        # Set core affinity for performance
        set_core_affinity()

        # Initialize controller instruments
        self.controller.initialize_instruments(self.trial_run)

        print("\n✓ All instruments initialized successfully")

    def run_experiment(self, use_sweep: bool = True):
        """
        Run the experiment.

        Args:
            use_sweep: If True, use ParameterSweep system (recommended)
                      If False, use legacy single-parameter scan

        Returns:
            data_array: NumPy array with experiment results
        """
        print("\n" + "="*70)
        print("Running Experiment...")
        print("="*70)

        if use_sweep:
            return self._run_with_parameter_sweep()
        else:
            return self._run_legacy_scan()

    def _run_with_parameter_sweep(self):
        """
        Run experiment using new ParameterSweep system.

        This is the recommended method - enables multi-parameter sweeps
        and function pointer optimization.
        """
        # Get instruments
        instruments = self.controller.get_instruments_dict()

        # Create parameter sweep
        sequence = self.params_dict['seq']['sequence']
        sweep = ParameterSweep(instruments, sequence)

        # Add sweep parameters from config
        scan_params = self.params_dict.get('scan', {})
        param_names = scan_params.get('names', [])
        param_values = scan_params.get('values', [])

        if not param_names:
            # Legacy single-parameter format
            param_names = ['frequency']  # Default
            param_values = [self.params_dict['mw'].get('freq', np.array([2.87e9]))]

        # Add sweeps (automatically registers fast setters!)
        for name, values in zip(param_names, param_values):
            sweep.add_sweep(name, values)

        # Set measurement function
        def measure_point(instr_dict, current_values):
            """Measurement function for each parameter point."""
            return self.controller._acquire_single_point(current_values)

        sweep.set_measurement_function(measure_point)

        # Run sweep
        print("\nStarting parameter sweep...")
        results = sweep.run()

        # Convert results to numpy array
        data_array = self._convert_results_to_array(results)

        return data_array

    def _run_legacy_scan(self):
        """
        Legacy single-parameter scan (for backward compatibility).

        This method replicates the old mainControl_diode.py behavior.
        """
        print("\nUsing legacy single-parameter scan...")

        # Extract scan parameters
        scan_params = self.params_dict.get('scan', {})
        Nruns = scan_params.get('Nruns', 1)
        Nscanpts = scan_params.get('Nscanpts', [101])[0]

        # Get parameter values
        param_values = self.params_dict['mw'].get('freq', np.linspace(2.82e9, 2.92e9, Nscanpts))

        # Initialize data array
        if self.mode == 'diode':
            Nsamples = self.params_dict['seq']['Nsamples']
            data_array = np.zeros((Nruns, Nscanpts, 2))  # 2 channels
        else:
            # Camera mode
            data_array = []  # Will be filled with images

        # Scan loop
        for run in range(Nruns):
            print(f"\nRun {run+1}/{Nruns}")

            for i, value in enumerate(param_values):
                # Set parameter (e.g., frequency)
                current_values = {'frequency': value}

                # Acquire data at this point
                result = self.controller._acquire_single_point(current_values)

                # Store result
                if self.mode == 'diode':
                    data_array[run, i, 0] = result.get('signal', 0)
                    data_array[run, i, 1] = result.get('reference', 0)
                else:
                    data_array.append(result.get('images', []))

                # Progress
                if (i+1) % 10 == 0:
                    print(f"  Progress: {i+1}/{Nscanpts}")

        return data_array

    def _convert_results_to_array(self, results: list) -> np.ndarray:
        """
        Convert ParameterSweep results to numpy array.

        Args:
            results: List of result dicts from sweep.run()

        Returns:
            NumPy array matching existing format
        """
        Nruns = self.params_dict['scan'].get('Nruns', 1)
        Nscanpts = len(results)

        if self.mode == 'diode':
            # Diode: (Nruns, Nscanpts, 2)
            data_array = np.zeros((Nruns, Nscanpts, 2))

            for i, result in enumerate(results):
                data_array[0, i, 0] = result.get('signal', 0)
                data_array[0, i, 1] = result.get('reference', 0)

        else:
            # Camera: list of images
            data_array = [r.get('images', []) for r in results]

        return data_array

    def save_data(self, data_array, folder_number: str = None):
        """
        Save experiment data and parameters.

        Args:
            data_array: NumPy array or list with experiment results
            folder_number: Optional folder number (auto-generated if None)
        """
        print("\n" + "="*70)
        print("Saving Results...")
        print("="*70)

        # Get sequence name for folder
        sequence = self.params_dict['seq']['sequence']

        # Auto-generate folder number if not provided
        if folder_number is None:
            folder_number = self.save_manager.get_next_folder_number(sequence)

        # Save data and parameters (SEPARATE files!)
        self.save_manager.save_experiment(
            data=data_array,
            params=self.params_dict,
            folder_number=folder_number,
            sequence=sequence,
            param_format='yaml'
        )

        print(f"\n✓ Results saved successfully")
        print(f"  Folder: Saved_Data/{time.strftime('%Y-%m-%d')}/{sequence}_{folder_number}/")
        print(f"  Data: data_{folder_number}.npy")
        print(f"  Params: params_{folder_number}.yaml (SEPARATE!)")

    def cleanup(self):
        """Close all instruments and cleanup."""
        print("\n" + "="*70)
        print("Cleaning up...")
        print("="*70)

        if self.controller:
            self.controller.cleanup()

        print("\n✓ Cleanup complete")

    def run_full_experiment(self, save_data: bool = True):
        """
        Run complete experiment workflow.

        This is the main entry point for running experiments.

        Args:
            save_data: Whether to save results to disk
        """
        try:
            # Initialize instruments
            self.initialize()

            # Run experiment
            data_array = self.run_experiment(use_sweep=True)

            # Save data
            if save_data:
                self.save_data(data_array)

            print("\n" + "="*70)
            print(" "*20 + "Experiment Complete!")
            print("="*70)

            return data_array

        except Exception as e:
            print(f"\n❌ Error during experiment: {e}")
            logging.exception(e)
            raise

        finally:
            # Always cleanup
            self.cleanup()


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Unified NV Experiment Control',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mainControl.py --config esr_config --mode diode
  python mainControl.py --config rabi_config --mode camera
  python mainControl.py --config esr_config  # Auto-detect mode
  python mainControl.py --config esr_config --trial y n  # Simulated SG
        """
    )

    parser.add_argument(
        '--config', '-c',
        type=str,
        default='esr_config',
        help='Config file name (without .py extension)'
    )

    parser.add_argument(
        '--mode', '-m',
        type=str,
        choices=['diode', 'camera', 'auto'],
        default='auto',
        help='Acquisition mode (auto-detect if not specified)'
    )

    parser.add_argument(
        '--trial', '-t',
        nargs=2,
        metavar=('SG', 'PB/CAM'),
        default=['n', 'n'],
        help='Trial run modes: y=simulation, n=real hardware'
    )

    parser.add_argument(
        '--no-save',
        action='store_true',
        help='Do not save results to disk'
    )

    parser.add_argument(
        '--legacy',
        action='store_true',
        help='Use legacy single-parameter scan instead of ParameterSweep'
    )

    return parser.parse_args()


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """Main entry point for unified experiment control."""

    # Parse command line arguments
    args = parse_arguments()

    # Setup
    setup_plotting()

    # Create runner
    runner = UnifiedExperimentRunner(
        config_file=args.config,
        mode=args.mode,
        trial_run=args.trial
    )

    # Run experiment
    data_array = runner.run_full_experiment(save_data=not args.no_save)

    return data_array


# ============================================================================
# BACKWARD COMPATIBILITY MODE
# ============================================================================

# For legacy code that imports from mainControl_diode.py or mainControl_camera.py
# These can be set by external code to maintain compatibility

global trial_run, expCfgFile, params

# Default values for backward compatibility
trial_run = DEFAULT_TRIAL_RUN
expCfgFile = 'esr_config'
expCfg = None
params = None

def initialize_instr(sequence):
    """
    Legacy function for backward compatibility.

    This function maintains the old API from mainControl_diode.py
    """
    print("⚠ Using legacy initialize_instr() - consider migrating to UnifiedExperimentRunner")

    global expCfg, params

    # Load config if not already loaded
    if expCfg is None:
        expCfg = import_module(expCfgFile)
        params = expCfg.params_dict

    # Create controller
    mode = detect_acquisition_mode(params)

    if mode == 'diode':
        controller = DiodeExperimentController(params)
    else:
        controller = CameraExperimentController(params, camera_mode='level')

    # Initialize
    controller.initialize_instruments(trial_run)

    # Return instruments in old format
    instruments = controller.get_instruments_dict()
    sg = instruments.get('sg')
    ao_task = instruments.get('ao_task')

    return sg, ao_task


def close_all(sg=None, ao_task=None, ai_task=None):
    """
    Legacy function for backward compatibility.

    This function maintains the old API for closing instruments.
    """
    print("⚠ Using legacy close_all() - consider migrating to UnifiedExperimentRunner")

    if sg is not None:
        try:
            sg.uninit_sg()
            print('✔ SG closed...')
        except Exception as e:
            print(f'⚠ Error closing SG: {e}')

    if ao_task is not None:
        try:
            ao_task.stop()
            print('✔ AO task stopped...')
        except Exception as e:
            print(f'⚠ Error closing AO: {e}')

    if ai_task is not None:
        try:
            ai_task.stop()
            print('✔ AI task stopped...')
        except Exception as e:
            print(f'⚠ Error closing AI: {e}')


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    """
    Main entry point when running as script.

    Usage:
        python mainControl.py --config esr_config --mode diode
        python mainControl.py --config rabi_config --mode camera
    """
    main()
