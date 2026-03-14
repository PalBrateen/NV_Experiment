"""
Core parameter system for NV Experiment control.

This module provides a flexible, performance-optimized framework for controlling
multiple instruments in physics experiments with multi-parameter sweeps.

Key Features:
- Automatic parameter-to-instrument association via ParameterRegistry
- Function pointer optimization for inner loop performance
- Flexible parameter sweep ordering
- N-dimensional nested sweeps
- YAML-based configuration support

Date: December 2024
"""
#%%
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional, Callable, Tuple
import numpy as np
import yaml
import json
from datetime import datetime
from pathlib import Path


# ============================================================================
# CORE REGISTRY SYSTEM
# ============================================================================

class ParameterRegistry:
    """
    Global registry mapping parameter names to instrument types.

    This enables automatic instrument lookup when adding parameter sweeps,
    allowing users to write: add_sweep('frequency', values)
    instead of: add_sweep('signal_generator', 'frequency', values)
    """

    _registry: Dict[str, str] = {}

    @classmethod
    def register(cls, parameter: str, instrument_type: str):
        """Register a parameter with its instrument type."""
        cls._registry[parameter] = instrument_type

    @classmethod
    def get_instrument_for_parameter(cls, parameter: str) -> Optional[str]:
        """Get the instrument type for a given parameter."""
        return cls._registry.get(parameter)

    @classmethod
    def get_all_parameters(cls) -> Dict[str, str]:
        """Get a copy of all registered parameters."""
        return cls._registry.copy()

    @classmethod
    def clear(cls):
        """Clear the registry (useful for testing)."""
        cls._registry.clear()


class Parameter:
    """
    Represents an instrument parameter with validation.

    Attributes:
        name: Parameter name
        value: Current value
        valid_range: Optional tuple (min, max) for validation
        units: Physical units (e.g., "Hz", "dBm", "s")
    """

    def __init__(self, name: str, initial_value: Any,
                 valid_range: Optional[Tuple] = None, units: str = ""):
        self.name = name
        self.value = initial_value
        self.valid_range = valid_range
        self.units = units

    def validate(self, value: Any) -> bool:
        """
        Validate a value against the parameter's valid range.

        Returns:
            True if valid, False otherwise
        """
        if self.valid_range is None:
            return True
        return self.valid_range[0] <= value <= self.valid_range[1]

    def __repr__(self):
        return f"Parameter({self.name}={self.value} {self.units}, range={self.valid_range})"


# ============================================================================
# INSTRUMENT BASE CLASS
# ============================================================================

class Instrument(ABC):
    """
    Abstract base class for all instruments.

    Subclasses must implement:
    - _register_parameters(): Register all instrument parameters
    - _update_instrument(): Send commands to actual hardware

    Provides two parameter setting methods:
    - change_parameter(): Safe with validation (use for setup/outer loops)
    - _set_<param>_direct(): Fast, no validation (implement for inner loops)
    """

    def __init__(self, name: str):
        """
        Initialize instrument.

        Args:
            name: Unique identifier for this instrument instance
        """
        self.name = name
        self._parameters: Dict[str, Parameter] = {}
        self._register_parameters()

    @abstractmethod
    def _register_parameters(self):
        """
        Register all parameters for this instrument.

        Implement in subclasses using:
            self.register_parameter(name, initial_value, valid_range, units)
        """
        pass

    def register_parameter(self, param_name: str, initial_value: Any,
                          valid_range: Optional[Tuple] = None, units: str = ""):
        """
        Register a parameter with automatic registry update.

        This method:
        1. Creates a Parameter object
        2. Registers it with the global ParameterRegistry
        3. Sets it as an instance attribute

        Args:
            param_name: Name of the parameter
            initial_value: Initial/default value
            valid_range: Optional (min, max) tuple for validation
            units: Physical units string
        """
        self._parameters[param_name] = Parameter(param_name, initial_value, valid_range, units)
        ParameterRegistry.register(param_name, self.__class__.__name__)
        setattr(self, param_name, initial_value)

    def change_parameter(self, parameter_name: str, new_value: Any):
        """
        Safe parameter setting with validation.

        Use this for:
        - Setup/initialization
        - Outer loop parameters
        - User-facing parameter changes

        Args:
            parameter_name: Name of parameter to change
            new_value: New value to set

        Raises:
            AttributeError: If parameter not registered
            ValueError: If value outside valid range
        """
        if parameter_name not in self._parameters:
            raise AttributeError(f"{parameter_name} not valid for {self.name}")

        param = self._parameters[parameter_name]
        if not param.validate(new_value):
            raise ValueError(f"Value {new_value} outside range {param.valid_range}")

        param.value = new_value
        setattr(self, parameter_name, new_value)
        self._update_instrument(parameter_name, new_value)

    @abstractmethod
    def _update_instrument(self, parameter_name: str, new_value: Any):
        """
        Update actual hardware with new parameter value.

        Implement in subclasses with SCPI/VISA/API commands.

        Args:
            parameter_name: Name of parameter that changed
            new_value: New value to send to hardware
        """
        pass

    def get_parameters(self) -> Dict[str, Any]:
        """Get current values of all parameters."""
        return {name: param.value for name, param in self._parameters.items()}

    def get_parameter_info(self, param_name: str) -> Parameter:
        """Get Parameter object with full info (value, range, units)."""
        return self._parameters.get(param_name)

    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}', params={list(self._parameters.keys())})"


# ============================================================================
# PARAMETER SWEEP SYSTEM WITH FUNCTION POINTER OPTIMIZATION
# ============================================================================

class ParameterSweep:
    """
    Automatic multi-parameter sweep with smart setters.

    - First added parameter = innermost loop (fastest changing)
    - Each instrument handles its own update logic via _set_X_direct() methods
    - Falls back to change_parameter() if no direct setter exists

    Performance: <10μs Python overhead per inner loop iteration.
    """

    def __init__(self, instruments: Dict[str, Instrument]):
        """
        Initialize parameter sweep.

        Args:
            instruments: Dict of instrument instances {name: instrument}
        """
        self.instruments = instruments
        self.sweep_params: List[Tuple[str, str, np.ndarray]] = []  # (inst_name, param, values)
        self._setters: Dict[str, Callable] = {}  # {param_name: setter_function}

    def add_sweep(self, parameter: str, values, position: Optional[int] = None):
        """
        Add parameter to sweep. First added = innermost loop.

        Automatically:
        - Finds instrument via ParameterRegistry
        - Registers fast setter if available, else uses change_parameter()

        Args:
            parameter: Parameter name (e.g., 'frequency', 'pulse_duration')
            values: Array/list of values to sweep
            position: Optional position (None = append)

        Raises:
            ValueError: If parameter not registered or instrument not found
        """
        # Find instrument type from registry
        inst_type = ParameterRegistry.get_instrument_for_parameter(parameter)
        if not inst_type:
            raise ValueError(f"Unknown parameter '{parameter}'. Not registered with any instrument.")

        # Find instrument instance
        instrument = None
        for inst in self.instruments.values():
            if inst.__class__.__name__ == inst_type:
                instrument = inst
                break

        if not instrument:
            raise ValueError(f"No {inst_type} instance in instruments dict.")

        # Add to sweep_params
        entry = (instrument.name, parameter, np.array(values))
        if position is None:
            self.sweep_params.append(entry)
        else:
            self.sweep_params.insert(position, entry)

        # Register setter (fast direct method or fallback)
        direct_method = f'_set_{parameter}_direct'
        if hasattr(instrument, direct_method):
            self._setters[parameter] = getattr(instrument, direct_method)
            print(f"  ✓ {parameter}: using fast setter")
        else:
            # Fallback to safe method
            self._setters[parameter] = lambda v, i=instrument, p=parameter: i.change_parameter(p, v)
            print(f"  ✓ {parameter}: using change_parameter()")

    def get_total_points(self) -> int:
        """Get total number of measurement points."""
        if not self.sweep_params:
            return 0
        return int(np.prod([len(vals) for _, _, vals in self.sweep_params]))

    def print_sweep_info(self):
        """Print sweep configuration summary."""
        print("\n" + "=" * 50)
        print("Sweep Configuration")
        print("=" * 50)
        for i, (inst, param, vals) in enumerate(self.sweep_params):
            level = "INNER" if i == 0 else f"Level {i}"
            print(f"  {param}: {len(vals)} points ({level})")
        print(f"\nTotal points: {self.get_total_points():,}")
        print("=" * 50 + "\n")

    def run(self, measure_func: Callable) -> List[Dict]:
        """
        Execute the sweep.

        At each parameter combination:
        1. Setters are called (instruments handle reprogramming internally)
        2. measure_func is called to acquire data

        Args:
            measure_func: Callable(instruments, current_params) -> dict

        Returns:
            List of measurement result dicts
        """
        print("Inside Run...")
        
        results = []
        total = self.get_total_points()
        count = [0]  # Use list for mutable closure

        def recursive_sweep(depth: int, current_params: dict):
            # Map depth to sweep_params index (reversed: first added = innermost)
            actual_index = len(self.sweep_params) - 1 - depth

            if actual_index < 0:
                # All parameters set - measure!
                count[0] += 1
                if count[0] % 10 == 0 or count[0] == total:
                    print(f"  Progress: {count[0]}/{total}")

                result = measure_func(self.instruments, current_params.copy())
                if result is not None:
                    results.append(result)
                return

            # Get this level's sweep info
            inst_name, param, values = self.sweep_params[actual_index]
            setter = self._setters[param]

            # Iterate over values
            for value in values:
                setter(value)  # Set parameter (may trigger PB reprogram)
                current_params[param] = value
                recursive_sweep(depth + 1, current_params)

        print(f"\nStarting sweep: {total} points")
        recursive_sweep(0, {})
        print(f"✓ Sweep complete: {len(results)} results\n")

        return results


# ============================================================================
# EXPERIMENT CONFIGURATION
# ============================================================================

class ExperimentConfig:
    """
    Configuration for experiment sequences.

    Supports both programmatic creation and loading from YAML files.
    """

    def __init__(self, sequence_type: str, instruments: Dict[str, Instrument],
                 config_path: Optional[str] = None):
        """
        Initialize experiment configuration.

        Args:
            sequence_type: Experiment type ('esr', 'rabi', 'hahn_echo', etc.)
            instruments: Dict of instrument instances
            config_path: Optional path to YAML config file
        """
        self.sequence_type = sequence_type
        self.base_params = {
            'num_runs': 1,
            'samples_per_point': 1,
            'save_dir': './data',
            'contrast_mode': 'ratio_signal_over_reference'
        }
        self.parameter_values = {}
        self.sweep = ParameterSweep(instruments)

        if config_path:
            self.load_config(config_path)
        else:
            self._setup_sequence()

    def _setup_sequence(self):
        """Setup default parameters based on sequence type."""
        if self.sequence_type == 'esr':
            self._setup_esr()
        elif self.sequence_type == 'rabi':
            self._setup_rabi()
        elif self.sequence_type == 'hahn_echo':
            self._setup_hahn_echo()
        else:
            print(f"ℹ Unknown sequence type: {self.sequence_type}. Using default settings.")

    def _setup_esr(self):
        """ESR sequence defaults - CUSTOMIZE THIS."""
        self.base_params.update({'num_runs': 30, 'samples_per_point': 8})
        self.parameter_values.update({'frequency': 2.87e9, 'power': 8})

        # Default frequency sweep
        freq_points = np.arange(2.82e9, 2.92e9, 1e6)
        self.sweep.add_sweep("frequency", freq_points)

    def _setup_rabi(self):
        """Rabi sequence defaults - CUSTOMIZE THIS."""
        self.base_params.update({'num_runs': 1, 'samples_per_point': 2})
        self.parameter_values.update({'frequency': 2.87e9, 'power': 8})

        # Default pulse duration sweep
        duration_points = np.arange(10e-9, 500e-9, 2e-9)
        self.sweep.add_sweep("pulse_duration", duration_points)

    def _setup_hahn_echo(self):
        """Hahn Echo sequence defaults - CUSTOMIZE THIS."""
        self.base_params.update({'num_runs': 20, 'samples_per_point': 4})
        self.parameter_values.update({'frequency': 2.87e9, 'power': 8})

        # Default tau sweep
        tau_points = np.linspace(100e-9, 10e-6, 100)
        self.sweep.add_sweep("tau", tau_points)

    @classmethod
    def from_legacy_dict(cls, params_dict: Dict, instruments: Dict[str, Instrument]) -> 'ExperimentConfig':
        """
        Create ExperimentConfig from legacy params_dict format.

        This enables backward compatibility with existing config files.

        Args:
            params_dict: Legacy parameter dictionary
            instruments: Dict of instrument instances

        Returns:
            ExperimentConfig instance
        """
        sequence = params_dict.get('seq', {}).get('sequence', 'unknown')
        config = cls(sequence, instruments, config_path=None)

        # Import base parameters
        scan_params = params_dict.get('scan', {})
        config.base_params['num_runs'] = scan_params.get('Nruns', 1)

        # Import parameter values
        if 'mw' in params_dict:
            config.parameter_values.update(params_dict['mw'])

        # Import sweep parameters
        if 'names' in scan_params and 'values' in scan_params:
            param_names = scan_params['names']
            param_values = scan_params['values']

            for name, values in zip(param_names, param_values):
                config.sweep.add_sweep(name, values)

        return config

    def save_config(self, path: str):
        """
        Save configuration to YAML file.

        Args:
            path: Output file path
        """
        config = {
            'sequence_type': self.sequence_type,
            'base_params': self.base_params,
            'parameter_values': self.parameter_values,
            'sweeps': [
                {
                    'parameter': p,
                    'values': v.tolist() if isinstance(v, np.ndarray) else v
                }
                for _, p, v in self.sweep.sweep_params
            ]
        }

        with open(path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, indent=4)

        print(f"✓ Configuration saved to {path}")

    def load_config(self, path: str):
        """
        Load configuration from YAML file.

        Args:
            path: Input file path
        """
        with open(path, 'r') as f:
            config = yaml.safe_load(f)

        self.sequence_type = config['sequence_type']
        self.base_params = config['base_params']
        self.parameter_values = config['parameter_values']

        # Clear and reload sweeps
        self.sweep.sweep_params = []
        for s in config.get('sweeps', []):
            self.sweep.add_sweep(s['parameter'], np.array(s['values']))

        print(f"✓ Configuration loaded from {path}")


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def run_experiment(config: ExperimentConfig) -> List[Dict]:
    """
    Run experiment using configuration.

    This function:
    1. Sets non-swept parameters to their configured values
    2. Sets swept parameters to their first values
    3. Executes the sweep
    4. Returns results

    Args:
        config: ExperimentConfig instance

    Returns:
        List of measurement result dictionaries
    """

    # Set non-swept parameters
    for param, value in config.parameter_values.items():
        is_swept = any(param == sp for _, sp, _ in config.sweep.sweep_params)
        if not is_swept:
            inst_type = ParameterRegistry.get_instrument_for_parameter(param)
            if inst_type:
                for inst in config.sweep.instruments.values():
                    if inst.__class__.__name__ == inst_type:
                        inst.change_parameter(param, value)
                        break

    # Set swept parameters to first values
    for inst_name, param, values in config.sweep.sweep_params:
        config.sweep.instruments[inst_name].change_parameter(param, values[0])

    # Define placeholder measurement function
    def measure_point(instruments, current_params):
        """Placeholder measurement function - override in config."""
        return {**current_params, 'signal': 0.0, 'timestamp': 0}

    # Run sweep with measurement function
    results = config.sweep.run(measure_point)

    return results


if __name__ == "__main__":
    print("\n" + "="*60)
    print("Parameter System Module Loaded Successfully")
    print("="*60)
    print("\nAvailable classes:")
    print("  - ParameterRegistry: Global parameter→instrument mapping")
    print("  - Parameter: Parameter with validation")
    print("  - Instrument: Abstract base class for instruments")
    print("  - ParameterSweep: N-dimensional sweep executor (function pointer optimized)")
    print("  - ExperimentConfig: Configuration management")
    print("\nUtility functions:")
    print("  - run_experiment(config): Execute experiment from config")
    print("="*60 + "\n")
