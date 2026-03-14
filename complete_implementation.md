Use this implementation to replace/update the instrument control system in my codebase. Integrate it with the existing esr_config.py and rabi_config.py files.

Key Features:

*   Automatic parameter-to-instrument association via ParameterRegistry
*   Sequence-specific optimizations for inner loop performance
*   Flexible parameter sweep ordering
*   YAML-based configuration files
*   Clean API: add_sweep('frequency', values)

Date: December 2024
```python
from abc import ABC, abstractmethod from typing import Dict, List, Any, Optional, Callable
import numpy as np import yaml from datetime import datetime

# ============================================================================
# =============================================================================

# CORE REGISTRY SYSTEM
# ====================

# ============================================================================
# =============================================================================

class ParameterRegistry: """Global registry mapping parameter names to instrument types."""

    _registry: Dict[str, str] = {}

    @classmethod
    def register(cls, parameter: str, instrument_type: str):
        cls._registry[parameter] = instrument_type
    
    @classmethod
    def get_instrument_for_parameter(cls, parameter: str) -> Optional[str]:
        return cls._registry.get(parameter)
    
    @classmethod
    def get_all_parameters(cls) -> Dict[str, str]:
        return cls._registry.copy()

class Parameter: """Represents an instrument parameter with validation."""
    def __init__(self, name: str, initial_value: Any, valid_range: Optional[tuple] = None, units: str = ""):
        self.name = name self.value = initial_value self.valid_range = valid_range self.units = units

    def validate(self, value: Any) -> bool:
        if self.valid_range is None:
            return True
        return self.valid_range[0] <= value <= self.valid_range[1]

class Instrument(ABC): """Base class for all instruments."""

    def __init__(self, name: str):
        self.name = name
        self._parameters: Dict[str, Parameter] = {}
        self._register_parameters()
    
    @abstractmethod
    def _register_parameters(self):
        """Register all parameters. Implement in subclasses."""
        pass
    
    def register_parameter(self, param_name: str, initial_value: Any, valid_range: Optional[tuple] = None, units: str = ""):
        self._parameters[param_name] = Parameter(param_name, initial_value, valid_range, units)
        ParameterRegistry.register(param_name, self.__class__.__name__)
        setattr(self, param_name, initial_value)
    
    def change_parameter(self, parameter_name: str, new_value: Any):
        """Safe parameter setting with validation (use for setup/outer loops)."""
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
        """Update actual hardware. Implement in subclasses."""
        pass
    
    def get_parameters(self) -> Dict[str, Any]:
        return {name: param.value for name, param in self._parameters.items()}

# ============================================================================
# =============================================================================

# EXAMPLE INSTRUMENTS - REPLACE WITH YOUR ACTUAL INSTRUMENT CLASSES
# =================================================================

# ============================================================================
# =============================================================================

class SignalGenerator(Instrument):
    def _register_parameters(self):
        self.register_parameter("frequency", 1e9, (0, 20e9), "Hz")
        self.register_parameter("power", 0, (-130, 25), "dBm")
        self.register_parameter("phase", 0, (0, 360), "degrees")
        self.register_parameter("pulse_duration", 100e-9, (0, 1), "seconds")

    def _update_instrument(self, parameter_name: str, new_value: Any):
        # Replace with actual SCPI/VISA commands
        print(f"[{self.name}] {parameter_name} = {new_value}")
    
    def _set_frequency_direct(self, value: float):
        """Direct frequency control for inner loops - CUSTOMIZE THIS"""
        # Example: self.device.write(f":FREQ {value}")
        pass
    
    def _set_pulse_duration_direct(self, value: float):
        """Direct pulse duration control - CUSTOMIZE THIS"""
        pass

# ============================================================================
# =============================================================================

# PARAMETER SWEEP SYSTEM
# ======================

# ============================================================================
# =============================================================================

class ParameterSweep: """Handles multi-parameter sweeps with nested loops."""

    def __init__(self, instruments: Dict[str, Instrument], sequence_type: str):
        self.instruments = instruments
        self.sequence_type = sequence_type
        self.sweep_params: List[tuple] = []
        self.measurement_function: Optional[Callable] = None
    
    def add_sweep(self, parameter: str, values: List[Any], position: Optional[int] = None):
        """Add parameter sweep. Instrument looked up automatically."""
        instrument_type = ParameterRegistry.get_instrument_for_parameter(parameter)
        if not instrument_type:
            raise ValueError(f"Unknown parameter: {parameter}")
            
        instrument = None
        for inst in self.instruments.values():
            if inst.__class__.__name__ == instrument_type:
                instrument = inst
                break
                
        if not instrument:
            raise ValueError(f"No {instrument_type} found")
        
        sweep_entry = (instrument.name, parameter, np.array(values))
        if position is None:
            self.sweep_params.append(sweep_entry)
        else:
            self.sweep_params.insert(position, sweep_entry)
    
    def reorder_sweeps(self, parameter_order: List[str]):
        """Reorder sweeps by parameter names (outermost to innermost)."""
        new_sweep_params = []
        for param_name in parameter_order:
            for sweep_entry in self.sweep_params:
                if sweep_entry[1] == param_name:
                    new_sweep_params.append(sweep_entry)
                    break
        self.sweep_params = new_sweep_params
    
    def print_sweep_order(self):
        """Print sweep order for verification."""
        print("Sweep order (outermost → innermost):")
        for i, (_, param, vals) in enumerate(self.sweep_params):
            print(f"  {i+1}. {param}: {len(vals)} points")
        if self.sweep_params:
            total = np.prod([len(v) for _, _, v in self.sweep_params])
            print(f"\nTotal measurements: {total}")
    
    def set_measurement_function(self, func: Callable):
        self.measurement_function = func
    
    def run(self) -> List[Dict]:
        """Execute the parameter sweep."""
        def set_parameter(inst_name: str, param: str, value: Any, depth: int):
            """Set parameter with depth-aware optimization."""
            instrument = self.instruments[inst_name]
            is_innermost = (depth == len(self.sweep_params) - 1)
            
            if is_innermost:
                # CUSTOMIZE THESE FOR YOUR SEQUENCES
                if self.sequence_type == 'esr' and param == 'frequency':
                    instrument._set_frequency_direct(value)
                    return
                elif self.sequence_type == 'rabi' and param == 'pulse_duration':
                    instrument._set_pulse_duration_direct(value)
                    return
            
            instrument.change_parameter(param, value)
        
        def recursive_sweep(index: int = 0, current_values: Optional[Dict] = None):
            if current_values is None:
                current_values = {}
            
            if index == len(self.sweep_params):
                if self.measurement_function:
                    result = self.measurement_function(self.instruments, current_values)
                    return [result] if result is not None else []
                return []
            
            results = []
            inst_name, param, values = self.sweep_params[index]
            
            for value in values:
                set_parameter(inst_name, param, value, loop_depth=index)
                current_values[param] = value
                results.extend(recursive_sweep(index + 1, current_values.copy()))
            
            return results
        
        return recursive_sweep()

# ============================================================================
# =============================================================================

# EXPERIMENT CONFIGURATION
# ========================

# ============================================================================
# =============================================================================

class ExperimentConfig:
    """Configuration for experiment sequences."""

    def __init__(self, sequence_type: str, instruments: Dict[str, Instrument],
                 config_path: Optional[str] = None):
        self.sequence_type = sequence_type
        self.base_params = {
            'num_runs': 1,
            'samples_per_point': 1,
            'save_dir': './data',
            'contrast_mode': 'ratio_signal_over_reference'
        }
        self.parameter_values = {}
        self.sweep = ParameterSweep(instruments, sequence_type)
        
        if config_path:
            self.load_config(config_path)
        else:
            self._setup_sequence()
    
    def _setup_sequence(self):
        if self.sequence_type == 'esr':
            self._setup_esr()
        elif self.sequence_type == 'rabi':
            self._setup_rabi()
        else:
            raise ValueError(f"Unknown sequence: {self.sequence_type}")
    
    def _setup_esr(self):
        """ESR sequence defaults - CUSTOMIZE THIS"""
        self.base_params.update({'num_runs': 30, 'samples_per_point': 8})
        self.parameter_values.update({'frequency': 2.87e9, 'power': 8})
        
        freq_points = np.arange(2.82e9, 2.92e9, 1e6)
        self.sweep.add_sweep("frequency", freq_points)
    
    def _setup_rabi(self):
        """Rabi sequence defaults - CUSTOMIZE THIS"""
        self.base_params.update({'num_runs': 1, 'samples_per_point': 2})
        self.parameter_values.update({'frequency': 2.87e9, 'power': 8})
        
        duration_points = np.arange(10e-9, 500e-9, 2e-9)
        self.sweep.add_sweep("pulse_duration", duration_points)
    
    def save_config(self, path: str):
        """Save to YAML."""
        config = {
            'sequence_type': self.sequence_type,
            'base_params': self.base_params,
            'parameter_values': self.parameter_values,
            'sweeps': [{'parameter': p, 'values': v.tolist()} 
                      for _, p, v in self.sweep.sweep_params]
        }
        with open(path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
    
    def load_config(self, path: str):
        """Load from YAML."""
        with open(path, 'r') as f:
            config = yaml.safe_load(f)
        self.sequence_type = config['sequence_type']
        self.base_params = config['base_params']
        self.parameter_values = config['parameter_values']
        self.sweep.sweep_params = []
        for s in config.get('sweeps', []):
            self.sweep.add_sweep(s['parameter'], np.array(s['values']))

def run_experiment(config: ExperimentConfig) -> List[Dict]:
"""Run experiment using configuration."""
    
    # Set non-swept parameters
    for param, value in config.parameter_values.items():
        is_swept = any(param == sp for _, sp, _ in config.sweep.sweep_params)
        if not is_swept:
            inst_type = ParameterRegistry.get_instrument_for_parameter(param)
            for inst in config.sweep.instruments.values():
                if inst.name == inst_type:
                    inst.change_parameter(param, value)
                    break

    # Set swept parameters to first values
    for inst_name, param, values in config.sweep.sweep_params:
        config.sweep.instruments[inst_name].change_parameter(param, values[0])
    
    print(f"\nRunning {config.sequence_type} experiment...")
    config.sweep.print_sweep_order()
    
    results = config.sweep.run()
    print(f"\nComplete: {len(results)} measurements")
    return results

```