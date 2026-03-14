Now let me create the documentation file:

Overview
--------

This system provides a flexible, performance-optimized framework for controlling multiple instruments in physics experiments (ESR, Rabi, Hahn Echo, etc.) with multi-parameter sweeps.

Key Design Decisions
--------------------

### 1\. Parameter Registry

**Why**: Automatically associates parameters with instruments

*   Users write: `add_sweep('frequency', values)`
*   Instead of: `add_sweep('signal_generator', 'frequency', values)`
*   Registry maintains mapping: `'frequency' → 'SignalGenerator'`

### 2\. Performance Optimization

**Two-tier parameter setting**:

*   `change_parameter()`: Safe with validation (use for setup/outer loops)
*   Direct methods `_set_<param>_direct()`: Fast, no validation (use for inner loops)

**Timing comparison**:

*   Python overhead (if checks, dict lookups): ~nanoseconds
*   Hardware communication: ~milliseconds
*   **Conclusion**: Python overhead is negligible (<0.01%)

### 3\. Sweep Ordering

**Critical for performance**: Last added = innermost = fastest changing

```python
sweep.add_sweep("power", [5, 8, 10])      # Outer: changes 3 times
sweep.add_sweep("frequency", freq_points)  # Inner: changes 303 times
```

Integration with Existing Code
------------------------------

### Step 1: Update Your Instrument Classes

Replace your current instrument classes with ones inheriting from `Instrument`:

```python
class YourSignalGenerator(Instrument):
    def _register_parameters(self):
        # Register all parameters
        self.register_parameter("frequency", 2.87e9, (0, 20e9), "Hz")
        self.register_parameter("power", 0, (-130, 25), "dBm")
        # Add all your parameters
    
    def _update_instrument(self, parameter_name: str, new_value: Any):
        # Your existing SCPI/VISA commands
        if parameter_name == "frequency":
            self.device.write(f":FREQ {new_value}")
        elif parameter_name == "power":
            self.device.write(f":POW {new_value}")
    
    # Add direct methods for inner loop parameters
    def _set_frequency_direct(self, value: float):
        """Fast frequency setting for ESR inner loop"""
        self.device.write(f":FREQ {value}")
```

### Step 2: Update Config Files

Convert your `esr_config.py` and `rabi_config.py` to use `ExperimentConfig`:

```python
# esr_config.py
def create_esr_config(instruments):
    config = ExperimentConfig('esr', instruments)
    
    # Customize defaults
    config.base_params['num_runs'] = 30
    config.parameter_values['power'] = 8
    
    # Customize sweep
    config.sweep.sweep_params = []  # Clear default
    freq_points = np.linspace(2.82e9, 2.92e9, 101)
    config.sweep.add_sweep("frequency", freq_points)
    
    return config
```
### Step 3: Update Measurement Function

```python
def measure(instruments, current_values):
    """
    Your measurement code.
    
    Args:
        instruments: Dict of instrument instances
        current_values: Dict like {'frequency': 2.87e9, 'power': 8}
    
    Returns:
        Dict with measurement results
    """
    # Get instruments
    sig_gen = instruments['sg1']
    camera = instruments['camera']
    
    # Trigger your sequence
    # ... your existing measurement code ...
    
    # Return results
    return {
        **current_values,  # Include all parameters
        'signal': signal_value,
        'contrast': contrast,
        'timestamp': datetime.now().isoformat()
    }
```

### Step 4: Run Experiments

```python
# Create instruments
sig_gen = YourSignalGenerator("sg1")
laser = YourLaser("laser1")
instruments = {"sg1": sig_gen, "laser1": laser}

# Create config
config = create_esr_config(instruments)

# Set measurement function
config.sweep.set_measurement_function(measure)

# Run
results = run_experiment(config)

# Save results
config.save_config("esr_experiment.yaml")
```

Multi-Parameter Sweeps
----------------------

### Order Control

```python

# Method 1: Add in order (first = outer, last = inner)
sweep.add_sweep("power", [5, 8, 10])        # Outermost
sweep.add_sweep("phase", [0, 90, 180])      # Middle
sweep.add_sweep("frequency", freq_points)   # Innermost (optimized)

# Method 2: Specify position
sweep.add_sweep("frequency", freq_points, position=2)  # Innermost

# Method 3: Reorder after adding
sweep.reorder_sweeps(["power", "phase", "frequency"])

# Verify order
sweep.print_sweep_order()
```

### Performance Best Practices

1.  **Put fastest-changing parameter innermost**
2.  **Implement direct methods for innermost parameters**
3.  **Use validation only for outer loops**

```python

# Good: frequency changes often, put it innermost
sweep.add_sweep("power", [5, 8, 10])         # 3 changes
sweep.add_sweep("frequency", freq_points)    # 303 changes

# Bad: opposite order
sweep.add_sweep("frequency", freq_points)    # 101 changes
sweep.add_sweep("power", [5, 8, 10])         # 1111 changes (wasteful!)
```

Adding New Sequence Types
-------------------------

### 1\. Add setup method to ExperimentConfig

```python
def _setup_hahn_echo(self):
    """Hahn Echo defaults"""
    self.base_params.update({'num_runs': 20})
    self.parameter_values.update({
        'frequency': 2.87e9,
        'power': 8
    })
    
    # Default tau sweep
    tau_points = np.linspace(100e-9, 10e-6, 100)
    self.sweep.add_sweep("tau", tau_points)
```

### 2\. Add optimization in ParameterSweep.run()

In the `set_parameter` function inside `run()`:

```python
elif self.sequence_type == 'hahn_echo' and param == 'tau':
    instrument._set_tau_direct(value)
    return
```

### 3\. Implement direct method in instrument

```python
class PulseGenerator(Instrument):
    def _set_tau_direct(self, value: float):
        """Fast tau setting for Hahn Echo"""
        self.device.set_delay(value)
```

YAML Configuration Files
------------------------

Example `esr_experiment.yaml`:

```yaml

    sequence_type: esr
    base_params:
      num_runs: 30
      samples_per_point: 8
      save_dir: ./data
    parameter_values:
      frequency: 2870000000
      power: 8
    sweeps:
      - parameter: frequency
        values: [2820000000, 2821000000, ...]
```
    
## Troubleshooting

### "No instrument found for parameter X"
- Parameter not registered in any instrument's `_register_parameters()`
- Check: `ParameterRegistry.get_all_parameters()`

### "No X instance found"
- Instrument not passed to ExperimentConfig
- Check: `instruments` dict contains the right instrument types

### Slow inner loop
- Implement `_set_<param>_direct()` method
- Add optimization in `ParameterSweep.run()`
- Verify sweep order with `print_sweep_order()`

## Migration Checklist

- [ ] Update instrument classes to inherit from `Instrument`
- [ ] Register all parameters in `_register_parameters()`
- [ ] Implement `_update_instrument()` with SCPI/VISA commands
- [ ] Add direct methods for performance-critical parameters
- [ ] Convert config files to use `ExperimentConfig` if required
- [ ] Update measurement function signature
- [ ] Test with small sweeps first
- [ ] Verify sweep order is correct
- [ ] Benchmark inner loop performance


