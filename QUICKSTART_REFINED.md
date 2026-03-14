# Quick Start Guide (Refined System)

**Date**: December 16, 2024

## TL;DR

The refined system is **simpler** and features **auto-reprogramming**. Here's everything you need to know:

---

## Run an Experiment (3 Steps)

### Step 1: Create Config

```python
# esr_config_new.py
from spinapi import ns, us, ms
from SGcontrol import GHz
import numpy as np

params_dict = {
    'seq': {
        'sequence': 'esr_seq',
        't_AOM': 8000 * ms,
        'Nsamples': 8,
        'PBchannels': {...}
    },
    'mw': {
        'power': 8,
    },
    'scan': {
        'Nruns': 30,
    },
    'sweep': {  # NEW FORMAT!
        'frequency': np.linspace(2.82*GHz, 2.92*GHz, 101)
    }
}
```

### Step 2: Initialize Instruments

```python
from mainControl import run_sweep_experiment
from SGcontrol import SignalGenerator
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogInputTask
from esr_config_new import params_dict

instruments = {
    'sg': SignalGenerator(name="sg1"),
    'pb': PulseBlaster(params_dict, name="pb1"),
    'ai_task': AnalogInputTask(..., name="ai_task")
}
```

### Step 3: Run Experiment

```python
results = run_sweep_experiment(params_dict, instruments, mode='diode')
```

**Done!** That's it. The PB auto-reprograms, data is acquired, results returned.

---

## What's New?

### 1. Auto-Reprogramming 🚀

**PulseBlaster now auto-reprograms** when sweep parameters change!

```python
# Rabi sweep: PB reprograms for EACH pulse_duration!
'sweep': {
    'pulse_duration': np.arange(10*ns, 500*ns, 2*ns)
}

# You don't need to manually reprogram - it's automatic!
```

### 2. Simpler Config Format

**Before**:
```python
'scan': {
    'names': ['frequency'],
    'values': [freq_values],
    'Nscanpts': [101],
}
```

**After**:
```python
'sweep': {
    'frequency': freq_values
}
```

### 3. One Function for Everything

```python
# Works for ESR, Rabi, T2, Hahn Echo, etc.
results = run_sweep_experiment(params_dict, instruments, mode='diode')
```

---

## Example Experiments

### ESR (Frequency Sweep)

```python
# esr_config_new.py
params_dict = {
    'seq': {'sequence': 'esr_seq', ...},
    'sweep': {
        'frequency': np.linspace(2.82*GHz, 2.92*GHz, 101)
    }
}

results = run_sweep_experiment(params_dict, instruments, 'diode')
```

### Rabi (Pulse Duration Sweep)

```python
# rabi_config_new.py
params_dict = {
    'seq': {'sequence': 'rabi_seq', ...},
    'sweep': {
        'pulse_duration': np.arange(10*ns, 500*ns, 2*ns)
        # PB auto-reprograms for each value!
    }
}

results = run_sweep_experiment(params_dict, instruments, 'diode')
```

### T2/Hahn Echo (Tau Sweep)

```python
# t2_config.py
params_dict = {
    'seq': {'sequence': 'hahn_echo_seq', ...},
    'sweep': {
        'tau': np.linspace(1*us, 100*us, 50)
        # PB auto-reprograms for each tau!
    }
}

results = run_sweep_experiment(params_dict, instruments, 'diode')
```

### 2D Sweep (Power × Frequency)

```python
params_dict = {
    'seq': {'sequence': 'esr_seq', ...},
    'sweep': {
        'power': np.array([5, 8, 10]),           # Outer loop
        'frequency': np.linspace(2.82e9, 2.92e9, 51),  # Inner loop
    }
}

results = run_sweep_experiment(params_dict, instruments, 'diode')
# Total: 3 × 51 = 153 measurements
```

---

## API Reference

### run_sweep_experiment()

```python
run_sweep_experiment(params_dict, instruments, mode='diode')
```

**Args**:
- `params_dict`: Config dict with `'sweep'` key
- `instruments`: Dict of instrument instances
- `mode`: `'diode'` or `'camera'`

**Returns**:
- List of measurement result dicts

**What it does**:
1. Creates `ParameterSweep`
2. Adds all parameters from `params_dict['sweep']`
3. Sets PB instrument type
4. Runs sweep with appropriate measurement function
5. Returns results

### diode_measurement()

```python
diode_measurement(instruments, current_params)
```

**Generic measurement function for diode/APD experiments.**

**Args**:
- `instruments`: Dict with `'pb'`, `'ai_task'`, `'_config'`
- `current_params`: Dict of current parameter values

**Returns**:
- Dict with `signal`, `reference`, `timestamp`, and all param values

### camera_measurement()

```python
camera_measurement(instruments, current_params)
```

**Generic measurement function for camera experiments.**

**Args**:
- `instruments`: Dict with `'pb'`, `'camera'`, `'_config'`
- `current_params`: Dict of current parameter values

**Returns**:
- Dict with `frames`, `mean_intensity`, `timestamp`, and all param values

---

## Config File Structure

```python
params_dict = {
    # Sequence parameters
    'seq': {
        'sequence': 'esr_seq',      # Sequence name
        't_AOM': 8000 * ms,         # Timing parameters
        'Nsamples': 8,              # Number of samples
        'PBchannels': {...},        # PB channel mapping
    },

    # MW parameters
    'mw': {
        'power': 8,                 # Fixed power (if not swept)
        'freq': 2.87e9,             # Fixed frequency (if not swept)
    },

    # Scan parameters
    'scan': {
        'Nruns': 30,                # Number of runs
    },

    # PB parameters
    'pb': {
        'clk_cyc': 1e3/500e6,       # Clock cycle time
        'channels': {...},
    },

    # SWEEP PARAMETERS (NEW!)
    'sweep': {
        # Add parameters to sweep
        # First added = innermost loop (fastest changing)
        'frequency': np.linspace(2.82*GHz, 2.92*GHz, 101),
        # 'power': np.array([5, 8, 10]),  # Optional 2nd parameter
    },

    # Plotting/saving
    'plot': {...},
    'save': {...},
}
```

---

## Smart Setters (Auto-Reprogram)

The following parameters **auto-reprogram** the PulseBlaster:

| Parameter | Use Case | Fast Setter |
|-----------|----------|-------------|
| `pulse_duration` | Rabi oscillations | `_set_pulse_duration_direct()` |
| `tau` | T2, Hahn Echo | `_set_tau_direct()` |
| `t_AOM` | AOM timing optimization | `_set_t_AOM_direct()` |
| `ro_delay` | Readout delay optimization | `_set_ro_delay_direct()` |

**How it works**:
```python
# When you add sweep:
sweep.add_sweep('pulse_duration', values)

# Each value triggers:
1. Parameter update: pb.parameter_dict['seq']['pulse_duration'] = value
2. Auto-reprogram: pb._reprogram_sequence()
3. Ready for measurement!
```

---

## Advanced Usage

### Custom Measurement Function

```python
def my_custom_measurement(instruments, current_params):
    """Custom measurement with special processing."""
    pb = instruments['pb']
    ai_task = instruments['ai_task']

    # Your custom logic here
    pb.start_sequence()
    data = ai_task.read_daq(1000)

    # Custom processing
    processed = np.fft.fft(data)

    return {
        **current_params,
        'raw_data': data,
        'fft': processed,
        'timestamp': time.time()
    }

# Use it:
results = sweep.run(my_custom_measurement)
```

### Manual Sweep Control

```python
from parameter_system import ParameterSweep

# Create sweep
sweep = ParameterSweep(instruments)

# Add parameters (first = innermost)
sweep.add_sweep('frequency', freq_values)
sweep.add_sweep('power', power_values)

# Print info
sweep.print_sweep_info()

# Run with custom function
results = sweep.run(my_custom_measurement)
```

---

## Troubleshooting

### "Unknown parameter 'X'"

**Problem**: Parameter not registered with any instrument.

**Solution**: Check that the instrument class has `register_parameter()` call:

```python
class MyInstrument(Instrument):
    def _register_parameters(self):
        self.register_parameter("my_param", 0, (0, 100), "units")
```

### "No fast setter for parameter"

**Problem**: Parameter uses safe method instead of fast setter.

**Solution**: This is OK! It just means parameter uses `change_parameter()` instead of direct setter. Add direct setter if needed:

```python
def _set_my_param_direct(self, value: float):
    """Fast setter for my_param."""
    self.my_param = value
    # Update hardware directly
```

### PB not reprogramming

**Problem**: PB sequence not updating when parameters change.

**Solution**: Make sure you're using sweep parameters that have smart setters:
- `pulse_duration`
- `tau`
- `t_AOM`
- `ro_delay`

Or add your own in PBcontrol.py:
```python
def _set_MY_PARAM_direct(self, value: float):
    self.parameter_dict['seq']['MY_PARAM'] = value
    self._reprogram_sequence()
```

---

## Files to Check

| File | Purpose |
|------|---------|
| [esr_config_new.py](esr_config_new.py) | Example ESR config |
| [rabi_config_new.py](rabi_config_new.py) | Example Rabi config |
| [mainControl.py](mainControl.py) | Generic measurement functions |
| [parameter_system.py](parameter_system.py) | Simplified ParameterSweep |
| [PBcontrol.py](PBcontrol.py) | Smart setters with auto-reprogram |
| [REFINEMENTS_COMPLETE.md](REFINEMENTS_COMPLETE.md) | Complete documentation |

---

## Summary

**3 steps to run any experiment**:

1. **Create config** with `'sweep': {...}`
2. **Initialize instruments**
3. **Call** `run_sweep_experiment()`

**Key features**:
- ✅ Auto-reprogramming PulseBlaster
- ✅ Generic measurement functions
- ✅ Simpler config format
- ✅ One function for all experiments

**Ready to use!** 🚀
