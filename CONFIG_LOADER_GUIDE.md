# Configuration Loader Integration Guide

**Date**: December 16, 2024

## Overview

The config loader system now supports **both Python and YAML configuration files** seamlessly. The system auto-detects the format and loads configurations appropriately.

---

## Quick Start

### Load Any Config (Auto-Detect Format)

```python
from mainControl import load_experiment_config

# Auto-detects .yaml, .yml, or .py
params = load_experiment_config('configs/esr_config')

# Or specify full path with extension
params = load_experiment_config('configs/esr_config.yaml')
params = load_experiment_config('esr_config_new.py')
```

### Run Experiment from YAML

```python
from mainControl import run_sweep_experiment, load_experiment_config

# Load YAML config
params_dict = load_experiment_config('configs/esr_config.yaml')

# Initialize instruments
instruments = {...}

# Run experiment
results = run_sweep_experiment(params_dict, instruments, mode='diode')
```

---

## Configuration Formats

### 1. YAML Format (Recommended for GUI)

**File**: `configs/esr_config.yaml`

```yaml
experiment:
  name: ESR
  description: Electron Spin Resonance frequency sweep
  sequence: esr_seq
  mode: diode

signal_generator:
  power_dBm: 8
  modulation: pulse

pulse_blaster:
  t_AOM_us: 8000
  clk_MHz: 100

daq:
  Nsamples: 8
  sampling_rate_Hz: 2000000

acquisition:
  Nruns: 30

sweep:
  - parameter: frequency
    instrument: signal_generator
    start: 2.82e9
    stop: 2.92e9
    step: 1.0e6
    units: Hz

output:
  filename_prefix: ESR
```

**Benefits**:
- ✅ Human-readable
- ✅ Easy to edit (no Python knowledge required)
- ✅ GUI-friendly (easy to generate from UI)
- ✅ Version control friendly
- ✅ Cross-platform

**Limitations**:
- ❌ No computed values
- ❌ No helper functions
- ❌ Must pre-calculate sweep arrays

### 2. Python Format (Recommended for Complex Configs)

**File**: `esr_config_new.py`

```python
from spinapi import ns, us, ms
from SGcontrol import GHz
import numpy as np

params_dict = {
    'seq': {
        'sequence': 'esr_seq',
        't_AOM': 8000 * ms,
        'Nsamples': 8,
    },
    'mw': {
        'power': 8,
    },
    'scan': {
        'Nruns': 30,
    },
    'sweep': {
        'frequency': np.linspace(2.82*GHz, 2.92*GHz, 101)  # Computed!
    }
}
```

**Benefits**:
- ✅ Computed values (`np.linspace`, formulas)
- ✅ Helper functions (custom sweep generation)
- ✅ Unit imports (`GHz`, `MHz`, `ns`)
- ✅ Conditional logic
- ✅ Full Python flexibility

**Limitations**:
- ❌ Requires Python knowledge
- ❌ Less GUI-friendly

---

## File Structure

### Created Files

```
NV_Experiment/
├── configs/                          # YAML config directory
│   ├── esr_config.yaml              # ESR YAML config
│   └── rabi_config.yaml             # Rabi YAML config
│
├── config_loader.py                 # YAML loader (updated)
├── mainControl.py                   # Updated with load_experiment_config()
├── example_yaml_config.py           # Complete YAML usage examples
│
├── esr_config_new.py                # Python config (new format)
├── rabi_config_new.py               # Python config (new format)
└── CONFIG_LOADER_GUIDE.md           # This file
```

---

## Usage Examples

### Example 1: Load YAML Config

```python
from mainControl import load_experiment_config

# Method 1: Auto-detect format
params = load_experiment_config('configs/esr_config')  # Finds .yaml

# Method 2: Explicit path
params = load_experiment_config('configs/esr_config.yaml')

# Verify
print(f"Loaded: {params['experiment']['name']}")
print(f"Sequence: {params['seq']['sequence']}")
print(f"Sweep: {list(params['sweep'].keys())}")
```

### Example 2: Run Experiment from YAML

```python
from mainControl import run_sweep_experiment, load_experiment_config
from SGcontrol import SignalGenerator
from PBcontrol import PulseBlaster
from DAQcontrol import AnalogInputTask

# Load YAML config
params_dict = load_experiment_config('configs/rabi_config.yaml')

# Initialize instruments
instruments = {
    'sg': SignalGenerator(name="sg1"),
    'pb': PulseBlaster(params_dict, name="pb1"),
    'ai_task': AnalogInputTask(..., name="ai_task")
}

# Run experiment - PB auto-reprograms!
results = run_sweep_experiment(params_dict, instruments, mode='diode')

print(f"✓ Complete: {len(results)} measurements")
```

### Example 3: Load Python Config

```python
# Same function works for Python configs!
params = load_experiment_config('esr_config_new')  # Finds .py
params = load_experiment_config('esr_config_new.py')  # Explicit

# Or old way (still works)
from esr_config_new import params_dict
```

### Example 4: Mixed Workflow (YAML + Python)

```python
# Load base YAML config
params = load_experiment_config('configs/esr_config.yaml')

# Modify in Python (e.g., based on previous results)
import numpy as np
params['sweep']['frequency'] = np.linspace(2.86e9, 2.88e9, 51)
params['scan']['Nruns'] = 50

# Run experiment with modified config
results = run_sweep_experiment(params, instruments, mode='diode')
```

### Example 5: Create YAML Config Programmatically

```python
from config_loader import create_esr_config, save_config

# Create config with custom parameters
config = create_esr_config(
    freq_start=2.85e9,
    freq_stop=2.89e9,
    freq_step=0.5e6,
    power_dBm=10,
    Nsamples=16
)

# Save to YAML
save_config(config, 'configs/my_custom_esr.yaml')

# Load and use
params = load_experiment_config('configs/my_custom_esr.yaml')
```

---

## YAML Config Structure

### Complete Structure

```yaml
# Experiment metadata
experiment:
  name: <string>
  description: <string>
  sequence: <string>  # e.g., 'esr_seq', 'rabi_seq'
  mode: <string>      # 'diode' or 'camera'

# Signal generator parameters
signal_generator:
  power_dBm: <float>
  frequency_Hz: <float>  # Optional, omit if sweeping
  modulation: <string>   # 'pulse', 'cw', etc.

# PulseBlaster timing parameters
pulse_blaster:
  t_AOM_us: <float>       # AOM pulse duration (microseconds)
  ro_delay_ns: <float>    # Readout delay (nanoseconds)
  AOM_lag_ns: <float>     # AOM lag (nanoseconds)
  MW_lag_ns: <float>      # MW lag (nanoseconds)
  clk_MHz: <float>        # Clock frequency (MHz)

# DAQ parameters
daq:
  Nsamples: <int>
  sampling_rate_Hz: <float>
  voltage_range_V: [<float>, <float>]

# Camera parameters (optional, for camera mode)
camera:
  exposure_ms: <float>
  trigger_mode: <string>  # 'level', 'sync', 'timeseries'
  roi: [<int>, <int>, <int>, <int>]  # [x, y, width, height]

# Acquisition settings
acquisition:
  Nruns: <int>
  contrast_mode: <string>
  randomize: <bool>

# Sweep configuration (ORDER PRESERVED!)
sweep:
  - parameter: <string>     # Parameter name
    instrument: <string>    # Instrument name
    start: <float>          # Start value
    stop: <float>           # Stop value
    step: <float>           # Step size
    units: <string>         # Units for reference

  # Add more sweep parameters for multi-parameter sweeps
  - parameter: <string>
    # ...

# Output settings
output:
  save_path: <string>
  filename_prefix: <string>
  auto_save_interval: <int>

  plot:
    x_label: <string>
    x_scale: <float>
    y_label: <string>
    live_update: <bool>
```

### Minimal YAML Config

```yaml
experiment:
  name: ESR
  sequence: esr_seq
  mode: diode

signal_generator:
  power_dBm: 8

pulse_blaster:
  t_AOM_us: 8000

daq:
  Nsamples: 8

acquisition:
  Nruns: 30

sweep:
  - parameter: frequency
    start: 2.82e9
    stop: 2.92e9
    step: 1e6
```

---

## Unit Conversions

The config loader automatically converts units to SI base units:

| YAML Field | Units | Converted To | Example |
|------------|-------|--------------|---------|
| `t_AOM_us` | μs | seconds (s) | 8000 μs → 0.008 s |
| `ro_delay_ns` | ns | seconds (s) | 1800 ns → 1.8e-6 s |
| `AOM_lag_ns` | ns | seconds (s) | 800 ns → 8e-7 s |
| `MW_lag_ns` | ns | seconds (s) | 150 ns → 1.5e-7 s |
| `clk_MHz` | MHz | MHz (kept) | 100 MHz → 100 |
| `power_dBm` | dBm | dBm (kept) | 8 dBm → 8 |
| `frequency_Hz` | Hz | Hz (kept) | 2.87e9 Hz → 2.87e9 |

---

## Integration with Existing System

### How It Works

1. **`load_experiment_config()`** in [mainControl.py](mainControl.py:111-176):
   - Auto-detects file format (.yaml, .yml, or .py)
   - Calls appropriate loader
   - Returns standardized `params_dict`

2. **`load_config()`** in [config_loader.py](config_loader.py:14-27):
   - Loads YAML file
   - Calls `process_config()` to convert to params_dict format

3. **`process_config()`** in [config_loader.py](config_loader.py:30-145):
   - Converts YAML structure to params_dict
   - Handles unit conversions
   - Generates sweep arrays from start/stop/step

4. **`run_sweep_experiment()`** in [mainControl.py](mainControl.py:168-204):
   - Works with params_dict from any source
   - Doesn't care about config format!

### Backward Compatibility

✅ **All existing code works unchanged:**

```python
# Old way (still works)
from esr_config_new import params_dict
results = run_sweep_experiment(params_dict, instruments, 'diode')

# New way (also works)
params_dict = load_experiment_config('configs/esr_config.yaml')
results = run_sweep_experiment(params_dict, instruments, 'diode')
```

---

## Command Line Usage

```bash
# With YAML config
python mainControl.py --config configs/esr_config.yaml

# With Python config (auto-detect)
python mainControl.py --config esr_config_new

# Explicit Python config
python mainControl.py --config esr_config_new.py
```

---

## Creating YAML Configs

### Method 1: Manual Creation

Create a `.yaml` file following the structure above.

### Method 2: Programmatic Creation

```python
from config_loader import create_esr_config, create_rabi_config, save_config

# ESR config
esr_config = create_esr_config(
    freq_start=2.82e9,
    freq_stop=2.92e9,
    freq_step=1e6,
    power_dBm=8,
    Nsamples=8,
    Nruns=30
)
save_config(esr_config, 'configs/my_esr.yaml')

# Rabi config
rabi_config = create_rabi_config(
    duration_start=10e-9,
    duration_stop=500e-9,
    duration_step=2e-9,
    mw_freq=2.87e9,
    power_dBm=8
)
save_config(rabi_config, 'configs/my_rabi.yaml')
```

### Method 3: From GUI

GUI can generate YAML configs easily:

```python
def gui_save_config(gui_params):
    """Save GUI parameters to YAML."""
    config = {
        'experiment': {
            'name': gui_params['exp_name'],
            'sequence': gui_params['sequence'],
            'mode': gui_params['mode']
        },
        'sweep': [
            {
                'parameter': gui_params['sweep_param'],
                'start': gui_params['start_val'],
                'stop': gui_params['stop_val'],
                'step': gui_params['step_val']
            }
        ],
        # ... more fields
    }

    save_config(config, f"configs/{gui_params['exp_name']}.yaml")
```

---

## Best Practices

### When to Use YAML

✅ **Use YAML when:**
- GUI is generating configs
- Sharing configs with non-Python users
- Simple, straightforward experiments
- Version controlling experiment parameters
- Need human-readable, editable configs

### When to Use Python

✅ **Use Python when:**
- Need computed sweep values
- Using helper functions for sweep generation
- Complex conditional logic
- Advanced users comfortable with Python
- Need unit imports (GHz, MHz, ns, etc.)

### Hybrid Approach

**Best of both worlds:**

1. **Save base config as YAML** (GUI-generated)
2. **Load and modify in Python** (for advanced tweaks)
3. **Run experiment** with modified params

```python
# Load YAML base
params = load_experiment_config('configs/base_esr.yaml')

# Tweak in Python
params['sweep']['frequency'] = custom_frequency_array()
params['scan']['Nruns'] = 100

# Run
results = run_sweep_experiment(params, instruments, 'diode')
```

---

## Troubleshooting

### "Config file not found"

**Problem**: `FileNotFoundError`

**Solution**: Check file path and extension
```python
# Try with full path
params = load_experiment_config('configs/esr_config.yaml')

# Or use Path
from pathlib import Path
config_path = Path('configs') / 'esr_config.yaml'
params = load_experiment_config(str(config_path))
```

### "Unknown config format"

**Problem**: Unrecognized extension

**Solution**: Use .py, .yaml, or .yml
```python
# Not supported
params = load_experiment_config('config.txt')  # ❌

# Supported
params = load_experiment_config('config.yaml')  # ✅
params = load_experiment_config('config.py')    # ✅
```

### YAML parsing errors

**Problem**: Invalid YAML syntax

**Solution**: Validate YAML online or with yamllint
```bash
# Install yamllint
pip install yamllint

# Validate file
yamllint configs/esr_config.yaml
```

---

## Summary

**Configuration loader integration complete!**

### Features

- ✅ **Dual format support**: Python (.py) and YAML (.yaml/.yml)
- ✅ **Auto-detection**: `load_experiment_config()` detects format automatically
- ✅ **Unit conversion**: YAML values converted to SI units
- ✅ **Sweep generation**: Auto-generates arrays from start/stop/step
- ✅ **GUI-friendly**: Easy for GUIs to generate YAML configs
- ✅ **Backward compatible**: All existing code works unchanged

### Files Created/Updated

| File | Status | Purpose |
|------|--------|---------|
| [configs/esr_config.yaml](configs/esr_config.yaml) | ✅ Created | YAML ESR config |
| [configs/rabi_config.yaml](configs/rabi_config.yaml) | ✅ Created | YAML Rabi config |
| [config_loader.py](config_loader.py) | ✅ Updated | Fixed unit conversions, added pb params |
| [mainControl.py](mainControl.py) | ✅ Updated | Added `load_experiment_config()` |
| [example_yaml_config.py](example_yaml_config.py) | ✅ Created | Complete usage examples |
| [CONFIG_LOADER_GUIDE.md](CONFIG_LOADER_GUIDE.md) | ✅ Created | This guide |

### Next Steps

1. **Test YAML configs** with hardware
2. **GUI integration** - generate YAML from GUI inputs
3. **Extend YAML schemas** - add validation
4. **Create more example configs** - T2, Hahn Echo, etc.

**The system is ready for both Python and YAML configs!** 🎉
