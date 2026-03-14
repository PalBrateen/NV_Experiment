## Phase 3 & 4 Complete: Experiment Controllers and Data Management

## Summary

Phases 3 and 4 are complete! The experiment controller infrastructure and data management system have been successfully implemented, providing a unified class-based approach for both diode and camera-based experiments.

## Files Created

### 1. [save_manager.py](save_manager.py) - Phase 4
**Complete data management system** with:
- ✅ NumPy array storage (.npy) - backward compatible
- ✅ **Separate** YAML parameter files - human-readable, view without loading data
- ✅ JSON support as alternative format
- ✅ Automatic folder creation: `Saved_Data/YYYY-MM-DD/sequence_XXX/`
- ✅ NumPy type conversion for serialization
- ✅ Load/save convenience methods
- ✅ Experiment folder management utilities

**File Structure**:
```
Saved_Data/
└── 2024-12-15/
    ├── esr_001/
    │   ├── data_001.npy      # Raw data array
    │   └── params_001.yaml   # Parameters (SEPARATE file!)
    └── rabi_002/
        ├── data_002.npy
        └── params_002.yaml
```

**Key Features**:
- **Separate parameter files**: View parameters without loading large data arrays
- **Automatic type conversion**: NumPy → Python types for YAML/JSON
- **Smart folder management**: Date-based organization, auto-numbering
- **Convenience methods**: `save_experiment()`, `load_experiment()`, `get_next_folder_number()`

### 2. [experiment_controller.py](experiment_controller.py) - Phase 3
**Unified controller infrastructure** with:
- ✅ `ExperimentController` (ABC) - Base class for all experiments
- ✅ `DiodeExperimentController` - For APD/diode measurements
- ✅ `CameraExperimentController` - For camera-based measurements
- ✅ Common interface for both acquisition modes
- ✅ Instrument initialization and cleanup
- ✅ Integration with SaveManager

**Benefits**:
- Common interface → easier to maintain
- Future unification straightforward → just merge subclasses
- Helps now → cleaner code, reusable components
- Helps later → easy to add unified GUI

---

## Usage Examples

### SaveManager Usage

#### Basic Save/Load
```python
from save_manager import SaveManager
import numpy as np

# Create save manager
save_mgr = SaveManager()  # Uses default: ../Saved_Data

# Save data and parameters separately
data_array = np.random.rand(10, 101, 1000)  # (Nruns, Nscanpts, samples)
params_dict = {'mw': {'freq': 2.87e9, 'power': 8}, 'scan': {...}}

save_mgr.save_numpy(data_array, "001", "esr")
save_mgr.save_params_yaml(params_dict, "001", "esr")

# Or save both together
save_mgr.save_experiment(data_array, params_dict, "001", "esr", param_format='yaml')

# Load data and parameters
data, params = save_mgr.load_experiment("Saved_Data/2024-12-15/esr_001/data_001.npy")
# Automatically finds and loads params_001.yaml
```

#### Smart Folder Management
```python
# Get next available folder number
next_num = save_mgr.get_next_folder_number("esr")
# If esr_001 and esr_002 exist → returns "003"

# List all experiments for a date
experiments = save_mgr.list_experiments("2024-12-15")
# Returns: [Path('esr_001'), Path('esr_002'), Path('rabi_001'), ...]
```

#### Separate Parameter Files
```python
# Key feature: Parameters saved separately!
# ✅ Can view params without loading 10 GB data file
# ✅ Human-readable YAML format
# ✅ Version control friendly

# Load only parameters (fast!)
params = save_mgr.load_params("Saved_Data/2024-12-15/esr_001/params_001.yaml")
print(f"Experiment was: {params['seq']['sequence']}")
print(f"Frequency range: {params['scan']['values'][0]}")

# Load data only when needed (later)
data = save_mgr.load_numpy("Saved_Data/2024-12-15/esr_001/data_001.npy")
```

---

### ExperimentController Usage

#### Diode/APD Experiments
```python
from experiment_controller import DiodeExperimentController
from importlib import import_module

# Load configuration
expCfg = import_module('esr_config')
params_dict = expCfg.params_dict

# Create controller
controller = DiodeExperimentController(params_dict)

# Initialize instruments
trial_run = ['n', 'n']  # Both real hardware
controller.initialize_instruments(trial_run)

# Get instruments
instruments = controller.get_instruments_dict()
sg = instruments['sg']
pb = instruments['pb']
ai_task = instruments['ai_task']

# Use instruments...
# (Acquisition logic to be filled from mainControl_diode.py)

# Cleanup when done
controller.cleanup()
```

#### Camera Experiments
```python
from experiment_controller import CameraExperimentController

# Create controller
controller = CameraExperimentController(params_dict, camera_mode='level')

# Initialize instruments
trial_run = ['n', 'y']  # Real SG, simulated camera
controller.initialize_instruments(trial_run)

# Get instruments
instruments = controller.get_instruments_dict()
camera = instruments['camera']
sg = instruments['sg']
pb = instruments['pb']

# Use instruments...
# (Acquisition logic to be filled from mainControl_camera.py)

# Cleanup
controller.cleanup()
```

---

## Architecture Benefits

### 1. Unified Interface
Both diode and camera controllers share the same base class:
```python
# Common methods for both
controller.initialize_instruments(trial_run)
controller.get_instruments_dict()
controller.cleanup()
controller._acquire_single_point(param_values)  # Abstract
```

### 2. Mode-Specific Implementation
Each subclass handles its specific requirements:
- **DiodeExperimentController**: AI task, optional AO task
- **CameraExperimentController**: Camera, AO task, trigger modes

### 3. Future Unification
When ready, easy to merge:
```python
# Future unified controller (example)
class UnifiedExperimentController(ExperimentController):
    def __init__(self, params_dict, mode='diode'):
        self.mode = mode
        # Choose instruments based on mode

    def initialize_instruments(self, trial_run):
        if self.mode == 'diode':
            # Initialize diode instruments
        elif self.mode == 'camera':
            # Initialize camera instruments
```

---

## Integration with Parameter System

The controllers are designed to work seamlessly with the parameter system:

```python
from parameter_system import ParameterSweep
from experiment_controller import DiodeExperimentController

# Create controller
controller = DiodeExperimentController(params_dict)
controller.initialize_instruments(['n', 'n'])

# Get instruments dict
instruments = controller.get_instruments_dict()

# Create parameter sweep
sweep = ParameterSweep(instruments, 'esr')
sweep.add_sweep("frequency", freq_values)

# Set measurement function using controller
def measure(instruments, current_values):
    return controller._acquire_single_point(current_values)

sweep.set_measurement_function(measure)

# Run sweep
results = sweep.run()

# Save with SaveManager
save_mgr = controller.save_manager
folder_num = save_mgr.get_next_folder_number("esr")

# Convert results to numpy array (implementation needed)
data_array = ...  # Convert list of dicts to numpy array

save_mgr.save_experiment(data_array, params_dict, folder_num, "esr")
```

---

## Data Format: Multi-Parameter Sweeps

SaveManager supports multi-parameter sweeps with flattened arrays:

### Structure
```python
# 2D sweep: 51 freq × 3 power = 153 combinations
data_array.shape = (Nruns, 153, Nsamples*Nchannels)

# Parameters stored in separate YAML
params_dict = {
    'scan': {
        'names': ['frequency', 'power'],
        'values': [freq_array, power_array],
        'Nscanpts': [51, 3],
        'total_scanpts': 153,
        'param_indices': {  # Mapping flat index → (i_freq, i_power)
            0: (0, 0),
            1: (0, 1),
            ...
            152: (50, 2)
        }
    },
    'mw': {...},
    'seq': {...}
}
```

### Reconstruction for Plotting
```python
# Load experiment
data, params = save_mgr.load_experiment("path/to/data_001.npy")

# Extract sweep info
scan_info = params['scan']
param_shapes = scan_info['Nscanpts']  # [51, 3]

# Reshape for 2D plotting
data_2d = data[0, :, :].mean(axis=1)  # Average over samples
data_reshaped = data_2d.reshape(param_shapes)  # (51, 3)

# Now can plot as 2D
import matplotlib.pyplot as plt
plt.imshow(data_reshaped, aspect='auto')
plt.xlabel('Power index')
plt.ylabel('Frequency index')
```

---

## Backward Compatibility

### SaveManager
- ✅ Same NumPy .npy format as before
- ✅ Same folder structure
- ✅ Added: Separate parameter files (new feature)
- ✅ No breaking changes to existing save/load code

### ExperimentController
- ✅ Uses same instruments (SignalGenerator, PulseBlaster, DAQ)
- ✅ Same parameter dictionaries (params_dict)
- ✅ Same trial_run configuration
- ✅ Wraps existing functionality, doesn't replace it

---

## Next Steps: Integration

### Remaining Work

1. **Fill acquisition logic in DiodeExperimentController._acquire_single_point()**
   - Extract from mainControl_diode.py's acquire_data() and process_data()
   - Implement:
     - Sequence programming
     - AI task triggering
     - Data acquisition
     - Processing (contrast, signal extraction)

2. **Fill acquisition logic in CameraExperimentController._acquire_single_point()**
   - Extract from mainControl_camera.py
   - Implement:
     - Camera trigger setup
     - Image acquisition
     - Image processing

3. **Update mainControl_diode.py to use DiodeExperimentController**
   - Keep existing code working
   - Add option to use new controller
   - Gradual migration

4. **Update mainControl_camera.py to use CameraExperimentController**
   - Similar migration approach

---

## Example: Complete ESR Experiment Flow

```python
from experiment_controller import DiodeExperimentController
from parameter_system import ParameterSweep
from save_manager import SaveManager
from importlib import import_module
import numpy as np

# 1. Load configuration
expCfg = import_module('esr_config')
params_dict = expCfg.params_dict

# 2. Create controller
controller = DiodeExperimentController(params_dict)
controller.initialize_instruments(trial_run=['n', 'n'])

# 3. Create parameter sweep
instruments = controller.get_instruments_dict()
sweep = ParameterSweep(instruments, 'esr')

# 4. Add frequency sweep
freq_values = np.linspace(2.82e9, 2.92e9, 101)
sweep.add_sweep("frequency", freq_values)
# Output: ✓ Registered fast setter for frequency

# 5. Set measurement function
def measure(instruments, current_values):
    return controller._acquire_single_point(current_values)

sweep.set_measurement_function(measure)

# 6. Run sweep
print("\nRunning ESR sweep...")
results = sweep.run()

# 7. Convert to numpy array (simplified)
# TODO: Implement proper conversion from list of dicts
data_array = np.array([r['signal'] for r in results])

# 8. Save data and parameters
save_mgr = controller.save_manager
folder_num = save_mgr.get_next_folder_number("esr")
save_mgr.save_experiment(data_array, params_dict, folder_num, "esr")

# 9. Cleanup
controller.cleanup()

print(f"\n✔ Experiment complete!")
print(f"Data saved to: Saved_Data/{time.strftime('%Y-%m-%d')}/esr_{folder_num}/")
```

---

## Architecture Validation

Phases 3 & 4 follow the approved architecture:

### Phase 3 Checklist
- ✅ Base `ExperimentController` class created
- ✅ `DiodeExperimentController` implemented
- ✅ `CameraExperimentController` implemented
- ✅ Common interface for both modes
- ✅ Abstract methods for mode-specific logic
- ✅ Instrument initialization and cleanup
- ✅ Integration with SaveManager

### Phase 4 Checklist
- ✅ SaveManager with NumPy array storage
- ✅ Separate YAML parameter files
- ✅ JSON support as alternative
- ✅ Automatic folder creation
- ✅ NumPy type conversion
- ✅ Convenience methods (save_experiment, load_experiment)
- ✅ Folder management utilities

---

## Summary

**Phases 3 & 4 Status**: ✅ COMPLETE

### Files Created
1. **save_manager.py** (289 lines)
   - Complete data management system
   - NumPy + YAML/JSON support
   - Separate parameter files
   - Smart folder management

2. **experiment_controller.py** (343 lines)
   - Base ExperimentController class
   - DiodeExperimentController for APD/diode
   - CameraExperimentController for camera
   - Common interface, mode-specific implementation

### Key Achievements
- ✅ Unified controller infrastructure
- ✅ Separate parameter files (human-readable, version control friendly)
- ✅ Backward compatible data format
- ✅ Future unification enabled
- ✅ Integration with parameter system ready
- ✅ Clean separation of concerns

### Benefits
- **Now**: Cleaner code, reusable components, better organization
- **Later**: Easy to add unified GUI, straightforward to merge controllers
- **Always**: Separate parameter files make experiments more manageable

**Ready for integration with existing mainControl files!**
