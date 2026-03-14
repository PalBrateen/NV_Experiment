# NV Experiment Refactoring: Multi-Parameter System with Function Pointer Optimization

## Overview

Refactor the NV_Experiment codebase to implement a flexible multi-parameter sweep system using **function pointer optimization** for maximum performance, following the architecture defined in `ARCHITECTURE_GUIDE.md`, `complete_implementation.md`, and `Function pointer approach.md`.

**Note**: Qt GUI development is deferred to a later phase after core system is validated.

---

## Key Design Decisions (User Feedback Addressed)

### 1. Main Control Files: Unified Class Structure

**Problem**: `mainControl_diode.py` (voltage/count) and `mainControl_camera.py` (image) are separate entry points

**Solution**: Create base `ExperimentController` class with two subclasses:
- `DiodeExperimentController` - for APD/diode measurements
- `CameraExperimentController` - for image-based measurements

**Benefits**:
- ✅ Common interface for both modes
- ✅ Future unification is straightforward (just merge subclasses)
- ✅ Helps now: cleaner code, reusable components
- ✅ Helps later: easy to add unified GUI for both modes

### 2. Data Storage: NumPy Arrays (Not xarray)

**Decision**: Keep NumPy `.npy` files for data storage

**Rationale**:
- ✅ Backward compatible with existing code
- ✅ No new dependencies
- ✅ Fast and efficient
- ✅ xarray can be added later for post-processing/analysis if needed

**Structure**: Flattened multi-dimensional sweeps
```python
# 2D sweep: 51 freq × 3 power = 153 combinations
data_array.shape = (Nruns, 153, Nsamples*Nchannels)

# Metadata in params_dict saves mapping
params['scan']['param_indices'] = {0: (0,0), 1: (0,1), ..., 152: (50,2)}
```

### 3. Config Files: Python .py (Not YAML)

**Decision**: Keep `.py` config files, save params as YAML after experiment

**Rationale - Why .py configs**:
- ✅ Computed ranges: `freq_values = np.linspace(2.82*GHz, 2.92*GHz, 101)`
- ✅ Helper functions: `_nonlin_freq_sweep(f)`
- ✅ Unit imports: `from SGcontrol import GHz`
- ✅ Conditional logic for different experiments
- ✅ Users already familiar with this format

**YAML limitation**: Static only, cannot compute ranges

**Workflow**:
1. User edits `esr_config.py` (Python) → defines `params_dict`
2. Ranges auto-computed on module import
3. Experiment runs using `params_dict`
4. **SaveManager saves params_dict → YAML** (separate from data file!)

**Benefits**:
- ✅ Config: Full Python flexibility
- ✅ Save: Human-readable YAML with all computed values
- ✅ Separate files: Can view params without loading data
- ✅ Format: `data_001.npy` + `params_001.yaml`

---

## Current State

- **Procedural control**: mainControl_diode.py (717 lines) with global variables
- **Single-parameter scanning**: Only 1D sweeps supported
- **Instrument classes**: SignalGenerator, PulseBlaster, DAQcontrol (partially OOP)
- **Config files**: Dictionary-based (esr_config.py, rabi_config.py, t2_config.py)
- **Data storage**: 3D numpy arrays (Nruns, Nscanpts, samples)
- **No parameter abstraction**: Direct instrument method calls

## Target Architecture

- **Parameter Registry**: Automatic parameter → instrument mapping
- **N-dimensional sweeps**: Arbitrary multi-parameter combinations
- **Function Pointer Optimization**: Pre-register fast setters in `_fast_setters` dict (see "Function pointer approach.md")
- **Performance**: <10μs Python overhead per inner loop iteration
- **Flexible data storage**: xarray for labeled multi-dimensional data
- **YAML configs**: Human-readable configuration files

---

## Implementation Plan

### PHASE 1: Core Parameter System (Week 1-2)

#### 1.1 Create parameter_system.py

**File**: `d:\Brateen\NV_Experiment\parameter_system.py` (NEW)

**Implementation**:
```python
# Core classes from complete_implementation.md:
- ParameterRegistry: Global mapping (parameter_name → instrument_type)
- Parameter: Validation, ranges, units
- Instrument (ABC): Base class for all instruments
- ParameterSweep: N-dimensional nested loop executor with FUNCTION POINTER optimization
- ExperimentConfig: YAML-based configuration
```

**Key Features**:
- `ParameterRegistry.register(param, instrument_type)` - Automatic mapping
- `Instrument.register_parameter()` - Register parameters in subclasses
- `Instrument.change_parameter()` - Safe setting with validation
- `ParameterSweep.add_sweep(param, values)` - Automatic instrument lookup + fast setter registration
- `ParameterSweep.run()` - Returns xarray.Dataset

**CRITICAL: Function Pointer Approach for Performance**

Per "Function pointer approach.md", use function pointers instead of if-checks:

```python
class ParameterSweep:
    def __init__(self, instruments, sequence_type):
        self.instruments = instruments
        self.sequence_type = sequence_type
        self.sweep_params = []
        self._fast_setters = {}  # Dictionary of function pointers

    def add_sweep(self, parameter, values):
        # Find instrument and add sweep
        instrument = self._find_instrument_for_parameter(parameter)
        self.sweep_params.append((instrument.name, parameter, values))

        # Register fast setter if it exists (FUNCTION POINTER)
        direct_method = f'_set_{parameter}_direct'
        if hasattr(instrument, direct_method):
            fast_func = getattr(instrument, direct_method)
            self._fast_setters[(instrument.name, parameter)] = fast_func
            print(f"✓ Registered fast setter for {parameter}")

    def run(self):
        def set_parameter(inst_name, param, value, depth):
            is_innermost = (depth == len(self.sweep_params) - 1)

            if is_innermost:
                # Use function pointer (1 check + 1 dict lookup)
                fast_setter = self._fast_setters.get((inst_name, param))
                if fast_setter:
                    fast_setter(value)  # Direct call to stored function
                    return

            # Fallback to safe method with validation
            self.instruments[inst_name].change_parameter(param, value)
```

**Advantages over if-checks**:
- Old: 3+ if-checks every iteration (sequence_type? param name?)
- New: 1 check (innermost?) + 1 dict lookup → call function
- Cleaner, faster, more maintainable

**Data Structure Decision**: Use **NumPy arrays** (xarray deferred for later)
```python
# Multi-parameter data: flattened sweep dimension
# Example: 2D sweep (51 freq × 3 power = 153 total)
data_array.shape = (Nruns, 153, Nsamples*Nchannels)

# Parameter mapping stored in params_dict['scan']
params['scan']['param_indices'] = {0: (0,0), 1: (0,1), ..., 152: (50,2)}
```

**Testing**:
- Unit tests for Parameter validation
- Test _fast_setters dictionary population
- Mock instruments for ParameterSweep
- Test 1D, 2D, 3D sweeps with different orderings
- Benchmark function pointer vs if-check approach

---

### PHASE 2: Instrument Refactoring (Week 2-3)

#### 2.1 Refactor SignalGenerator

**File**: `d:\Brateen\NV_Experiment\SGcontrol.py`

**Changes**:
```python
class SignalGenerator(Instrument):
    def __init__(self, name: str = "sg1"):
        super().__init__(name, "SignalGenerator")
        # Existing VISA initialization
        self.init_sg(sg_addr)

    def _register_parameters(self):
        self.register_parameter("frequency", 2.87e9, (0, 20e9), "Hz")
        self.register_parameter("power", -2, (-130, 25), "dBm")
        self.register_parameter("phase", 0, (0, 360), "degrees")

    def _update_instrument(self, parameter_name, new_value):
        if parameter_name == "frequency":
            self.sg.write(f'FREQ{new_value}Hz')
        elif parameter_name == "power":
            self.sg.write(f'AMPR{new_value}dBm')

    # Performance optimization for ESR inner loop
    def _set_frequency_direct(self, value: float):
        """Fast frequency setting - no validation"""
        self.freq = value
        self.sg.write(f'FREQ{value}Hz')
```

**Backward Compatibility**: Keep existing methods, add deprecation warnings

#### 2.2 Refactor PulseBlaster

**File**: `d:\Brateen\NV_Experiment\PBcontrol.py`

**Changes**:
```python
class PulseBlaster(Instrument):
    def __init__(self, name: str = "pb1", parameter_dict: dict = None):
        super().__init__(name, "PulseBlaster")
        self.parameter_dict = parameter_dict or {}
        self.configure()

    def _register_parameters(self):
        self.register_parameter("pulse_duration", 100e-9, (0, 1), "s")
        self.register_parameter("tau", 1e-6, (0, 1e-3), "s")
        self.register_parameter("interval", 500e-9, (0, 1e-3), "s")

    def _update_instrument(self, parameter_name, new_value):
        # PB requires sequence reprogramming
        # Handled by sequence controller
        pass

    def _set_pulse_duration_direct(self, value: float):
        """Fast pulse duration for Rabi"""
        self.pulse_duration = value
        # Reprogram PB sequence
```

#### 2.3 Refactor DAQcontrol

**File**: `d:\Brateen\NV_Experiment\DAQcontrol.py`

**Changes**: Dual inheritance from Instrument and DAQTask
```python
class AnalogInputTask(Instrument, DAQTask):
    def __init__(self, name: str = "ai_task", dev="P6363", ...):
        Instrument.__init__(self, name, "AnalogInputTask")
        DAQTask.__init__(self)
        # Existing initialization

    def _register_parameters(self):
        self.register_parameter("sampling_rate", 2e6, (1, 2e6), "Sa/s")
        self.register_parameter("voltage_range_min", -10, (-10, 10), "V")
        self.register_parameter("voltage_range_max", 10, (-10, 10), "V")
```

---

### PHASE 3: Main Control Refactoring (Week 3-4)

#### 3.1 Create Base ExperimentController Class

**File**: `d:\Brateen\NV_Experiment\experiment_controller.py` (NEW)

**Purpose**: Base class for both diode and camera-based experiments

```python
from abc import ABC, abstractmethod

class ExperimentController(ABC):
    """Base class for experiment control - voltage/count and image-based"""

    def __init__(self, config: ExperimentConfig, acquisition_mode: str):
        self.config = config
        self.acquisition_mode = acquisition_mode  # 'diode' or 'camera'
        self.instruments = {}
        self.save_manager = SaveManager()

    @abstractmethod
    def initialize_instruments(self):
        """Initialize instruments - implemented by subclasses"""
        pass

    @abstractmethod
    def _acquire_single_point(self, param_values: dict) -> dict:
        """Acquire data at one parameter point - implemented by subclasses"""
        pass

    def run_experiment(self) -> np.ndarray:
        """Execute parameter sweep - returns numpy array"""
        def measure(instruments, current_values):
            return self._acquire_single_point(current_values)

        self.config.sweep.set_measurement_function(measure)
        results = self.config.sweep.run()  # Returns list of dicts

        # Convert to numpy array for compatibility
        return self._results_to_numpy(results)

    def _results_to_numpy(self, results: list) -> np.ndarray:
        """Convert results list to numpy array matching current format"""
        # Reconstruct (Nruns, Nscanpts, samples) structure
        pass

    def cleanup(self):
        """Close all instruments"""
        pass


class DiodeExperimentController(ExperimentController):
    """Controller for voltage/count-based experiments (mainControl_diode.py)"""

    def __init__(self, config: ExperimentConfig):
        super().__init__(config, 'diode')

    def initialize_instruments(self):
        """Initialize SG, PB, DAQ for diode experiments"""
        self.instruments['sg'] = SignalGenerator() if trial_run[0]=='n' else SignalGenerator_sim()
        self.instruments['pb'] = PulseBlaster(self.config.parameter_values)
        self.instruments['ai_task'] = AnalogInputTask(...)
        self.instruments['ao_task'] = AnalogOutputTask(...) if needed

    def _acquire_single_point(self, param_values: dict) -> dict:
        """Acquire diode/APD counts"""
        # Existing acquire_data() logic from mainControl_diode.py
        [cts, scan_time] = acquire_data(...)
        processed = process_data(...)
        return {**param_values, 'signal': ..., 'contrast': ..., 'scan_time': ...}


class CameraExperimentController(ExperimentController):
    """Controller for image-based experiments (mainControl_camera.py)"""

    def __init__(self, config: ExperimentConfig, camera_mode: str = 'level'):
        super().__init__(config, 'camera')
        self.camera_mode = camera_mode  # 'level', 'sync', 'timeseries'

    def initialize_instruments(self):
        """Initialize SG, PB, camera, DAQ AO for camera experiments"""
        self.instruments['sg'] = SignalGenerator() if trial_run[0]=='n' else SignalGenerator_sim()
        self.instruments['pb'] = PulseBlaster()
        self.instruments['ao_task'] = AnalogOutputTask()
        self.instruments['camera'] = CameraWorker(simulate=trial_run[1]=='y')
        # Configure camera trigger mode
        if 'level' in self.camera_mode:
            self.instruments['camera'].triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.LEVEL
        elif 'sync' in self.camera_mode:
            self.instruments['camera'].triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.SYNCREADOUT

    def _acquire_single_point(self, param_values: dict) -> dict:
        """Acquire camera image(s)"""
        # Existing acquire_data() logic from mainControl_camera.py
        images = acquire_images(...)
        processed = process_images(...)
        return {**param_values, 'images': ..., 'mean_signal': ..., 'contrast': ...}
```

**Benefits**:
- Common interface for both acquisition modes
- Future unification is straightforward
- Shared functionality in base class
- Each subclass handles mode-specific logic

#### 3.2 Convert mainControl_diode.py to Class

**File**: `d:\Brateen\NV_Experiment\mainControl_diode.py`

**Approach**: Convert to use DiodeExperimentController

```python
from experiment_controller import DiodeExperimentController
from esr_config import create_esr_config

if __name__ == "__main__":
    # Load config (still .py file)
    expCfgFile = 'esr_config'
    expCfg = import_module(expCfgFile)

    # Create config from dict
    config = ExperimentConfig.from_legacy_dict(expCfg.params_dict, {})

    # Create controller
    controller = DiodeExperimentController(config)
    controller.initialize_instruments()

    # Display parameters dialog
    if display_parameters == 'yes':
        results = controller.run_experiment()  # Returns numpy array

        # Save (numpy + YAML)
        if savefile_yn == 'yes':
            controller.save_manager.save_numpy(results, params)
            controller.save_manager.save_params_yaml(params)

    controller.cleanup()
```

#### 3.3 Convert mainControl_camera.py to Class

**File**: `d:\Brateen\NV_Experiment\mainControl_camera.py`

**Similar conversion** using CameraExperimentController

---

### PHASE 4: Data Management (Week 4)

#### 4.1 Implement SaveManager with NumPy Arrays

**File**: `d:\Brateen\NV_Experiment\save_manager.py` (NEW)

**Decision**: Keep NumPy arrays for now (xarray can be added later for post-processing)

**Features**:
- NumPy (.npy) for data storage - maintains current format
- YAML (.yaml) for parameters - separate file, human-readable
- JSON (.json) for parameters - alternative format
- Automatic folder creation with date/sequence structure

```python
class SaveManager:
    def __init__(self, base_dir: Path = None):
        self.base_dir = base_dir or Path.cwd().parent / "Saved_Data"

    def save_numpy(self, data: np.ndarray, folder_number: str, sequence: str):
        """Save numpy data array"""
        file_dir = self._create_folder(sequence, folder_number)
        np.save(file_dir / f"data_{folder_number}.npy", data, allow_pickle=False)
        print(f"✔ Data saved to {file_dir / f'data_{folder_number}.npy'}")

    def save_params_yaml(self, params: dict, folder_number: str, sequence: str):
        """Save parameters as YAML (separate file)"""
        file_dir = self._create_folder(sequence, folder_number)
        yaml_file = file_dir / f"params_{folder_number}.yaml"

        # Convert numpy types to Python types for serialization
        params_serializable = self._convert_numpy_to_python(params)

        with open(yaml_file, 'w') as f:
            yaml.dump(params_serializable, f, indent=4, default_flow_style=False)
        print(f"✔ Parameters saved to {yaml_file}")

    def save_params_json(self, params: dict, folder_number: str, sequence: str):
        """Save parameters as JSON (alternative to YAML)"""
        file_dir = self._create_folder(sequence, folder_number)
        json_file = file_dir / f"params_{folder_number}.json"

        params_serializable = self._convert_numpy_to_python(params)

        with open(json_file, 'w') as f:
            json.dump(params_serializable, f, indent=4)
        print(f"✔ Parameters saved to {json_file}")

    def _convert_numpy_to_python(self, obj):
        """Recursively convert numpy types to Python types"""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.generic):
            return obj.item()
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_to_python(value)
                    for key, value in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_numpy_to_python(item) for item in obj]
        elif isinstance(obj, Path):
            return str(obj)
        else:
            return obj

    def _create_folder(self, sequence: str, folder_number: str) -> Path:
        """Create save directory with date/sequence structure"""
        date_dir = self.base_dir / time.strftime("%Y-%m-%d", time.localtime())
        file_dir = date_dir / f"{sequence}_{folder_number}"

        if not date_dir.exists():
            date_dir.mkdir(parents=True, exist_ok=True)
        if not file_dir.exists():
            file_dir.mkdir(exist_ok=True)

        return file_dir

    def load_numpy(self, filepath: Path) -> np.ndarray:
        """Load numpy data"""
        return np.load(filepath, allow_pickle=False)

    def load_params(self, filepath: Path) -> dict:
        """Load parameters from YAML or JSON"""
        if filepath.suffix == '.yaml':
            with open(filepath, 'r') as f:
                return yaml.safe_load(f)
        elif filepath.suffix == '.json':
            with open(filepath, 'r') as f:
                return json.load(f)
```

#### 4.2 NumPy Data Structure for Multi-Parameter Sweeps

**Current**: 3D array `(Nruns, Nscanpts, samples)`

**Multi-parameter**: N-D array with flattened sweep dimensions

```python
# Example: 2D sweep (frequency × power)
# freq_points = 51, power_points = 11
# Flattened sweep dimension: 51 × 11 = 561

data_array.shape = (Nruns, 561, Nsamples*Nchannels)

# Store parameter mapping in params dict
params['scan'] = {
    'names': ['frequency', 'power'],
    'values': [freq_array, power_array],
    'Nscanpts': [51, 11],  # Per-parameter counts
    'total_scanpts': 561,   # Total combinations
    'param_indices': {      # Mapping flat index -> (i_freq, i_power)
        0: (0, 0),
        1: (0, 1),
        ...
        560: (50, 10)
    }
}
```

**Reconstruction for plotting**:
```python
def reshape_for_plotting(data_1d, param_shapes):
    """Reshape (561,) → (51, 11) for 2D plotting"""
    return data_1d.reshape(param_shapes)
```

**Benefits**:
- Backward compatible with existing code (still numpy)
- No new dependencies (xarray not needed yet)
- Metadata in separate YAML file (human-readable)
- Easy to load: `data = np.load(...)`, `params = yaml.safe_load(...)`

---

### PHASE 5: Config File Strategy (Week 5)

#### 5.1 Decision: Keep Python .py Config Files

**Rationale**: Python configs allow:
- ✅ Computed ranges (`np.linspace`, custom functions)
- ✅ Helper functions (e.g., `_nonlin_freq_sweep`)
- ✅ Unit imports (`from SGcontrol import GHz`)
- ✅ Conditional logic
- ✅ Familiar to users

**YAML limitation**: Static only, no computation

#### 5.2 Enhanced params_dict Structure

**File**: `d:\Brateen\NV_Experiment\esr_config.py` (example)

```python
from SGcontrol import GHz
import numpy as np

# Computed frequency range (auto-calculated on load!)
freq_values = np.linspace(2.82*GHz, 2.92*GHz, 101)

t_AOM = 1e-3
MW_power = -2
Nsamples = 300000

# Enhanced params_dict for multi-parameter support
params_dict = {
    'scan': {
        'names': ['frequency'],          # List (multi-param ready)
        'values': [freq_values],         # List of arrays
        'Nscanpts': [101],               # Per-parameter
        'total_scanpts': 101,            # Total combinations
        'Nruns': 1,
    },
    'mw': {'power': MW_power, 'freq': freq_values},
    'seq': {
        'Nsamples': Nsamples,
        'sequence': 'esr_seq',
        'args': [t_AOM],
    },
    # ... rest of params
}
```

**Multi-parameter example** (frequency × power):
```python
freq_values = np.linspace(2.82*GHz, 2.92*GHz, 51)
power_values = np.array([5, 8, 10])

params_dict['scan'] = {
    'names': ['frequency', 'power'],     # TWO params
    'values': [freq_values, power_values],
    'Nscanpts': [51, 3],
    'total_scanpts': 153,  # 51 × 3
}
```

#### 5.3 Workflow: .py Config → params_dict → YAML Save

1. User edits `esr_config.py` (Python)
2. Ranges computed on import
3. Experiment uses `params_dict`
4. SaveManager saves params_dict → YAML (after experiment)

**Benefits**:
- Config: Full Python flexibility
- Save: Human-readable YAML with computed values

---

### PHASE 6: Performance Optimization & Validation (Week 6)

#### 6.1 Verify Function Pointer Performance

**Already implemented in Phase 1** - Function pointer approach in `ParameterSweep.add_sweep()` and `run()`

**Performance target**: <1% Python overhead vs hardware communication time

**Verification**:
- Benchmark function pointer lookup vs if-check cascade
- Measure overhead: function call + dict lookup should be <10 microseconds
- Profile with `cProfile` to ensure no bottlenecks

---

### PHASE 7: Testing & Validation (Week 7)

#### 7.1 Unit Tests

**File**: `d:\Brateen\NV_Experiment\tests\test_parameter_system.py` (NEW)

- Parameter validation tests
- ParameterRegistry tests
- Instrument mock tests
- ParameterSweep execution tests

#### 7.2 Integration Tests

**File**: `d:\Brateen\NV_Experiment\tests\test_esr_experiment.py` (NEW)

- 1D ESR sweep test
- 2D multi-parameter sweep test
- Data structure validation
- Backward compatibility tests

#### 7.3 Performance Benchmarks

**File**: `d:\Brateen\NV_Experiment\tests\benchmark_inner_loop.py` (NEW)

- Compare change_parameter vs direct methods
- Measure overhead: target <1ms per point
- Profile sweep execution

---

## Critical Files Summary

### New Files to Create
1. `parameter_system.py` - Core abstraction layer with function pointer optimization
2. `experiment_controller.py` - Base ExperimentController + DiodeExperimentController + CameraExperimentController
3. `save_manager.py` - Unified save/load system (NumPy + YAML)
4. `tests/test_parameter_system.py` - Unit tests
5. `tests/test_esr_experiment.py` - Integration tests
6. `tests/benchmark_inner_loop.py` - Performance tests (function pointer vs if-checks)

### Files to Modify
1. `SGcontrol.py` - Inherit from Instrument, add _register_parameters(), add _set_frequency_direct()
2. `PBcontrol.py` - Inherit from Instrument, add parameter abstraction
3. `DAQcontrol.py` - Dual inheritance (Instrument + DAQTask)
4. `mainControl_diode.py` - Convert to use DiodeExperimentController
5. `mainControl_camera.py` - Convert to use CameraExperimentController
6. `esr_config.py` - Enhance params_dict for multi-parameter support
7. `rabi_config.py` - Enhance params_dict for multi-parameter support
8. `t2_config.py` - Enhance params_dict for multi-parameter support

### Files to Keep Unchanged (Initially)
- `sequencecontrol.py` - Sequence logic remains same
- `control_daq_sequences.py` - Pulse sequences unchanged
- `connectionConfig.py` - Hardware definitions unchanged
- `Camcontrol.py` - Camera control unchanged (optional integration later)

---

## Migration Strategy

### Backward Compatibility
1. **Dual mode**: Old procedural code still works during transition
2. **Config compatibility**: Support both dict and ExperimentConfig
3. **Data format**: SaveManager supports both .npy and .nc
4. **Gradual migration**: 6-month deprecation timeline

### Deprecation Timeline
- **Month 1-2**: Phases 1-3 (core system, no breaking changes)
- **Month 3**: Phases 4-5 (config updates, data management)
- **Month 4**: Phases 6-7 (optimization, testing)
- **Month 5**: Soft deprecation warnings
- **Month 6**: Hard deprecation

---

## Risk Mitigation

1. **Hardware communication breakage**
   - Mitigation: Extensive simulation mode testing
   - Fallback: Keep old code path available

2. **Performance regression**
   - Mitigation: Benchmark at every phase
   - Target: <1% overhead

3. **Data loss during migration**
   - Mitigation: Dual save (old + new format)
   - Validation: Checksum verification

4. **GUI thread safety**
   - Mitigation: QThread for experiments
   - Testing: Stress tests with rapid start/stop

---

## Success Criteria

✅ Multi-parameter sweeps work (2D, 3D, arbitrary N-D)
✅ Performance: Function pointer approach - inner loop <10μs overhead per point
✅ Data: NumPy arrays with flattened multi-dimensional sweeps
✅ Params: Separate YAML save files (human-readable, view without loading data)
✅ Config: Python .py files with computed ranges
✅ Main control: Both diode and camera modes use ExperimentController classes
✅ Backward compatibility: Existing experiments still run
✅ Tests: >80% code coverage
✅ Function pointer optimization: Fast setters registered automatically in add_sweep()

---

## Next Steps

1. **User confirmation**: Approve this plan
2. **Phase 1 execution**: Create parameter_system.py
3. **Incremental testing**: Test each phase before proceeding
4. **User feedback**: Review at end of each phase
