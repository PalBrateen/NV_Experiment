Implementation Guide for Claude Code**
--------------------------------------------

Here's a structured guide:

### **PART A: Update `parameter_system.py`**

**A1. Simplify `ParameterSweep` class:**

```python
class ParameterSweep:
    """
    Automatic multi-parameter sweep with smart setters.
    
    - First added parameter = innermost loop (fastest changing)
    - Each instrument handles its own update logic via _set_X_direct() methods
    - Falls back to change_parameter() if no direct setter exists
    """
    
    def __init__(self, instruments: Dict[str, Instrument]):
        """
        Args:
            instruments: Dict of {name: Instrument} instances
        """
        self.instruments = instruments
        self.sweep_params = []      # [(inst_name, param_name, values), ...]
        self._setters = {}          # {param_name: setter_function}
    
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
```
***
### **PART B: Update `PBcontrol.py`**

**B1. Add smart setters that auto-reprogram:**

Add these methods to the `PulseBlaster` class:

```python
class PulseBlaster(Instrument):
    
    def __init__(self, parameter_dict, name: str = "pb1"):
        # Store reference to parameter_dict (will be modified by setters)
        self.parameter_dict = parameter_dict
        self._instr_type = 'diode'  # Set by experiment controller
        
        # ... rest of existing __init__ ...
        
        super().__init__(name)
    
    def _register_parameters(self):
        """Register PB timing parameters."""
        seq = self.parameter_dict.get('seq', {})
        self.register_parameter("pulse_duration", seq.get('pulse_duration', 100e-9), (0, 1), "s")
        self.register_parameter("tau", seq.get('tau', 1e-6), (0, 1e-3), "s")
        self.register_parameter("t_AOM", seq.get('t_AOM', 5e-6), (0, 1), "s")
        self.register_parameter("ro_delay", seq.get('ro_delay', 1e-6), (0, 1e-3), "s")
    
    def set_instr_type(self, instr_type: str):
        """Set instrument type ('diode' or 'cam_levelm', etc.)"""
        self._instr_type = instr_type
    
    def _get_sequence_name(self) -> str:
        """Get sequence name from parameter_dict."""
        return self.parameter_dict.get('seq', {}).get('sequence', '')
    
    def _build_sequence_args(self) -> list:
        """Build sequence args from current parameter_dict."""
        seq = self.parameter_dict.get('seq', {})
        sequence = self._get_sequence_name()
        
        if 'rabi' in sequence:
            return [seq.get('t_AOM'), seq.get('ro_delay'), 
                    seq.get('AOM_lag'), seq.get('MW_lag')]
        elif 'esr' in sequence:
            return [seq.get('t_AOM')]
        elif 't1' in sequence:
            return [seq.get('t_AOM'), seq.get('ro_delay'), seq.get('AOM_lag')]
        elif 't2' in sequence or 'hahn' in sequence:
            return [seq.get('t_AOM'), seq.get('ro_delay'), 
                    seq.get('AOM_lag'), seq.get('MW_lag'), seq.get('tau')]
        else:
            # Default
            return [seq.get('t_AOM')]
    
    def _reprogram_sequence(self):
        """Reprogram PB with current parameter_dict values."""
        sequence = self._get_sequence_name()
        seq_args = self._build_sequence_args()
        
        # Add PB channels
        pb_channels = self.parameter_dict.get('seq', {}).get('PBchannels', {})
        seq_args.append(pb_channels)
        
        # Program PB
        _, instr_list = PulseBlaster.PB_program(self._instr_type, sequence, seq_args)
        
        # Load to hardware
        if self._instr_type == 'diode':
            self.run_sequence_for_diode(instr_list)
        # Add other instr_types as needed
    
    # ================================================================
    # SMART SETTERS - Update value + auto-reprogram
    # ================================================================
    
    def _set_pulse_duration_direct(self, value: float):
        """Set pulse duration and reprogram PB. Used for Rabi."""
        self.parameter_dict['seq']['pulse_duration'] = value
        self._reprogram_sequence()
    
    def _set_tau_direct(self, value: float):
        """Set tau and reprogram PB. Used for T2/Hahn Echo."""
        self.parameter_dict['seq']['tau'] = value
        self._reprogram_sequence()
    
    def _set_t_AOM_direct(self, value: float):
        """Set AOM duration and reprogram PB."""
        self.parameter_dict['seq']['t_AOM'] = value
        self._reprogram_sequence()
    
    def _set_ro_delay_direct(self, value: float):
        """Set readout delay and reprogram PB."""
        self.parameter_dict['seq']['ro_delay'] = value
        self._reprogram_sequence()
    
    def _update_instrument(self, parameter_name: str, new_value):
        """Fallback update method (used by change_parameter)."""
        if parameter_name in ['pulse_duration', 'tau', 't_AOM', 'ro_delay']:
            self.parameter_dict['seq'][parameter_name] = new_value
            self._reprogram_sequence()
```
***
### **PART C: Update `SGcontrol.py`**

**Already has direct setters!** Just verify these exist:

```python

def _set_frequency_direct(self, value: float):
    """Fast frequency setting."""
    self.freq = value
    self.sg.write(f'FREQ{value}Hz')

def _set_power_direct(self, value: float):
    """Fast power setting."""
    self.rfamp = value
    self.sg.write(f'AMPR{value}dBm')
```
* * *

### **PART D: Update `mainControl.py`**

**D1. Add generic measurement functions:**

```python
def diode_measurement(instruments, current_params):
    """
    Generic diode/APD measurement.
    Parameters already set by sweep (including PB reprogrammed if needed).
    """
    pb = instruments['pb']
    ai_task = instruments['ai_task']
    Nsamples = instruments['_config']['seq']['Nsamples']
    
    pb.start_sequence()
    data = ai_task.read_daq(Nsamples)
    
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
```

**D2. Add universal experiment runner:**

```python
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
```
* * *

### **PART E: Config File Format**

**New config format with sweep definitions:**

```python
# esr_config.py
params_dict = {
    'seq': {
        'sequence': 'esr_seq',
        't_AOM': 8000 * ms,
        'Nsamples': 8,
        'PBchannels': {...}
    },
    'mw': {
        'power': 8,
        # Note: frequency is swept, not fixed
    },
    'scan': {
        'Nruns': 30,
    },
    # NEW: Sweep definitions (order matters: first = innermost)
    'sweep': {
        'frequency': np.linspace(2.82e9, 2.92e9, 101)
    }
}
```

```python
# rabi_config.py  
params_dict = {
    'seq': {
        'sequence': 'rabi_seq',
        't_AOM': 20 * us,
        'ro_delay': 1800 * ns,
        'AOM_lag': 800 * ns,
        'MW_lag': 150 * ns,
        'Nsamples': 2,
        'PBchannels': {...}
    },
    'mw': {
        'power': 8,
        'freq': 2.87e9,
    },
    'scan': {
        'Nruns': 1,
    },
    'sweep': {
        'pulse_duration': np.arange(10e-9, 500e-9, 2e-9)
    }
}
```
* * *

### **Summary: What Changes Where**

| File | Changes |
| --- | --- |
| `parameter_system.py` | Simplify `ParameterSweep` with `run(measure_func)` |
| `PBcontrol.py` | Add `_set_X_direct()` methods that auto-reprogram |
| `SGcontrol.py` | Already has direct setters (verify) |
| `mainControl.py` | Add `diode_measurement()`, `camera_measurement()`, `run_sweep_experiment()` |
| `*_config.py` | Add `'sweep': {...}` section |
