# Complete example of Function pointer approach
```python
class SignalGenerator(Instrument):
    def _register_parameters(self):
        self.register_parameter("frequency", 1e9, (0, 20e9), "Hz")
        self.register_parameter("power", 0, (-130, 25), "dBm")
    
    def _update_instrument(self, param, value):
        # Slow safe method with overhead
        print(f"Slow method: setting {param} to {value}")
        time.sleep(0.01)  # Simulate overhead
    
    def _set_frequency_direct(self, value):
        # Fast direct method
        print(f"FAST: frequency = {value}")
        # Direct hardware command with no validation

class ParameterSweep:
    def __init__(self, instruments, sequence_type):
        self.instruments = instruments
        self.sequence_type = sequence_type
        self.sweep_params = []
        self._fast_setters = {}  # The magic dictionary!
    
    def add_sweep(self, parameter, values):
        # Find instrument
        instrument_type = ParameterRegistry.get_instrument_for_parameter(parameter)
        instrument = None
        for inst in self.instruments.values():
            if inst.__class__.__name__ == instrument_type:
                instrument = inst
                break
        
        # Store sweep info
        self.sweep_params.append((instrument.name, parameter, values))
        
        # Register fast setter if it exists
        direct_method = f'_set_{parameter}_direct'
        if hasattr(instrument, direct_method):
            # Get the function and store it
            fast_func = getattr(instrument, direct_method)
            self._fast_setters[(instrument.name, parameter)] = fast_func
            print(f"✓ Registered fast setter for {parameter}")
    
    def run(self):
        def set_parameter(inst_name, param, value, depth):
            is_innermost = (depth == len(self.sweep_params) - 1)
            
            if is_innermost:
                # Try to use fast setter
                fast_setter = self._fast_setters.get((inst_name, param))
                if fast_setter:
                    fast_setter(value)  # Call stored function
                    return
            
            # Fallback to slow method
            self.instruments[inst_name].change_parameter(param, value)
        
        # ... recursive sweep code ...
```

# Usage example
```python
sig_gen = SignalGenerator("sg1")
instruments = {"sg1": sig_gen}

sweep = ParameterSweep(instruments, 'esr')
sweep.add_sweep("frequency", [2.8e9, 2.85e9, 2.9e9])
# Output: "✓ Registered fast setter for frequency"

# Now when sweep runs:
sweep.run()
# Inner loop calls: sig_gen._set_frequency_direct(2.8e9)  <- FAST
#                   sig_gen._set_frequency_direct(2.85e9) <- FAST
#                   sig_gen._set_frequency_direct(2.9e9)  <- FAST
```

## Visual Comparison

**Old way** (if-checks):
```
For each value in inner loop:
    ┌─ Check: is this innermost? → YES
    ├─ Check: sequence_type == 'esr'? → YES
    ├─ Check: param == 'frequency'? → YES
    └─ Call: _set_frequency_direct(value)
    
    (3 checks every iteration!)
```

**New way** (function pointer):
```
Setup (once):
    Store: _fast_setters[('sg1','frequency')] = <function>

For each value in inner loop:
    ┌─ Check: is this innermost? → YES
    ├─ Lookup: _fast_setters[('sg1','frequency')] → <function>
    └─ Call: <function>(value)
    
    (1 check + 1 dictionary lookup every iteration!)