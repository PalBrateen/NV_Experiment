# Plan: Fix mainControl.py Logic Issue

## Problem Analysis

The user identified that `run_sweep_experiment()` in mainControl.py has a logic issue:

1. **Issue**: For ESR experiments, PulseBlaster is never programmed initially
   - ESR sweeps `frequency` (SignalGenerator parameter)
   - No PB parameters (`pulse_duration`, `tau`, etc.) in the sweep
   - Current code only programs PB when PB parameters change during sweep
   - Result: PB never gets programmed → measurement fails

2. **Root Cause**: Missing initial PB programming before sweep starts

## Solution

Add initial PB programming in `run_sweep_experiment()` before the sweep begins.

### Logic Flow

**Current (broken)**:
```
1. Create ParameterSweep
2. Add sweep parameters (e.g., 'frequency')
3. Run sweep
   - For each frequency: set frequency (SG)
   - Call measurement → PB.start_sequence() ← ERROR: PB was never programmed!
```

**Fixed**:
```
1. Create ParameterSweep
2. Check if PB parameters are in sweep
3. If NOT (like ESR), program PB once now
4. Add sweep parameters
5. Run sweep
   - For each value: set parameter
   - Call measurement → PB.start_sequence() ← Now works!
```

## Implementation

### File to Modify

- `d:\Brateen\NV_Experiment\mainControl.py` (lines ~237-273)

### Changes

In `run_sweep_experiment()`, after setting PB instrument type, add:

```python
# Initialize PB sequence BEFORE sweep (important for ESR where PB params aren't swept)
if 'pb' in instruments:
    pb = instruments['pb']
    # Check if any PB parameters are being swept
    sweep_params = params_dict.get('sweep', {}).keys()
    pb_params = ['pulse_duration', 'tau', 't_AOM', 'ro_delay']
    pb_is_swept = any(param in pb_params for param in sweep_params)

    # If no PB parameters are swept, program it once now
    if not pb_is_swept:
        print("  ℹ Programming PulseBlaster (no PB params in sweep)")
        pb._reprogram_sequence()
```

### Why This Works

1. **ESR Case** (frequency sweep):
   - No PB params in sweep → PB programmed once before sweep
   - SG frequency changes during sweep (fast setter)
   - Measurement works!

2. **Rabi Case** (pulse_duration sweep):
   - PB param in sweep → Skip initial programming
   - PB reprograms automatically on each pulse_duration change (smart setter)
   - Measurement works!

3. **Multi-param Case** (power + frequency):
   - No PB params in sweep → PB programmed once
   - Both params change during sweep
   - Measurement works!

## Testing

After fix, test:

```python
# Test 1: ESR (no PB params in sweep)
params = load_experiment_config('configs/esr_config.yaml')
results = run_sweep_experiment(params, instruments, 'diode')

# Test 2: Rabi (PB param in sweep)
params = load_experiment_config('configs/rabi_config.yaml')
results = run_sweep_experiment(params, instruments, 'diode')
```

## Summary

**Fix**: Add initial PB programming check before sweep starts in `run_sweep_experiment()`

**Lines to change**: mainControl.py, after line 256 (after `set_instr_type`)

**Impact**: Fixes ESR and other experiments where PB parameters aren't swept
