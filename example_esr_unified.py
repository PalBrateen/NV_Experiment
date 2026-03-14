"""
End-to-End Example: ESR Experiment with Unified System

This example demonstrates the complete workflow using the new refactored system:
1. Load configuration from esr_config.py
2. Create DiodeExperimentController
3. Initialize instruments
4. Create ParameterSweep with function pointer optimization
5. Run experiment
6. Save data and parameters (separate files)

Date: 15 December 2024
"""

from importlib import import_module
import numpy as np
from pathlib import Path

# Import new unified system
from experiment_controller import DiodeExperimentController
from parameter_system import ParameterSweep
from save_manager import SaveManager
from SGcontrol import GHz

# ============================================================================
# STEP 1: Load Configuration
# ============================================================================

print("\n" + "="*70)
print(" "*20 + "ESR Experiment Example")
print("="*70)

# Load config from esr_config.py (still .py file for flexibility!)
expCfgFile = 'esr_config'
expCfg = import_module(expCfgFile)
params_dict = expCfg.params_dict

print(f"\n✓ Configuration loaded: {expCfgFile}")
print(f"  Sequence: {params_dict['seq']['sequence']}")
print(f"  Scan points: {params_dict['scan']['Nscanpts'][0]}")
print(f"  Runs: {params_dict['scan']['Nruns']}")

# ============================================================================
# STEP 2: Create Experiment Controller
# ============================================================================

print("\n" + "="*70)
print("Creating Experiment Controller...")
print("="*70)

# Create controller for diode-based experiments
controller = DiodeExperimentController(params_dict)

# ============================================================================
# STEP 3: Initialize Instruments
# ============================================================================

print("\nInitializing instruments...")

# trial_run: [sg_mode, pb_mode] where 'n'=real, 'y'=simulation
trial_run = ['y', 'n']  # Simulated SG, real PB for this example

controller.initialize_instruments(trial_run)

# Get instruments dictionary
instruments = controller.get_instruments_dict()
sg = instruments['sg']
pb = instruments['pb']
ai_task = instruments['ai_task']

print("\n✓ Instruments initialized successfully")

# ============================================================================
# STEP 4: Create Parameter Sweep with Function Pointer Optimization
# ============================================================================

print("\n" + "="*70)
print("Setting up Parameter Sweep...")
print("="*70)

# Create sweep manager
sweep = ParameterSweep(instruments, 'esr')

# Add frequency sweep - THIS AUTOMATICALLY:
# 1. Finds SignalGenerator instrument
# 2. Registers fast setter (_set_frequency_direct) as function pointer
# 3. Stores in _fast_setters dict for <10μs overhead
freq_values = params_dict['mw']['freq']
sweep.add_sweep("frequency", freq_values)
# Output: ✓ Registered fast setter for frequency

print(f"\n✓ Sweep configured:")
print(f"  Parameter: frequency")
print(f"  Range: {freq_values[0]/GHz:.3f} - {freq_values[-1]/GHz:.3f} GHz")
print(f"  Points: {len(freq_values)}")
print(f"  Using: Function pointer optimization (fast!)")

# ============================================================================
# STEP 5: Define Measurement Function
# ============================================================================

def measure_esr(instruments_dict, current_values):
    """
    Measurement function called at each parameter point.

    This is where acquisition logic from mainControl_diode.py goes.
    For this example, we'll use a simplified version.

    Args:
        instruments_dict: Dictionary of initialized instruments
        current_values: Dict of current parameter values (e.g., {'frequency': 2.87e9})

    Returns:
        Dict with measurement results (must include all param_values)
    """
    # Extract instruments
    pb = instruments_dict['pb1']
    ai_task = instruments_dict['ai_task']

    # Get DAQ parameters
    daq_params = controller.params_dict.get('daq', {})
    Nsamples = daq_params.get('Nsamples', 300000)

    # ========================================================================
    # ACQUISITION LOGIC (simplified for example)
    # ========================================================================

    # 1. Start PulseBlaster sequence
    pb.start_sequence()

    # 2. Acquire data from AI task
    # NOTE: In real implementation, this would be:
    # counts = ai_task.read_daq(Nsamples)
    # For this example, simulate data
    counts = np.random.rand(Nsamples) * 1000  # Simulated counts

    # 3. Process data
    # Signal region (reference) vs background
    sig_region = counts[0:Nsamples//2]
    ref_region = counts[Nsamples//2:]

    signal = np.mean(sig_region)
    reference = np.mean(ref_region)
    contrast = (signal - reference) / reference if reference != 0 else 0

    # 4. Return results (must include current_values!)
    return {
        **current_values,  # Include all swept parameters
        'signal': signal,
        'reference': reference,
        'contrast': contrast,
        'raw_counts': counts  # Optional: include raw data
    }

# Register measurement function
sweep.set_measurement_function(measure_esr)
print("\n✓ Measurement function registered")

# ============================================================================
# STEP 6: Run Experiment
# ============================================================================

print("\n" + "="*70)
print("Running ESR Sweep...")
print("="*70)

# Run the sweep - uses function pointer optimization for inner loop!
results = sweep.run()

print("\n✓ Sweep completed successfully")
print(f"  Total measurements: {len(results)}")

# ============================================================================
# STEP 7: Process Results
# ============================================================================

print("\n" + "="*70)
print("Processing Results...")
print("="*70)

# Extract data arrays
frequencies = np.array([r['frequency'] for r in results])
signals = np.array([r['signal'] for r in results])
contrasts = np.array([r['contrast'] for r in results])

# Organize into standard format: (Nruns, Nscanpts, samples)
# For this example, single run:
Nruns = params_dict['scan']['Nruns']
Nscanpts = len(freq_values)
Nsamples = params_dict['seq']['Nsamples']

# Create data array matching existing format
data_array = np.zeros((Nruns, Nscanpts, 2))  # 2 channels: signal, reference
for i, result in enumerate(results):
    data_array[0, i, 0] = result['signal']
    data_array[0, i, 1] = result['reference']

print(f"\n✓ Data processed:")
print(f"  Shape: {data_array.shape}")
print(f"  Mean signal: {np.mean(signals):.2f}")
print(f"  Mean contrast: {np.mean(contrasts)*100:.2f}%")

# ============================================================================
# STEP 8: Save Data and Parameters (SEPARATE FILES!)
# ============================================================================

print("\n" + "="*70)
print("Saving Results...")
print("="*70)

save_mgr = controller.save_manager

# Get next available folder number
folder_num = save_mgr.get_next_folder_number("esr")

# Save data and parameters SEPARATELY
# This is the key feature: params can be viewed without loading data!
save_mgr.save_experiment(
    data=data_array,
    params=params_dict,
    folder_number=folder_num,
    sequence="esr",
    param_format='yaml'  # Human-readable YAML
)

print(f"\n✓ Results saved:")
print(f"  Folder: Saved_Data/{save_mgr.base_dir.name}/{time.strftime('%Y-%m-%d')}/esr_{folder_num}/")
print(f"  Data file: data_{folder_num}.npy")
print(f"  Params file: params_{folder_num}.yaml (SEPARATE!)")

# ============================================================================
# STEP 9: Load and Verify (Example)
# ============================================================================

print("\n" + "="*70)
print("Load Example (verification)...")
print("="*70)

# Example: Load only parameters (fast! no data loading)
param_file = save_mgr.base_dir / time.strftime('%Y-%m-%d') / f"esr_{folder_num}" / f"params_{folder_num}.yaml"
loaded_params = save_mgr.load_params(param_file)

print(f"\n✓ Parameters loaded separately:")
print(f"  Sequence: {loaded_params['seq']['sequence']}")
print(f"  Frequency range: {loaded_params['mw']['freq'][0]/GHz:.3f} - {loaded_params['mw']['freq'][-1]/GHz:.3f} GHz")

# Later: Load data when needed
data_file = save_mgr.base_dir / time.strftime('%Y-%m-%d') / f"esr_{folder_num}" / f"data_{folder_num}.npy"
loaded_data, loaded_params_full = save_mgr.load_experiment(data_file)

print(f"\n✓ Data loaded:")
print(f"  Shape: {loaded_data.shape}")

# ============================================================================
# STEP 10: Cleanup
# ============================================================================

print("\n" + "="*70)
print("Cleaning up...")
print("="*70)

controller.cleanup()

print("\n✓ Experiment complete!")
print("="*70)

# ============================================================================
# SUMMARY OF KEY FEATURES DEMONSTRATED
# ============================================================================

print("\n" + "="*70)
print(" "*15 + "Key Features Demonstrated")
print("="*70)
print("""
1. ✅ Unified class structure (DiodeExperimentController)
2. ✅ Automatic instrument initialization
3. ✅ Function pointer optimization for inner loop (<10μs overhead)
4. ✅ Parameter sweep with automatic instrument lookup
5. ✅ Separate parameter files (view without loading data!)
6. ✅ NumPy data format (backward compatible)
7. ✅ YAML parameter format (human-readable)
8. ✅ Automatic folder creation and numbering
9. ✅ Clean separation of concerns
10. ✅ Easy to extend for multi-parameter sweeps

Next steps:
- Fill acquisition logic from mainControl_diode.py
- Unify mainControl_diode.py and mainControl_camera.py
- Add GUI for experiment control
""")
print("="*70 + "\n")


if __name__ == "__main__":
    print("Run this file to execute the ESR experiment example!")
    print("Note: Requires esr_config.py and all instrument modules")
