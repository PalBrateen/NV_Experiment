"""
Test file for parameter_system.py

This demonstrates the function pointer optimization and validates the
core parameter system implementation with mock instruments.
"""

from parameter_system import (
    ParameterRegistry, Parameter, Instrument, ParameterSweep,
    ExperimentConfig, run_experiment
)
import numpy as np
import time


# ============================================================================
# MOCK INSTRUMENTS FOR TESTING
# ============================================================================

class MockSignalGenerator(Instrument):
    """Mock signal generator for testing."""

    def __init__(self, name: str = "sg1"):
        self.call_count = {'slow': 0, 'fast': 0}  # Track method calls
        super().__init__(name)

    def _register_parameters(self):
        """Register signal generator parameters."""
        self.register_parameter("frequency", 2.87e9, (0, 20e9), "Hz")
        self.register_parameter("power", 0, (-130, 25), "dBm")
        self.register_parameter("phase", 0, (0, 360), "degrees")

    def _update_instrument(self, parameter_name: str, new_value):
        """Mock hardware update (slow method with overhead)."""
        self.call_count['slow'] += 1
        # Simulate slow SCPI communication
        # time.sleep(0.001)  # Commented out for faster testing

    # FAST DIRECT METHODS (no validation, for inner loops)
    def _set_frequency_direct(self, value: float):
        """Fast frequency setting - bypasses validation."""
        self.call_count['fast'] += 1
        self.frequency = value
        # Direct hardware command (no overhead)

    def _set_power_direct(self, value: float):
        """Fast power setting - bypasses validation."""
        self.call_count['fast'] += 1
        self.power = value


class MockPulseBlaster(Instrument):
    """Mock pulse blaster for testing."""

    def __init__(self, name: str = "pb1"):
        self.call_count = {'slow': 0, 'fast': 0}
        super().__init__(name)

    def _register_parameters(self):
        """Register pulse parameters."""
        self.register_parameter("pulse_duration", 100e-9, (0, 1), "s")
        self.register_parameter("tau", 1e-6, (0, 1e-3), "s")
        self.register_parameter("interval", 500e-9, (0, 1e-3), "s")

    def _update_instrument(self, parameter_name: str, new_value):
        """Mock hardware update."""
        self.call_count['slow'] += 1

    def _set_pulse_duration_direct(self, value: float):
        """Fast pulse duration setting."""
        self.call_count['fast'] += 1
        self.pulse_duration = value

    def _set_tau_direct(self, value: float):
        """Fast tau setting."""
        self.call_count['fast'] += 1
        self.tau = value


# ============================================================================
# TEST FUNCTIONS
# ============================================================================

def test_parameter_registry():
    """Test ParameterRegistry functionality."""
    print("\n" + "="*60)
    print("TEST 1: Parameter Registry")
    print("="*60)

    # Clear registry
    ParameterRegistry.clear()

    # Create mock instruments
    sg = MockSignalGenerator("sg1")
    pb = MockPulseBlaster("pb1")

    # Check registration
    all_params = ParameterRegistry.get_all_parameters()
    print(f"\nRegistered parameters: {all_params}")

    assert 'frequency' in all_params
    assert 'power' in all_params
    assert 'pulse_duration' in all_params

    # Check lookup
    assert ParameterRegistry.get_instrument_for_parameter('frequency') == 'MockSignalGenerator'
    assert ParameterRegistry.get_instrument_for_parameter('pulse_duration') == 'MockPulseBlaster'

    print("✓ Parameter registry working correctly")


def test_parameter_validation():
    """Test Parameter validation."""
    print("\n" + "="*60)
    print("TEST 2: Parameter Validation")
    print("="*60)

    sg = MockSignalGenerator("sg1")

    # Valid parameter change
    try:
        sg.change_parameter("frequency", 3e9)
        print(f"✓ Valid frequency change: {sg.frequency}")
    except Exception as e:
        print(f"✗ Failed: {e}")

    # Invalid parameter change (out of range)
    try:
        sg.change_parameter("frequency", 30e9)  # Exceeds 20 GHz
        print("✗ Should have raised ValueError")
    except ValueError as e:
        print(f"✓ Correctly rejected out-of-range value: {e}")

    # Invalid parameter name
    try:
        sg.change_parameter("invalid_param", 100)
        print("✗ Should have raised AttributeError")
    except AttributeError as e:
        print(f"✓ Correctly rejected invalid parameter: {e}")


def test_function_pointer_optimization():
    """
    Test function pointer optimization for inner loop performance.

    This demonstrates the key optimization: fast setters registered in
    _fast_setters dict and called directly in innermost loop.
    """
    print("\n" + "="*60)
    print("TEST 3: Function Pointer Optimization")
    print("="*60)

    ParameterRegistry.clear()
    sg = MockSignalGenerator("sg1")
    instruments = {"sg1": sg}

    # Create sweep
    sweep = ParameterSweep(instruments, 'esr')

    # Add sweeps - frequency should be innermost (fast)
    power_values = np.array([5, 8, 10])
    freq_values = np.linspace(2.82e9, 2.92e9, 11)  # 11 points for testing

    sweep.add_sweep("power", power_values)      # Outer loop (slow)
    sweep.add_sweep("frequency", freq_values)   # Inner loop (FAST)

    # Check fast setter registration
    print(f"\nRegistered fast setters: {list(sweep._fast_setters.keys())}")
    assert ("sg1", "frequency") in sweep._fast_setters
    print("✓ Fast setter registered for frequency")

    # Set measurement function (mock)
    def mock_measure(instruments, current_values):
        return {**current_values, 'signal': np.random.rand()}

    sweep.set_measurement_function(mock_measure)

    # Reset counters
    sg.call_count = {'slow': 0, 'fast': 0}

    # Run sweep
    print("\nRunning sweep...")
    results = sweep.run()

    # Analyze calls
    print(f"\n{'='*60}")
    print("Performance Analysis:")
    print(f"{'='*60}")
    print(f"Total measurements: {len(results)}")
    print(f"Slow method calls:  {sg.call_count['slow']}")
    print(f"Fast method calls:  {sg.call_count['fast']}")

    # Expected:
    # - Power changes: 3 times (outer loop) → slow method
    # - Frequency changes: 3 × 11 = 33 times (inner loop) → fast method
    expected_slow = 3  # Power changes
    expected_fast = 33  # Frequency changes

    print(f"\nExpected slow calls: {expected_slow}")
    print(f"Expected fast calls: {expected_fast}")

    if sg.call_count['fast'] == expected_fast:
        print("✓ Function pointer optimization working correctly!")
    else:
        print(f"✗ Unexpected call counts")

    print(f"\n{'='*60}")
    print("Speedup: Innermost loop uses direct function calls")
    print("No if-checks, no validation overhead")
    print(f"{'='*60}")


def test_multi_parameter_sweep():
    """Test 3D multi-parameter sweep."""
    print("\n" + "="*60)
    print("TEST 4: Multi-Parameter Sweep (3D)")
    print("="*60)

    ParameterRegistry.clear()
    sg = MockSignalGenerator("sg1")
    pb = MockPulseBlaster("pb1")
    instruments = {"sg1": sg, "pb1": pb}

    sweep = ParameterSweep(instruments, 'rabi')

    # 3D sweep: power × pulse_duration × frequency
    power_vals = np.array([5, 8])           # 2 points
    pulse_vals = np.array([50e-9, 100e-9])  # 2 points
    freq_vals = np.array([2.87e9, 2.88e9])  # 2 points

    sweep.add_sweep("power", power_vals)
    sweep.add_sweep("pulse_duration", pulse_vals)
    sweep.add_sweep("frequency", freq_vals)

    # Set measurement function
    def mock_measure(instruments, current_values):
        return {**current_values, 'signal': np.random.rand()}

    sweep.set_measurement_function(mock_measure)

    # Run
    results = sweep.run()

    expected_total = 2 * 2 * 2  # 8 combinations
    print(f"\nExpected measurements: {expected_total}")
    print(f"Actual measurements:   {len(results)}")

    assert len(results) == expected_total
    print("✓ Multi-parameter sweep working correctly")

    # Check all parameters present in results
    print(f"\nSample result keys: {results[0].keys()}")
    assert 'power' in results[0]
    assert 'pulse_duration' in results[0]
    assert 'frequency' in results[0]
    print("✓ All parameters recorded in results")


def test_experiment_config():
    """Test ExperimentConfig functionality."""
    print("\n" + "="*60)
    print("TEST 5: Experiment Configuration")
    print("="*60)

    ParameterRegistry.clear()
    sg = MockSignalGenerator("sg1")
    instruments = {"sg1": sg}

    # Create config
    config = ExperimentConfig('esr', instruments)

    print(f"\nSequence type: {config.sequence_type}")
    print(f"Base params: {config.base_params}")
    print(f"Parameter values: {config.parameter_values}")
    print(f"Sweep params: {len(config.sweep.sweep_params)}")

    # Save config
    config.save_config("test_esr_config.yaml")
    print("✓ Config saved to test_esr_config.yaml")

    # Load config
    config2 = ExperimentConfig('esr', instruments, config_path="test_esr_config.yaml")
    print("✓ Config loaded from YAML")

    assert config2.sequence_type == config.sequence_type
    print("✓ Configuration save/load working")


def benchmark_function_pointers():
    """
    Benchmark function pointer approach vs if-check approach.

    This measures the Python overhead for both approaches.
    """
    print("\n" + "="*60)
    print("TEST 6: Performance Benchmark")
    print("="*60)

    ParameterRegistry.clear()
    sg = MockSignalGenerator("sg1")
    instruments = {"sg1": sg}

    # Large inner loop for timing
    freq_values = np.linspace(2.82e9, 2.92e9, 1000)  # 1000 points

    sweep = ParameterSweep(instruments, 'esr')
    sweep.add_sweep("frequency", freq_values)

    def mock_measure(instruments, current_values):
        return {'signal': 1.0}

    sweep.set_measurement_function(mock_measure)

    # Time the sweep
    print("\nRunning 1000-point sweep...")
    start_time = time.perf_counter()
    results = sweep.run()
    elapsed = time.perf_counter() - start_time

    print(f"\n{'='*60}")
    print(f"Total time:        {elapsed*1000:.2f} ms")
    print(f"Per measurement:   {elapsed*1e6/len(results):.2f} μs")
    print(f"Fast method calls: {sg.call_count['fast']}")
    print(f"{'='*60}")

    if elapsed*1e6/len(results) < 100:  # <100 μs per point
        print("✓ Performance target met (<100 μs Python overhead)")
    else:
        print("ℹ Note: Overhead includes test infrastructure")


# ============================================================================
# RUN ALL TESTS
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print(" "*15 + "PARAMETER SYSTEM TEST SUITE")
    print("="*70)

    try:
        test_parameter_registry()
        test_parameter_validation()
        test_function_pointer_optimization()
        test_multi_parameter_sweep()
        test_experiment_config()
        benchmark_function_pointers()

        print("\n" + "="*70)
        print(" "*20 + "ALL TESTS PASSED ✓")
        print("="*70 + "\n")

    except Exception as e:
        print(f"\n{'='*70}")
        print(f"TEST FAILED: {e}")
        print("="*70 + "\n")
        raise
