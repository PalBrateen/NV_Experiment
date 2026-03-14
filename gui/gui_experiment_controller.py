"""
GUI Experiment Controller Wrapper

Thread-safe wrapper around DiodeExperimentController and CameraExperimentController
that provides Qt signals for GUI integration.

Key Features:
- Runs acquisition in separate QThread
- Emits signals for progress, data points, errors
- Supports pause/resume/stop operations
- Integrates with resource manager
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, Any, Optional, List
import numpy as np
import logging
import time
from PySide6.QtCore import QObject, Signal, QThread, QMutex, QMutexLocker

from experiment_controller import DiodeExperimentController, CameraExperimentController
from parameter_system import ParameterSweep
from gui.resource_manager import get_instrument_manager, InstrumentState, InstrumentType


class AcquisitionWorker(QObject):
    """
    Worker thread for running experiment acquisition.

    Signals:
        progress_updated: (current_point, total_points)
        data_point_acquired: (point_index, param_values_dict, signal_value, metadata_dict)
        run_completed: (run_number, total_runs)
        experiment_finished: ()
        error_occurred: (error_message)
    """

    # Signals
    progress_updated = Signal(int, int)  # current, total
    data_point_acquired = Signal(int, dict, float, dict)  # index, params, signal, metadata
    run_completed = Signal(int, int)  # run_number, total_runs
    experiment_finished = Signal()
    error_occurred = Signal(str)

    def __init__(
        self,
        controller: Any,  # DiodeExperimentController or CameraExperimentController
        sweep: ParameterSweep,
        num_runs: int
    ):
        """
        Initialize acquisition worker.

        Args:
            controller: Experiment controller instance
            sweep: Parameter sweep configuration
            num_runs: Number of averaging runs
        """
        super().__init__()
        self.controller = controller
        self.sweep = sweep
        self.num_runs = num_runs

        # Control flags
        self._running = False
        self._paused = False
        self._stop_requested = False
        self._mutex = QMutex()

    def run_acquisition(self):
        """Main acquisition loop (runs in worker thread)."""
        try:
            with QMutexLocker(self._mutex):
                self._running = True
                self._paused = False
                self._stop_requested = False

            total_points = len(self.sweep.param_combinations)

            # Run multiple averaging runs
            for run_num in range(1, self.num_runs + 1):
                if self._stop_requested:
                    break

                # Acquire all sweep points
                for point_idx, param_values in enumerate(self.sweep.param_combinations):
                    # Check for pause
                    while self._paused and not self._stop_requested:
                        time.sleep(0.1)

                    if self._stop_requested:
                        break

                    # Acquire single point
                    try:
                        result = self.controller._acquire_single_point(param_values)

                        # Extract signal value (voltage, counts, or intensity)
                        signal_value = result.get('signal', 0.0)

                        # Build metadata
                        metadata = {
                            'run': run_num,
                            'total_runs': self.num_runs,
                            'timestamp': time.time()
                        }

                        # Emit data point
                        self.data_point_acquired.emit(
                            point_idx,
                            param_values,
                            signal_value,
                            metadata
                        )

                        # Update progress
                        current_point = (run_num - 1) * total_points + point_idx + 1
                        total = self.num_runs * total_points
                        self.progress_updated.emit(current_point, total)

                    except Exception as e:
                        error_msg = f"Error acquiring point {point_idx}: {str(e)}"
                        logging.exception(error_msg)
                        self.error_occurred.emit(error_msg)

                # Run completed
                if not self._stop_requested:
                    self.run_completed.emit(run_num, self.num_runs)

            # Experiment finished
            if not self._stop_requested:
                self.experiment_finished.emit()

        except Exception as e:
            error_msg = f"Fatal error in acquisition: {str(e)}"
            logging.exception(error_msg)
            self.error_occurred.emit(error_msg)

        finally:
            with QMutexLocker(self._mutex):
                self._running = False
                self._paused = False

    def pause(self):
        """Pause acquisition."""
        with QMutexLocker(self._mutex):
            self._paused = True

    def resume(self):
        """Resume acquisition."""
        with QMutexLocker(self._mutex):
            self._paused = False

    def stop(self):
        """Stop acquisition."""
        with QMutexLocker(self._mutex):
            self._stop_requested = True
            self._paused = False

    def is_running(self) -> bool:
        """Check if acquisition is running."""
        with QMutexLocker(self._mutex):
            return self._running

    def is_paused(self) -> bool:
        """Check if acquisition is paused."""
        with QMutexLocker(self._mutex):
            return self._paused


class GuiExperimentController(QObject):
    """
    GUI-friendly wrapper for experiment controllers.

    Manages acquisition thread, resource claiming, and provides
    signals for GUI updates.

    Signals:
        progress_updated: (current_point, total_points)
        data_point_acquired: (point_index, param_values, signal, metadata)
        run_completed: (run_number, total_runs)
        experiment_finished: ()
        error_occurred: (error_message)
        state_changed: (new_state: str)
    """

    # Signals
    progress_updated = Signal(int, int)
    data_point_acquired = Signal(int, dict, float, dict)
    run_completed = Signal(int, int)
    experiment_finished = Signal()
    error_occurred = Signal(str)
    state_changed = Signal(str)  # 'idle', 'running', 'paused', 'finished', 'error'

    def __init__(self, owner_id: str, acquisition_mode: str = 'diode'):
        """
        Initialize GUI experiment controller.

        Args:
            owner_id: Identifier for resource management (e.g., "ESR Window #1")
            acquisition_mode: 'diode' or 'camera'
        """
        super().__init__()

        self.owner_id = owner_id
        self.acquisition_mode = acquisition_mode

        # Controllers and resources
        self.controller: Optional[Any] = None
        self.sweep: Optional[ParameterSweep] = None
        self.params_dict: Optional[Dict] = None

        # Threading
        self.worker: Optional[AcquisitionWorker] = None
        self.worker_thread: Optional[QThread] = None

        # Resource management
        self.resource_manager = get_instrument_manager()
        self.claimed_instruments: List[str] = []

        # State
        self.state = 'idle'

        # Data storage
        self.acquired_data: List[Dict] = []

    def configure(
        self,
        params_dict: Dict[str, Any],
        trial_run: List[str] = ['n', 'n']
    ) -> tuple[bool, Optional[str]]:
        """
        Configure experiment with parameters.

        Args:
            params_dict: Parameter dictionary from config
            trial_run: [sg_mode, pb_mode] simulation flags

        Returns:
            (success, error_message)
        """
        try:
            self.params_dict = params_dict

            # Create controller
            if self.acquisition_mode == 'diode':
                self.controller = DiodeExperimentController(params_dict)
            elif self.acquisition_mode == 'camera':
                self.controller = CameraExperimentController(params_dict)
            else:
                return False, f"Unknown acquisition mode: {self.acquisition_mode}"

            # Initialize instruments (will be done after claiming resources)
            # Store trial_run for later
            self.trial_run = trial_run

            # Create parameter sweep
            sweep_config = params_dict.get('sweep', {})
            if sweep_config:
                self.sweep = ParameterSweep(params_dict, self.controller.get_instruments_dict())
                for param_name, param_values in sweep_config.items():
                    self.sweep.add_sweep(param_name, param_values)
            else:
                return False, "No sweep configuration found in params_dict"

            self._set_state('idle')
            return True, None

        except Exception as e:
            error_msg = f"Configuration error: {str(e)}"
            logging.exception(error_msg)
            return False, error_msg

    def claim_resources(self) -> tuple[bool, Optional[str]]:
        """
        Claim instruments from resource manager.

        Returns:
            (success, error_message)
        """
        # Determine required instruments based on mode
        if self.acquisition_mode == 'diode':
            required = ['signal_generator', 'pulseblaster', 'analog_input']
        elif self.acquisition_mode == 'camera':
            required = ['signal_generator', 'pulseblaster', 'camera']
        else:
            return False, f"Unknown mode: {self.acquisition_mode}"

        # Attempt to claim
        success, error_msg = self.resource_manager.claim_instruments(
            self.owner_id,
            required,
            force=False
        )

        if success:
            self.claimed_instruments = required

        return success, error_msg

    def initialize_instruments(self) -> tuple[bool, Optional[str]]:
        """
        Initialize instruments (must be called after claim_resources).

        Returns:
            (success, error_message)
        """
        if not self.claimed_instruments:
            return False, "No instruments claimed. Call claim_resources() first."

        try:
            self.controller.initialize_instruments(trial_run=self.trial_run)

            # Update sweep with instrument instances
            self.sweep.instruments = self.controller.get_instruments_dict()

            # Update resource manager states to IN_USE
            for inst_name in self.claimed_instruments:
                self.resource_manager.set_instrument_state(
                    inst_name,
                    InstrumentState.FREE  # Ready but not yet in use
                )

            return True, None

        except Exception as e:
            error_msg = f"Instrument initialization error: {str(e)}"
            logging.exception(error_msg)
            return False, error_msg

    def start_acquisition(self) -> tuple[bool, Optional[str]]:
        """
        Start acquisition in worker thread.

        Returns:
            (success, error_message)
        """
        if not self.controller or not self.sweep:
            return False, "Not configured. Call configure() first."

        if not self.claimed_instruments:
            return False, "No instruments claimed. Call claim_resources() first."

        if self.worker_thread and self.worker_thread.isRunning():
            return False, "Acquisition already running"

        try:
            # Clear previous data
            self.acquired_data.clear()

            # Create worker and thread
            num_runs = self.params_dict.get('scan', {}).get('Nruns', 1)
            self.worker = AcquisitionWorker(self.controller, self.sweep, num_runs)
            self.worker_thread = QThread()

            # Move worker to thread
            self.worker.moveToThread(self.worker_thread)

            # Connect signals
            self.worker.progress_updated.connect(self._on_progress_updated)
            self.worker.data_point_acquired.connect(self._on_data_point_acquired)
            self.worker.run_completed.connect(self._on_run_completed)
            self.worker.experiment_finished.connect(self._on_experiment_finished)
            self.worker.error_occurred.connect(self._on_error_occurred)

            # Connect thread start to worker
            self.worker_thread.started.connect(self.worker.run_acquisition)

            # Update instrument states
            for inst_name in self.claimed_instruments:
                self.resource_manager.set_instrument_state(
                    inst_name,
                    InstrumentState.IN_USE
                )

            # Start thread
            self.worker_thread.start()
            self._set_state('running')

            return True, None

        except Exception as e:
            error_msg = f"Failed to start acquisition: {str(e)}"
            logging.exception(error_msg)
            return False, error_msg

    def pause_acquisition(self):
        """Pause acquisition."""
        if self.worker:
            self.worker.pause()
            self._set_state('paused')

    def resume_acquisition(self):
        """Resume acquisition."""
        if self.worker:
            self.worker.resume()
            self._set_state('running')

    def stop_acquisition(self):
        """Stop acquisition."""
        if self.worker:
            self.worker.stop()

        if self.worker_thread:
            self.worker_thread.quit()
            self.worker_thread.wait(5000)  # Wait up to 5 seconds

        self._set_state('idle')

    def release_resources(self):
        """Release claimed instruments."""
        if self.claimed_instruments:
            self.resource_manager.release_instruments(
                self.owner_id,
                self.claimed_instruments
            )
            self.claimed_instruments.clear()

    def cleanup(self):
        """Clean up controller and release resources."""
        self.stop_acquisition()

        if self.controller:
            try:
                self.controller.cleanup()
            except Exception as e:
                logging.exception(f"Error during cleanup: {e}")

        self.release_resources()

    def get_acquired_data(self) -> List[Dict]:
        """Get all acquired data points."""
        return self.acquired_data.copy()

    # Internal signal handlers
    def _on_progress_updated(self, current: int, total: int):
        """Handle progress update from worker."""
        self.progress_updated.emit(current, total)

    def _on_data_point_acquired(self, point_idx: int, params: dict, signal: float, metadata: dict):
        """Handle data point acquisition from worker."""
        # Store data
        data_point = {
            'index': point_idx,
            'parameters': params,
            'signal': signal,
            'metadata': metadata
        }
        self.acquired_data.append(data_point)

        # Forward signal
        self.data_point_acquired.emit(point_idx, params, signal, metadata)

    def _on_run_completed(self, run_num: int, total_runs: int):
        """Handle run completion from worker."""
        self.run_completed.emit(run_num, total_runs)

    def _on_experiment_finished(self):
        """Handle experiment completion from worker."""
        # Update instrument states back to CLAIMED (not in use)
        for inst_name in self.claimed_instruments:
            self.resource_manager.set_instrument_state(
                inst_name,
                InstrumentState.CLAIMED
            )

        self._set_state('finished')
        self.experiment_finished.emit()

    def _on_error_occurred(self, error_msg: str):
        """Handle error from worker."""
        self._set_state('error')
        self.error_occurred.emit(error_msg)

    def _set_state(self, new_state: str):
        """Update internal state and emit signal."""
        if self.state != new_state:
            self.state = new_state
            self.state_changed.emit(new_state)
