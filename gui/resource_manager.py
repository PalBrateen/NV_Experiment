"""
Resource Manager for NV Experiment GUI

Centralized instrument resource management system that ensures only one
acquisition can run at a time. Handles claim/release of instruments by
different GUI components (main window, experiment windows).

Thread-safe implementation using QMutex.
"""

from enum import Enum
from typing import Optional, List, Dict, Any
from PySide6.QtCore import QObject, Signal, QMutex, QMutexLocker


class InstrumentState(Enum):
    """Instrument availability states."""
    FREE = "free"              # Available for use
    CLAIMED = "claimed"        # Reserved but not actively in use
    IN_USE = "in_use"         # Currently being used for acquisition
    ERROR = "error"           # Instrument in error state
    DISCONNECTED = "disconnected"  # Not connected


class InstrumentType(Enum):
    """Supported instrument types."""
    SIGNAL_GENERATOR = "signal_generator"
    PULSEBLASTER = "pulseblaster"
    ANALOG_INPUT = "analog_input"
    ANALOG_OUTPUT = "analog_output"
    CAMERA = "camera"


class InstrumentManager(QObject):
    """
    Singleton class for managing instrument resources across GUI components.

    Ensures thread-safe claim/release of instruments and prevents conflicts
    when multiple experiment windows are open.

    Signals:
        instrument_claimed: Emitted when an instrument is claimed
        instrument_released: Emitted when an instrument is released
        instrument_state_changed: Emitted when instrument state changes
        conflict_detected: Emitted when a claim attempt fails due to conflict
    """

    # Signals
    instrument_claimed = Signal(str, str)  # instrument_name, owner
    instrument_released = Signal(str, str)  # instrument_name, previous_owner
    instrument_state_changed = Signal(str, str)  # instrument_name, new_state
    conflict_detected = Signal(str, str, str)  # instrument_name, current_owner, requester

    _instance: Optional['InstrumentManager'] = None
    _mutex = QMutex()

    def __new__(cls):
        """Singleton pattern - only one instance allowed."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize the instrument manager."""
        # Only initialize once (singleton)
        if hasattr(self, '_initialized'):
            return

        super().__init__()
        self._initialized = True

        # Instrument registry: {instrument_name: {instance, state, owner, type}}
        self._instruments: Dict[str, Dict[str, Any]] = {}

        # Owner registry: {owner_id: [list of claimed instruments]}
        self._owners: Dict[str, List[str]] = {}

        # Thread safety
        self._lock = QMutex()

    def register_instrument(
        self,
        name: str,
        instrument_instance: Any,
        instrument_type: InstrumentType,
        initial_state: InstrumentState = InstrumentState.DISCONNECTED
    ) -> None:
        """
        Register an instrument with the manager.

        Args:
            name: Unique identifier for the instrument
            instrument_instance: The actual instrument object
            instrument_type: Type of instrument (from InstrumentType enum)
            initial_state: Initial state (default: DISCONNECTED)
        """
        with QMutexLocker(self._lock):
            self._instruments[name] = {
                'instance': instrument_instance,
                'type': instrument_type,
                'state': initial_state,
                'owner': None,
                'error_message': None
            }

            self.instrument_state_changed.emit(name, initial_state.value)

    def unregister_instrument(self, name: str) -> None:
        """
        Unregister an instrument from the manager.

        Args:
            name: Instrument identifier
        """
        with QMutexLocker(self._lock):
            if name in self._instruments:
                # Force release if claimed
                if self._instruments[name]['owner']:
                    self._release_instrument_unsafe(name)

                del self._instruments[name]

    def claim_instruments(
        self,
        requester: str,
        instrument_names: List[str],
        force: bool = False
    ) -> tuple[bool, Optional[str]]:
        """
        Attempt to claim multiple instruments atomically.

        All instruments must be available, otherwise none are claimed.

        Args:
            requester: Identifier of the requesting component (e.g., "ESR Window #1")
            instrument_names: List of instrument names to claim
            force: If True, forcibly release instruments from current owner

        Returns:
            (success, error_message): Success status and error message if failed
        """
        with QMutexLocker(self._lock):
            # Validation
            for name in instrument_names:
                if name not in self._instruments:
                    return False, f"Instrument '{name}' not registered"

            # Check availability
            conflicts = []
            for name in instrument_names:
                inst = self._instruments[name]
                if inst['state'] in (InstrumentState.CLAIMED, InstrumentState.IN_USE):
                    if inst['owner'] != requester:
                        conflicts.append((name, inst['owner']))
                elif inst['state'] == InstrumentState.DISCONNECTED:
                    return False, f"Instrument '{name}' is disconnected"
                elif inst['state'] == InstrumentState.ERROR:
                    return False, f"Instrument '{name}' is in error state: {inst['error_message']}"

            # Handle conflicts
            if conflicts:
                if not force:
                    conflict_msg = ", ".join([f"{name} (owner: {owner})" for name, owner in conflicts])
                    for name, owner in conflicts:
                        self.conflict_detected.emit(name, owner, requester)
                    return False, f"Instruments in use: {conflict_msg}"
                else:
                    # Force release
                    for name, owner in conflicts:
                        self._release_instrument_unsafe(name)

            # Claim all instruments
            for name in instrument_names:
                self._instruments[name]['state'] = InstrumentState.CLAIMED
                self._instruments[name]['owner'] = requester
                self.instrument_claimed.emit(name, requester)
                self.instrument_state_changed.emit(name, InstrumentState.CLAIMED.value)

            # Track owner
            if requester not in self._owners:
                self._owners[requester] = []
            self._owners[requester].extend(instrument_names)

            return True, None

    def release_instruments(self, requester: str, instrument_names: Optional[List[str]] = None) -> bool:
        """
        Release instruments claimed by a requester.

        Args:
            requester: Identifier of the component releasing instruments
            instrument_names: Specific instruments to release, or None for all

        Returns:
            Success status
        """
        with QMutexLocker(self._lock):
            if requester not in self._owners:
                return True  # Nothing to release

            if instrument_names is None:
                # Release all instruments owned by requester
                instrument_names = self._owners[requester].copy()

            for name in instrument_names:
                if name in self._instruments and self._instruments[name]['owner'] == requester:
                    self._release_instrument_unsafe(name)

            # Clean up owner registry
            if requester in self._owners:
                self._owners[requester] = [
                    n for n in self._owners[requester] if n not in instrument_names
                ]
                if not self._owners[requester]:
                    del self._owners[requester]

            return True

    def _release_instrument_unsafe(self, name: str) -> None:
        """
        Release a single instrument (internal use, assumes lock is held).

        Args:
            name: Instrument identifier
        """
        if name not in self._instruments:
            return

        inst = self._instruments[name]
        previous_owner = inst['owner']

        # Set state back to FREE if not disconnected/error
        if inst['state'] not in (InstrumentState.DISCONNECTED, InstrumentState.ERROR):
            inst['state'] = InstrumentState.FREE

        inst['owner'] = None

        if previous_owner:
            self.instrument_released.emit(name, previous_owner)
            self.instrument_state_changed.emit(name, inst['state'].value)

    def set_instrument_state(
        self,
        name: str,
        state: InstrumentState,
        error_message: Optional[str] = None
    ) -> bool:
        """
        Update instrument state (e.g., CONNECTED, IN_USE, ERROR).

        Args:
            name: Instrument identifier
            state: New state
            error_message: Optional error message if state is ERROR

        Returns:
            Success status
        """
        with QMutexLocker(self._lock):
            if name not in self._instruments:
                return False

            self._instruments[name]['state'] = state
            if error_message:
                self._instruments[name]['error_message'] = error_message

            self.instrument_state_changed.emit(name, state.value)
            return True

    def get_instrument_status(self, name: str) -> Optional[Dict[str, Any]]:
        """
        Get current status of an instrument.

        Args:
            name: Instrument identifier

        Returns:
            Dictionary with 'state', 'owner', 'type', 'instance', 'error_message'
            or None if instrument not found
        """
        with QMutexLocker(self._lock):
            if name not in self._instruments:
                return None

            inst = self._instruments[name]
            return {
                'state': inst['state'],
                'owner': inst['owner'],
                'type': inst['type'],
                'instance': inst['instance'],
                'error_message': inst['error_message']
            }

    def get_instrument_instance(self, name: str) -> Optional[Any]:
        """
        Get the actual instrument instance object.

        Args:
            name: Instrument identifier

        Returns:
            Instrument instance or None if not found
        """
        with QMutexLocker(self._lock):
            if name not in self._instruments:
                return None
            return self._instruments[name]['instance']

    def get_all_instruments(self) -> Dict[str, Dict[str, Any]]:
        """
        Get status of all registered instruments.

        Returns:
            Dictionary mapping instrument names to their status
        """
        with QMutexLocker(self._lock):
            return {
                name: {
                    'state': inst['state'],
                    'owner': inst['owner'],
                    'type': inst['type'],
                    'error_message': inst['error_message']
                }
                for name, inst in self._instruments.items()
            }

    def is_available(self, name: str) -> bool:
        """
        Check if an instrument is available for claiming.

        Args:
            name: Instrument identifier

        Returns:
            True if instrument is FREE, False otherwise
        """
        with QMutexLocker(self._lock):
            if name not in self._instruments:
                return False
            return self._instruments[name]['state'] == InstrumentState.FREE

    def get_owner(self, name: str) -> Optional[str]:
        """
        Get the current owner of an instrument.

        Args:
            name: Instrument identifier

        Returns:
            Owner identifier or None if not claimed
        """
        with QMutexLocker(self._lock):
            if name not in self._instruments:
                return None
            return self._instruments[name]['owner']

    def cleanup(self) -> None:
        """Release all instruments and reset manager."""
        with QMutexLocker(self._lock):
            # Release all claimed instruments
            for name, inst in self._instruments.items():
                if inst['owner']:
                    self._release_instrument_unsafe(name)

            self._owners.clear()


# Singleton instance getter
def get_instrument_manager() -> InstrumentManager:
    """Get the singleton InstrumentManager instance."""
    return InstrumentManager()
