"""
Unit Conversion Utility for NV Experiment GUI

Provides conversion between different units for frequency, time, power, voltage, etc.
Used by parameter editor widgets for unit-aware input fields.

Features:
- Frequency: Hz, kHz, MHz, GHz
- Time: ns, us, ms, s
- Power: dBm, W, mW
- Voltage: V, mV, uV
- Magnetic field: G, mT, T
"""

from typing import Dict, List, Tuple, Optional
from enum import Enum


class UnitType(Enum):
    """Types of physical units."""
    FREQUENCY = "frequency"
    TIME = "time"
    POWER = "power"
    VOLTAGE = "voltage"
    MAGNETIC_FIELD = "magnetic_field"
    TEMPERATURE = "temperature"
    DIMENSIONLESS = "dimensionless"


class UnitConverter:
    """
    Converts values between different units.

    Usage:
        converter = UnitConverter()
        value_hz = converter.to_base(2.87, 'GHz', UnitType.FREQUENCY)  # 2.87e9
        value_ghz = converter.from_base(2.87e9, 'GHz', UnitType.FREQUENCY)  # 2.87
    """

    # Conversion factors to base units
    # Base units: Hz, s, W, V, T, K
    CONVERSION_FACTORS = {
        UnitType.FREQUENCY: {
            'Hz': 1.0,
            'kHz': 1e3,
            'MHz': 1e6,
            'GHz': 1e9,
        },
        UnitType.TIME: {
            'ns': 1e-9,
            'us': 1e-6,
            'μs': 1e-6,  # Alternative symbol
            'ms': 1e-3,
            's': 1.0,
        },
        UnitType.POWER: {
            'W': 1.0,
            'mW': 1e-3,
            'uW': 1e-6,
            'μW': 1e-6,
            'nW': 1e-9,
        },
        UnitType.VOLTAGE: {
            'V': 1.0,
            'mV': 1e-3,
            'uV': 1e-6,
            'μV': 1e-6,
        },
        UnitType.MAGNETIC_FIELD: {
            'G': 1e-4,  # Gauss to Tesla
            'mT': 1e-3,
            'T': 1.0,
        },
        UnitType.TEMPERATURE: {
            'K': 1.0,
            'C': 1.0,  # Handled separately (offset)
        },
        UnitType.DIMENSIONLESS: {
            '': 1.0,
            'counts': 1.0,
            '%': 0.01,
        }
    }

    def __init__(self):
        """Initialize unit converter."""
        pass

    def to_base(self, value: float, unit: str, unit_type: UnitType) -> float:
        """
        Convert value to base unit.

        Args:
            value: Value in given unit
            unit: Unit string (e.g., 'GHz', 'us')
            unit_type: Type of unit

        Returns:
            Value in base unit

        Example:
            >>> converter.to_base(2.87, 'GHz', UnitType.FREQUENCY)
            2870000000.0
        """
        if unit_type not in self.CONVERSION_FACTORS:
            raise ValueError(f"Unknown unit type: {unit_type}")

        factors = self.CONVERSION_FACTORS[unit_type]

        if unit not in factors:
            raise ValueError(f"Unknown unit '{unit}' for type {unit_type}")

        # Special handling for temperature
        if unit_type == UnitType.TEMPERATURE and unit == 'C':
            return value + 273.15

        return value * factors[unit]

    def from_base(self, value: float, unit: str, unit_type: UnitType) -> float:
        """
        Convert value from base unit to target unit.

        Args:
            value: Value in base unit
            unit: Target unit string
            unit_type: Type of unit

        Returns:
            Value in target unit

        Example:
            >>> converter.from_base(2.87e9, 'GHz', UnitType.FREQUENCY)
            2.87
        """
        if unit_type not in self.CONVERSION_FACTORS:
            raise ValueError(f"Unknown unit type: {unit_type}")

        factors = self.CONVERSION_FACTORS[unit_type]

        if unit not in factors:
            raise ValueError(f"Unknown unit '{unit}' for type {unit_type}")

        # Special handling for temperature
        if unit_type == UnitType.TEMPERATURE and unit == 'C':
            return value - 273.15

        return value / factors[unit]

    def convert(self, value: float, from_unit: str, to_unit: str,
                unit_type: UnitType) -> float:
        """
        Convert value from one unit to another.

        Args:
            value: Value in source unit
            from_unit: Source unit
            to_unit: Target unit
            unit_type: Type of unit

        Returns:
            Value in target unit

        Example:
            >>> converter.convert(2870, 'MHz', 'GHz', UnitType.FREQUENCY)
            2.87
        """
        base_value = self.to_base(value, from_unit, unit_type)
        return self.from_base(base_value, to_unit, unit_type)

    def get_available_units(self, unit_type: UnitType) -> List[str]:
        """
        Get list of available units for a given type.

        Args:
            unit_type: Type of unit

        Returns:
            List of unit strings
        """
        if unit_type not in self.CONVERSION_FACTORS:
            return []

        return list(self.CONVERSION_FACTORS[unit_type].keys())

    def get_preferred_unit(self, value: float, unit_type: UnitType,
                          base_unit: str = None) -> Tuple[float, str]:
        """
        Get the most appropriate unit for displaying a value.
        Chooses unit that gives value in range [0.1, 1000].

        Args:
            value: Value in base unit
            unit_type: Type of unit
            base_unit: Override base unit (if None, uses default)

        Returns:
            (converted_value, unit_string)

        Example:
            >>> converter.get_preferred_unit(2.87e9, UnitType.FREQUENCY)
            (2.87, 'GHz')
        """
        if unit_type not in self.CONVERSION_FACTORS:
            return value, ""

        factors = self.CONVERSION_FACTORS[unit_type]

        # Sort units by scale (largest to smallest)
        sorted_units = sorted(factors.items(), key=lambda x: x[1], reverse=True)

        # Find first unit that gives value >= 0.1
        for unit, factor in sorted_units:
            converted = value / factor
            if abs(converted) >= 0.1 or factor == sorted_units[-1][1]:
                return converted, unit

        # Fallback to smallest unit
        return value / sorted_units[-1][1], sorted_units[-1][0]

    def format_value(self, value: float, unit: str, unit_type: UnitType,
                    decimals: int = 3) -> str:
        """
        Format value with unit for display.

        Args:
            value: Value to format
            unit: Unit string
            unit_type: Type of unit
            decimals: Number of decimal places

        Returns:
            Formatted string (e.g., "2.870 GHz")
        """
        return f"{value:.{decimals}f} {unit}"

    def parse_value(self, text: str, unit_type: UnitType) -> Tuple[Optional[float], Optional[str]]:
        """
        Parse value and unit from string.

        Args:
            text: String to parse (e.g., "2.87 GHz" or "2870 MHz")
            unit_type: Expected unit type

        Returns:
            (value, unit) or (None, None) if parsing fails

        Example:
            >>> converter.parse_value("2.87 GHz", UnitType.FREQUENCY)
            (2.87, 'GHz')
        """
        text = text.strip()

        # Get available units for this type
        available_units = self.get_available_units(unit_type)

        # Try to find unit at end of string
        for unit in available_units:
            if text.endswith(unit):
                value_str = text[:-len(unit)].strip()
                try:
                    value = float(value_str)
                    return value, unit
                except ValueError:
                    continue

        # Try parsing as just a number (assume base unit)
        try:
            value = float(text)
            # Return with first available unit (typically base unit)
            if available_units:
                return value, available_units[0]
            return value, ""
        except ValueError:
            return None, None


# Power conversion (dBm <-> W)
def dbm_to_watts(dbm: float) -> float:
    """
    Convert power from dBm to Watts.

    Args:
        dbm: Power in dBm

    Returns:
        Power in Watts
    """
    return 10 ** ((dbm - 30) / 10)


def watts_to_dbm(watts: float) -> float:
    """
    Convert power from Watts to dBm.

    Args:
        watts: Power in Watts

    Returns:
        Power in dBm
    """
    import math
    return 10 * math.log10(watts) + 30


# Example usage and testing
if __name__ == "__main__":
    converter = UnitConverter()

    print("=== Frequency Conversions ===")
    print(f"2.87 GHz = {converter.to_base(2.87, 'GHz', UnitType.FREQUENCY)} Hz")
    print(f"2.87e9 Hz = {converter.from_base(2.87e9, 'GHz', UnitType.FREQUENCY)} GHz")
    print(f"2870 MHz = {converter.convert(2870, 'MHz', 'GHz', UnitType.FREQUENCY)} GHz")

    print("\n=== Time Conversions ===")
    print(f"8000 ns = {converter.to_base(8000, 'ns', UnitType.TIME)} s")
    print(f"8000 ns = {converter.convert(8000, 'ns', 'us', UnitType.TIME)} us")

    print("\n=== Auto Unit Selection ===")
    value, unit = converter.get_preferred_unit(2.87e9, UnitType.FREQUENCY)
    print(f"2.87e9 Hz -> {value} {unit}")

    value, unit = converter.get_preferred_unit(8000e-9, UnitType.TIME)
    print(f"8000e-9 s -> {value} {unit}")

    print("\n=== Power Conversions ===")
    print(f"10 dBm = {dbm_to_watts(10):.6f} W")
    print(f"0.01 W = {watts_to_dbm(0.01):.2f} dBm")

    print("\n=== Parse Value ===")
    value, unit = converter.parse_value("2.87 GHz", UnitType.FREQUENCY)
    print(f"Parsed '2.87 GHz' -> {value} {unit}")

    print("\n=== Available Units ===")
    print(f"Frequency: {converter.get_available_units(UnitType.FREQUENCY)}")
    print(f"Time: {converter.get_available_units(UnitType.TIME)}")
