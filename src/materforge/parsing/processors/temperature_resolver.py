import logging
import numpy as np
import re
from typing import Dict, List, Union, Optional

from materforge.core.materials import Material
from materforge.parsing.io.data_handler import load_property_data
from materforge.parsing.config.yaml_keys import (
    MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
    SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
    INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY,
    FILE_PATH_KEY, DEPENDENCY_KEY, VALUE_KEY
)
from materforge.data.constants import PhysicalConstants, ProcessingConstants

logger = logging.getLogger(__name__)


class TemperatureResolver:
    """Handles processing of different temperature definition formats in YAML configurations."""
    # --- Class Constants ---
    ABSOLUTE_ZERO = PhysicalConstants.ABSOLUTE_ZERO
    EPSILON = ProcessingConstants.TEMPERATURE_EPSILON
    MIN_POINTS = ProcessingConstants.MIN_TEMPERATURE_POINTS
    # Temperature reference mapping - now extensible
    @classmethod
    def get_temperature_reference_map(cls) -> Dict[str, str]:
        """Get temperature reference mapping for all registered material types."""
        base_map = {
            MELTING_TEMPERATURE_KEY: "melting_temperature",
            BOILING_TEMPERATURE_KEY: "boiling_temperature",
            SOLIDUS_TEMPERATURE_KEY: "solidus_temperature",
            LIQUIDUS_TEMPERATURE_KEY: "liquidus_temperature",
            INITIAL_BOILING_TEMPERATURE_KEY: "initial_boiling_temperature",
            FINAL_BOILING_TEMPERATURE_KEY: "final_boiling_temperature",
        }
        # Add extensions for other material types
        extension_map = {}
        return {**base_map, **extension_map}

    # --- Main Public API ---
    @staticmethod
    def resolve_temperature_definition(temp_def: Union[List, str, int, float],
                                       n_values: Optional[int] = None,
                                       material: Optional[Material] = None) -> np.ndarray:
        """
        Process different temperature definition formats with optional material reference support.
        Args:
            temp_def: Temperature definition (list, string, numeric, or reference)
                - List: [300, 400, 500] (explicit temperatures)
                - String: "(300, 500, 50)" (range format)
                - String: "(300, 50)" (equidistant format, requires n_values)
                - Float: 500.0 (single temperature)
                - String: "melting_temperature + 50" (reference format, requires material)
            n_values: Number of values (required for equidistant format)
            material: Material object for temperature reference resolution (optional)
        Returns:
            np.ndarray: Processed temperature array
        Examples:
            # Direct numeric value
            resolve_temperature_definition(500.0) # Returns [500.0]
            # List of temperatures
            resolve_temperature_definition([300, 400, 500]) # Returns [300, 400, 500]
            # Equidistant format
            resolve_temperature_definition("(300, 50)", n_values=5) # Returns [300, 350, 400, 450, 500]
            # Range format
            resolve_temperature_definition("(300, 500, 50)") # Returns [300, 350, 400, 450, 500]
            # Temperature reference (requires material)
            resolve_temperature_definition("melting_temperature", material=material)
        """
        if isinstance(temp_def, list):
            return TemperatureResolver._resolve_list_format(temp_def, material)
        elif isinstance(temp_def, str):
            return TemperatureResolver._resolve_string_format(temp_def, n_values, material)
        elif isinstance(temp_def, (int, float)):
            temp_val = float(temp_def)
            if temp_val <= TemperatureResolver.ABSOLUTE_ZERO:
                raise ValueError(
                    f"Temperature must be above absolute zero ({TemperatureResolver.ABSOLUTE_ZERO}K), got {temp_val}K")
            return np.array([temp_val], dtype=float)
        else:
            raise ValueError(f"Unsupported temperature definition format: {type(temp_def)}")

    @staticmethod
    def extract_from_config(prop_config: dict, material: Material) -> np.ndarray:
        """
        Extract temperature array from property configuration.
        Args:
            prop_config: Property configuration dictionary
            material: Material object for reference resolution
        Returns:
            np.ndarray: Temperature array extracted from configuration
        """
        # Handle FILE properties (need to re-read temperature data from file)
        if FILE_PATH_KEY in prop_config:
            try:
                temp_array, _ = load_property_data(prop_config)
                return temp_array
            except Exception as e:
                raise ValueError(f"Failed to extract dependency array from file: {str(e)}") from e
        # Handle properties with explicit temperature definitions
        if DEPENDENCY_KEY in prop_config:
            temp_def = prop_config[DEPENDENCY_KEY]
            n_values = len(prop_config[VALUE_KEY]) if VALUE_KEY in prop_config else None
            return TemperatureResolver.resolve_temperature_definition(temp_def, n_values, material)
        raise ValueError("Cannot extract dependency array: no dependency information in config")

    # --- Temperature Reference Resolution ---
    @staticmethod
    def resolve_temperature_reference(temp_ref: Union[str, float, int], material: Material) -> float:
        """
        Consolidated temperature reference resolution.
        Handles numeric values, simple references, and arithmetic expressions.
        Args:
            temp_ref: Temperature reference (string, float, or int)
            material: Material object for reference resolution
        Returns:
            float: Resolved temperature value
        Examples:
            >>> resolver = TemperatureResolver()
            >>> resolver.resolve_temperature_reference(500.0, material)
            500.0
            >>> resolver.resolve_temperature_reference("melting_temperature + 50", material)
            1083.5  # Assuming copper with melting point 1033.5K
        """
        # Handle direct numeric values (int/float)
        if isinstance(temp_ref, (int, float)):
            result = float(temp_ref)
            if result <= TemperatureResolver.ABSOLUTE_ZERO:
                raise ValueError(f"Temperature must be above absolute zero, got {result}K")
            return result
        # Handle string-based definitions
        if isinstance(temp_ref, str):
            # Try direct numeric conversion first
            try:
                result = float(temp_ref)
                if result <= TemperatureResolver.ABSOLUTE_ZERO:
                    raise ValueError(f"Temperature must be above absolute zero, got {result}K")
                return result
            except ValueError:
                pass  # Not a numeric string, continue with reference resolution
            # Handle arithmetic expressions with temperature references
            if '+' in temp_ref or '-' in temp_ref:
                # match = re.match(r'(\w+)\s*([+-])\s*(\d+(?:\.\d+)?)', temp_ref.strip())
                match = re.match(ProcessingConstants.TEMP_ARITHMETIC_REGEX, temp_ref.strip())
                if match:
                    base_temp_name, operator, offset = match.groups()
                    base_temp = TemperatureResolver.get_temperature_value(base_temp_name, material)
                    offset_val = float(offset)
                    result = base_temp + offset_val if operator == '+' else base_temp - offset_val
                    return result
            # Direct temperature reference
            result = TemperatureResolver.get_temperature_value(temp_ref, material)
            return result
        raise ValueError(f"Unsupported temperature reference type: {type(temp_ref)} for value {temp_ref}")

    @staticmethod
    def get_temperature_value(temp_ref: Union[str, float, int], material: Material) -> float:
        """
        Helper function to get temperature value from material or direct numeric input.
        Args:
            temp_ref: Temperature reference (string, float, or int)
            material: Material object for reference resolution
        Returns:
            float: Temperature value
        """
        # Handle direct numeric values
        if isinstance(temp_ref, (int, float)):
            result = float(temp_ref)
            if result <= TemperatureResolver.ABSOLUTE_ZERO:
                raise ValueError(f"Temperature must be above absolute zero, got {result}K")
            return result
        # Handle string references
        if isinstance(temp_ref, str):
            # Try numeric conversion first
            try:
                result = float(temp_ref)
                if result <= TemperatureResolver.ABSOLUTE_ZERO:
                    raise ValueError(f"Temperature must be above absolute zero, got {result}K")
                return result
            except ValueError:
                pass  # Not numeric, try material reference
            # Dynamic material reference lookup
            temp_map = TemperatureResolver.get_temperature_reference_map()
            if temp_ref in temp_map:
                attr_name = temp_map[temp_ref]
                if hasattr(material, attr_name):
                    return float(getattr(material, attr_name))
                else:
                    raise ValueError(f"Material does not have attribute '{attr_name}' for reference '{temp_ref}'")
            else:
                raise ValueError(f"Unknown temperature reference: '{temp_ref}'")
        raise ValueError(f"Unsupported temperature value type: {type(temp_ref)}")

    # --- Private Processing Methods ---
    @staticmethod
    def _resolve_list_format(temp_list: List[Union[int, float, str]],
                             material: Optional[Material] = None) -> np.ndarray:
        """
        Process explicit temperature list with optional material reference support.
        Args:
            temp_list: List of temperature values or references
            material: Material object for reference resolution (optional)
        Returns:
            np.ndarray: Processed temperature array
        """
        try:
            temp_array = []
            for temp_item in temp_list:
                if isinstance(temp_item, str) and material is not None:
                    # Handle temperature references like "solidus_temperature", "melting_temperature + 50"
                    temp_array.append(TemperatureResolver.resolve_temperature_reference(temp_item, material))
                else:
                    temp_array.append(float(temp_item))
            temp_array = np.array(temp_array)
            # Validate all temperatures are above absolute zero
            if np.any(temp_array <= TemperatureResolver.ABSOLUTE_ZERO):
                invalid_temps = temp_array[temp_array <= TemperatureResolver.ABSOLUTE_ZERO]
                raise ValueError(f"Temperature must be above absolute zero ({TemperatureResolver.ABSOLUTE_ZERO}K), "
                                 f"got {invalid_temps}")
            return temp_array
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid temperature list: {temp_list} \n -> {str(e)}") from e

    @staticmethod
    def _resolve_string_format(temp_str: str,
                               n_values: Optional[int] = None,
                               material: Optional[Material] = None) -> np.ndarray:
        """
        Process string-based temperature definitions.
        Args:
            temp_str: String temperature definition
            n_values: Number of values (required for equidistant format)
            material: Material object for reference resolution (optional)
        Returns:
            np.ndarray: Processed temperature array
        """
        if TemperatureResolver._is_temperature_reference(temp_str):
            return TemperatureResolver._process_simple_reference(temp_str, material)
        return TemperatureResolver._process_temperature_range_format(temp_str, n_values)

    @staticmethod
    def _is_temperature_reference(temp_str: str) -> bool:
        """
        Check if the temperature string is a simple reference (not parenthesized format).
        Args:
            temp_str: Temperature string to check
        Returns:
            bool: True if it's a simple reference, False if parenthesized format
        """
        return not (temp_str.startswith('(') and temp_str.endswith(')'))

    @staticmethod
    def _process_simple_reference(temp_str: str, material: Optional[Material] = None) -> np.ndarray:
        """
        Process simple temperature reference.
        Args:
            temp_str: Simple temperature reference string
            material: Material object for reference resolution
        Returns:
            np.ndarray: Single-element temperature array
        """
        if material is not None:
            return np.array([TemperatureResolver.resolve_temperature_reference(temp_str, material)])
        else:
            raise ValueError(f"String temperature definition must be enclosed in parentheses"
                             f"or require material for reference: {temp_str}")

    @staticmethod
    def _process_temperature_range_format(temp_str: str, n_values: Optional[int] = None) -> np.ndarray:
        """
        Process parenthesized temperature format.
        Args:
            temp_str: Parenthesized temperature string
            n_values: Number of values for equidistant format
        Returns:
            np.ndarray: Processed temperature array
        """
        try:
            content = temp_str.strip('()')
            values = [x.strip() for x in content.split(',')]
            if len(values) == 2:
                # Format: (start, increment/decrement) - requires n_values
                return TemperatureResolver._resolve_equidistant_format(values, n_values)
            elif len(values) == 3:
                # Format: (start, stop, step/points)
                return TemperatureResolver._resolve_range_format(values)
            else:
                raise ValueError(f"Temperature string must have 2 or 3 comma-separated values, got {len(values)}")
        except Exception as e:
            raise ValueError(f"Invalid temperature string format: {temp_str} \n -> {str(e)}") from e

    @staticmethod
    def _resolve_equidistant_format(values: List[str], n_values: Optional[int]) -> np.ndarray:
        """
        Process equidistant temperature format: (start, increment).
        Args:
            values: List containing [start, increment] as strings
            n_values: Number of values to generate
        Returns:
            np.ndarray: Generated temperature array
        """
        if n_values is None:
            raise ValueError(
                "Number of values required for equidistant temperature format (start, increment/decrement)")
        if n_values < TemperatureResolver.MIN_POINTS:
            raise ValueError(f"Number of values must be at least {TemperatureResolver.MIN_POINTS}, got {n_values}")
        try:
            start, increment = float(values[0]), float(values[1])
            if abs(increment) <= TemperatureResolver.EPSILON:
                raise ValueError("Temperature increment/decrement cannot be zero")
            if start <= TemperatureResolver.ABSOLUTE_ZERO:
                raise ValueError(f"Start temperature must be above absolute zero"
                                 f"({TemperatureResolver.ABSOLUTE_ZERO}K), got {start}K")
            # Generate temperature array
            temp_array = np.array([start + i * increment for i in range(n_values)])
            # Validate all temperatures are above absolute zero
            if np.any(temp_array <= TemperatureResolver.ABSOLUTE_ZERO):
                invalid_temps = temp_array[temp_array <= TemperatureResolver.ABSOLUTE_ZERO]
                raise ValueError(f"Generated temperatures must be above absolute zero, got {invalid_temps}")
            return temp_array
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid equidistant temperature format: ({values[0]}, {values[1]}) \n -> {str(e)}") from e

    # --- Range Format Methods ---
    @staticmethod
    def _resolve_range_format(values: List[str]) -> np.ndarray:
        """
        Process range temperature format: (start, stop, step/points).
        Args:
            values: List containing [start, stop, step_or_points] as strings
        Returns:
            np.ndarray: Generated temperature array
        """
        try:
            start, stop = float(values[0]), float(values[1])
            # Validate basic temperature constraints
            TemperatureResolver._validate_range_temperatures(start, stop)
            # Determine format type and delegate to appropriate method
            third_param_str = values[2].strip()
            format_type = TemperatureResolver._determine_format_type(third_param_str)
            if format_type == "points":
                return TemperatureResolver._resolve_points_format(values)
            elif format_type == "step":
                return TemperatureResolver._resolve_step_format(values)
            else:
                raise ValueError(f"Unknown format type: {format_type}")
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid range temperature format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _determine_format_type(third_param: str) -> str:
        """
        Determine if third parameter is step size or number of points.
        Args:
            third_param: Third parameter as string
        Returns:
            str: Either "step" or "points"
        """
        try:
            third_param_value = float(third_param)
            # Check if string represents an integer (no decimal point)
            is_integer_format = (
                    '.' not in third_param and
                    'e' not in third_param.lower() and
                    third_param_value == int(third_param_value) and
                    0 < third_param_value <= 1000
            )
            # Determine if it's number of points or step size
            if is_integer_format:
                return "points"
            else:
                return "step"
        except (ValueError, TypeError):
            raise ValueError(f"Third parameter must be numeric, got: {third_param}")

    @staticmethod
    def _resolve_step_format(values: List[str]) -> np.ndarray:
        """
        Handle (start, end, step) format.
        Args:
            values: List containing [start, stop, step] as strings
        Returns:
            np.ndarray: Generated temperature array using step size
        """
        try:
            start, stop, step = float(values[0]), float(values[1]), float(values[2])
            # Validate step size
            if abs(step) <= TemperatureResolver.EPSILON:
                raise ValueError("Temperature step cannot be zero")
            # Validate step direction matches range direction
            if (start < stop and step <= 0) or (start > stop and step >= 0):
                raise ValueError("Step sign must match range direction")
            # Validate step size is reasonable for the range
            if abs(step) > abs(stop - start):
                raise ValueError(f"Absolute value of step ({abs(step)}) is too large for the range. "
                                 f"It should be <= {abs(stop - start)}")
            # Generate temperature array using step size
            temp_array = np.arange(start, stop + step / 2, step)
            return temp_array
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid step format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _resolve_points_format(values: List[str]) -> np.ndarray:
        """
        Handle (start, end, num_points) format.
        Args:
            values: List containing [start, stop, num_points] as strings
        Returns:
            np.ndarray: Generated temperature array with specified number of points
        """
        try:
            start, stop = float(values[0]), float(values[1])
            n_points = int(float(values[2]))
            # Validate number of points
            if n_points < TemperatureResolver.MIN_POINTS:
                raise ValueError(f"Number of points must be at least {TemperatureResolver.MIN_POINTS}, got {n_points}")
            if n_points > 10000:  # Reasonable upper limit to prevent memory issues
                raise ValueError(f"Number of points ({n_points}) is too large. Maximum allowed is 10000.")
            # Generate temperature array using linspace
            temp_array = np.linspace(start, stop, n_points)
            return temp_array
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid points format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _validate_range_temperatures(start: float, stop: float) -> None:
        """
        Validate start and stop temperatures for range formats.
        Args:
            start: Start temperature
            stop: Stop temperature
        Raises:
            ValueError: If temperatures are invalid
        """
        if start <= TemperatureResolver.ABSOLUTE_ZERO or stop <= TemperatureResolver.ABSOLUTE_ZERO:
            raise ValueError(f"Temperatures must be above absolute zero ({TemperatureResolver.ABSOLUTE_ZERO}K), "
                             f"got start={start}K, stop={stop}K")
        if abs(start - stop) <= TemperatureResolver.EPSILON:
            raise ValueError(f"Start and stop temperatures must be different, got start={start}K, stop={stop}K")

    # --- Validation Methods ---
    @staticmethod
    def validate_temperature_array(temp_array: np.ndarray, context: str = "") -> None:
        """
        Validate a temperature array for common issues.
        Args:
            temp_array: Temperature array to validate
            context: Context string for error messages
        """
        if len(temp_array) == 0:
            raise ValueError(f"Temperature array is empty{' for ' + context if context else ''}")
        if len(temp_array) < TemperatureResolver.MIN_POINTS:
            raise ValueError(f"Temperature array must have at least {TemperatureResolver.MIN_POINTS} points, "
                             f"got {len(temp_array)}{' for ' + context if context else ''}")
        if np.any(temp_array <= TemperatureResolver.ABSOLUTE_ZERO):
            invalid_temps = temp_array[temp_array <= TemperatureResolver.ABSOLUTE_ZERO]
            raise ValueError(f"All temperatures must be above absolute zero ({TemperatureResolver.ABSOLUTE_ZERO}K), "
                             f"got {invalid_temps}{' for ' + context if context else ''}")
        if not np.all(np.isfinite(temp_array)):
            raise ValueError(f"Temperature array contains non-finite values{' for ' + context if context else ''}")
