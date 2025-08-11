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


class DependencyResolver:
    """Handles processing of different dependency definition formats in YAML configurations."""
    # --- Class Constants ---
    ABSOLUTE_ZERO = PhysicalConstants.ABSOLUTE_ZERO
    EPSILON = ProcessingConstants.TEMPERATURE_EPSILON
    MIN_POINTS = ProcessingConstants.MIN_TEMPERATURE_POINTS
    # Dependency reference mapping - now extensible
    @classmethod
    def get_dependency_reference_map(cls) -> Dict[str, str]:
        """Get dependency reference mapping for all registered material types."""
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
    def resolve_dependency_definition(dep_def: Union[List, str, int, float],
                                      n_values: Optional[int] = None,
                                      material: Optional[Material] = None) -> np.ndarray:
        """
        Process different dependency definition formats with optional material reference support.
        Args:
            dep_def: Dependency definition (list, string, numeric, or reference)
                - List: [300, 400, 500] (explicit values)
                - String: "(300, 500, 50)" (range format)
                - String: "(300, 50)" (equidistant format, requires n_values)
                - Float: 500.0 (single value)
                - String: "melting_temperature + 50" (reference format, requires material)
            n_values: Number of values (required for equidistant format)
            material: Material object for temperature reference resolution (optional)
        Returns:
            np.ndarray: Processed dependency array
        Examples:
            # Direct numeric value
            resolve_dependency_definition(500.0) # Returns [500.0]
            # List of values
            resolve_dependency_definition([300, 400, 500]) # Returns [300, 400, 500]
            # Equidistant format
            resolve_dependency_definition("(300, 50)", n_values=5) # Returns [300, 350, 400, 450, 500]
            # Range format
            resolve_dependency_definition("(300, 500, 50)") # Returns [300, 350, 400, 450, 500]
            # Temperature reference (requires material)
            resolve_dependency_definition("melting_temperature", material=material)
        """
        if isinstance(dep_def, list):
            return DependencyResolver._resolve_list_format(dep_def, material)
        elif isinstance(dep_def, str):
            return DependencyResolver._resolve_string_format(dep_def, n_values, material)
        elif isinstance(dep_def, (int, float)):
            dep_val = float(dep_def)
            if dep_val <= DependencyResolver.ABSOLUTE_ZERO:
                raise ValueError(
                    f"Dependency must be above absolute zero ({DependencyResolver.ABSOLUTE_ZERO}K), got {dep_val}K")
            return np.array([dep_val], dtype=float)
        else:
            raise ValueError(f"Unsupported dependency definition format: {type(dep_def)}")

    @staticmethod
    def extract_from_config(prop_config: dict, material: Material) -> np.ndarray:
        """
        Extract dependency array from property configuration.
        Args:
            prop_config: Property configuration dictionary
            material: Material object for reference resolution
        Returns:
            np.ndarray: Dependency array extracted from configuration
        """
        # Handle FILE properties (need to re-read dependency data from file)
        if FILE_PATH_KEY in prop_config:
            try:
                data_array, _ = load_property_data(prop_config)
                return data_array
            except Exception as e:
                raise ValueError(f"Failed to extract dependency array from file: {str(e)}") from e
        # Handle properties with explicit dependency definitions
        if DEPENDENCY_KEY in prop_config:
            dep_def = prop_config[DEPENDENCY_KEY]
            n_values = len(prop_config[VALUE_KEY]) if VALUE_KEY in prop_config else None
            return DependencyResolver.resolve_dependency_definition(dep_def, n_values, material)
        raise ValueError("Cannot extract dependency array: no dependency information in config")

    # --- Dependency Reference Resolution ---
    @staticmethod
    def resolve_dependency_reference(dep_ref: Union[str, float, int], material: Material) -> float:
        """
        Consolidated dependency reference resolution.
        Handles numeric values, simple references, and arithmetic expressions.
        Args:
            dep_ref: Dependency reference (string, float, or int)
            material: Material object for reference resolution
        Returns:
            float: Resolved dependency value
        Examples:
            >>> resolver = DependencyResolver()
            >>> resolver.resolve_dependency_reference(500.0, material)
            500.0
            >>> resolver.resolve_dependency_reference("melting_temperature + 50", material)
            1083.5  # Assuming copper with melting point 1033.5K
        """
        # Handle direct numeric values (int/float)
        if isinstance(dep_ref, (int, float)):
            result = float(dep_ref)
            if result <= DependencyResolver.ABSOLUTE_ZERO:
                raise ValueError(f"Temperature must be above absolute zero, got {result}K")
            return result
        # Handle string-based definitions
        if isinstance(dep_ref, str):
            # Try direct numeric conversion first
            try:
                result = float(dep_ref)
                if result <= DependencyResolver.ABSOLUTE_ZERO:
                    raise ValueError(f"Temperature must be above absolute zero, got {result}K")
                return result
            except ValueError:
                pass  # Not a numeric string, continue with reference resolution
            # Handle arithmetic expressions with dependency references
            if '+' in dep_ref or '-' in dep_ref:
                # match = re.match(r'(\w+)\s*([+-])\s*(\d+(?:\.\d+)?)', dep_ref.strip())
                match = re.match(ProcessingConstants.TEMP_ARITHMETIC_REGEX, dep_ref.strip())
                if match:
                    base_dep_name, operator, offset = match.groups()
                    base_dep = DependencyResolver.get_dependency_value(base_dep_name, material)
                    offset_val = float(offset)
                    result = base_dep + offset_val if operator == '+' else base_dep - offset_val
                    return result
            # Direct dependency reference
            result = DependencyResolver.get_dependency_value(dep_ref, material)
            return result
        raise ValueError(f"Unsupported temperature reference type: {type(dep_ref)} for value {dep_ref}")

    @staticmethod
    def get_dependency_value(dep_ref: Union[str, float, int], material: Material) -> float:
        """
        Helper function to get dependency value from material or direct numeric input.
        Args:
            dep_ref: Dependency reference (string, float, or int)
            material: Material object for reference resolution
        Returns:
            float: Dependency value
        """
        # Handle direct numeric values
        if isinstance(dep_ref, (int, float)):
            result = float(dep_ref)
            if result <= DependencyResolver.ABSOLUTE_ZERO:
                raise ValueError(f"Temperature must be above absolute zero, got {result}K")
            return result
        # Handle string references
        if isinstance(dep_ref, str):
            # Try numeric conversion first
            try:
                result = float(dep_ref)
                if result <= DependencyResolver.ABSOLUTE_ZERO:
                    raise ValueError(f"Dependency must be above absolute zero, got {result}K")
                return result
            except ValueError:
                pass  # Not numeric, try material reference
            # Dynamic material reference lookup
            dep_map = DependencyResolver.get_dependency_reference_map()
            if dep_ref in dep_map:
                attr_name = dep_map[dep_ref]
                if hasattr(material, attr_name):
                    return float(getattr(material, attr_name))
                else:
                    raise ValueError(f"Material does not have attribute '{attr_name}' for reference '{dep_ref}'")
            else:
                raise ValueError(f"Unknown dependency reference: '{dep_ref}'")
        raise ValueError(f"Unsupported dependency value type: {type(dep_ref)}")

    # --- Private Processing Methods ---
    @staticmethod
    def _resolve_list_format(dep_list: List[Union[int, float, str]],
                             material: Optional[Material] = None) -> np.ndarray:
        """
        Process explicit dependency list with optional material reference support.
        Args:
            dep_list: List of dependency values or references
            material: Material object for reference resolution (optional)
        Returns:
            np.ndarray: Processed dependency array
        """
        try:
            dep_array = []
            for dep_item in dep_list:
                if isinstance(dep_item, str) and material is not None:
                    # Handle dependency references like "solidus_temperature", "melting_temperature + 50"
                    dep_array.append(DependencyResolver.resolve_dependency_reference(dep_item, material))
                else:
                    dep_array.append(float(dep_item))
            dep_array = np.array(dep_array)
            # Validate all dependencies are above absolute zero
            if np.any(dep_array <= DependencyResolver.ABSOLUTE_ZERO):
                invalid_deps = dep_array[dep_array <= DependencyResolver.ABSOLUTE_ZERO]
                raise ValueError(f"Dependency must be above absolute zero ({DependencyResolver.ABSOLUTE_ZERO}K), "
                                 f"got {invalid_deps}")
            return dep_array
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid dependency list: {dep_list} \n -> {str(e)}") from e

    @staticmethod
    def _resolve_string_format(dep_str: str,
                               n_values: Optional[int] = None,
                               material: Optional[Material] = None) -> np.ndarray:
        """
        Process string-based dependency definitions.
        Args:
            dep_str: String temperature definition
            n_values: Number of values (required for equidistant format)
            material: Material object for reference resolution (optional)
        Returns:
            np.ndarray: Processed dependency array
        """
        if DependencyResolver._is_dependency_reference(dep_str):
            return DependencyResolver._process_simple_reference(dep_str, material)
        return DependencyResolver._process_dependency_range_format(dep_str, n_values)

    @staticmethod
    def _is_dependency_reference(dep_str: str) -> bool:
        """
        Check if the dependency string is a simple reference (not parenthesized format).
        Args:
            dep_str: Dependency string to check
        Returns:
            bool: True if it's a simple reference, False if parenthesized format
        """
        return not (dep_str.startswith('(') and dep_str.endswith(')'))

    @staticmethod
    def _process_simple_reference(dep_str: str, material: Optional[Material] = None) -> np.ndarray:
        """
        Process simple dependency reference.
        Args:
            dep_str: Simple dependency reference string
            material: Material object for reference resolution
        Returns:
            np.ndarray: Single-element dependency array
        """
        if material is not None:
            return np.array([DependencyResolver.resolve_dependency_reference(dep_str, material)])
        else:
            raise ValueError(f"String temperature definition must be enclosed in parentheses"
                             f"or require material for reference: {dep_str}")

    @staticmethod
    def _process_dependency_range_format(dep_str: str, n_values: Optional[int] = None) -> np.ndarray:
        """
        Process parenthesized dependency format.
        Args:
            dep_str: Parenthesized temperature string
            n_values: Number of values for equidistant format
        Returns:
            np.ndarray: Processed temperature array
        """
        try:
            content = dep_str.strip('()')
            values = [x.strip() for x in content.split(',')]
            if len(values) == 2:
                # Format: (start, increment/decrement) - requires n_values
                return DependencyResolver._resolve_equidistant_format(values, n_values)
            elif len(values) == 3:
                # Format: (start, stop, step/points)
                return DependencyResolver._resolve_range_format(values)
            else:
                raise ValueError(f"Dependency string must have 2 or 3 comma-separated values, got {len(values)}")
        except Exception as e:
            raise ValueError(f"Invalid dependency string format: {dep_str} \n -> {str(e)}") from e

    @staticmethod
    def _resolve_equidistant_format(values: List[str], n_values: Optional[int]) -> np.ndarray:
        """
        Process equidistant dependency format: (start, increment).
        Args:
            values: List containing [start, increment] as strings
            n_values: Number of values to generate
        Returns:
            np.ndarray: Generated dependency array
        """
        if n_values is None:
            raise ValueError(
                "Number of values required for equidistant dependency format (start, increment/decrement)")
        if n_values < DependencyResolver.MIN_POINTS:
            raise ValueError(f"Number of values must be at least {DependencyResolver.MIN_POINTS}, got {n_values}")
        try:
            start, increment = float(values[0]), float(values[1])
            if abs(increment) <= DependencyResolver.EPSILON:
                raise ValueError("Temperature increment/decrement cannot be zero")
            if start <= DependencyResolver.ABSOLUTE_ZERO:
                raise ValueError(f"Start temperature must be above absolute zero"
                                 f"({DependencyResolver.ABSOLUTE_ZERO}K), got {start}K")
            # Generate dependency array
            dep_array = np.array([start + i * increment for i in range(n_values)])
            # Validate all dependencies are above absolute zero
            if np.any(dep_array <= DependencyResolver.ABSOLUTE_ZERO):
                invalid_deps = dep_array[dep_array <= DependencyResolver.ABSOLUTE_ZERO]
                raise ValueError(f"Generated dependencies must be above absolute zero, got {invalid_deps}")
            return dep_array
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid equidistant dependency format: ({values[0]}, {values[1]}) \n -> {str(e)}") from e

    # --- Range Format Methods ---
    @staticmethod
    def _resolve_range_format(values: List[str]) -> np.ndarray:
        """
        Process range dependency format: (start, stop, step/points).
        Args:
            values: List containing [start, stop, step_or_points] as strings
        Returns:
            np.ndarray: Generated dependency array
        """
        try:
            start, stop = float(values[0]), float(values[1])
            # Validate basic dependency constraints
            DependencyResolver._validate_range_dependencies(start, stop)
            # Determine format type and delegate to appropriate method
            third_param_str = values[2].strip()
            format_type = DependencyResolver._determine_format_type(third_param_str)
            if format_type == "points":
                return DependencyResolver._resolve_points_format(values)
            elif format_type == "step":
                return DependencyResolver._resolve_step_format(values)
            else:
                raise ValueError(f"Unknown format type: {format_type}")
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid range dependency format: ({', '.join(values)}) \n -> {str(e)}") from e

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
            np.ndarray: Generated dependency array using step size
        """
        try:
            start, stop, step = float(values[0]), float(values[1]), float(values[2])
            # Validate step size
            if abs(step) <= DependencyResolver.EPSILON:
                raise ValueError("Temperature step cannot be zero")
            # Validate step direction matches range direction
            if (start < stop and step <= 0) or (start > stop and step >= 0):
                raise ValueError("Step sign must match range direction")
            # Validate step size is reasonable for the range
            if abs(step) > abs(stop - start):
                raise ValueError(f"Absolute value of step ({abs(step)}) is too large for the range. "
                                 f"It should be <= {abs(stop - start)}")
            # Generate dependency array using step size
            dep_array = np.arange(start, stop + step / 2, step)
            return dep_array
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid step format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _resolve_points_format(values: List[str]) -> np.ndarray:
        """
        Handle (start, end, num_points) format.
        Args:
            values: List containing [start, stop, num_points] as strings
        Returns:
            np.ndarray: Generated dependency array with specified number of points
        """
        try:
            start, stop = float(values[0]), float(values[1])
            n_points = int(float(values[2]))
            # Validate number of points
            if n_points < DependencyResolver.MIN_POINTS:
                raise ValueError(f"Number of points must be at least {DependencyResolver.MIN_POINTS}, got {n_points}")
            if n_points > 10000:  # Reasonable upper limit to prevent memory issues
                raise ValueError(f"Number of points ({n_points}) is too large. Maximum allowed is 10000.")
            # Generate dependency array using linspace
            dep_array = np.linspace(start, stop, n_points)
            return dep_array
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid points format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _validate_range_dependencies(start: float, stop: float) -> None:
        """
        Validate start and stop dependencies for range formats.
        Args:
            start: Start dependency
            stop: Stop dependency
        Raises:
            ValueError: If dependencies are invalid
        """
        if start <= DependencyResolver.ABSOLUTE_ZERO or stop <= DependencyResolver.ABSOLUTE_ZERO:
            raise ValueError(f"Dependencies must be above absolute zero ({DependencyResolver.ABSOLUTE_ZERO}K), "
                             f"got start={start}K, stop={stop}K")
        if abs(start - stop) <= DependencyResolver.EPSILON:
            raise ValueError(f"Start and stop dependencies must be different, got start={start}K, stop={stop}K")

    # --- Validation Methods ---
    @staticmethod
    def validate_dependency_array(dep_array: np.ndarray, context: str = "") -> None:
        """
        Validate a dependency array for common issues.
        Args:
            dep_array: Dependency array to validate
            context: Context string for error messages
        """
        if len(dep_array) == 0:
            raise ValueError(f"Dependency array is empty{' for ' + context if context else ''}")
        if len(dep_array) < DependencyResolver.MIN_POINTS:
            raise ValueError(f"Dependency array must have at least {DependencyResolver.MIN_POINTS} points, "
                             f"got {len(dep_array)}{' for ' + context if context else ''}")
        if np.any(dep_array <= DependencyResolver.ABSOLUTE_ZERO):
            invalid_deps = dep_array[dep_array <= DependencyResolver.ABSOLUTE_ZERO]
            raise ValueError(f"All dependencies must be above absolute zero ({DependencyResolver.ABSOLUTE_ZERO}K), "
                             f"got {invalid_deps}{' for ' + context if context else ''}")
        if not np.all(np.isfinite(dep_array)):
            raise ValueError(f"Dependency array contains non-finite values{' for ' + context if context else ''}")
