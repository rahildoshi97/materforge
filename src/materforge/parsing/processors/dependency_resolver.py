import logging
import numpy as np
import re
from typing import Dict, List, Union, Optional, Any

from materforge.core.materials import Material
from materforge.parsing.io.data_handler import load_property_data
from materforge.parsing.config.yaml_keys import (
    MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
    SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
    INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY,
    FILE_PATH_KEY, TEMPERATURE_KEY, DEPENDENCY_KEY, VALUE_KEY,
    DEPENDENCY_COLUMN_KEY, TEMPERATURE_COLUMN_KEY, PROPERTY_COLUMN_KEY,
    INDEPENDENT_VARIABLES_KEY, DEPENDENCIES_KEY, RANGES_KEY, COLUMNS_KEY,
    MAX_DEPENDENCIES, SUPPORTED_DEPENDENCY_NAMES
)
from materforge.data.constants import PhysicalConstants, ProcessingConstants

logger = logging.getLogger(__name__)

class DependencyResolver:
    """Handles multi-dependency resolution and validation for the new YAML format."""

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.independent_vars = self._extract_independent_variables()
        self.is_multi_dependency = INDEPENDENT_VARIABLES_KEY in config

    def _extract_independent_variables(self) -> Dict[str, str]:
        """Extract independent variable mapping from config."""
        if INDEPENDENT_VARIABLES_KEY in self.config:
            mapping = self.config[INDEPENDENT_VARIABLES_KEY]
            if not isinstance(mapping, dict):
                raise ValueError("independent_variables must be a dictionary")
            return mapping
        else:
            # Default single dependency (backward compatibility)
            return {'temperature': 'T'}

    def validate_dependency_structure(self, prop_name: str, prop_config: Dict[str, Any]) -> None:
        """Validate multi-dependency configuration structure."""
        if not isinstance(prop_config, dict):
            return  # Not a dependency-based property

        # Check for new multi-dependency format
        if DEPENDENCIES_KEY in prop_config:
            self._validate_multi_dependency_format(prop_name, prop_config)
        elif DEPENDENCY_KEY in prop_config:
            # Legacy single dependency format - still valid
            pass
        else:
            # No dependency information - could be constant or file-based
            pass

    def _validate_multi_dependency_format(self, prop_name: str, prop_config: Dict[str, Any]) -> None:
        """Validate new multi-dependency format."""
        dependencies = prop_config[DEPENDENCIES_KEY]

        if not isinstance(dependencies, list):
            raise ValueError(f"Property '{prop_name}': dependencies must be a list")

        if len(dependencies) > MAX_DEPENDENCIES:
            raise ValueError(f"Property '{prop_name}': maximum {MAX_DEPENDENCIES} dependencies allowed, got {len(dependencies)}")

        # Validate dependency names
        for dep in dependencies:
            if dep not in self.independent_vars:
                available = list(self.independent_vars.keys())
                raise ValueError(f"Property '{prop_name}': dependency '{dep}' not found in independent_variables. Available: {available}")

            if dep not in SUPPORTED_DEPENDENCY_NAMES:
                logger.warning(f"Dependency '{dep}' is not in standard supported list: {SUPPORTED_DEPENDENCY_NAMES}")

        # Validate ranges configuration
        if RANGES_KEY in prop_config:
            ranges_config = prop_config[RANGES_KEY]
            if not isinstance(ranges_config, dict):
                raise ValueError(f"Property '{prop_name}': ranges must be a dictionary")

            # Check that all dependencies have ranges
            for dep in dependencies:
                if dep not in ranges_config:
                    raise ValueError(f"Property '{prop_name}': missing range for dependency '{dep}'")

        # Validate bounds configuration
        if 'bounds' in prop_config:
            bounds_config = prop_config['bounds']
            if isinstance(bounds_config, dict):
                # Multi-dependency bounds format
                for dep in dependencies:
                    if dep not in bounds_config:
                        raise ValueError(f"Property '{prop_name}': missing bounds for dependency '{dep}'")
            elif isinstance(bounds_config, list) and len(dependencies) == 1:
                # Single dependency with list bounds (backward compatible)
                pass
            else:
                raise ValueError(f"Property '{prop_name}': bounds format invalid for multi-dependency")

    def resolve_dependency_ranges(self, prop_config: Dict[str, Any], material: Material) -> Dict[str, np.ndarray]:
        """Resolve ranges for all dependencies in a property."""
        if DEPENDENCIES_KEY not in prop_config:
            # Legacy single dependency format
            return self._resolve_legacy_format(prop_config, material)

        dependencies = prop_config[DEPENDENCIES_KEY]
        ranges_config = prop_config.get(RANGES_KEY, {})

        resolved_ranges = {}
        for dep_name in dependencies:
            if dep_name not in ranges_config:
                raise ValueError(f"Missing range definition for dependency '{dep_name}'")

            range_def = ranges_config[dep_name]
            # Use existing temperature resolver logic (works for any dependency type)
            resolved_range = self.resolve_dependency_definition(range_def, material=material, dependency_type=dep_name)
            resolved_ranges[dep_name] = resolved_range

            logger.debug(f"Resolved dependency '{dep_name}': {len(resolved_range)} points, range [{resolved_range[0]:.1f}, {resolved_range[-1]:.1f}]")

        return resolved_ranges

    def _resolve_legacy_format(self, prop_config: Dict[str, Any], material: Material) -> Dict[str, np.ndarray]:
        """Handle legacy single dependency format."""
        if DEPENDENCY_KEY in prop_config:
            temp_def = prop_config[DEPENDENCY_KEY]
            n_values = len(prop_config.get(VALUE_KEY, [])) if VALUE_KEY in prop_config else None
            resolved_range = self.resolve_dependency_definition(temp_def, n_values, material, "temperature")
            return {'temperature': resolved_range}
        else:
            # No explicit dependency (file-based or other)
            return {}

    def get_dependency_symbols(self, dependencies: List[str], symbol_mapping: Dict[str, Any]) -> Dict[str, Any]:
        """Get symbol mapping for specific dependencies."""
        dep_symbols = {}
        for dep_name in dependencies:
            if dep_name not in self.independent_vars:
                raise ValueError(f"Dependency '{dep_name}' not found in independent_variables")

            yaml_symbol = self.independent_vars[dep_name]
            if dep_name not in symbol_mapping:
                raise ValueError(f"No symbol provided for dependency '{dep_name}'")

            dep_symbols[yaml_symbol] = symbol_mapping[dep_name]

        return dep_symbols

    # === Core Resolution Methods ===

    @staticmethod
    def resolve_dependency_definition(dep_def: Union[List, str, int, float],
                                      n_values: Optional[int] = None,
                                      material: Optional[Material] = None,
                                      dependency_type: str = "temperature") -> np.ndarray:
        """Process different dependency definition formats with optional material reference support."""
        logger.debug("Resolving %s dependency: %s", dependency_type, dep_def)

        if isinstance(dep_def, list):
            return DependencyResolver._resolve_list_format(dep_def, material, dependency_type)
        elif isinstance(dep_def, str):
            return DependencyResolver._resolve_string_format(dep_def, n_values, material, dependency_type)
        elif isinstance(dep_def, (int, float)):
            dep_val = float(dep_def)
            DependencyResolver._validate_dependency_value(dep_val, dependency_type)
            return np.array([dep_val], dtype=float)
        else:
            raise ValueError(f"Unsupported {dependency_type} definition format: {type(dep_def)}")

    @staticmethod
    def _validate_dependency_value(value: float, dependency_type: str) -> None:
        """Validate dependency value based on type."""
        if dependency_type == "temperature" and value <= PhysicalConstants.ABSOLUTE_ZERO:
            raise ValueError(f"Temperature must be above absolute zero ({PhysicalConstants.ABSOLUTE_ZERO}K), got {value}K")
        elif dependency_type == "strain_rate" and value < 0:
            raise ValueError(f"Strain rate must be non-negative, got {value}")
        elif dependency_type == "concentration" and (value < 0 or value > 1):
            raise ValueError(f"Concentration must be between 0 and 1, got {value}")
        elif dependency_type == "pressure" and value < 0:
            raise ValueError(f"Pressure must be non-negative, got {value}")

    @staticmethod
    def extract_from_config(prop_config: dict, material: Material, dependency_type: str = "temperature") -> np.ndarray:
        """Extract dependency array from property configuration."""
        # Handle FILE properties
        if FILE_PATH_KEY in prop_config:
            try:
                # For multi-dependency files, we need to handle column mapping
                if COLUMNS_KEY in prop_config:
                    columns = prop_config[COLUMNS_KEY]
                    if dependency_type in columns:
                        # Create temporary config for load_property_data
                        temp_config = {
                            **prop_config,
                            DEPENDENCY_COLUMN_KEY: columns[dependency_type],
                            PROPERTY_COLUMN_KEY: columns.get(PROPERTY_COLUMN_KEY, columns.get("property"))
                        }
                        temp_array, _ = load_property_data(temp_config)
                        return temp_array

                # Fallback to legacy format
                temp_array, _ = load_property_data(prop_config)
                return temp_array
            except Exception as e:
                raise ValueError(f"Failed to extract {dependency_type} array from file: {str(e)}") from e

        # Handle properties with explicit dependency definitions
        dependency_key = DEPENDENCY_KEY if DEPENDENCY_KEY in prop_config else TEMPERATURE_KEY
        if dependency_key in prop_config:
            dep_def = prop_config[dependency_key]
            n_values = len(prop_config[VALUE_KEY]) if VALUE_KEY in prop_config else None
            return DependencyResolver.resolve_dependency_definition(dep_def, n_values, material, dependency_type)

        raise ValueError(f"Cannot extract {dependency_type} array: no dependency information in config")

    @staticmethod
    def resolve_dependency_reference(dep_ref: Union[str, float, int], material: Material,
                                     dependency_type: str = "temperature") -> float:
        """Consolidated dependency reference resolution."""
        logger.debug("Resolving %s reference: %s", dependency_type, dep_ref)

        # Handle direct numeric values
        if isinstance(dep_ref, (int, float)):
            result = float(dep_ref)
            DependencyResolver._validate_dependency_value(result, dependency_type)
            return result

        # Handle string-based definitions
        if isinstance(dep_ref, str):
            # Try direct numeric conversion first
            try:
                result = float(dep_ref)
                DependencyResolver._validate_dependency_value(result, dependency_type)
                return result
            except ValueError:
                pass

            # Handle arithmetic expressions with dependency references
            if '+' in dep_ref or '-' in dep_ref:
                match = re.match(ProcessingConstants.TEMP_ARITHMETIC_REGEX, dep_ref.strip())
                if match:
                    base_dep_name, operator, offset = match.groups()
                    base_dep = DependencyResolver.get_dependency_value(base_dep_name, material, dependency_type)
                    offset_val = float(offset)
                    result = base_dep + offset_val if operator == '+' else base_dep - offset_val
                    return result

            # Direct dependency reference
            result = DependencyResolver.get_dependency_value(dep_ref, material, dependency_type)
            return result

        raise ValueError(f"Unsupported {dependency_type} reference type: {type(dep_ref)} for value {dep_ref}")

    @staticmethod
    def get_dependency_value(dep_ref: Union[str, float, int], material: Material,
                             dependency_type: str = "temperature") -> float:
        """Helper function to get dependency value from material or direct numeric input."""
        # Handle direct numeric values
        if isinstance(dep_ref, (int, float)):
            result = float(dep_ref)
            DependencyResolver._validate_dependency_value(result, dependency_type)
            return result

        # Handle string references
        if isinstance(dep_ref, str):
            # Try numeric conversion first
            try:
                result = float(dep_ref)
                DependencyResolver._validate_dependency_value(result, dependency_type)
                return result
            except ValueError:
                pass

            # Dynamic material reference lookup (currently only supports temperature references)
            if dependency_type == "temperature":
                dep_map = DependencyResolver.get_dependency_reference_map()
                if dep_ref in dep_map:
                    attr_name = dep_map[dep_ref]
                    if hasattr(material, attr_name):
                        return float(getattr(material, attr_name))
                    else:
                        raise ValueError(f"Material does not have attribute '{attr_name}' for reference '{dep_ref}'")
                else:
                    raise ValueError(f"Unknown temperature reference: '{dep_ref}'")
            else:
                raise ValueError(f"Reference resolution not yet supported for {dependency_type}: '{dep_ref}'")

        raise ValueError(f"Unsupported {dependency_type} value type: {type(dep_ref)}")

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
        extension_map = {}
        return {**base_map, **extension_map}

    # === Private Helper Methods ===

    @staticmethod
    def _resolve_list_format(dep_list: List[Union[int, float, str]], material: Optional[Material] = None,
                             dependency_type: str = "temperature") -> np.ndarray:
        """Process explicit dependency list with optional material reference support."""
        try:
            dep_array = []
            for dep_item in dep_list:
                if isinstance(dep_item, str) and material is not None:
                    dep_array.append(DependencyResolver.resolve_dependency_reference(dep_item, material, dependency_type))
                else:
                    val = float(dep_item)
                    DependencyResolver._validate_dependency_value(val, dependency_type)
                    dep_array.append(val)

            dep_array = np.array(dep_array)
            return dep_array

        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid {dependency_type} list: {dep_list} \n -> {str(e)}") from e

    @staticmethod
    def _resolve_string_format(dep_str: str, n_values: Optional[int] = None, material: Optional[Material] = None,
                               dependency_type: str = "temperature") -> np.ndarray:
        """Process string-based dependency definitions."""
        if DependencyResolver._is_dependency_reference(dep_str):
            return DependencyResolver._process_simple_reference(dep_str, material, dependency_type)
        return DependencyResolver._process_dependency_range_format(dep_str, n_values, dependency_type)

    @staticmethod
    def _is_dependency_reference(dep_str: str) -> bool:
        """Check if the dependency string is a simple reference (not parenthesized format)."""
        return not (dep_str.startswith('(') and dep_str.endswith(')'))

    @staticmethod
    def _process_simple_reference(dep_str: str, material: Optional[Material] = None,
                                  dependency_type: str = "temperature") -> np.ndarray:
        """Process simple dependency reference."""
        if material is not None:
            return np.array([DependencyResolver.resolve_dependency_reference(dep_str, material, dependency_type)])
        else:
            raise ValueError(f"String {dependency_type} definition must be enclosed in parentheses "
                             f"or require material for reference: {dep_str}")

    @staticmethod
    def _process_dependency_range_format(dep_str: str, n_values: Optional[int] = None,
                                         dependency_type: str = "temperature") -> np.ndarray:
        """Process parenthesized dependency format."""
        try:
            content = dep_str.strip('()')
            values = [x.strip() for x in content.split(',')]

            if len(values) == 2:
                return DependencyResolver._resolve_equidistant_format(values, n_values, dependency_type)
            elif len(values) == 3:
                return DependencyResolver._resolve_range_format(values, dependency_type)
            else:
                raise ValueError(f"{dependency_type.capitalize()} string must have 2 or 3 comma-separated values, got {len(values)}")

        except Exception as e:
            raise ValueError(f"Invalid {dependency_type} string format: {dep_str} \n -> {str(e)}") from e

    @staticmethod
    def _resolve_equidistant_format(values: List[str], n_values: Optional[int],
                                    dependency_type: str = "temperature") -> np.ndarray:
        """Process equidistant dependency format: (start, increment)."""
        if n_values is None:
            raise ValueError(f"Number of values required for equidistant {dependency_type} format (start, increment/decrement)")

        if n_values < ProcessingConstants.MIN_TEMPERATURE_POINTS:
            raise ValueError(f"Number of values must be at least {ProcessingConstants.MIN_TEMPERATURE_POINTS}, got {n_values}")

        try:
            start, increment = float(values[0]), float(values[1])

            if abs(increment) <= ProcessingConstants.TEMPERATURE_EPSILON:
                raise ValueError(f"{dependency_type.capitalize()} increment/decrement cannot be zero")

            DependencyResolver._validate_dependency_value(start, dependency_type)

            # Generate dependency array
            dep_array = np.array([start + i * increment for i in range(n_values)])

            # Validate all values
            for val in dep_array:
                DependencyResolver._validate_dependency_value(val, dependency_type)

            return dep_array

        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid equidistant {dependency_type} format: ({values[0]}, {values[1]}) \n -> {str(e)}") from e

    @staticmethod
    def _resolve_range_format(values: List[str], dependency_type: str = "temperature") -> np.ndarray:
        """Process range dependency format: (start, stop, step/points)."""
        try:
            start, stop = float(values[0]), float(values[1])

            # Validate basic dependency constraints
            DependencyResolver._validate_dependency_value(start, dependency_type)
            DependencyResolver._validate_dependency_value(stop, dependency_type)

            if abs(start - stop) <= ProcessingConstants.TEMPERATURE_EPSILON:
                raise ValueError(f"Start and stop {dependency_type} must be different, got start={start}, stop={stop}")

            # Determine format type and delegate
            third_param_str = values[2].strip()
            format_type = DependencyResolver._determine_format_type(third_param_str)

            if format_type == "points":
                return DependencyResolver._resolve_points_format(values, dependency_type)
            elif format_type == "step":
                return DependencyResolver._resolve_step_format(values, dependency_type)
            else:
                raise ValueError(f"Unknown format type: {format_type}")

        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid range {dependency_type} format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _determine_format_type(third_param: str) -> str:
        """Determine if third parameter is step size or number of points."""
        try:
            third_param_value = float(third_param)
            is_integer_format = (
                    '.' not in third_param and
                    'e' not in third_param.lower() and
                    third_param_value == int(third_param_value) and
                    0 < third_param_value <= 1000
            )
            return "points" if is_integer_format else "step"
        except (ValueError, TypeError):
            raise ValueError(f"Third parameter must be numeric, got: {third_param}")

    @staticmethod
    def _resolve_step_format(values: List[str], dependency_type: str = "temperature") -> np.ndarray:
        """Handle (start, end, step) format."""
        try:
            start, stop, step = float(values[0]), float(values[1]), float(values[2])

            if abs(step) <= ProcessingConstants.TEMPERATURE_EPSILON:
                raise ValueError(f"{dependency_type.capitalize()} step cannot be zero")

            if (start < stop and step <= 0) or (start > stop and step >= 0):
                raise ValueError("Step sign must match range direction")

            if abs(step) > abs(stop - start):
                raise ValueError(f"Absolute value of step ({abs(step)}) is too large for the range. "
                                 f"It should be <= {abs(stop - start)}")

            dep_array = np.arange(start, stop + step / 2, step)
            return dep_array

        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid step format: ({', '.join(values)}) \n -> {str(e)}") from e

    @staticmethod
    def _resolve_points_format(values: List[str], dependency_type: str = "temperature") -> np.ndarray:
        """Handle (start, end, num_points) format."""
        try:
            start, stop = float(values[0]), float(values[1])
            n_points = int(float(values[2]))

            if n_points < ProcessingConstants.MIN_TEMPERATURE_POINTS:
                raise ValueError(f"Number of points must be at least {ProcessingConstants.MIN_TEMPERATURE_POINTS}, got {n_points}")

            if n_points > 10000:
                raise ValueError(f"Number of points ({n_points}) is too large. Maximum allowed is 10000.")

            dep_array = np.linspace(start, stop, n_points)
            return dep_array

        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid points format: ({', '.join(values)}) \n -> {str(e)}") from e
