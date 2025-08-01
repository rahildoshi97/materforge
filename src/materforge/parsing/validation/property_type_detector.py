import logging
import re
from enum import auto, Enum
from typing import Any, Dict, Set
import sympy as sp

from materforge.parsing.config.yaml_keys import (
    FILE_PATH_KEY, COLUMNS_KEY, BOUNDS_KEY, REGRESSION_KEY, DEPENDENCIES_KEY,
    EQUATION_KEY, CONSTANT_KEY, EXTRAPOLATE_KEY, SIMPLIFY_KEY, DEGREE_KEY,
    SEGMENTS_KEY, PRE_KEY, POST_KEY, MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
    SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY, INITIAL_BOILING_TEMPERATURE_KEY,
    FINAL_BOILING_TEMPERATURE_KEY, VALUE_KEY, RANGES_KEY, MAX_STEP_FUNCTION_DEPENDENCIES,
    MAX_DEPENDENCIES_KEY
)
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)

# --- Enum ---

class PropertyType(Enum):
    CONSTANT_VALUE = auto()
    STEP_FUNCTION = auto()
    FILE_IMPORT = auto()
    TABULAR_DATA = auto()
    PIECEWISE_EQUATION = auto()
    COMPUTED_PROPERTY = auto()
    INVALID = auto()


class PropertyTypeDetector:
    """Utility class for detecting and validating property types from configuration values."""

    # Detection rules (order matters - more specific patterns first)
    DETECTION_RULES = [
        # File import (unique key check)
        (lambda c: DEPENDENCIES_KEY in c and FILE_PATH_KEY in c, PropertyType.FILE_IMPORT),
        # Step function (dependencies + value with 2 elements + single dependency)
        (lambda c: (DEPENDENCIES_KEY in c and VALUE_KEY in c and
                    PropertyTypeDetector._is_step_function(c)), PropertyType.STEP_FUNCTION),
        # Tabular data (dependencies + value)
        (lambda c: DEPENDENCIES_KEY in c and VALUE_KEY in c, PropertyType.TABULAR_DATA),
        # Piecewise equation (dependencies + equation as list)
        (lambda c: (DEPENDENCIES_KEY in c and EQUATION_KEY in c and
                    isinstance(c.get(EQUATION_KEY), list)), PropertyType.PIECEWISE_EQUATION),
        # Computed property (dependencies + equation as string)
        (lambda c: (DEPENDENCIES_KEY in c and EQUATION_KEY in c and
                    isinstance(c.get(EQUATION_KEY), str)), PropertyType.COMPUTED_PROPERTY),
    ]

    @staticmethod
    def determine_property_type(prop_name: str, config: Any) -> PropertyType:
        """Determines the property type using a declarative, rule-based approach."""
        logger.debug(f"Determining property type for '{prop_name}'")
        if PropertyTypeDetector._is_constant_format(config):
            return PropertyType.CONSTANT_VALUE
        if not isinstance(config, dict):
            raise ValueError(f"Property '{prop_name}' has an invalid format. "
                             f"Expected a dictionary or a numeric constant, but got {type(config).__name__}.")
        for detector, prop_type in PropertyTypeDetector.DETECTION_RULES:
            if detector(config):
                logger.debug(f"Detected property '{prop_name}' as type: {prop_type.name}")
                return prop_type
        present_keys = sorted(config.keys())
        raise ValueError(f"Property '{prop_name}' doesn't match any known configuration pattern. "
                         f"Present keys: {present_keys}.")

    @staticmethod
    def _is_multi_dependency_step_function(config: Dict[str, Any]) -> bool:
        """Check if config matches multi-dependency step function format."""
        if not (DEPENDENCIES_KEY in config and VALUE_KEY in config and RANGES_KEY in config):
            return False

        val_list = config.get(VALUE_KEY)
        is_two_values = isinstance(val_list, list) and len(val_list) == 2

        # Check if ranges contain single values (not arrays) for each dependency
        ranges = config.get(RANGES_KEY, {})
        dependencies = config.get(DEPENDENCIES_KEY, [])

        if len(dependencies) != 1:  # Step functions should have exactly one dependency
            return False

        dep_name = dependencies[0]
        if dep_name not in ranges:
            return False

        range_def = ranges[dep_name]
        is_single_range = not isinstance(range_def, list) or (isinstance(range_def, list) and len(range_def) == 1)

        return is_two_values and is_single_range

    @staticmethod
    def _is_multi_dependency_tabular_data(config: Dict[str, Any]) -> bool:
        """Check if config matches multi-dependency tabular data format."""
        return (DEPENDENCIES_KEY in config and VALUE_KEY in config and RANGES_KEY in config
                and EQUATION_KEY not in config)

    @staticmethod
    def _is_multi_dependency_piecewise_equation(config: Dict[str, Any]) -> bool:
        """Check if config matches multi-dependency piecewise equation format."""
        return (DEPENDENCIES_KEY in config and EQUATION_KEY in config and RANGES_KEY in config
                and isinstance(config.get(EQUATION_KEY), list))

    @staticmethod
    def _is_multi_dependency_computed_property(config: Dict[str, Any]) -> bool:
        """Check if config matches multi-dependency computed property format."""
        return (DEPENDENCIES_KEY in config and EQUATION_KEY in config and RANGES_KEY in config
                and isinstance(config.get(EQUATION_KEY), str))

    # --- High-Level Detectors (for DETECTION_RULES) ---
    @staticmethod
    def _is_constant_format(val: Any) -> bool:
        """Checks if the value has the format of a numeric constant."""
        if isinstance(val, int):
            raise ValueError(f"must be defined as a float, not an integer. Use decimal format like '{val}.0'")
        return isinstance(val, float) or (isinstance(val, str) and ('.' in val or 'e' in val.lower()))

    @staticmethod
    def _is_step_function(config: Dict[str, Any]) -> bool:
        """Check if config matches legacy step function format."""
        val_list = config.get(VALUE_KEY)
        dependencies = config.get(DEPENDENCIES_KEY, [])
        is_two_values = isinstance(val_list, list) and len(val_list) == 2
        is_single_dependency = isinstance(dependencies, list) and len(dependencies) == 1
        return is_two_values and is_single_dependency

    @staticmethod
    def validate_property_config(prop_name: str, config: Any, prop_type: PropertyType) -> None:
        """Performs strict validation based on the detected property type."""
        logger.debug(f"Validating property '{prop_name}' for type: {prop_type.name}")
        validator_map = {
            PropertyType.CONSTANT_VALUE: PropertyTypeDetector._validate_constant_value,
            PropertyType.STEP_FUNCTION: PropertyTypeDetector._validate_step_function,
            PropertyType.FILE_IMPORT: PropertyTypeDetector._validate_file_import,
            PropertyType.TABULAR_DATA: PropertyTypeDetector._validate_tabular_data,
            PropertyType.PIECEWISE_EQUATION: PropertyTypeDetector._validate_piecewise_equation,
            PropertyType.COMPUTED_PROPERTY: PropertyTypeDetector._validate_computed_property,
        }
        validator = validator_map.get(prop_type)
        if validator:
            try:
                validator(prop_name, config)
            except Exception as e:
                raise ValueError(
                    f"Invalid configuration for '{prop_name}' (expected type {prop_type.name}): {str(e)}") from e
        else:
            raise NotImplementedError(f"No validation implemented for property type: {prop_type.name}")

    @staticmethod
    def _validate_constant_value(prop_name: str, val: Any) -> None:
        try:
            float(val)
        except (ValueError, TypeError):
            raise ValueError(f"'{prop_name}' could not be converted to a float. Invalid value: '{val}'")

    @staticmethod
    def _validate_step_function(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCIES_KEY, RANGES_KEY, VALUE_KEY, BOUNDS_KEY}
        optional = {""}
        PropertyTypeDetector._check_keys(config, required, optional, "STEP_FUNCTION")
        # Validate dependencies
        dependencies = config[DEPENDENCIES_KEY]
        if not isinstance(dependencies, list) or len(dependencies) != MAX_STEP_FUNCTION_DEPENDENCIES:
            raise ValueError("Step functions must have exactly one dependency")
        dependency = dependencies[0]  # Assuming single dependency for step functions
        if not isinstance(dependency, str):
            raise ValueError(f"Dependency for step function '{prop_name}' must be a string, "
                             f"got {type(dependency).__name__}")
        # Validate ranges
        ranges = config[RANGES_KEY]
        if not isinstance(ranges, dict):
            raise ValueError("'ranges' must be a dictionary mapping dependencies to their ranges")
        if dependency not in ranges:
            raise ValueError(f"Missing range for dependency '{dependency}'")
        # Validate value
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list) or len(val_list) != 2:
            raise ValueError(f"'value' for a step function must be a list of exactly two numbers, got {val_list}")
        try:
            float(val_list[0])
            float(val_list[1])
        except (ValueError, TypeError):
            raise ValueError(f"step function values must be numeric, got {val_list}")
        # Validate bounds structure
        bounds = config[BOUNDS_KEY]
        if not isinstance(bounds, dict):
            raise ValueError("'bounds' must be a dictionary mapping dependencies to bound pairs")
        if dependency not in bounds:
            raise ValueError(f"Missing bounds for dependency '{dependency}'")
        PropertyTypeDetector._check_bounds(bounds[dependency])

    @staticmethod
    def _validate_file_import(prop_name: str, config: Dict[str, Any]) -> None:
        required = {FILE_PATH_KEY, DEPENDENCIES_KEY, COLUMNS_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "FILE_IMPORT")
        # Validate dependencies
        dependencies = config[DEPENDENCIES_KEY]
        if not isinstance(dependencies, list) or len(dependencies) != MAX_DEPENDENCIES_KEY:
            raise ValueError(f"File import supports at most {MAX_DEPENDENCIES_KEY} dependency")
        dependency = dependencies[0]
        if not isinstance(dependency, str):
            raise ValueError(f"Dependency for file import '{prop_name}' must be a string, "
                             f"got {type(dependency).__name__}")
        # Validate bounds and columns structure
        bounds = config[BOUNDS_KEY]
        columns = config[COLUMNS_KEY]
        if not isinstance(bounds, dict) or dependency not in bounds:
            raise ValueError(f"Missing bounds for dependency '{dependency}'")
        if not isinstance(columns, dict):
            raise ValueError("'columns' must be a dictionary")
        PropertyTypeDetector._check_bounds(bounds[dependency])
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])

    @staticmethod
    def _validate_tabular_data(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCIES_KEY, RANGES_KEY, VALUE_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "TABULAR_DATA")
        # Validate dependencies
        dependencies = config[DEPENDENCIES_KEY]
        if not isinstance(dependencies, list) or len(dependencies) != MAX_DEPENDENCIES_KEY:
            raise ValueError(f"Tabular data supports at most {MAX_DEPENDENCIES_KEY} dependency")
        dependency = dependencies[0]
        if not isinstance(dependency, str):
            raise ValueError(f"Dependency for tabular data '{prop_name}' must be a string, "
                             f"got {type(dependency).__name__}")
        # Validate bounds and ranges structure
        bounds = config[BOUNDS_KEY]
        ranges = config[RANGES_KEY]
        if not isinstance(bounds, dict) or dependency not in bounds:
            raise ValueError(f"Missing bounds for dependency '{dependency}'")
        if not isinstance(ranges, dict) or dependency not in ranges:
            raise ValueError(f"Missing range for dependency '{dependency}'")
        PropertyTypeDetector._check_bounds(bounds[dependency])
        # Validate value
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list):
            raise ValueError("'value' for tabular data must be a list")
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])
        if isinstance(ranges[dependency], list) and len(ranges[dependency]) != len(val_list):
            raise ValueError(f"temperature list (length {len(dependencies)}) and value list (length {len(val_list)}) "
                             f"must have the same length")

    @staticmethod
    def _validate_piecewise_equation(prop_name: str, config: Dict[str, Any]) -> None:
        # Check for multi-dependency format
        required = {DEPENDENCIES_KEY, RANGES_KEY, EQUATION_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "PIECEWISE_EQUATION")
        # Validate dependencies
        dependencies = config[DEPENDENCIES_KEY]
        if not isinstance(dependencies, list) or len(dependencies) != MAX_DEPENDENCIES_KEY:
            raise ValueError(f"Piecewise equations support at most {MAX_DEPENDENCIES_KEY} dependency")
        dependency = dependencies[0]
        if not isinstance(dependency, str):
            raise ValueError(f"Dependency for piecewise equation '{prop_name}' must be a string, "
                             f"got {type(dependency).__name__}")
        # Validate bounds and ranges structure
        bounds = config[BOUNDS_KEY]
        ranges = config[RANGES_KEY]
        if not isinstance(bounds, dict) or dependency not in bounds:
            raise ValueError(f"Missing bounds for dependency '{dependency}'")
        if not isinstance(ranges, dict) or dependency not in ranges:
            raise ValueError(f"Missing range for dependency '{dependency}'")
        PropertyTypeDetector._check_bounds(bounds[dependency])
        # Validate equation
        if not isinstance(config[EQUATION_KEY], list):
            raise ValueError("'equation' for a piecewise equation must be a list of strings")
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])

    @staticmethod
    def _validate_computed_property(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCIES_KEY, RANGES_KEY, EQUATION_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "COMPUTED_PROPERTY")
        # Validate dependencies
        dependencies = config[DEPENDENCIES_KEY]
        if not isinstance(dependencies, list) or len(dependencies) != MAX_DEPENDENCIES_KEY:
            raise ValueError(f"Computed properties support at most {MAX_DEPENDENCIES_KEY} dependencies")
        dependency = dependencies[0]
        # Validate bounds and ranges structure
        bounds = config[BOUNDS_KEY]
        ranges = config[RANGES_KEY]
        if not isinstance(bounds, dict) or dependency not in bounds:
            raise ValueError(f"Missing bounds for dependency '{dependency}'")
        if not isinstance(ranges, dict) or dependency not in ranges:
            raise ValueError(f"Missing range for dependency '{dependency}'")
        PropertyTypeDetector._check_bounds(bounds[dependency])
        # Validate equation
        if not isinstance(config[EQUATION_KEY], str):
            raise ValueError("'equation' for a computed property must be a string")
        try:
            sp.sympify(config[EQUATION_KEY])
        except (sp.SympifyError, TypeError) as e:
            raise ValueError(f"invalid mathematical expression in 'equation': {str(e)}")
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])

    @staticmethod
    def _check_keys(value: Dict[str, Any], required: Set[str], optional: Set[str], context: str) -> None:
        keys = set(value.keys())
        missing = required - keys
        if missing:
            raise ValueError(f"missing required keys for {context} property: {sorted(list(missing))}")
        extra = keys - required - optional
        if extra:
            raise ValueError(f"found unexpected keys for {context} property: {sorted(list(extra))}")

    @staticmethod
    def _check_bounds(bounds: Any) -> None:
        if not isinstance(bounds, list) or len(bounds) != 2:
            raise ValueError("'bounds' must be a list of exactly two elements")
        valid = {CONSTANT_KEY, EXTRAPOLATE_KEY}
        if bounds[0] not in valid or bounds[1] not in valid:
            raise ValueError(f"bound types must be one of {valid}, got {bounds}")

    @staticmethod
    def _check_regression(reg: Dict[str, Any]) -> None:
        PropertyTypeDetector._check_keys(reg, {SIMPLIFY_KEY, DEGREE_KEY, SEGMENTS_KEY}, set(), "regression")
        if reg[SIMPLIFY_KEY] not in {PRE_KEY, POST_KEY}:
            raise ValueError(f"regression 'simplify' must be '{PRE_KEY}' or '{POST_KEY}'")
        if not isinstance(reg[DEGREE_KEY], int) or reg[DEGREE_KEY] < 1:
            raise ValueError("regression 'degree' must be a positive integer")
        if not isinstance(reg[SEGMENTS_KEY], int) or reg[SEGMENTS_KEY] < 1:
            raise ValueError("regression 'segments' must be a positive integer")
