import logging
import re
from enum import auto, Enum
from typing import Any, Dict, Set
import sympy as sp
from materforge.parsing.config.yaml_keys import (
    FILE_PATH_KEY, COLUMNS_KEY, PROPERTY_KEY, BOUNDS_KEY, REGRESSION_KEY,
    EQUATION_KEY, CONSTANT_KEY, EXTRAPOLATE_KEY, SIMPLIFY_KEY,
    DEGREE_KEY, SEGMENTS_KEY, PRE_KEY, POST_KEY, VALUE_KEY,
    DEPENDENCIES_KEY, RANGES_KEY, MAX_STEP_FUNCTION_DEPENDENCIES,
    STEP_FUNCTION_LIMITATION_MESSAGE
)
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)

class PropertyType(Enum):
    CONSTANT_VALUE = auto()
    STEP_FUNCTION = auto()
    FILE_IMPORT = auto()
    TABULAR_DATA = auto()
    PIECEWISE_EQUATION = auto()
    COMPUTED_PROPERTY = auto()
    INVALID = auto()

class PropertyTypeDetector:
    """Utility class for detecting and validating property types from multi-dependency configurations."""

    # Detection rules for multi-dependency format
    DETECTION_RULES = [
        # File import (unique key combination)
        (lambda c: DEPENDENCIES_KEY in c and FILE_PATH_KEY in c, PropertyType.FILE_IMPORT),

        # Multi-dependency format checks
        (lambda c: DEPENDENCIES_KEY in c and VALUE_KEY in c and PropertyTypeDetector._is_step_function(c),
         PropertyType.STEP_FUNCTION),
        (lambda c: DEPENDENCIES_KEY in c and VALUE_KEY in c, PropertyType.TABULAR_DATA),
        (lambda c: DEPENDENCIES_KEY in c and EQUATION_KEY in c and isinstance(c.get(EQUATION_KEY), list),
         PropertyType.PIECEWISE_EQUATION),
        (lambda c: DEPENDENCIES_KEY in c and EQUATION_KEY in c and isinstance(c.get(EQUATION_KEY), str),
         PropertyType.COMPUTED_PROPERTY),
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
    def _is_constant_format(val: Any) -> bool:
        """Checks if the value has the format of a numeric constant."""
        if isinstance(val, int):
            raise ValueError(f"must be defined as a float, not an integer. Use decimal format like '{val}.0'")
        return isinstance(val, float) or (isinstance(val, str) and ('.' in val or 'e' in val.lower()))

    @staticmethod
    def _is_step_function(config: Dict[str, Any]) -> bool:
        """Check for step function in multi-dependency format."""
        val_list = config.get(VALUE_KEY)
        dependencies = config.get(DEPENDENCIES_KEY, [])
        ranges = config.get(RANGES_KEY, {})

        # Must have exactly 2 values and single dependency
        is_two_values = isinstance(val_list, list) and len(val_list) == 2
        is_single_dep = len(dependencies) == 1

        if is_two_values and is_single_dep:
            dep_name = dependencies[0]
            if dep_name in ranges:
                # Check if range is a single value (not a list/tuple)
                range_def = ranges[dep_name]
                return not isinstance(range_def, (list, tuple))

        return False

    @staticmethod
    def _validate_step_function(prop_name: str, config: Dict[str, Any]) -> None:
        """Validate step function - ENFORCES single dependency limitation."""
        required = {DEPENDENCIES_KEY, RANGES_KEY, VALUE_KEY}
        optional = {BOUNDS_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "STEP_FUNCTION")

        dependencies = config[DEPENDENCIES_KEY]

        # **KEY VALIDATION: Enforce single dependency for step functions**
        if len(dependencies) != MAX_STEP_FUNCTION_DEPENDENCIES:
            raise ValueError(
                f"Step functions must have exactly one dependency, got {len(dependencies)}. "
                f"{STEP_FUNCTION_LIMITATION_MESSAGE}"
            )

        # Validate bounds format
        bounds_config = config.get(BOUNDS_KEY, {})
        if isinstance(bounds_config, dict):
            # Multi-dependency bounds format
            dep_name = dependencies[0]
            if dep_name not in bounds_config:
                raise ValueError(f"Missing bounds for dependency '{dep_name}'")
            PropertyTypeDetector._check_bounds(bounds_config[dep_name])
        elif isinstance(bounds_config, list):
            # Legacy bounds format fallback
            PropertyTypeDetector._check_bounds(bounds_config)

        # Validate value array
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list) or len(val_list) != 2:
            raise ValueError(f"'value' for a step function must be a list of exactly two numbers, got {val_list}")

        try:
            float(val_list[0])
            float(val_list[1])
        except (ValueError, TypeError):
            raise ValueError(f"step function values must be numeric, got {val_list}")

    @staticmethod
    def _validate_tabular_data(prop_name: str, config: Dict[str, Any]) -> None:
        """Validate tabular data for multi-dependency format."""
        required = {DEPENDENCIES_KEY, RANGES_KEY, VALUE_KEY}
        optional = {BOUNDS_KEY, REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "TABULAR_DATA")

        # Validate bounds
        bounds_config = config.get(BOUNDS_KEY, {})
        if isinstance(bounds_config, dict):
            dependencies = config[DEPENDENCIES_KEY]
            for dep in dependencies:
                if dep in bounds_config:
                    PropertyTypeDetector._check_bounds(bounds_config[dep])

        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])

        dependencies = config[DEPENDENCIES_KEY]
        val_list = config[VALUE_KEY]
        ranges_config = config.get(RANGES_KEY, {})

        if not isinstance(val_list, list):
            raise ValueError("'value' for tabular data must be a list.")

        # For single dependency, check if ranges and values match
        if len(dependencies) == 1:
            dep_name = dependencies[0]
            if dep_name in ranges_config:
                range_def = ranges_config[dep_name]
                if isinstance(range_def, list) and len(range_def) != len(val_list):
                    raise ValueError(f"Range list (length {len(range_def)}) and value list (length {len(val_list)}) "
                                     f"must have the same length for single dependency")

    @staticmethod
    def _validate_piecewise_equation(prop_name: str, config: Dict[str, Any]) -> None:
        """Validate piecewise equation for multi-dependency format."""
        required = {DEPENDENCIES_KEY, RANGES_KEY, EQUATION_KEY}
        optional = {BOUNDS_KEY, REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "PIECEWISE_EQUATION")

        # Validate bounds if present
        if BOUNDS_KEY in config:
            bounds_config = config[BOUNDS_KEY]
            if isinstance(bounds_config, dict):
                dependencies = config[DEPENDENCIES_KEY]
                for dep in dependencies:
                    if dep in bounds_config:
                        PropertyTypeDetector._check_bounds(bounds_config[dep])
            elif isinstance(bounds_config, list):
                PropertyTypeDetector._check_bounds(bounds_config)

        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])

        if not isinstance(config[EQUATION_KEY], list):
            raise ValueError("'equation' for a piecewise equation must be a list of strings")

    @staticmethod
    def _validate_computed_property(prop_name: str, config: Dict[str, Any]) -> None:
        """Validate computed property for multi-dependency format."""
        required = {DEPENDENCIES_KEY, RANGES_KEY, EQUATION_KEY}
        optional = {BOUNDS_KEY, REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "COMPUTED_PROPERTY")

        # Validate bounds if present
        if BOUNDS_KEY in config:
            bounds_config = config[BOUNDS_KEY]
            if isinstance(bounds_config, dict):
                dependencies = config[DEPENDENCIES_KEY]
                for dep in dependencies:
                    if dep in bounds_config:
                        PropertyTypeDetector._check_bounds(bounds_config[dep])
            elif isinstance(bounds_config, list):
                PropertyTypeDetector._check_bounds(bounds_config)

        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])

        if not isinstance(config[EQUATION_KEY], str):
            raise ValueError("'equation' for a computed property must be a string")

        try:
            sp.sympify(config[EQUATION_KEY])
        except (sp.SympifyError, TypeError) as e:
            raise ValueError(f"invalid mathematical expression in 'equation': {str(e)}")

    @staticmethod
    def _validate_file_import(prop_name: str, config: Dict[str, Any]) -> None:
        """Validate file import - enhanced for multi-dependency."""
        required = {DEPENDENCIES_KEY, FILE_PATH_KEY, COLUMNS_KEY}
        optional = {BOUNDS_KEY, REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "FILE_IMPORT")

        columns = config[COLUMNS_KEY]
        if not isinstance(columns, dict):
            raise ValueError("'columns' must be a dictionary for multi-dependency files")

        if PROPERTY_KEY not in columns:
            raise ValueError(f"'{COLUMNS_KEY}' must contain '{PROPERTY_KEY}' key")

        # Validate dependencies if present
        if DEPENDENCIES_KEY in config:
            dependencies = config[DEPENDENCIES_KEY]
            for dep in dependencies:
                if dep not in columns:
                    raise ValueError(f"Column mapping missing for dependency '{dep}'")

    @staticmethod
    def _validate_constant_value(prop_name: str, val: Any) -> None:
        try:
            float(val)
        except (ValueError, TypeError):
            raise ValueError(f"'{prop_name}' could not be converted to a float. Invalid value: '{val}'")

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
