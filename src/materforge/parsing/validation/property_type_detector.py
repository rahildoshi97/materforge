import logging
import re
from enum import auto, Enum
from typing import Any, Dict, Set

import sympy as sp

from materforge.parsing.config.yaml_keys import (
    FILE_PATH_KEY, DEPENDENCY_COLUMN_KEY, PROPERTY_COLUMN_KEY, BOUNDS_KEY,
    REGRESSION_KEY, DEPENDENCY_KEY, EQUATION_KEY, CONSTANT_KEY,
    EXTRAPOLATE_KEY, SIMPLIFY_KEY, DEGREE_KEY, SEGMENTS_KEY, PRE_KEY, POST_KEY,
    MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY, SOLIDUS_TEMPERATURE_KEY,
    LIQUIDUS_TEMPERATURE_KEY, INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY, VALUE_KEY
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


# --- Main Class ---
class PropertyTypeDetector:
    """Utility class for detecting and validating property types from configuration values."""

    # --- DETECTION RULES ---
    # The order is crucial: more specific patterns must come before general ones.
    DETECTION_RULES = [
        # Unique key checks first (most efficient)
        (lambda c: FILE_PATH_KEY in c, PropertyType.FILE_IMPORT),
        # Patterns sharing keys (order matters)
        (lambda c: DEPENDENCY_KEY in c and VALUE_KEY in c and PropertyTypeDetector._is_step_function(c),
         PropertyType.STEP_FUNCTION),
        (lambda c: DEPENDENCY_KEY in c and VALUE_KEY in c, PropertyType.TABULAR_DATA),
        (lambda c: DEPENDENCY_KEY in c and EQUATION_KEY in c and isinstance(c.get(EQUATION_KEY), list),
         PropertyType.PIECEWISE_EQUATION),
        (lambda c: DEPENDENCY_KEY in c and EQUATION_KEY in c and isinstance(c.get(EQUATION_KEY), str),
         PropertyType.COMPUTED_PROPERTY),
    ]

    # --- Main Public API ---
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

    # --- High-Level Detectors (for DETECTION_RULES) ---
    @staticmethod
    def _is_constant_format(val: Any) -> bool:
        """Checks if the value has the format of a numeric constant."""
        if isinstance(val, int):
            raise ValueError(f"must be defined as a float, not an integer. Use decimal format like '{val}.0'")
        return isinstance(val, float) or (isinstance(val, str) and ('.' in val or 'e' in val.lower()))

    @staticmethod
    def _is_step_function(config: Dict[str, Any]) -> bool:
        """
        A quick, non-validating check if a config looks like a step function.
        A step function has a list of 2 values AND a single temperature point (not a list).
        """
        val_list = config.get(VALUE_KEY)
        temp_def = config.get(DEPENDENCY_KEY)
        is_two_values = isinstance(val_list, list) and len(val_list) == 2
        is_single_temp = not isinstance(temp_def, list)  # Must be a string or number
        return is_two_values and is_single_temp

    # --- Strict Validators (called by the parser) ---
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

    # --- Strict Validators (called by validate_property_config) ---
    @staticmethod
    def _validate_constant_value(prop_name: str, val: Any) -> None:
        try:
            float(val)
        except (ValueError, TypeError):
            raise ValueError(f"'{prop_name}' could not be converted to a float. Invalid value: '{val}'")

    @staticmethod
    def _validate_step_function(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCY_KEY, VALUE_KEY}
        optional = {BOUNDS_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "STEP_FUNCTION")
        if BOUNDS_KEY in config:
            PropertyTypeDetector._check_bounds(config[BOUNDS_KEY])
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list) or len(val_list) != 2:
            raise ValueError(f"'value' for a step function must be a list of exactly two numbers, got {val_list}")
        try:
            float(val_list[0])
            float(val_list[1])
        except (ValueError, TypeError):
            raise ValueError(f"step function values must be numeric, got {val_list}")
        temp_def = config[DEPENDENCY_KEY]
        if isinstance(temp_def, str):  # Check if it's a valid arithmetic expression
            match = re.match(ProcessingConstants.TEMP_ARITHMETIC_REGEX, temp_def.strip())
            if match:  # If it matches, check if the base reference is valid
                base_ref = match.group(1)
                valid_refs = {
                    MELTING_TEMPERATURE_KEY, SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
                    INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY
                }
                if base_ref not in valid_refs:
                    raise ValueError(f"invalid base temperature reference '{base_ref}' in expression '{temp_def}'. "
                                     f"Allowed base references are: {sorted(list(valid_refs))}")
            else:  # If not arithmetic, it must be an exact reference
                valid_refs = {
                    MELTING_TEMPERATURE_KEY, SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
                    INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY
                }
                if temp_def not in valid_refs:
                    raise ValueError(f"invalid temperature reference '{temp_def}'. "
                                     f"Must be a numeric value, a valid transition name, or an arithmetic expression "
                                     f"(e.g., 'melting_temperature + 10').")
        elif not isinstance(temp_def, (int, float, str)):
            raise ValueError(f"'temperature' must be a numeric value or a valid transition reference, got '{temp_def}'")

    @staticmethod
    def _validate_file_import(prop_name: str, config: Dict[str, Any]) -> None:
        required = {FILE_PATH_KEY, DEPENDENCY_COLUMN_KEY, PROPERTY_COLUMN_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "FILE_IMPORT")
        PropertyTypeDetector._check_bounds(config[BOUNDS_KEY])
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])



    @staticmethod
    def _validate_tabular_data(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCY_KEY, VALUE_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "TABULAR_DATA")
        PropertyTypeDetector._check_bounds(config[BOUNDS_KEY])
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])
        temp_def = config[DEPENDENCY_KEY]
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list):
            raise ValueError("'value' for a key-val property must be a list.")
        if isinstance(temp_def, list) and len(temp_def) != len(val_list):
            raise ValueError(f"temperature list (length {len(temp_def)}) and value list (length {len(val_list)}) "
                             f"must have the same length")

    @staticmethod
    def _validate_piecewise_equation(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCY_KEY, EQUATION_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "PIECEWISE_EQUATION")
        PropertyTypeDetector._check_bounds(config[BOUNDS_KEY])
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])
        if not isinstance(config[EQUATION_KEY], list):
            raise ValueError("'equation' for a piecewise equation must be a list of strings")

    @staticmethod
    def _validate_computed_property(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCY_KEY, EQUATION_KEY, BOUNDS_KEY}
        optional = {REGRESSION_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "COMPUTED_PROPERTY")
        PropertyTypeDetector._check_bounds(config[BOUNDS_KEY])
        if REGRESSION_KEY in config:
            PropertyTypeDetector._check_regression(config[REGRESSION_KEY])
        if not isinstance(config[EQUATION_KEY], str):
            raise ValueError("'equation' for a computed property must be a string")
        try:
            sp.sympify(config[EQUATION_KEY])
        except (sp.SympifyError, TypeError) as e:
            raise ValueError(f"invalid mathematical expression in 'equation': {str(e)}")

    # --- Low-Level Validation Helpers ---
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
