# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
import re
from enum import auto, Enum
from typing import Any, Dict, Set
import sympy as sp

from materforge.parsing.config.yaml_keys import (
    FILE_PATH_KEY, DEPENDENCY_COLUMN_KEY, PROPERTY_COLUMN_KEY, BOUNDS_KEY,
    REGRESSION_KEY, DEPENDENCY_KEY, EQUATION_KEY, CONSTANT_KEY,
    LINEAR_KEY, SIMPLIFY_KEY, DEGREE_KEY, SEGMENTS_KEY, PRE_KEY, POST_KEY, VALUE_KEY,
)
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)

# Valid Python identifier - accepts any user-defined scalar reference.
_IDENTIFIER_RE = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*$')


class PropertyType(Enum):
    CONSTANT_VALUE = auto()
    STEP_FUNCTION = auto()
    FILE_IMPORT = auto()
    TABULAR_DATA = auto()
    PIECEWISE_EQUATION = auto()
    COMPUTED_PROPERTY = auto()
    INVALID = auto()


class PropertyTypeDetector:
    """Detects and validates property types from YAML config values."""

    # Detection rules - order is critical: most specific first.
    DETECTION_RULES = [
        (lambda c: FILE_PATH_KEY in c, PropertyType.FILE_IMPORT),
        (lambda c: DEPENDENCY_KEY in c and VALUE_KEY in c and PropertyTypeDetector._is_step_function(c),
         PropertyType.STEP_FUNCTION),
        (lambda c: DEPENDENCY_KEY in c and VALUE_KEY in c, PropertyType.TABULAR_DATA),
        (lambda c: DEPENDENCY_KEY in c and EQUATION_KEY in c and isinstance(c.get(EQUATION_KEY), list),
         PropertyType.PIECEWISE_EQUATION),
        (lambda c: DEPENDENCY_KEY in c and EQUATION_KEY in c and isinstance(c.get(EQUATION_KEY), str),
         PropertyType.COMPUTED_PROPERTY),
    ]
    # Validator dispatch table - built once at class definition time.
    _VALIDATOR_MAP = None  # populated lazily via _get_validator_map()

    # --- Public API ---
    @staticmethod
    def determine_property_type(prop_name: str, config: Any) -> PropertyType:
        """Determines the PropertyType for a config entry using rule-based detection.

        Args:
            prop_name: Name of the property (used in error messages only).
            config:    Raw config value from the YAML properties block.
        Returns:
            The detected PropertyType.
        Raises:
            ValueError: If config is not numeric, not a dict, or matches no rule.
        """
        logger.debug("Determining property type for '%s'", prop_name)
        if PropertyTypeDetector._is_constant_format(config):
            return PropertyType.CONSTANT_VALUE
        if not isinstance(config, dict):
            raise ValueError(f"Property '{prop_name}' has an invalid format. "
                f"Expected a dictionary or a numeric constant, got {type(config).__name__}.")
        for detector, prop_type in PropertyTypeDetector.DETECTION_RULES:
            if detector(config):
                logger.debug("Detected '%s' as %s", prop_name, prop_type.name)
                return prop_type
        raise ValueError(f"Property '{prop_name}' does not match any known configuration pattern. "
            f"Present keys: {sorted(config.keys())}.")

    @staticmethod
    def validate_property_config(prop_name: str, config: Any,
                                  prop_type: PropertyType) -> None:
        """Performs strict structural validation for the detected property type.

        Args:
            prop_name: Name of the property.
            config:    Raw config value from the YAML properties block.
            prop_type: The already-detected PropertyType.
        Raises:
            ValueError:          If config structure is invalid for the given type.
            NotImplementedError: If no validator exists for prop_type.
        """
        logger.debug("Validating '%s' as %s", prop_name, prop_type.name)
        validator = PropertyTypeDetector._get_validator(prop_type)
        if validator is None:
            raise NotImplementedError(f"No validator registered for property type: {prop_type.name}")
        try:
            validator(prop_name, config)
        except Exception as e:
            raise ValueError(
                f"Invalid configuration for '{prop_name}' (type {prop_type.name}): {str(e)}") from e

    # --- High-level detectors (used in DETECTION_RULES) ---
    @staticmethod
    def _is_constant_format(val: Any) -> bool:
        """Returns True if val is a plain numeric constant (float or float-like string).

        Raises:
            ValueError: If val is a bare integer (must use float notation).
        """
        if isinstance(val, int):
            raise ValueError(
                f"must be defined as a float, not an integer. "
                f"Use decimal format like '{val}.0'")
        return isinstance(val, float) or (
            isinstance(val, str) and ('.' in val or 'e' in val.lower()))

    @staticmethod
    def _is_step_function(config: Dict[str, Any]) -> bool:
        """Returns True if config looks like a step function (non-validating).

        A step function has exactly 2 values and a single (non-list) dependency.
        """
        val_list = config.get(VALUE_KEY)
        dep_def = config.get(DEPENDENCY_KEY)
        return (isinstance(val_list, list)
                and len(val_list) == 2
                and not isinstance(dep_def, list))

    # --- Strict validators ---
    @staticmethod
    def _validate_constant_value(prop_name: str, val: Any) -> None:
        try:
            float(val)
        except (ValueError, TypeError):
            raise ValueError(f"could not be converted to float. Invalid value: '{val}'")

    @staticmethod
    def _validate_step_function(prop_name: str, config: Dict[str, Any]) -> None:
        required = {DEPENDENCY_KEY, VALUE_KEY}
        optional = {BOUNDS_KEY}
        PropertyTypeDetector._check_keys(config, required, optional, "STEP_FUNCTION")
        if BOUNDS_KEY in config:
            PropertyTypeDetector._check_bounds(config[BOUNDS_KEY])
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list) or len(val_list) != 2:
            raise ValueError(f"'value' must be a list of exactly 2 numbers, got {val_list}")
        try:
            float(val_list[0])
            float(val_list[1])
        except (ValueError, TypeError):
            raise ValueError(f"step function values must be numeric, got {val_list}")
        dep_def = config[DEPENDENCY_KEY]
        if isinstance(dep_def, str):
            stripped = dep_def.strip()
            is_identifier = bool(_IDENTIFIER_RE.match(stripped))
            is_arithmetic = bool(re.match(ProcessingConstants.PROPERTY_ARITHMETIC_REGEX, stripped))
            if not is_identifier and not is_arithmetic:
                raise ValueError(f"dependency '{dep_def}' must be a valid property name or arithmetic expression "
                    f"(e.g. 'solidus_temp', 'solidus_temp + 10', 'my_ref - 1')")
        elif not isinstance(dep_def, (int, float)):
            raise ValueError(f"'dependency' must be a numeric value or a property name reference, got '{dep_def}'")

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
        dep_def = config[DEPENDENCY_KEY]
        val_list = config[VALUE_KEY]
        if not isinstance(val_list, list):
            raise ValueError("'value' for a tabular property must be a list")
        if isinstance(dep_def, list) and len(dep_def) != len(val_list):
            raise ValueError(f"dependency list (length {len(dep_def)}) and value list (length {len(val_list)}) "
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
            raise ValueError("'equation' for a piecewise property must be a list of strings")

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

    # --- Validator dispatch ---
    @staticmethod
    def _get_validator(prop_type: PropertyType):
        """Returns the validator function for the given PropertyType, or None."""
        return {
            PropertyType.CONSTANT_VALUE:    PropertyTypeDetector._validate_constant_value,
            PropertyType.STEP_FUNCTION:     PropertyTypeDetector._validate_step_function,
            PropertyType.FILE_IMPORT:       PropertyTypeDetector._validate_file_import,
            PropertyType.TABULAR_DATA:      PropertyTypeDetector._validate_tabular_data,
            PropertyType.PIECEWISE_EQUATION: PropertyTypeDetector._validate_piecewise_equation,
            PropertyType.COMPUTED_PROPERTY: PropertyTypeDetector._validate_computed_property,
        }.get(prop_type)

    # --- Low-level validation helpers ---
    @staticmethod
    def _check_keys(value: Dict[str, Any], required: Set[str], optional: Set[str], context: str) -> None:
        keys = set(value.keys())
        missing = required - keys
        if missing:
            raise ValueError(f"missing required keys for {context}: {sorted(missing)}")
        extra = keys - required - optional
        if extra:
            raise ValueError(f"unexpected keys for {context}: {sorted(extra)}")

    @staticmethod
    def _check_bounds(bounds: Any) -> None:
        if not isinstance(bounds, list) or len(bounds) != 2:
            raise ValueError("'bounds' must be a list of exactly two elements")
        valid = {CONSTANT_KEY, LINEAR_KEY}
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
