# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
import numpy as np
import re
from typing import List, Union, Optional
from materforge.core.materials import Material
from materforge.parsing.io.data_handler import load_property_data
from materforge.parsing.config.yaml_keys import FILE_PATH_KEY, DEPENDENCY_KEY, VALUE_KEY
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


class DependencyResolver:
    """Handles processing of different dependency definition formats in YAML configurations."""
    EPSILON = ProcessingConstants.DEPENDENCY_EPSILON
    MIN_POINTS = ProcessingConstants.MIN_DEPENDENCY_POINTS

    # --- Public API ---
    @staticmethod
    def resolve_dependency_definition(dep_def: Union[List, str, int, float],
                                      n_values: Optional[int] = None,
                                      material: Optional[Material] = None) -> np.ndarray:
        """Processes different dependency definition formats.

        Args:
            dep_def: Dependency definition. Accepted formats:
                - list:      explicit values, e.g. [300, 400, 500]
                - str tuple: range with step, e.g. "(300, 3000, 50.0)"
                - str tuple: range with points, e.g. "(300, 3000, 71)"
                - str tuple: equidistant, e.g. "(300, 50)" (requires n_values)
                - float/int: single value, e.g. 500.0
                - str:       property reference or arithmetic, e.g. "solidus_temp + 5"
            n_values: Required for the 2-element equidistant tuple format.
            material: Required for string property-name references.
        Returns:
            Processed dependency array.
        """
        if isinstance(dep_def, list):
            return DependencyResolver._resolve_list_format(dep_def, material)
        if isinstance(dep_def, str):
            return DependencyResolver._resolve_string_format(dep_def, n_values, material)
        if isinstance(dep_def, (int, float)):
            return np.array([float(dep_def)], dtype=float)
        raise ValueError(f"Unsupported dependency definition format: {type(dep_def)}")

    @staticmethod
    def extract_from_config(prop_config: dict, material: Material) -> np.ndarray:
        """Extracts the dependency array from a property configuration dict.

        Args:
            prop_config: Property configuration dictionary.
            material:    Material instance for reference resolution.
        Returns:
            Dependency array extracted from the configuration.
        """
        if FILE_PATH_KEY in prop_config:
            try:
                data_array, _ = load_property_data(prop_config)
                return data_array
            except Exception as e:
                raise ValueError(f"Failed to extract dependency array from file: {str(e)}") from e
        if DEPENDENCY_KEY in prop_config:
            dep_def = prop_config[DEPENDENCY_KEY]
            n_values = (len(prop_config[VALUE_KEY]) if VALUE_KEY in prop_config else None)
            return DependencyResolver.resolve_dependency_definition(dep_def, n_values, material)
        raise ValueError("Cannot extract dependency array: no dependency information in config")

    @staticmethod
    def resolve_dependency_reference(dep_ref: Union[str, float, int], material: Material) -> float:
        """Resolves a dependency reference to a concrete float value.

        Handles numeric values, plain property-name references, and arithmetic
        expressions such as "solidus_temp + 5".

        Args:
            dep_ref: Numeric value, a property name, or an arithmetic expression.
            material: Material instance for property-name lookups.
        Returns:
            Resolved float value.
        """
        if isinstance(dep_ref, (int, float)):
            return float(dep_ref)
        if isinstance(dep_ref, str):
            # Numeric string
            try:
                return float(dep_ref)
            except ValueError:
                pass
            # Arithmetic expression: "some_ref + 50" or "some_ref - 10"
            if '+' in dep_ref or '-' in dep_ref:
                match = re.match(ProcessingConstants.PROPERTY_ARITHMETIC_REGEX, dep_ref.strip())
                if match:
                    base_name, operator, offset = match.groups()
                    base = DependencyResolver.get_dependency_value(base_name, material)
                    offset_val = float(offset)
                    return base + offset_val if operator == '+' else base - offset_val
            # Plain property-name reference
            return DependencyResolver.get_dependency_value(dep_ref, material)
        raise ValueError(f"Unsupported dependency reference type: {type(dep_ref)} for value '{dep_ref}'")

    @staticmethod
    def get_dependency_value(dep_ref: Union[str, float, int], material: Material) -> float:
        """Resolves a scalar dependency value from a numeric input or a dynamic
        property on the material.

        Any scalar constant assigned to the material (e.g. solidus_temp,
        my_custom_ref) is resolvable here - no hardcoded name whitelist.

        Args:
            dep_ref: Numeric value, numeric string, or property name.
            material: Material instance to look up named references on.
        Returns:
            Resolved float value.
        Raises:
            ValueError: If dep_ref is not numeric and not found on the material,
                or if the found property is not a scalar constant.
        """
        if isinstance(dep_ref, (int, float)):
            return float(dep_ref)
        if isinstance(dep_ref, str):
            try:
                return float(dep_ref)
            except ValueError:
                pass
            if dep_ref in material.property_names():
                val = getattr(material, dep_ref)
                try:
                    return float(val)
                except (TypeError, ValueError):
                    raise ValueError(f"Dependency reference '{dep_ref}' exists on the material "
                        f"but is not a scalar constant (got {type(val).__name__}). "
                        f"Only CONSTANT_VALUE properties can be used as references.")
            raise ValueError(f"Unknown dependency reference: '{dep_ref}'. "
                f"Available scalar properties: {sorted(material.property_names())}")
        raise ValueError(f"Unsupported dependency value type: {type(dep_ref)}")

    # --- Private helpers ---
    @staticmethod
    def _resolve_list_format(dep_list: List[Union[int, float, str]],
                             material: Optional[Material] = None) -> np.ndarray:
        """Processes an explicit dependency list, resolving string references."""
        try:
            dep_array = []
            for item in dep_list:
                if isinstance(item, str) and material is not None:
                    dep_array.append(DependencyResolver.resolve_dependency_reference(item, material))
                else:
                    dep_array.append(float(item))
            return np.array(dep_array, dtype=float)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid dependency list: {dep_list}\n -> {str(e)}") from e

    @staticmethod
    def _resolve_string_format(dep_str: str,
                               n_values: Optional[int] = None,
                               material: Optional[Material] = None) -> np.ndarray:
        """Dispatches string dependency definitions to the correct sub-handler."""
        if not (dep_str.startswith('(') and dep_str.endswith(')')):
            if material is not None:
                return np.array([DependencyResolver.resolve_dependency_reference(dep_str, material)])
            raise ValueError(f"String dependency '{dep_str}' must be parenthesised or "
                f"a material reference (requires material argument)")
        return DependencyResolver._process_dependency_range_format(dep_str, n_values)

    @staticmethod
    def _process_dependency_range_format(dep_str: str, n_values: Optional[int] = None) -> np.ndarray:
        """Processes parenthesised range/equidistant dependency formats."""
        try:
            content = dep_str.strip('()')
            values = [x.strip() for x in content.split(',')]
            if len(values) == 2:
                return DependencyResolver._resolve_equidistant_format(values, n_values)
            if len(values) == 3:
                return DependencyResolver._resolve_range_format(values)
            raise ValueError(f"Dependency string must have 2 or 3 comma-separated values, got {len(values)}")
        except Exception as e:
            raise ValueError(f"Invalid dependency string format: {dep_str}\n -> {str(e)}") from e

    @staticmethod
    def _resolve_equidistant_format(values: List[str], n_values: Optional[int]) -> np.ndarray:
        """Processes (start, increment) equidistant format."""
        if n_values is None:
            raise ValueError("Number of values required for equidistant format (start, increment)")
        if n_values < DependencyResolver.MIN_POINTS:
            raise ValueError(f"Number of values must be at least {DependencyResolver.MIN_POINTS}, got {n_values}")
        try:
            start, increment = float(values[0]), float(values[1])
            if abs(increment) <= DependencyResolver.EPSILON:
                raise ValueError("Increment cannot be zero")
            return np.array([start + i * increment for i in range(n_values)], dtype=float)
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"Invalid equidistant format: ({values[0]}, {values[1]})\n -> {str(e)}") from e

    @staticmethod
    def _resolve_range_format(values: List[str]) -> np.ndarray:
        """Processes (start, stop, step/points) range format."""
        try:
            start, stop = float(values[0]), float(values[1])
            DependencyResolver._validate_range_endpoints(start, stop)
            fmt = DependencyResolver._determine_format_type(values[2].strip())
            if fmt == "points":
                return DependencyResolver._resolve_points_format(values)
            return DependencyResolver._resolve_step_format(values)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid range format: ({', '.join(values)})\n -> {str(e)}") from e

    @staticmethod
    def _determine_format_type(third_param: str) -> str:
        """Returns 'points' or 'step' based on the third range parameter."""
        try:
            val = float(third_param)
            is_integer_format = (
                '.' not in third_param
                and 'e' not in third_param.lower()
                and val == int(val)
                and 0 < val <= 1000
            )
            return "points" if is_integer_format else "step"
        except (ValueError, TypeError):
            raise ValueError(f"Third parameter must be numeric, got: {third_param}")

    @staticmethod
    def _resolve_step_format(values: List[str]) -> np.ndarray:
        """Processes (start, stop, step) format."""
        try:
            start, stop, step = float(values[0]), float(values[1]), float(values[2])
            if abs(step) <= DependencyResolver.EPSILON:
                raise ValueError("Step cannot be zero")
            if (start < stop and step <= 0) or (start > stop and step >= 0):
                raise ValueError("Step sign must match range direction")
            if abs(step) > abs(stop - start):
                raise ValueError(f"Step ({abs(step)}) is too large for range ({abs(stop - start)})")
            return np.arange(start, stop + step / 2, step)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid step format: ({', '.join(values)})\n -> {str(e)}") from e

    @staticmethod
    def _resolve_points_format(values: List[str]) -> np.ndarray:
        """Processes (start, stop, num_points) format."""
        try:
            start, stop = float(values[0]), float(values[1])
            n_points = int(float(values[2]))
            if n_points < DependencyResolver.MIN_POINTS:
                raise ValueError(f"Number of points must be at least {DependencyResolver.MIN_POINTS}, got {n_points}")
            if n_points > 10000:
                raise ValueError(f"Number of points ({n_points}) exceeds maximum (10000)")
            return np.linspace(start, stop, n_points)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid points format: ({', '.join(values)})\n -> {str(e)}") from e

    @staticmethod
    def _validate_range_endpoints(start: float, stop: float) -> None:
        """Validates that start and stop are distinct."""
        if abs(start - stop) <= DependencyResolver.EPSILON:
            raise ValueError(f"Start and stop must be different, got start={start}, stop={stop}")

    @staticmethod
    def validate_dependency_array(dep_array: np.ndarray, context: str = "") -> None:
        """Validates a dependency array for common issues.

        Args:
            dep_array: Array to validate.
            context:   Optional context string for error messages.
        """
        ctx = f" for {context}" if context else ""
        if len(dep_array) == 0:
            raise ValueError(f"Dependency array is empty{ctx}")
        if len(dep_array) < DependencyResolver.MIN_POINTS:
            raise ValueError(f"Dependency array must have at least {DependencyResolver.MIN_POINTS} points, "
                f"got {len(dep_array)}{ctx}")
        if not np.all(np.isfinite(dep_array)):
            raise ValueError(f"Dependency array contains non-finite values{ctx}")
