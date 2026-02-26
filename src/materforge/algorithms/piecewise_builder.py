# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import Dict, List, Union
import numpy as np
import sympy as sp
from materforge.algorithms.interpolation import ensure_ascending_order
from materforge.algorithms.regression_processor import RegressionProcessor
from materforge.data.constants import ProcessingConstants
from materforge.parsing.config.yaml_keys import CONSTANT_KEY, EXTRAPOLATE_KEY, BOUNDS_KEY, PRE_KEY
from materforge.parsing.utils.utilities import ensure_sympy_compatible

logger = logging.getLogger(__name__)

# The placeholder symbol used in YAML equation strings.
_YAML_PLACEHOLDER = sp.Symbol('T')

class PiecewiseBuilder:
    """Centralised piecewise function creation with different build strategies."""

    @staticmethod
    def build_from_data(dep_array: np.ndarray, prop_array: np.ndarray,
                        dependency: sp.Symbol, config: Dict, prop_name: str) -> sp.Piecewise:
        """Builds a piecewise function from raw dependency–value arrays.

        Handles ascending/descending input, optional pre-regression, and
        symbol substitution when the caller uses a non-placeholder symbol.

        Args:
            dep_array:   Dependency axis data.
            prop_array:  Corresponding property values.
            dependency:  SymPy symbol for the resulting expression.
            config:      Property configuration dict (bounds, regression).
            prop_name:   Property name used in log and error messages.
        Returns:
            Symbolic piecewise function expressed in terms of dependency.
        Raises:
            ValueError: On null/empty/mismatched arrays or build failure.
        """
        logger.info("Building piecewise function for property: %s", prop_name)
        if dep_array is None or prop_array is None:
            raise ValueError(f"dep_array and prop_array cannot be None for '{prop_name}'")
        if len(dep_array) != len(prop_array):
            raise ValueError(f"Array length mismatch for '{prop_name}': "
                f"dep_array={len(dep_array)}, prop_array={len(prop_array)}")
        if len(dep_array) == 0:
            raise ValueError(f"Empty data arrays provided for '{prop_name}'")
        try:
            dep_array, prop_array = ensure_ascending_order(dep_array, prop_array)
            lower_bound_type, upper_bound_type = config[BOUNDS_KEY]
            has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
                config, prop_name, len(dep_array))
            if has_regression:
                logger.info("Regression enabled for %r: type=%s, degree=%d, segments=%d",
                            prop_name, simplify_type, degree, segments)
            if has_regression and simplify_type == PRE_KEY:
                pw_result = PiecewiseBuilder._build_with_regression(
                    dep_array, prop_array, _YAML_PLACEHOLDER,
                    lower_bound_type, upper_bound_type, degree, segments)
            else:
                pw_result = PiecewiseBuilder._build_without_regression(
                    dep_array, prop_array, _YAML_PLACEHOLDER,
                    lower_bound_type, upper_bound_type)
            if dependency != _YAML_PLACEHOLDER:
                logger.debug("Substituting %s -> %s for property %r",
                             _YAML_PLACEHOLDER, dependency, prop_name)
                pw_result = pw_result.subs(_YAML_PLACEHOLDER, dependency)
            logger.info("Successfully built piecewise function for property: %s", prop_name)
            return pw_result
        except Exception as e:
            raise ValueError(f"Failed building piecewise from data for '{prop_name}': {str(e)}") from e

    @staticmethod
    def build_from_formulas(dep_points: np.ndarray, equations: List[Union[str, sp.Expr]],
                            dependency: sp.Symbol, lower_bound_type: str = CONSTANT_KEY,
                            upper_bound_type: str = CONSTANT_KEY) -> sp.Piecewise:
        """Builds a piecewise function from symbolic equations and breakpoints.

        Args:
            dep_points:       Breakpoints defining segment boundaries (n+1 for n equations).
            equations:        Symbolic expressions as strings or SymPy Expr objects.
            dependency:       SymPy symbol used in the equations.
            lower_bound_type: Boundary behaviour below dep_points[0].
            upper_bound_type: Boundary behaviour above dep_points[-1].
        Returns:
            Symbolic piecewise function.
        Raises:
            ValueError: On mismatched counts, invalid symbols, or parse failures.
        """
        logger.info("Building piecewise from %d formulas and %d breakpoints",
                    len(equations), len(dep_points))
        if len(dep_points) < 2:
            raise ValueError("At least 2 breakpoints required for piecewise equations")
        if len(equations) != len(dep_points) - 1:
            raise ValueError(
                f"Number of equations ({len(equations)}) must be one less than breakpoints "
                f"({len(dep_points)})")
        try:
            dep_points = np.asarray(dep_points, dtype=float)
            parsed_equations = []
            for i, eqn_str in enumerate(equations):
                try:
                    expr = sp.sympify(eqn_str)
                    invalid = [str(s) for s in expr.free_symbols if s != dependency]
                    if invalid:
                        raise ValueError(
                            f"Invalid symbols {invalid} in equation '{eqn_str}'. Only '{dependency}' is allowed.")
                    parsed_equations.append(expr)
                except Exception as e:
                    raise ValueError(f"Failed to parse equation {i + 1}: '{eqn_str}': {e}") from e
            # Special case: single expression with extrapolation at both ends
            if (len(parsed_equations) == 1
                    and lower_bound_type == EXTRAPOLATE_KEY
                    and upper_bound_type == EXTRAPOLATE_KEY):
                logger.warning(
                    "Single expression with extrapolation at both ends for dependency %s. "
                    "Consider simplifying to a direct equation.", dependency
                )
                return sp.Piecewise((parsed_equations[0], dependency >= -sp.oo))
            conditions = []
            if lower_bound_type == CONSTANT_KEY:
                const_value = parsed_equations[0].subs(dependency, dep_points[0])
                conditions.append((const_value, dependency < dep_points[0]))
            for i, expr in enumerate(parsed_equations):
                if i == 0 and lower_bound_type == EXTRAPOLATE_KEY:
                    conditions.append((expr, dependency < dep_points[i + 1]))
                elif i == len(parsed_equations) - 1 and upper_bound_type == EXTRAPOLATE_KEY:
                    conditions.append((expr, dependency >= dep_points[i]))
                else:
                    conditions.append((expr, sp.And(dependency >= dep_points[i], dependency < dep_points[i + 1])))
            if upper_bound_type == CONSTANT_KEY:
                const_value = parsed_equations[-1].subs(dependency, dep_points[-1])
                conditions.append((const_value, dependency >= dep_points[-1]))
            logger.info("Successfully built piecewise from formulas with %d conditions", len(conditions))
            return sp.Piecewise(*conditions)
        except Exception as e:
            raise ValueError(f"Failed building piecewise from formulas: {str(e)}") from e

    @staticmethod
    def _build_without_regression(dep_array: np.ndarray, prop_array: np.ndarray,
                                   dependency: sp.Symbol, lower: str, upper: str) -> sp.Piecewise:
        """Builds a linear interpolation piecewise function from data arrays.

        Args:
            dep_array:  Dependency data points (must be sorted ascending).
            prop_array: Property values.
            dependency: SymPy symbol.
            lower:      Lower boundary type.
            upper:      Upper boundary type.
        Returns:
            Linear interpolation piecewise function.
        """
        logger.debug("Building linear interpolation piecewise: %d points, bounds=(%s, %s)",
                     len(dep_array), lower, upper)
        dep_array  = [ensure_sympy_compatible(x) for x in dep_array]
        prop_array = [ensure_sympy_compatible(x) for x in prop_array]
        conditions = []
        # Lower boundary
        if lower == CONSTANT_KEY:
            lower_expr = prop_array[0]
        else:
            if len(dep_array) > 1:
                slope = (prop_array[1] - prop_array[0]) / (dep_array[1] - dep_array[0])
                lower_expr = prop_array[0] + slope * (dependency - dep_array[0])
            else:
                lower_expr = prop_array[0]
        conditions.append((lower_expr, dependency < dep_array[0]))
        # Interior segments
        for i in range(len(dep_array) - 1):
            slope = (prop_array[i + 1] - prop_array[i]) / (dep_array[i + 1] - dep_array[i])
            expr = prop_array[i] + slope * (dependency - dep_array[i])
            conditions.append((expr, sp.And(dependency >= dep_array[i], dependency < dep_array[i + 1])))
        # Upper boundary
        if upper == CONSTANT_KEY:
            upper_expr = prop_array[-1]
        else:
            if len(dep_array) > 1:
                slope = (prop_array[-1] - prop_array[-2]) / (dep_array[-1] - dep_array[-2])
                upper_expr = prop_array[-1] + slope * (dependency - dep_array[-1])
            else:
                upper_expr = prop_array[-1]
        conditions.append((upper_expr, dependency >= dep_array[-1]))
        return sp.Piecewise(*conditions)

    @staticmethod
    def _build_with_regression(dep_array: np.ndarray, prop_array: np.ndarray,
                                dependency: sp.Symbol, lower: str, upper: str,
                                degree: int, segments: int) -> sp.Piecewise:
        """Builds a regression-based piecewise function.

        Delegates to RegressionProcessor, providing a unified build interface.

        Args:
            dep_array:  Dependency data points.
            prop_array: Property values.
            dependency: SymPy symbol.
            lower:      Lower boundary type.
            upper:      Upper boundary type.
            degree:     Polynomial degree for regression.
            segments:   Number of piecewise segments.
        Returns:
            Regression-based piecewise function.
        """
        logger.info("Building regression piecewise: degree=%d, segments=%d", degree, segments)
        return RegressionProcessor.process_regression(dep_array, prop_array, dependency, lower, upper, degree, segments,
            seed=ProcessingConstants.DEFAULT_REGRESSION_SEED)
