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


class PiecewiseBuilder:
    """Centralized piecewise function creation with different strategies."""

    @staticmethod
    def build_from_data(temp_array: np.ndarray, prop_array: np.ndarray,
                        T: sp.Symbol, config: Dict, prop_name: str) -> sp.Piecewise:
        """
        Main entry point for data-based piecewise creation.
        Args:
            temp_array: Temperature data points
            prop_array: Property values corresponding to temperatures
            T: Temperature symbol for the piecewise function
            config: Configuration dictionary containing bounds and regression settings
            prop_name: Name of the property (for logging and error messages)
        Returns:
            sp.Piecewise: Symbolic piecewise function
        """
        logger.info("Building piecewise function for property: %s", prop_name)
        logger.debug("Input data - temperature points: %d, property points: %d",
                     len(temp_array) if temp_array is not None else 0,
                     len(prop_array) if prop_array is not None else 0)
        # Validate input arrays
        if temp_array is None or prop_array is None:
            logger.error("Null arrays provided for property '%s'", prop_name)
            raise ValueError(f"Temperature and property arrays cannot be None for '{prop_name}'")
        if len(temp_array) != len(prop_array):
            logger.error("Array length mismatch for '%s': temp=%d, prop=%d",
                         prop_name, len(temp_array), len(prop_array))
            raise ValueError(f"Temperature and property arrays must have same length for '{prop_name}'")
        if len(temp_array) == 0:
            logger.error("Empty arrays provided for property '%s'", prop_name)
            raise ValueError(f"Empty data arrays provided for '{prop_name}'")
        try:
            # Ensure ascending order (handles both ascending and descending input)
            original_order = "ascending" if temp_array[0] < temp_array[-1] else "descending"
            temp_array, prop_array = ensure_ascending_order(temp_array, prop_array)
            logger.debug("Data order for '%s': %s (reordered if needed)", prop_name, original_order)
            # Extract configuration
            lower_bound_type, upper_bound_type = config[BOUNDS_KEY]
            logger.debug("Boundary types for '%s': lower=%s, upper=%s",
                         prop_name, lower_bound_type, upper_bound_type)
            T_standard = sp.Symbol('T')
            # Check for regression configuration
            has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
                config, prop_name, len(temp_array))
            if has_regression:
                logger.info("Regression enabled for '%s': type=%s, degree=%d, segments=%d",
                            prop_name, simplify_type, degree, segments)
            else:
                logger.debug("No regression configured for property '%s'", prop_name)
            # Create piecewise function based on regression settings
            if has_regression and simplify_type == PRE_KEY:
                logger.debug("Building piecewise with pre-regression for '%s'", prop_name)
                pw_result = PiecewiseBuilder._build_with_regression(
                    temp_array, prop_array, T_standard, lower_bound_type, upper_bound_type,
                    degree, segments)
            else:
                logger.debug("Building piecewise without regression for '%s'", prop_name)
                pw_result = PiecewiseBuilder._build_without_regression(
                    temp_array, prop_array, T_standard, lower_bound_type, upper_bound_type)
            # Handle symbol substitution if needed
            if isinstance(T, sp.Symbol) and str(T) != 'T':
                logger.debug("Substituting symbol T -> %s for property '%s'", T, prop_name)
                pw_result = pw_result.subs(T_standard, T)
            logger.info("Successfully built piecewise function for property: %s", prop_name)
            return pw_result
        except Exception as e:
            logger.error("Failed to build piecewise from data for '%s': %s", prop_name, e, exc_info=True)
            raise ValueError(f"Failed building piecewise from data for '{prop_name}': {str(e)}") from e

    @staticmethod
    def build_from_formulas(temp_points: np.ndarray, equations: List[Union[str, sp.Expr]],
                            T: sp.Symbol, lower_bound_type: str = CONSTANT_KEY,
                            upper_bound_type: str = CONSTANT_KEY) -> sp.Piecewise:
        """
        Create piecewise from symbolic equations.
        Args:
            temp_points: Temperature breakpoints
            equations: List of symbolic expressions (strings or SymPy expressions)
            T: Temperature symbol
            lower_bound_type: Boundary behavior below first breakpoint
            upper_bound_type: Boundary behavior above last breakpoint
        Returns:
            sp.Piecewise: Symbolic piecewise function
        """
        logger.info("Building piecewise function from %d formulas and %d breakpoints",
                    len(equations), len(temp_points))
        logger.debug("Temperature breakpoints: %s", temp_points.tolist() if len(temp_points) <= 10
                     else f"[{temp_points[0]}, ..., {temp_points[-1]}]")
        logger.debug("Boundary types: lower=%s, upper=%s", lower_bound_type, upper_bound_type)
        if len(equations) != len(temp_points) - 1:
            logger.error("Equation count mismatch: %d equations for %d breakpoints",
                         len(equations), len(temp_points))
            raise ValueError(f"Number of equations ({len(equations)}) must be one less than temperature/break points "
                             f"({len(temp_points)})")
        if len(temp_points) < 2:
            logger.error("Insufficient breakpoints: %d (minimum 2 required)", len(temp_points))
            raise ValueError("At least 2 temperature points required for piecewise equations")
        try:
            temp_points = np.asarray(temp_points, dtype=float)
            # Parse equations into SymPy expressions using the provided symbol T
            parsed_equations = []
            for i, eqn_str in enumerate(equations):
                try:
                    logger.debug("Parsing equation %d: %s", i + 1, eqn_str)
                    expr = sp.sympify(eqn_str)
                    # Validate that only T symbol is used
                    free_symbols = expr.free_symbols
                    invalid_symbols = [str(sym) for sym in free_symbols if str(sym) != str(T)]
                    if invalid_symbols:
                        logger.error("Invalid symbols in equation %d '%s': %s (only '%s' allowed)",
                                     i + 1, eqn_str, invalid_symbols, T)
                        raise ValueError(
                            f"Invalid symbols {invalid_symbols} in equation '{eqn_str}'. Only '{T}' is allowed.")
                    parsed_equations.append(expr)
                    logger.debug("Successfully parsed equation %d", i + 1)
                except Exception as e:
                    raise ValueError(f"Failed to parse equation {i + 1}: '{eqn_str}': {e}")
            # Special case: single expression with extrapolation at both ends
            if (len(parsed_equations) == 1 and
                    lower_bound_type == EXTRAPOLATE_KEY and
                    upper_bound_type == EXTRAPOLATE_KEY):
                logger.warning(
                    "Using a single expression with extrapolation at both ends. "
                    "Consider simplifying your YAML definition to use a direct equation."
                )
                # Return as Piecewise for consistency
                result = sp.Piecewise((parsed_equations[0], T >= -sp.oo))  # Always true, but explicit
                logger.debug("Created single-expression piecewise with universal extrapolation")
                return result
            # Build piecewise conditions for multiple equations or different boundary types
            conditions = []
            logger.debug("Building piecewise conditions for %d segments", len(parsed_equations))
            # Handle lower bound
            if lower_bound_type == CONSTANT_KEY:
                const_value = parsed_equations[0].subs(T, temp_points[0])
                conditions.append((const_value, T < temp_points[0]))
                logger.debug("Added lower constant boundary: value=%.3f at T<%.1f",
                             float(const_value), temp_points[0])
            # Handle intervals (including special cases for first and last)
            for i, expr in enumerate(parsed_equations):
                if i == 0 and lower_bound_type == EXTRAPOLATE_KEY:
                    # First segment with extrapolation
                    conditions.append((expr, T < temp_points[i + 1]))
                    logger.debug("Added first segment with extrapolation: T<%.1f", temp_points[i + 1])
                elif i == len(parsed_equations) - 1 and upper_bound_type == EXTRAPOLATE_KEY:
                    # Last segment with extrapolation
                    conditions.append((expr, T >= temp_points[i]))
                    logger.debug("Added last segment with extrapolation: T>=%.1f", temp_points[i])
                else:  # Regular interval
                    conditions.append((expr, sp.And(T >= temp_points[i], T < temp_points[i + 1])))
                    logger.debug("Added regular interval: %.1f<=T<%.1f", temp_points[i], temp_points[i + 1])
            # Handle upper bound
            if upper_bound_type == CONSTANT_KEY:
                const_value = parsed_equations[-1].subs(T, temp_points[-1])
                conditions.append((const_value, T >= temp_points[-1]))
                logger.debug("Added upper constant boundary: value=%.3f at T>=%.1f",
                             float(const_value), temp_points[-1])
            result = sp.Piecewise(*conditions)
            logger.info("Successfully built piecewise function from formulas with %d conditions", len(conditions))
            return result
        except Exception as e:
            logger.error("Failed to build piecewise from formulas: %s", e, exc_info=True)
            raise ValueError(f"Failed building piecewise from formulas: {str(e)}") from e

    @staticmethod
    def _build_without_regression(temp_array: np.ndarray, prop_array: np.ndarray,
                                  T: sp.Symbol, lower: str, upper: str) -> sp.Piecewise:
        """
        Create basic linear interpolation piecewise function.
        Args:
            temp_array: Temperature data points (must be sorted)
            prop_array: Property values
            T: Temperature symbol
            lower: Lower boundary behavior ('constant' or 'extrapolate')
            upper: Upper boundary behavior ('constant' or 'extrapolate')
        Returns:
            sp.Piecewise: Linear interpolation piecewise function
        """
        logger.debug("Building linear interpolation piecewise: %d data points, bounds=(%s,%s)",
                     len(temp_array), lower, upper)
        # Convert arrays to Python floats to ensure SymPy compatibility
        temp_array = [ensure_sympy_compatible(x) for x in temp_array]
        prop_array = [ensure_sympy_compatible(x) for x in prop_array]
        conditions = []
        # Handle lower bound (T < temp_array[0])
        if lower == CONSTANT_KEY:
            lower_expr = prop_array[0]
            logger.debug("Lower bound: constant value %.3f", lower_expr)
        else:  # extrapolate
            if len(temp_array) > 1:
                slope = (prop_array[1] - prop_array[0]) / (temp_array[1] - temp_array[0])
                lower_expr = prop_array[0] + slope * (T - temp_array[0])
                logger.debug("Lower bound: extrapolation with slope %.6f", slope)
            else:
                lower_expr = prop_array[0]
                logger.debug("Lower bound: single point, using constant %.3f", lower_expr)
        conditions.append((lower_expr, T < temp_array[0]))
        # Handle main interpolation segments
        logger.debug("Creating %d interpolation segments", len(temp_array) - 1)
        for i in range(len(temp_array) - 1):
            slope = (prop_array[i + 1] - prop_array[i]) / (temp_array[i + 1] - temp_array[i])
            expr = prop_array[i] + slope * (T - temp_array[i])
            condition = sp.And(T >= temp_array[i], T < temp_array[i + 1])
            conditions.append((expr, condition))
            # logger.debug("Segment %d: T∈[%.1f,%.1f), slope=%.6f",
            #              i + 1, temp_array[i], temp_array[i + 1], slope)  #TODO: Uncomment for detailed output
        # Handle upper bound (T >= temp_array[-1])
        if upper == CONSTANT_KEY:
            upper_expr = prop_array[-1]
            logger.debug("Upper bound: constant value %.3f", upper_expr)
        else:  # extrapolate
            if len(temp_array) > 1:
                slope = (prop_array[-1] - prop_array[-2]) / (temp_array[-1] - temp_array[-2])
                upper_expr = prop_array[-1] + slope * (T - temp_array[-1])
                logger.debug("Upper bound: extrapolation with slope %.6f", slope)
            else:
                upper_expr = prop_array[-1]
                logger.debug("Upper bound: single point, using constant %.3f", upper_expr)
        conditions.append((upper_expr, T >= temp_array[-1]))
        logger.debug("Created piecewise function with %d total conditions", len(conditions))
        return sp.Piecewise(*conditions)

    @staticmethod
    def _build_with_regression(temp_array: np.ndarray, prop_array: np.ndarray,
                               T: sp.Symbol, lower: str, upper: str,
                               degree: int, segments: int) -> sp.Piecewise:
        """
        Create piecewise with regression.
        This delegates to RegressionProcessor but provides a unified interface.
        Args:
            temp_array: Temperature data points
            prop_array: Property values
            T: Temperature symbol
            lower: Lower boundary behavior
            upper: Upper boundary behavior
            degree: Polynomial degree for regression
            segments: Number of segments for regression
        Returns:
            sp.Piecewise: Regression-based piecewise function
        """
        logger.info("Building piecewise with regression: degree=%d, segments=%d", degree, segments)
        logger.debug("Regression parameters: lower=%s, upper=%s, data_points=%d",
                     lower, upper, len(temp_array))
        try:
            result = RegressionProcessor.process_regression(temp_array, prop_array, T, lower, upper, degree, segments,
                                                            seed=ProcessingConstants.DEFAULT_REGRESSION_SEED)
            logger.info("Successfully created regression-based piecewise function")
            return result
        except Exception as e:
            logger.error("Regression processing failed: %s", e, exc_info=True)
            raise
