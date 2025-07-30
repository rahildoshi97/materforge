import logging
from typing import Tuple, Union

import pwlf
import sympy as sp

from materforge.data.constants import ProcessingConstants
from materforge.parsing.config.yaml_keys import (REGRESSION_KEY, SIMPLIFY_KEY, DEGREE_KEY, SEGMENTS_KEY, EXTRAPOLATE_KEY, CONSTANT_KEY)
from materforge.parsing.utils.utilities import ensure_sympy_compatible

logger = logging.getLogger(__name__)


class RegressionProcessor:
    """Handles all regression-related functionality."""

    @staticmethod
    def process_regression_params(prop_config: dict, prop_name: str, data_length: int) \
            -> Tuple[bool, Union[str, None], Union[int, None], Union[int, None]]:
        """Process regression parameters from configuration."""
        logger.debug("Processing regression parameters for property: %s", prop_name)
        has_regression = isinstance(prop_config, dict) and REGRESSION_KEY in prop_config
        if not has_regression:
            logger.debug("No regression configuration found for property: %s", prop_name)
            return False, None, None, None
        try:
            regression_config = prop_config[REGRESSION_KEY]
            simplify_type = regression_config[SIMPLIFY_KEY]
            degree = regression_config[DEGREE_KEY]
            segments = regression_config[SEGMENTS_KEY]
            logger.info("Regression config for '%s': type=%s, degree=%d, segments=%d",
                        prop_name, simplify_type, degree, segments)
            # Validation
            if segments >= data_length:
                logger.error("Too many segments for '%s': %d segments >= %d data points",
                             prop_name, segments, data_length)
                raise ValueError(
                    f"Number of segments ({segments}) must be less than number of data points ({data_length})")
            if segments > ProcessingConstants.MAX_REGRESSION_SEGMENTS:
                logger.error("Segments exceed maximum for '%s': %d > %d",
                             prop_name, segments, ProcessingConstants.MAX_REGRESSION_SEGMENTS)
                raise ValueError(
                    f"Number of segments ({segments}) is too high for {prop_name}. "
                    f"Please reduce it to {ProcessingConstants.MAX_REGRESSION_SEGMENTS} or less.")
            elif segments > ProcessingConstants.WARNING_REGRESSION_SEGMENTS:
                logger.warning("High segment count for '%s' (%d) may lead to overfitting",
                               prop_name, segments)
            if degree < 1:
                logger.error("Invalid degree for '%s': %d (must be >= 1)", prop_name, degree)
                raise ValueError(f"Regression degree must be at least 1, got {degree}")
            logger.debug("Regression parameters validated successfully for: %s", prop_name)
            return has_regression, simplify_type, degree, segments
        except KeyError as e:
            logger.error("Missing regression parameter for '%s': %s", prop_name, e)
            raise ValueError(f"Missing regression parameter for '{prop_name}': {str(e)}") from e
        except Exception as e:
            logger.error("Error processing regression parameters for '%s': %s", prop_name, e, exc_info=True)
            raise

    @staticmethod
    def process_regression(temp_array, prop_array, T, lower_bound_type, upper_bound_type,
                           degree, segments, seed=ProcessingConstants.DEFAULT_REGRESSION_SEED):
        """Centralized regression processing logic."""
        logger.info("Starting regression processing: degree=%d, segments=%d, seed=%d",
                    degree, segments, seed)
        logger.debug("Data range: T∈[%.1f, %.1f], prop∈[%.3e, %.3e]",
                     temp_array[0], temp_array[-1],
                     min(prop_array), max(prop_array))
        try:
            # Initialize piecewise linear fit
            logger.debug("Initializing PWLF with degree=%d, seed=%d", degree, seed)
            v_pwlf = pwlf.PiecewiseLinFit(temp_array, prop_array, degree=degree, seed=seed)
            # Perform fitting
            logger.debug("Fitting %d segments to %d data points", segments, len(temp_array))
            fit_result = v_pwlf.fit(n_segments=segments)
            if hasattr(fit_result, 'success') and not fit_result.success:
                logger.warning("PWLF fitting may not have converged optimally")
            # Log fit quality metrics if available
            if hasattr(v_pwlf, 'ssr'):
                logger.debug("Regression fit quality - SSR: %.6e", v_pwlf.ssr)
            # Log breakpoints
            if hasattr(v_pwlf, 'fit_breaks'):
                logger.debug("Fit breakpoints: %s", v_pwlf.fit_breaks.tolist())
            # Create symbolic conditions
            logger.debug("Creating symbolic conditions with bounds: lower=%s, upper=%s",
                         lower_bound_type, upper_bound_type)
            conditions = RegressionProcessor.get_symbolic_conditions(
                v_pwlf, T, lower_bound_type, upper_bound_type)
            result = sp.Piecewise(*conditions)
            logger.info("Successfully completed regression processing with %d conditions", len(conditions))
            return result
        except Exception as e:
            logger.error("Regression processing failed: %s", e, exc_info=True)
            raise ValueError(f"Regression processing failed: {str(e)}") from e

    @staticmethod
    def get_symbolic_conditions(pwlf_: pwlf.PiecewiseLinFit, x: sp.Symbol, lower_: str, upper_: str):
        """Create symbolic conditions for a piecewise function from a pwlf fit."""
        logger.debug("Creating symbolic conditions for %d segments", pwlf_.n_segments)
        fit_breaks = [ensure_sympy_compatible(b) for b in pwlf_.fit_breaks]
        conditions = []
        # Special case: single segment with extrapolation at both ends
        if pwlf_.n_segments == 1 and lower_ == EXTRAPOLATE_KEY and upper_ == EXTRAPOLATE_KEY:
            logger.debug("Single segment with full extrapolation")
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, 1, x)
            conditions.append((eqn, x >= -sp.oo))
            logger.debug("Created universal condition for single segment")
            return conditions
        # Handle lower bound for all cases
        if lower_ == CONSTANT_KEY:
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, 1, x)
            const_value = eqn.evalf(subs={x: fit_breaks[0]})
            conditions.append((const_value, x < fit_breaks[0]))
            logger.debug("Added lower constant boundary: value=%.3f at x<%.1f",
                         float(const_value), fit_breaks[0])
        # Process all segments
        for i in range(pwlf_.n_segments):
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, i + 1, x)
            # First segment with extrapolation
            if i == 0 and lower_ == EXTRAPOLATE_KEY:
                conditions.append((eqn, x < fit_breaks[i + 1]))
                logger.debug("Added first segment with extrapolation: x<%.1f", fit_breaks[i + 1])
            # Last segment with extrapolation
            elif i == pwlf_.n_segments - 1 and upper_ == EXTRAPOLATE_KEY:
                conditions.append((eqn, x >= fit_breaks[i]))
                logger.debug("Added last segment with extrapolation: x>=%.1f", fit_breaks[i])
            else:  # Regular intervals
                conditions.append((eqn, sp.And(x >= fit_breaks[i], x < fit_breaks[i + 1])))
                logger.debug("Added regular segment %d: %.1f<=x<%.1f",
                             i + 1, fit_breaks[i], fit_breaks[i + 1])
        # Handle upper bound
        if upper_ == CONSTANT_KEY:
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, pwlf_.n_segments, x)
            const_value = eqn.evalf(subs={x: fit_breaks[-1]})
            conditions.append((const_value, x >= fit_breaks[-1]))
            logger.debug("Added upper constant boundary: value=%.3f at x>=%.1f",
                         float(const_value), fit_breaks[-1])
        logger.debug("Created %d symbolic conditions total", len(conditions))
        return conditions

    @staticmethod
    def get_symbolic_eqn(pwlf_: pwlf.PiecewiseLinFit, segment_number: int, x: Union[float, sp.Symbol]):
        """Get symbolic equation for a specific segment."""
        logger.debug("Getting symbolic equation for segment %d (degree=%d)", segment_number, pwlf_.degree)
        if pwlf_.degree < 1:
            logger.error("Invalid degree: %d (must be >= 1)", pwlf_.degree)
            raise ValueError('Degree must be at least 1')
        if segment_number < 1 or segment_number > pwlf_.n_segments:
            logger.error("Invalid segment number: %d (valid range: 1-%d)", segment_number, pwlf_.n_segments)
            raise ValueError('segment_number not possible')
        try:
            # Convert NumPy types to Python floats for SymPy compatibility
            beta_values = [ensure_sympy_compatible(b) for b in pwlf_.beta]
            fit_breaks = [ensure_sympy_compatible(b) for b in pwlf_.fit_breaks]
            my_eqn = 0
            # Assemble degree = 1 first
            for line in range(segment_number):
                if line == 0:
                    my_eqn = beta_values[0] + (beta_values[1]) * (x - fit_breaks[0])
                    logger.debug("Base equation: %.6f + %.6f*(x - %.1f)",
                                 beta_values[0], beta_values[1], fit_breaks[0])
                else:
                    my_eqn += (beta_values[line + 1]) * (x - fit_breaks[line])
                    logger.debug("Added linear term %d: %.6f*(x - %.1f)",
                                 line, beta_values[line + 1], fit_breaks[line])
            # Assemble all other degrees
            if pwlf_.degree > 1:
                logger.debug("Adding higher-order terms (degree %d)", pwlf_.degree)
                for k in range(2, pwlf_.degree + 1):
                    for line in range(segment_number):
                        beta_index = pwlf_.n_segments * (k - 1) + line + 1
                        term = (beta_values[beta_index]) * (x - fit_breaks[line]) ** k
                        my_eqn += term
                        logger.debug("Added degree-%d term: %.6f*(x - %.1f)^%d",
                                     k, beta_values[beta_index], fit_breaks[line], k)
            # Only call simplify if x is symbolic
            if isinstance(x, (sp.Symbol, sp.Expr)):
                logger.debug("Simplifying symbolic equation for segment %d", segment_number)
                result = my_eqn.simplify()
            else:
                logger.debug("Returning numeric equation for segment %d", segment_number)
                result = my_eqn
            return result
        except Exception as e:
            logger.error("Error creating symbolic equation for segment %d: %s", segment_number, e, exc_info=True)
            raise ValueError(f"Failed to create equation for segment {segment_number}: {str(e)}") from e
