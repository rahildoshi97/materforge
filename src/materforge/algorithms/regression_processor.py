# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import List, Tuple, Union
import pwlf
import sympy as sp
from materforge.data.constants import ProcessingConstants
from materforge.parsing.config.yaml_keys import (
    REGRESSION_KEY, SIMPLIFY_KEY, DEGREE_KEY, SEGMENTS_KEY, LINEAR_KEY, CONSTANT_KEY,
)
from materforge.parsing.utils.utilities import ensure_sympy_compatible

logger = logging.getLogger(__name__)


class RegressionProcessor:
    """Handles all regression-related functionality."""

    @staticmethod
    def process_regression_params(
        prop_config: dict,
        prop_name: str,
        data_length: int,
    ) -> Tuple[bool, Union[str, None], Union[int, None], Union[int, None]]:
        """Extracts and validates regression parameters from a property config dict.

        Args:
            prop_config:  Full property configuration dictionary.
            prop_name:    Property name used in log and error messages.
            data_length:  Number of data points (used to validate segment count).
        Returns:
            Tuple of (has_regression, simplify_type, degree, segments).
            All values are None when has_regression is False.
        Raises:
            ValueError: On invalid or missing regression parameters.
        """
        has_regression = isinstance(prop_config, dict) and REGRESSION_KEY in prop_config
        if not has_regression:
            return False, None, None, None
        try:
            reg = prop_config[REGRESSION_KEY]
            simplify_type = reg[SIMPLIFY_KEY]
            degree        = reg[DEGREE_KEY]
            segments      = reg[SEGMENTS_KEY]
            logger.info("Regression config for %r: type=%s, degree=%d, segments=%d",
                        prop_name, simplify_type, degree, segments)
            if segments >= data_length:
                raise ValueError(
                    f"Number of segments ({segments}) must be less than "
                    f"number of data points ({data_length})")
            if segments > ProcessingConstants.MAX_REGRESSION_SEGMENTS:
                raise ValueError(
                    f"Number of segments ({segments}) exceeds the maximum ({ProcessingConstants.MAX_REGRESSION_SEGMENTS}) for {prop_name!r}.")
            if segments > ProcessingConstants.WARNING_REGRESSION_SEGMENTS:
                logger.warning("High segment count for %r (%d) - may lead to overfitting",
                               prop_name, segments)
            if degree < 1:
                raise ValueError(f"Regression degree must be at least 1, got {degree}")
            return True, simplify_type, degree, segments
        except KeyError as e:
            raise ValueError(f"Missing regression parameter for {prop_name!r}: {str(e)}") from e
        except Exception as e:
            logger.error("Error processing regression parameters for %r: %s", prop_name, e, exc_info=True)
            raise

    @staticmethod
    def process_regression(dep_array,
        prop_array, dependency: sp.Symbol,
        lower_bound_type: str,
        upper_bound_type: str, degree: int,
        segments: int, seed: int = ProcessingConstants.DEFAULT_REGRESSION_SEED,) -> sp.Piecewise:
        """Fits a piecewise polynomial regression and returns a SymPy Piecewise.

        Args:
            dep_array:        Dependency axis data.
            prop_array:       Corresponding property values.
            dependency:       SymPy symbol for the resulting expression.
            lower_bound_type: Boundary behaviour below dep_array[0].
            upper_bound_type: Boundary behaviour above dep_array[-1].
            degree:           Polynomial degree per segment.
            segments:         Number of piecewise segments.
            seed:             Random seed for PWLF optimisation.
        Returns:
            Symbolic piecewise function expressed in terms of dependency.
        Raises:
            ValueError: If PWLF fitting fails.
        """
        logger.info("Starting regression: degree=%d, segments=%d, seed=%d", degree, segments, seed)
        try:
            v_pwlf = pwlf.PiecewiseLinFit(dep_array, prop_array, degree=degree, seed=seed)
            fit_result = v_pwlf.fit(n_segments=segments)
            if hasattr(fit_result, "success") and not fit_result.success:
                logger.warning("PWLF fitting may not have converged optimally")
            if hasattr(v_pwlf, "ssr"):
                logger.debug("Regression SSR: %.6e", v_pwlf.ssr)
            if hasattr(v_pwlf, "fit_breaks"):
                logger.debug("Fit breakpoints: %s", v_pwlf.fit_breaks.tolist())
            conditions = RegressionProcessor.get_symbolic_conditions(
                v_pwlf, dependency, lower_bound_type, upper_bound_type)
            result = sp.Piecewise(*conditions)
            logger.info("Regression complete - %d conditions", len(conditions))
            return result
        except Exception as e:
            raise ValueError(f"Regression processing failed: {str(e)}") from e

    @staticmethod
    def get_symbolic_conditions(pwlf_: pwlf.PiecewiseLinFit,
        dependency: sp.Symbol, lower_: str, upper_: str) -> List[Tuple]:
        """Builds SymPy (expression, condition) pairs from a fitted PWLF model.

        Args:
            pwlf_:      Fitted PiecewiseLinFit instance.
            dependency: SymPy symbol used in the expressions.
            lower_:     Lower boundary type (CONSTANT_KEY or LINEAR_KEY).
            upper_:     Upper boundary type (CONSTANT_KEY or LINEAR_KEY).
        Returns:
            List of (expr, condition) tuples suitable for sp.Piecewise.
        """
        logger.debug("Building symbolic conditions for %d segments", pwlf_.n_segments)
        fit_breaks = [ensure_sympy_compatible(b) for b in pwlf_.fit_breaks]
        conditions = []
        # Special case: single segment with full extrapolation
        if pwlf_.n_segments == 1 and lower_ == LINEAR_KEY and upper_ == LINEAR_KEY:
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, 1, dependency)
            conditions.append((eqn, dependency >= -sp.oo))
            return conditions
        if lower_ == CONSTANT_KEY:
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, 1, dependency)
            const_val = eqn.evalf(subs={dependency: fit_breaks[0]})
            conditions.append((const_val, dependency < fit_breaks[0]))
        for i in range(pwlf_.n_segments):
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, i + 1, dependency)
            if i == 0 and lower_ == LINEAR_KEY:
                conditions.append((eqn, dependency < fit_breaks[i + 1]))
            elif i == pwlf_.n_segments - 1 and upper_ == LINEAR_KEY:
                conditions.append((eqn, dependency >= fit_breaks[i]))
            else:
                conditions.append(
                    (eqn, sp.And(dependency >= fit_breaks[i], dependency < fit_breaks[i + 1])))
        if upper_ == CONSTANT_KEY:
            eqn = RegressionProcessor.get_symbolic_eqn(pwlf_, pwlf_.n_segments, dependency)
            const_val = eqn.evalf(subs={dependency: fit_breaks[-1]})
            conditions.append((const_val, dependency >= fit_breaks[-1]))
        logger.debug("Built %d symbolic conditions", len(conditions))
        return conditions

    @staticmethod
    def get_symbolic_eqn(pwlf_: pwlf.PiecewiseLinFit,
        segment_number: int, dependency: Union[float, sp.Symbol]) -> sp.Expr:
        """Constructs the symbolic polynomial equation for a single PWLF segment.

        Args:
            pwlf_:          Fitted PiecewiseLinFit instance.
            segment_number: 1-based segment index.
            dependency:     SymPy symbol or numeric value.
        Returns:
            SymPy expression for the segment.
        Raises:
            ValueError: On invalid degree or out-of-range segment number.
        """
        if pwlf_.degree < 1:
            raise ValueError("Degree must be at least 1")
        if segment_number < 1 or segment_number > pwlf_.n_segments:
            raise ValueError(f"segment_number {segment_number} out of range (1-{pwlf_.n_segments})")
        try:
            beta       = [ensure_sympy_compatible(b) for b in pwlf_.beta]
            fit_breaks = [ensure_sympy_compatible(b) for b in pwlf_.fit_breaks]
            # Degree-1 base
            eqn = beta[0] + beta[1] * (dependency - fit_breaks[0])
            for line in range(1, segment_number):
                eqn += beta[line + 1] * (dependency - fit_breaks[line])
            # Higher-degree terms
            if pwlf_.degree > 1:
                for k in range(2, pwlf_.degree + 1):
                    for line in range(segment_number):
                        beta_index = pwlf_.n_segments * (k - 1) + line + 1
                        eqn += beta[beta_index] * (dependency - fit_breaks[line]) ** k
            if isinstance(dependency, (sp.Symbol, sp.Expr)):
                return eqn.simplify()
            return eqn
        except Exception as e:
            raise ValueError(f"Failed to create equation for segment {segment_number}: {str(e)}") from e
