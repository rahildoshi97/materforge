# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import List, Optional, Union
import sympy as sp

logger = logging.getLogger(__name__)


class PiecewiseInverter:
    """Creates inverse functions for linear piecewise functions (degree <= 1 only)."""

    def __init__(self, tolerance: float = 1e-12) -> None:
        self.tolerance = tolerance

    @staticmethod
    def create_inverse(piecewise_func: Union[sp.Piecewise, sp.Expr],
        input_symbol: Union[sp.Symbol, sp.Basic],
        output_symbol: Union[sp.Symbol, sp.Basic],
        tolerance: float = 1e-12) -> sp.Piecewise:
        """Creates the inverse of a linear piecewise function.

        Args:
            piecewise_func: Original piecewise function f(input_symbol).
            input_symbol:   Independent variable of the original function.
            output_symbol:  Symbol for the inverse function's argument.
            tolerance:      Numerical tolerance for inversion stability.
        Returns:
            Inverse piecewise function expressed in terms of output_symbol.
        Raises:
            ValueError: If any piece has degree > 1 or the function is non-monotonic.
        Example:
            >>> dep = sp.Symbol("T")
            >>> E   = sp.Symbol("E")
            >>> pw  = sp.Piecewise((2*dep + 100, dep < 500), (3*dep - 400, True))
            >>> inv = PiecewiseInverter.create_inverse(pw, dep, E)
        """
        logger.info("Creating inverse function: %s = f_inv(%s)", input_symbol, output_symbol)
        if not isinstance(piecewise_func, sp.Piecewise):
            raise ValueError(f"Expected Piecewise function, got {type(piecewise_func).__name__}")
        inverter = PiecewiseInverter(tolerance)
        return inverter._create_inverse_impl(piecewise_func, input_symbol, output_symbol)

    # --- Internal implementation ---
    def _create_inverse_impl(self, piecewise_func: sp.Piecewise,
        input_symbol: Union[sp.Symbol, sp.Basic],
        output_symbol: Union[sp.Symbol, sp.Basic]) -> sp.Piecewise:
        """Core inversion logic."""
        logger.info("Inverting piecewise with %d pieces", len(piecewise_func.args))
        self._validate_linear_only(piecewise_func, input_symbol)
        piece_boundaries: List[dict] = []
        for i, (expr, condition) in enumerate(piecewise_func.args):
            if condition is sp.true:
                last_boundary = piece_boundaries[-1]["boundary_val"] if piece_boundaries else None
                piece_boundaries.append({
                    "index":        i,
                    "expr":         expr,
                    "boundary_val": last_boundary,
                    "output_bound": float("inf"),
                    "is_final":     True,
                    "slope":        None,
                })
            else:
                try:
                    boundary_val = self._extract_boundary(condition, input_symbol)
                    output_bound = float(expr.subs(input_symbol, boundary_val))
                    slope = None
                    if sp.degree(expr, input_symbol) == 1:
                        slope = float(sp.Poly(expr, input_symbol).all_coeffs()[0])
                    piece_boundaries.append({
                        "index":        i,
                        "expr":         expr,
                        "boundary_val": boundary_val,
                        "output_bound": output_bound,
                        "is_final":     False,
                        "slope":        slope,
                    })
                except Exception as e:
                    raise ValueError(f"Error processing piece {i + 1}: {str(e)}") from e
        self._validate_monotonicity(piece_boundaries, input_symbol)
        inverse_conditions = []
        for i, piece in enumerate(piece_boundaries):
            inv_expr = self._invert_linear_expression(
                piece["expr"], input_symbol, output_symbol, piece["boundary_val"])
            if piece["is_final"]:
                inverse_conditions.append((inv_expr, True))
            else:
                slope = piece["slope"]
                ob    = piece["output_bound"]
                cond  = (output_symbol > ob) if (slope is not None and slope < 0) else (output_symbol < ob)
                inverse_conditions.append((inv_expr, cond))
        result = sp.Piecewise(*inverse_conditions)
        logger.info("Created inverse with %d conditions", len(inverse_conditions))
        return result

    def _validate_monotonicity(self, piece_boundaries: List[dict], input_symbol: sp.Symbol) -> None:
        """Raises ValueError if the piecewise function is non-monotonic."""
        non_final = [p for p in piece_boundaries if not p["is_final"]]
        if len(non_final) < 2:
            return
        slopes = []
        for piece in non_final:
            deg = sp.degree(piece["expr"], input_symbol)
            if deg == 0:
                slopes.append(0.0)
            elif deg == 1:
                slopes.append(float(sp.Poly(piece["expr"], input_symbol).all_coeffs()[0]))
        positive = [s for s in slopes if s >  self.tolerance]
        negative = [s for s in slopes if s < -self.tolerance]
        if positive and negative:
            raise ValueError("Piecewise function is not monotonic - mix of increasing and decreasing pieces.")
        logger.debug("Monotonicity check passed. Slopes: %s", slopes)

    @staticmethod
    def _validate_linear_only(piecewise_func: sp.Piecewise, input_symbol: sp.Symbol) -> None:
        """Raises ValueError if any piece has degree > 1."""
        for i, (expr, _) in enumerate(piecewise_func.args):
            try:
                deg = sp.degree(expr, input_symbol)
                if deg > 1:
                    raise ValueError(
                        f"Piece {i + 1} has degree {deg}. "
                        f"Only linear functions (degree <= 1) are supported.")
            except Exception as e:
                raise ValueError(f"Error validating piece {i + 1}: {str(e)}") from e

    @staticmethod
    def _extract_boundary(condition: sp.Basic, symbol: sp.Symbol) -> float:
        """Extracts the numeric boundary value from a condition such as dep < 300.0.

        Args:
            condition: A SymPy relational expression.
            symbol:    The independent variable symbol.
        Returns:
            Boundary value as a Python float.
        Raises:
            ValueError: If the boundary cannot be extracted.
        """
        try:
            if hasattr(condition, "rhs"):
                return float(condition.rhs)
            if hasattr(condition, "args") and len(condition.args) == 2:
                lhs, rhs = condition.args
                if lhs == symbol:
                    return float(rhs)
                if rhs == symbol:
                    return float(lhs)
            raise ValueError(f"Cannot extract boundary from condition: {condition}")
        except (ValueError, TypeError) as e:
            raise ValueError(f"Error extracting boundary from {condition}: {str(e)}") from e

    def _invert_linear_expression(self, expr: sp.Expr, input_symbol: sp.Symbol,
        output_symbol: sp.Symbol, boundary_val: Optional[float] = None) -> sp.Expr:
        """Inverts a linear expression: a*x + b = y  =>  x = (y - b) / a.

        For constant pieces (degree 0), returns boundary_val as the inverse.

        Args:
            expr:          Linear expression to invert.
            input_symbol:  Independent variable (x).
            output_symbol: Output variable (y).
            boundary_val:  Boundary of the piece (used for constant pieces).
        Returns:
            Inverted SymPy expression.
        Raises:
            ValueError: If the linear coefficient is near-zero or degree > 1.
        """
        try:
            deg = sp.degree(expr, input_symbol)
            if deg == 0:
                if boundary_val is not None and boundary_val != float("inf"):
                    return sp.sympify(boundary_val)
                return sp.sympify(float(expr))
            if deg == 1:
                coeffs = sp.Poly(expr, input_symbol).all_coeffs()
                a, b = float(coeffs[0]), float(coeffs[1])
                if abs(a) < self.tolerance:
                    raise ValueError(f"Linear coefficient {a:.2e} is too small for stable inversion "
                        f"(tolerance={self.tolerance:.2e})")
                return (output_symbol - b) / a
            raise ValueError(
                f"Expression has degree {deg}; only linear expressions are supported")
        except Exception as e:
            raise ValueError(f"Failed to invert expression {expr}: {str(e)}") from e
