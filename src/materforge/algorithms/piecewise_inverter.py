# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import sympy as sp
from typing import Union
import logging

logger = logging.getLogger(__name__)


class PiecewiseInverter:
    """
    Creates inverse functions for linear piecewise functions only.
    Simplified version that supports only degree 1 polynomial.
    """

    def __init__(self, tolerance: float = 1e-12):
        """Initialize the inverter with numerical tolerance."""
        self.tolerance = tolerance
        logger.debug("PiecewiseInverter initialized with tolerance: %.2e", tolerance)

    @staticmethod
    def create_inverse(piecewise_func: Union[sp.Piecewise, sp.Expr],
                       input_symbol: Union[sp.Symbol, sp.Basic],
                       output_symbol: Union[sp.Symbol, sp.Basic],
                       tolerance: float = 1e-12) -> sp.Piecewise:
        """
        Create the inverse of a linear piecewise function.

        This static method provides a convenient interface for creating inverse
        functions without requiring explicit instantiation of PiecewiseInverter.
        Args:
            piecewise_func: The original piecewise function E = f(T)
            input_symbol: Original input symbol (e.g., T)
            output_symbol: Output symbol for inverse function (e.g., E)
            tolerance: Numerical tolerance for inversion stability (default: 1e-12)
        Returns:
            Inverse piecewise function T = f_inv(E)
        Raises:
            ValueError: If any piece has degree > 1
        Examples:
            >>> T = sp.Symbol('T')
            >>> E = sp.Symbol('E')
            >>> piecewise_function = sp.Piecewise((2*T + 100, T < 500), (3*T - 400, True))
            >>> inverse = PiecewiseInverter.create_inverse(piecewise_function, T, E)
        """
        logger.info("Creating inverse function: %s = f_inv(%s)", input_symbol, output_symbol)
        logger.debug("Using tolerance: %.2e", tolerance)
        if not isinstance(piecewise_func, sp.Piecewise):
            logger.error("Input is not a piecewise function: %s", type(piecewise_func))
            raise ValueError(f"Expected Piecewise function, got {type(piecewise_func)}")
        inverter = PiecewiseInverter(tolerance)
        try:
            result = inverter._create_inverse_impl(piecewise_func, input_symbol, output_symbol)
            logger.info("Successfully created inverse function")
            return result
        except Exception as e:
            logger.error("Failed to create inverse function: %s", e, exc_info=True)
            raise

    def _create_inverse_impl(self, piecewise_func: Union[sp.Piecewise, sp.Expr],
                             input_symbol: Union[sp.Symbol, sp.Basic],
                             output_symbol: Union[sp.Symbol, sp.Basic]) -> sp.Piecewise:
        """
        Internal implementation for creating inverse functions.

        This method contains the core logic for inverting piecewise functions,
        separated from the public static interface for better maintainability.
        Args:
            piecewise_func: The original piecewise function
            input_symbol: Original input symbol
            output_symbol: Output symbol for inverse function
        Returns:
            Inverse piecewise function
        """
        num_pieces = len(piecewise_func.args)
        logger.info("Creating inverse for piecewise function with %d pieces", num_pieces)
        # Validate that all pieces are linear
        logger.debug("Validating linearity of all pieces")
        self._validate_linear_only(piecewise_func, input_symbol)
        # Process each piece
        inverse_conditions = []
        piece_boundaries = []  # Track piece information
        # First pass: collect all piece information
        for i, (expr, condition) in enumerate(piecewise_func.args):
            logger.debug("Processing piece %d: expr=%s, condition=%s", i + 1, expr, condition)
            if condition == True:  # Final piece (no condition)
                # For final piece, get boundary from the previous piece
                if piece_boundaries:
                    last_boundary_temp = piece_boundaries[-1]['boundary_temp']
                else:
                    last_boundary_temp = None
                piece_boundaries.append({
                    'index': i,
                    'expr': expr,
                    'condition': condition,
                    'boundary_temp': last_boundary_temp,  # Use previous piece's boundary
                    'boundary_energy': float('inf'),
                    'is_final': True,
                    'slope': None
                })
                logger.debug("Final piece: using boundary from previous piece T=%.3f", last_boundary_temp if last_boundary_temp else 0)
            else:
                # Extract temperature boundary
                try:
                    boundary_temp = self._extract_boundary(condition, input_symbol)
                    logger.debug("Extracted boundary: %.3f", boundary_temp)
                    # Calculate energy at boundary for domain condition
                    boundary_energy = float(expr.subs(input_symbol, boundary_temp))
                    # Get slope to determine condition direction
                    degree = sp.degree(expr, input_symbol)
                    slope = None
                    if degree == 1:
                        coeffs = sp.Poly(expr, input_symbol).all_coeffs()
                        slope = float(coeffs[0])
                    piece_boundaries.append({
                        'index': i,
                        'expr': expr,
                        'condition': condition,
                        'boundary_temp': boundary_temp,
                        'boundary_energy': boundary_energy,
                        'is_final': False,
                        'slope': slope
                    })
                    logger.debug("Boundary: T=%.3f -> E=%.3e (slope=%.6f)", boundary_temp, boundary_energy, slope if slope else 0)
                except Exception as e:
                    logger.error("Error processing piece %d: %s", i + 1, e)
                    raise ValueError(f"Error processing piece {i + 1}: {str(e)}") from e
        # Validate monotonicity
        self._validate_monotonicity(piece_boundaries, input_symbol)
        # Build inverse conditions
        inverse_conditions = []
        # Second pass: build inverse conditions with correct domain mapping
        for i, piece_info in enumerate(piece_boundaries):
            expr = piece_info['expr']
            boundary_temp = piece_info['boundary_temp']
            boundary_energy = piece_info['boundary_energy']
            is_final = piece_info['is_final']
            slope = piece_info['slope']
            # Always pass boundary_temp to _invert_linear_expression for correct handling
            inverse_expr = self._invert_linear_expression(expr, input_symbol, output_symbol, boundary_temp)
            if is_final:
                # Final piece: no upper bound
                inverse_conditions.append((inverse_expr, True))
                logger.debug("Added final piece with universal condition")
            else:
                # Determine condition direction based on slope
                if slope is not None and slope < 0:
                    # Decreasing function: use E > boundary
                    condition = output_symbol > boundary_energy
                    logger.debug("Added piece %d (decreasing): condition E > %.3e", i + 1, boundary_energy)
                else:
                    # Increasing or constant: use E < boundary
                    condition = output_symbol < boundary_energy
                    logger.debug("Added piece %d (increasing): condition E < %.3e", i + 1, boundary_energy)
                inverse_conditions.append((inverse_expr, condition))
        result = sp.Piecewise(*inverse_conditions)
        logger.info("Created inverse function with %d conditions", len(inverse_conditions))
        return result

    def _validate_monotonicity(self, piece_boundaries: list, input_symbol: sp.Symbol) -> None:
        """Validate that the piecewise function is monotonic."""
        logger.debug("Validating monotonicity")
        non_final_pieces = [p for p in piece_boundaries if not p['is_final']]
        if len(non_final_pieces) < 2:
            return
        # Get all slopes
        slopes = []
        for piece in non_final_pieces:
            degree = sp.degree(piece['expr'], input_symbol)
            if degree == 0:
                slopes.append(0)  # Constant
            elif degree == 1:
                coeffs = sp.Poly(piece['expr'], input_symbol).all_coeffs()
                slopes.append(float(coeffs[0]))
        # Check if all slopes have the same sign (all non-negative or all non-positive)
        positive_slopes = [s for s in slopes if s > self.tolerance]
        negative_slopes = [s for s in slopes if s < -self.tolerance]
        if positive_slopes and negative_slopes:
            logger.error("Non-monotonic function detected: mixed positive and negative slopes")
            raise ValueError("Piecewise function is not monotonic. Mix of increasing and decreasing pieces.")
        logger.debug("Monotonicity validation passed. Slopes: %s", slopes)

    @staticmethod
    def _validate_linear_only(piecewise_func: sp.Piecewise, input_symbol: sp.Symbol) -> None:
        """Validate that all pieces are linear (degree <= 1)."""
        logger.debug("Validating linearity for %d pieces", len(piecewise_func.args))
        for i, (expr, condition) in enumerate(piecewise_func.args):
            try:
                degree = sp.degree(expr, input_symbol)
                logger.debug("Piece %d degree: %d", i + 1, degree)
                if degree > 1:
                    logger.error("Non-linear piece found: piece %d has degree %d", i + 1, degree)
                    raise ValueError(
                        f"Piece {i + 1} has degree {degree}. Only linear functions (degree ≤ 1) are supported.")
            except Exception as e:
                logger.error("Error checking degree for piece %d: %s", i + 1, e)
                raise ValueError(f"Error validating piece {i + 1}: {str(e)}") from e
        logger.debug("All pieces validated as linear")

    @staticmethod
    def _extract_boundary(condition: sp.Basic, symbol: sp.Symbol) -> float:
        """Extract the boundary value from a condition like 'T < 300.0'."""
        logger.debug("Extracting boundary from condition: %s", condition)
        try:
            if hasattr(condition, 'rhs'):
                boundary = float(condition.rhs)
                logger.debug("Extracted boundary from rhs: %.3f", boundary)
                return boundary
            elif hasattr(condition, 'args') and len(condition.args) == 2:
                if condition.args[0] == symbol:
                    boundary = float(condition.args[1])
                    logger.debug("Extracted boundary from args[1]: %.3f", boundary)
                    return boundary
                elif condition.args[1] == symbol:
                    boundary = float(condition.args[0])
                    logger.debug("Extracted boundary from args[0]: %.3f", boundary)
                    return boundary
            logger.error("Cannot extract boundary from condition: %s", condition)
            raise ValueError(f"Cannot extract boundary from condition: {condition}")
        except (ValueError, TypeError) as e:
            logger.error("Error extracting boundary value: %s", e)
            raise ValueError(f"Error extracting boundary from {condition}: {str(e)}") from e

    def _invert_linear_expression(self, expr: sp.Expr, input_symbol: sp.Symbol,
                                  output_symbol: sp.Symbol, boundary_temp: float = None) -> sp.Expr:
        """
        Invert a linear expression: ax + b = y → x = (y - b) / a
        Args:
            expr: Linear expression to invert
            input_symbol: Input variable (x)
            output_symbol: Output variable (y)
            boundary_temp: Boundary temperature (for constant pieces)
        Returns:
            Inverted expression
        """
        logger.debug("Inverting expression: %s", expr)
        try:
            degree = sp.degree(expr, input_symbol)
            logger.debug("Expression degree: %d", degree)
            if degree == 0:
                # Constant function - return the boundary temperature (not the constant!)
                # When E = C (constant), the inverse should be T = boundary_temp
                if boundary_temp is not None and boundary_temp != float('inf'):
                    logger.debug("Constant expression: returning boundary temp %.3f", boundary_temp)
                    return sp.sympify(boundary_temp)
                else:
                    # Fallback: return the constant (should not happen for well-formed piecewise)
                    const_value = float(expr)
                    logger.debug("Constant expression (no boundary): returning %.3f", const_value)
                    return sp.sympify(const_value)
            elif degree == 1:
                # Linear function: ax + b = y → x = (y - b) / a
                coeffs = sp.Poly(expr, input_symbol).all_coeffs()
                a, b = float(coeffs[0]), float(coeffs[1])
                logger.debug("Linear coefficients: a=%.6f, b=%.6f", a, b)
                if abs(a) < self.tolerance:
                    logger.error("Linear coefficient too small for inversion: %.2e < %.2e", abs(a), self.tolerance)
                    raise ValueError("Linear coefficient is too small for stable inversion")
                result = (output_symbol - b) / a
                logger.debug("Inverted linear expression: (%s - %.6f) / %.6f", output_symbol, b, a)
                return result
            else:
                logger.error("Unsupported expression degree: %d", degree)
                raise ValueError(f"Expression has degree {degree}, only linear expressions are supported")
        except Exception as e:
            logger.error("Error inverting expression '%s': %s", expr, e, exc_info=True)
            raise ValueError(f"Failed to invert expression {expr}: {str(e)}") from e

    @staticmethod
    def create_energy_density_inverse(material, output_symbol_name: str = 'E') -> sp.Piecewise:
        """
        Create inverse function for energy density: T = f_inv(E)
        Args:
            material: Material object with energy_density property
            output_symbol_name: Symbol name for energy density (default: 'E')
        Returns:
            Inverse piecewise function
        Raises:
            ValueError: If energy density is not linear piecewise
        """
        logger.info("Creating energy density inverse for material: %s", getattr(material, 'name', 'Unknown'))
        logger.debug("Output symbol name: %s", output_symbol_name)
        if not hasattr(material, 'energy_density'):
            logger.error("Material missing energy_density property")
            raise ValueError("Material does not have energy_density property")
        energy_density_func = material.energy_density
        logger.debug("Energy density function type: %s", type(energy_density_func))
        if not isinstance(energy_density_func, sp.Piecewise):
            logger.error("Energy density is not piecewise: %s", type(energy_density_func))
            raise ValueError(f"Energy density must be a piecewise function, found: {type(energy_density_func)}")
        # Extract the temperature symbol
        symbols = energy_density_func.free_symbols
        logger.debug("Found %d free symbols: %s", len(symbols), symbols)
        if len(symbols) != 1:
            logger.error("Energy density function has %d symbols, expected 1: %s", len(symbols), symbols)
            raise ValueError(f"Energy density function must have exactly one symbol, found: {symbols}")
        temp_symbol = list(symbols)[0]
        logger.info("Using temperature symbol: %s (type: %s)", temp_symbol, type(temp_symbol))
        energy_symbol = sp.Symbol(output_symbol_name)
        logger.debug("Created energy symbol: %s", energy_symbol)
        # Create the inverter and generate inverse
        try:
            inverter = PiecewiseInverter()
            inverse_func = inverter.create_inverse(energy_density_func, temp_symbol, energy_symbol)
            logger.info("Successfully created inverse function: %s = f_inv(%s)", temp_symbol, output_symbol_name)
            return inverse_func
        except Exception as e:
            logger.error("Failed to create energy density inverse: %s", e, exc_info=True)
            raise ValueError(f"Failed to create energy density inverse: {str(e)}") from e
