# test_inverse_1.py
import sympy as sp
from typing import List, Tuple
import logging
from pathlib import Path
from pymatlib.parsing.api import create_material

logger = logging.getLogger(__name__)

T = sp.Symbol('T')
E = sp.Symbol('E')

current_file = Path(__file__)
yaml_path = current_file.parent.parent.parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"
ss316l = create_material(yaml_path=yaml_path, T=T, enable_plotting=True)

# Display the energy density function
print(f"Energy Density Function: {ss316l.energy_density}")
print("\n" + "=" * 80)

class PiecewiseInverter:
    """
    Creates inverse functions for piecewise polynomial functions.
    Supports linear (degree 1) and quadratic (degree 2) polynomials.
    """
    def __init__(self, tolerance: float = 1e-10):
        """
        Initialize the inverter.
        Args:
            tolerance: Numerical tolerance for comparisons
        """
        self.tolerance = tolerance

    def create_inverse(self, piecewise_func: sp.Piecewise,
                       input_symbol: sp.Symbol,
                       output_symbol: sp.Symbol) -> sp.Piecewise:
        """
        Create the inverse of a piecewise function.
        Args:
            piecewise_func: The original piecewise function E = f(T)
            input_symbol: Original input symbol (e.g., T)
            output_symbol: Output symbol for inverse function (e.g., E)
        Returns:
            Inverse piecewise function T = f_inv(E)
        """
        logger.info(f"Creating inverse for piecewise function with {len(piecewise_func.args)} pieces")
        # Extract pieces and analyze each segment
        pieces = []
        for i, (expr, condition) in enumerate(piecewise_func.args):
            logger.debug(f"Processing piece {i}: {expr} when {condition}")
            # Skip the final "True" condition for now
            if condition == True:
                continue
            # Extract the boundary from the condition
            boundary = self._extract_boundary(condition, input_symbol)
            # Determine the degree of the polynomial
            degree = sp.degree(expr, input_symbol)
            if degree > 2:
                raise ValueError(f"Unsupported polynomial degree {degree}. Only linear and quadratic are supported.")
            # Create inverse for this piece
            inverse_expr, domain_bounds = self._invert_polynomial(expr, input_symbol, output_symbol, boundary, degree)
            pieces.append((inverse_expr, domain_bounds, boundary))
        # Handle the final piece (True condition)
        if piecewise_func.args and piecewise_func.args[-1][1] == True:
            final_expr = piecewise_func.args[-1][0]
            degree = sp.degree(final_expr, input_symbol)
            if degree <= 2:
                # For the final piece, we need the starting boundary from the previous piece
                if pieces:
                    start_boundary = pieces[-1][2]
                else:
                    start_boundary = float('-inf')
                inverse_expr, domain_bounds = self._invert_polynomial(
                    final_expr, input_symbol, output_symbol, start_boundary, degree, is_final=True)
                pieces.append((inverse_expr, domain_bounds, float('inf')))
        # Build the inverse piecewise function
        return self._build_inverse_piecewise(pieces, output_symbol)

    @staticmethod
    def _extract_boundary(condition: sp.Basic, symbol: sp.Symbol) -> float:
        """Extract the boundary value from a condition like 'T < 300.0'."""
        if hasattr(condition, 'rhs'):
            return float(condition.rhs)
        elif hasattr(condition, 'args') and len(condition.args) == 2:
            # Handle cases like T < 300 or T <= 300
            if condition.args[0] == symbol:
                return float(condition.args[1])
            elif condition.args[1] == symbol:
                return float(condition.args[0])

        raise ValueError(f"Cannot extract boundary from condition: {condition}")

    def _invert_polynomial(self, expr: sp.Expr, input_symbol: sp.Symbol,
                           output_symbol: sp.Symbol, boundary: float,
                           degree: int, is_final: bool = False) -> Tuple[sp.Expr, Tuple[float, float]]:
        """Invert a polynomial expression with improved boundary handling."""
        if degree == 0:
            # Constant function - no inverse
            constant_val = float(expr)
            return output_symbol, (constant_val, constant_val)
        elif degree == 1:
            # Linear case (unchanged)
            coeffs = sp.Poly(expr, input_symbol).all_coeffs()
            a, b = float(coeffs[0]), float(coeffs[1])
            if abs(a) < self.tolerance:
                raise ValueError("Linear coefficient is too small for stable inversion")
            inverse_expr = (output_symbol - b) / a
            # Calculate domain bounds
            if not is_final:
                output_at_boundary = float(expr.subs(input_symbol, boundary))
                if a > 0:
                    domain_bounds = (float('-inf'), output_at_boundary)
                else:
                    domain_bounds = (output_at_boundary, float('inf'))
            else:
                # For final piece, domain extends to infinity
                if a > 0:
                    domain_bounds = (float(expr.subs(input_symbol, boundary)), float('inf'))
                else:
                    domain_bounds = (float('-inf'), float(expr.subs(input_symbol, boundary)))
            return inverse_expr, domain_bounds
        elif degree == 2:
            # Quadratic function: ax² + bx + c = y -> solve for x # IMPROVED QUADRATIC HANDLING
            coeffs = sp.Poly(expr, input_symbol).all_coeffs()
            a, b, c = float(coeffs[0]), float(coeffs[1]), float(coeffs[2])
            if abs(a) < self.tolerance:
                raise ValueError("Quadratic coefficient is too small for stable inversion")
            # Quadratic formula roots
            discriminant = b**2 - 4*a*(c - output_symbol)
            # We need to choose the correct root based on the monotonicity in the interval
            sqrt_discriminant = sp.sqrt(discriminant)
            root1 = (-b + sqrt_discriminant) / (2*a)
            root2 = (-b - sqrt_discriminant) / (2*a)
            # IMPROVED ROOT SELECTION
            if not is_final:
                # Test both roots at the boundary to see which gives the correct value
                boundary_energy = float(expr.subs(input_symbol, boundary))
                try:
                    # Evaluate both roots at the boundary energy
                    root1_at_boundary = float(root1.subs(output_symbol, boundary_energy))
                    root2_at_boundary = float(root2.subs(output_symbol, boundary_energy))
                    # Choose the root that gives a value closest to the boundary temperature
                    if abs(root1_at_boundary - boundary) < abs(root2_at_boundary - boundary):
                        inverse_expr = root1
                    else:
                        inverse_expr = root2
                except (ValueError, TypeError, ZeroDivisionError):
                    # Fallback to derivative-based selection
                    derivative = sp.diff(expr, input_symbol)
                    deriv_at_boundary = float(derivative.subs(input_symbol, boundary))
                    if deriv_at_boundary > 0:  # Increasing function
                        inverse_expr = root1 if a > 0 else root2
                    else:  # Decreasing function
                        inverse_expr = root2 if a > 0 else root1
                # Calculate domain bounds
                output_at_boundary = boundary_energy
                domain_bounds = (float('-inf'), output_at_boundary)
            else:
                # For final piece, use similar logic but test at a point beyond the boundary
                test_temp = boundary + 50.0  # Test point beyond boundary
                test_energy = float(expr.subs(input_symbol, test_temp))
                try:
                    root1_at_test = float(root1.subs(output_symbol, test_energy))
                    root2_at_test = float(root2.subs(output_symbol, test_energy))
                    if abs(root1_at_test - test_temp) < abs(root2_at_test - test_temp):
                        inverse_expr = root1
                    else:
                        inverse_expr = root2
                except (ValueError, TypeError, ZeroDivisionError):
                    inverse_expr = root1  # Default fallback
                boundary_energy = float(expr.subs(input_symbol, boundary))
                domain_bounds = (boundary_energy, float('inf'))
            return inverse_expr, domain_bounds
        else:
            raise ValueError(f"Unsupported polynomial degree: {degree}")

    @staticmethod
    def _build_inverse_piecewise(pieces: List[Tuple[sp.Expr, Tuple[float, float], float]],
                                 output_symbol: sp.Symbol) -> sp.Piecewise:
        """Build the final inverse piecewise function with precise boundary conditions."""
        conditions = []
        for i, (inverse_expr, (min_out, max_out), boundary) in enumerate(pieces):
            if i == len(pieces) - 1:  # Last piece
                conditions.append((inverse_expr, True))
            else:
                # Use the exact boundary energy value for the condition
                if max_out != float('inf') and max_out != float('-inf'):
                    # Use <= to be more inclusive at boundaries
                    conditions.append((inverse_expr, output_symbol <= max_out))
                else:
                    conditions.append((inverse_expr, output_symbol >= min_out))
        return sp.Piecewise(*conditions)

def create_energy_density_inverse(material, output_symbol_name: str = 'E') -> sp.Piecewise:
    """
    Create inverse function for energy density: T = f_inv(E)
    Args:
        material: Material object with energy_density property
        output_symbol_name: Symbol name for energy density (default: 'E')
    Returns:
        Inverse piecewise function
    """
    if not hasattr(material, 'energy_density'):
        raise ValueError("Material does not have energy_density property")
    energy_density_func = material.energy_density
    if not isinstance(energy_density_func, sp.Piecewise):
        raise ValueError("Energy density must be a piecewise function")
    # Extract the temperature symbol from the energy density function
    symbols = energy_density_func.free_symbols
    if len(symbols) != 1:
        raise ValueError(f"Energy density function must have exactly one symbol, found: {symbols}")
    temp_symbol = list(symbols)[0]
    energy_symbol = sp.Symbol(output_symbol_name)
    # Create the inverter and generate inverse
    inverter = PiecewiseInverter()
    inverse_func = inverter.create_inverse(energy_density_func, temp_symbol, energy_symbol)
    logger.info(f"Created inverse function: T = f_inv({output_symbol_name})")
    return inverse_func

# Test implementation
def test_with_real_material():
    """Test with actual SS316L material properties."""
    print("Testing with Real SS316L Material")
    print("=" * 50)
    # Create material with consistent symbol
    T = sp.Symbol('T')
    E = sp.Symbol('E')
    current_file = Path(__file__)
    yaml_path = current_file.parent.parent.parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"
    ss316l = create_material(yaml_path=yaml_path, T=T, enable_plotting=True)
    print(f"Energy Density Function: {ss316l.energy_density}")
    # Create inverter
    inverter = PiecewiseInverter()
    try:
        inverse_energy_density = inverter.create_inverse(ss316l.energy_density, T, E)
        print(f"Inverse function: {inverse_energy_density}")
        # Test comprehensive temperature range
        test_temperatures = [
            # Low temperature region
            0, 50, 150, 250, 299,
            # Boundary transitions
            300, 301, 350, 500,
            # Heating phase
            800, 1000, 1200, 1500, 1602.147, 1605, 1606.5, 1667, 1667.5,
            # Phase transition region
            1668, 1668.5, 1669, 1670, 1685, 1700, 1725, 1734, 1734.5, 1735,
            # High temperature liquid phase
            1800, 2000, 2500, 2999,
            # Ultra-high temperature
            3000, 3500, 4000
        ]
        print("\nRound-trip accuracy test:")
        max_error = 0.0
        for temp in test_temperatures:
            try:
                # Forward: T -> E
                energy = float(ss316l.energy_density.subs(T, temp))
                # Backward: E -> T
                recovered_temp = float(inverse_energy_density.subs(E, energy))
                error = abs(temp - recovered_temp)
                max_error = max(max_error, error)
                status = "✓" if error < 1e-6 else "⚠" if error < 1e-3 else "✗"
                print(f"{status} T={temp:6.1f}K -> E={energy:12.2e} -> T={recovered_temp:6.1f}K, Error={error:.2e}")
            except Exception as e:
                print(f"✗ Error at T={temp}K: {e}")
        print(f"\nMaximum error: {max_error:.2e}")
        if max_error < 1e-6:
            print("✓ EXCELLENT: All tests passed with high precision")
        elif max_error < 1e-3:
            print("✓ GOOD: All tests passed with acceptable precision")
        else:
            print("⚠ WARNING: Some tests have large errors")
    except Exception as e:
        print(f"✗ Failed to create inverse: {e}")

if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    # Run tests
    test_with_real_material()
