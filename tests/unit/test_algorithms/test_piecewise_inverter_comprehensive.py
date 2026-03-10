"""Comprehensive edge case tests for piecewise inverter."""

import pytest
import sympy as sp
from materforge.algorithms.piecewise_inverter import PiecewiseInverter
from materforge.core.materials import Material

def _make_material(name: str = "TestMaterial", **scalar_props) -> Material:
    """Construct a minimal Material with dynamic scalar properties."""
    mat = Material(name=name)
    for k, v in scalar_props.items():
        setattr(mat, k, v)
    return mat


class TestPiecewiseInverterComprehensive:
    """Comprehensive edge case tests for piecewise function inversion."""

    def test_create_inverse_non_piecewise_input(self):
        """Error raised for non-piecewise input."""
        T, E = sp.Symbol('T'), sp.Symbol('E')
        with pytest.raises(ValueError, match="Expected Piecewise function"):
            PiecewiseInverter.create_inverse(2*T + 100, T, E)

    def test_create_inverse_non_linear_piece(self):
        """Error raised for quadratic pieces."""
        T, E = sp.Symbol('T'), sp.Symbol('E')
        piecewise_func = sp.Piecewise((T**2 + 100, T < 500), (3*T - 400, True))
        with pytest.raises(ValueError, match="degree 2.*Only linear functions"):
            PiecewiseInverter.create_inverse(piecewise_func, T, E)

    def test_create_inverse_single_constant_piece(self):
        """Error raised when SymPy simplifies the expression to a Float."""
        T, E = sp.Symbol('T'), sp.Symbol('E')
        simplified_expr = sp.Piecewise((150.0, True))
        with pytest.raises(ValueError, match="Expected Piecewise function"):
            PiecewiseInverter.create_inverse(simplified_expr, T, E)

    def test_create_inverse_zero_slope_piece(self):
        """Error raised for constant (zero-slope) piecewise."""
        T, E = sp.Symbol('T'), sp.Symbol('E')
        simplified_expr = sp.Piecewise((100, T < 500), (100, True))
        with pytest.raises(ValueError, match="Expected Piecewise function"):
            PiecewiseInverter.create_inverse(simplified_expr, T, E)

    def test_create_inverse_very_small_slope(self):
        """Error raised for slopes below the numerical stability threshold."""
        T, E = sp.Symbol('T'), sp.Symbol('E')
        piecewise_func = sp.Piecewise((1e-15*T + 100, T < 500), (200, True))
        with pytest.raises(ValueError, match="too small for stable inversion"):
            PiecewiseInverter.create_inverse(piecewise_func, T, E, tolerance=1e-12)

    def test_create_inverse_complex_boundary_conditions(self):
        """Inversion succeeds for simple two-piece linear function."""
        T, E = sp.Symbol('T'), sp.Symbol('E')
        piecewise_func = sp.Piecewise((2*T + 100, T < 500), (3*T - 400, True))
        inverse = PiecewiseInverter.create_inverse(piecewise_func, T, E)
        assert isinstance(inverse, sp.Piecewise)

    def test_extract_boundary_edge_cases(self):
        """Boundary extraction handles both strict and non-strict inequalities."""
        T = sp.Symbol('T')
        assert PiecewiseInverter._extract_boundary(T < 500, T) == 500.0
        assert PiecewiseInverter._extract_boundary(T <= 500, T) == 500.0

    def test_extract_boundary_invalid_condition(self):
        """Error raised for non-relational boundary conditions."""
        T = sp.Symbol('T')
        with pytest.raises(ValueError, match="Cannot extract boundary"):
            PiecewiseInverter._extract_boundary(sp.sin(T), T)

    def test_invert_linear_expression_constant(self):
        """Constant expression returns its numeric value."""
        inverter = PiecewiseInverter(tolerance=1e-12)
        T, E = sp.Symbol('T'), sp.Symbol('E')
        result = inverter._invert_linear_expression(sp.Float(42), T, E)
        assert result == 42.0

    def test_invert_linear_expression_negative_slope(self):
        """Linear expression with negative slope inverts correctly."""
        inverter = PiecewiseInverter(tolerance=1e-12)
        T, E = sp.Symbol('T'), sp.Symbol('E')
        result = inverter._invert_linear_expression(-2*T + 100, T, E)
        assert result.equals((E - 100) / (-2))

    def test_invert_linear_expression_unsupported_degree(self):
        """Error raised for cubic (degree 3) expressions."""
        inverter = PiecewiseInverter()
        T, E = sp.Symbol('T'), sp.Symbol('E')
        with pytest.raises(ValueError, match="degree 3.*only linear expressions"):
            inverter._invert_linear_expression(T**3 + 2*T**2 + T + 1, T, E)
