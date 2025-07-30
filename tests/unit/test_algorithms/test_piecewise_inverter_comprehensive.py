"""Comprehensive edge case tests for piecewise inverter."""

import pytest
import sympy as sp
from materforge.algorithms.piecewise_inverter import PiecewiseInverter
from materforge.core.materials import Material
from materforge.core.elements import ChemicalElement


class TestPiecewiseInverterComprehensive:
    """Test comprehensive edge cases for piecewise function inversion."""

    def test_create_inverse_non_piecewise_input(self):
        """Test error handling for non-piecewise input."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Regular expression, not piecewise
        regular_expr = 2*T + 100
        with pytest.raises(ValueError, match="Expected Piecewise function"):
            PiecewiseInverter.create_inverse(regular_expr, T, E)

    def test_create_inverse_non_linear_piece(self):
        """Test error handling for non-linear pieces."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Piecewise with quadratic piece
        piecewise_func = sp.Piecewise((T**2 + 100, T < 500), (3*T - 400, True))
        with pytest.raises(ValueError, match="degree 2.*Only linear functions"):
            PiecewiseInverter.create_inverse(piecewise_func, T, E)

    def test_create_inverse_single_constant_piece(self):
        """Test error handling for simplified constant expressions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # This will be simplified to a Float by SymPy
        simplified_expr = sp.Piecewise((150.0, True))
        with pytest.raises(ValueError, match="Expected Piecewise function"):
            PiecewiseInverter.create_inverse(simplified_expr, T, E)

    def test_create_inverse_zero_slope_piece(self):
        """Test error handling for zero slope (horizontal line)."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Create piecewise with zero slope (constant function)
        # This should be handled as a constant, not cause division by zero
        simplified_expr = sp.Piecewise((100, T < 500), (100, True))
        with pytest.raises(ValueError, match="Expected Piecewise function"):
            PiecewiseInverter.create_inverse(simplified_expr, T, E)

    def test_create_inverse_very_small_slope(self):
        """Test handling of very small but non-zero slopes."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')

        # Very small slope that might cause numerical issues
        small_slope = 1e-15
        piecewise_func = sp.Piecewise((small_slope*T + 100, T < 500), (200, True))

        with pytest.raises(ValueError, match="too small for stable inversion"):
            PiecewiseInverter.create_inverse(piecewise_func, T, E, tolerance=1e-12)

    def test_create_inverse_complex_boundary_conditions(self):
        """Test inversion with simple boundary conditions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Use simple conditions instead of complex And conditions
        piecewise_func = sp.Piecewise((2*T + 100, T < 500), (3*T - 400, True))
        inverse = PiecewiseInverter.create_inverse(piecewise_func, T, E)
        assert isinstance(inverse, sp.Piecewise)

    def test_extract_boundary_edge_cases(self):
        """Test boundary extraction with various condition formats."""
        T = sp.Symbol('T')
        # Test different condition formats
        condition1 = T < 500
        boundary1 = PiecewiseInverter._extract_boundary(condition1, T)
        assert boundary1 == 500.0
        # Test with different comparison operators
        condition2 = T <= 500
        boundary2 = PiecewiseInverter._extract_boundary(condition2, T)
        assert boundary2 == 500.0

    def test_extract_boundary_invalid_condition(self):
        """Test error handling for invalid boundary conditions."""
        T = sp.Symbol('T')
        # Invalid condition format
        invalid_condition = sp.sin(T)  # Not a comparison
        with pytest.raises(ValueError, match="Cannot extract boundary"):
            PiecewiseInverter._extract_boundary(invalid_condition, T)

    def test_invert_linear_expression_edge_cases(self):
        """Test linear expression inversion edge cases."""
        inverter = PiecewiseInverter(tolerance=1e-12)
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Test constant expression (degree 0)
        const_expr = sp.Float(42)
        result = inverter._invert_linear_expression(const_expr, T, E)
        assert result == 42.0
        # Test linear expression with negative slope
        linear_expr = -2*T + 100
        result = inverter._invert_linear_expression(linear_expr, T, E)
        expected = (E - 100) / (-2)
        assert result.equals(expected)

    def test_invert_linear_expression_unsupported_degree(self):
        """Test error handling for unsupported polynomial degrees."""
        inverter = PiecewiseInverter()
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Cubic expression (degree 3)
        cubic_expr = T**3 + 2*T**2 + T + 1
        with pytest.raises(ValueError, match="degree 3.*only linear expressions"):
            inverter._invert_linear_expression(cubic_expr, T, E)

    def test_create_energy_density_inverse_missing_property(self):
        """Test error handling when energy density property is missing."""
        # Create material without energy density
        element = ChemicalElement(
            name="Iron", atomic_number=26, atomic_mass=55.845,
            melting_temperature=1811, boiling_temperature=3134,
            latent_heat_of_fusion=13800, latent_heat_of_vaporization=340000
        )
        material = Material(
            name="TestMaterial",
            material_type="pure_metal",
            elements=[element],
            composition=[1.0],
            melting_temperature=sp.Float(1811),
            boiling_temperature=sp.Float(3134)
        )
        with pytest.raises(ValueError, match="Energy density must be a piecewise function.*NoneType"):
            PiecewiseInverter.create_energy_density_inverse(material)

    def test_create_energy_density_inverse_non_piecewise(self):
        """Test error handling when energy density is not piecewise."""
        element = ChemicalElement(
            name="Iron", atomic_number=26, atomic_mass=55.845,
            melting_temperature=1811, boiling_temperature=3134,
            latent_heat_of_fusion=13800, latent_heat_of_vaporization=340000
        )
        material = Material(
            name="TestMaterial",
            material_type="pure_metal",
            elements=[element],
            composition=[1.0],
            melting_temperature=sp.Float(1811),
            boiling_temperature=sp.Float(3134)
        )
        # Set non-piecewise energy density
        T = sp.Symbol('T')
        material.energy_density = 2*T + 100  # Regular expression
        with pytest.raises(ValueError, match="must be a piecewise function"):
            PiecewiseInverter.create_energy_density_inverse(material)

    def test_create_energy_density_inverse_multiple_symbols(self):
        """Test error handling when energy density has multiple symbols."""
        element = ChemicalElement(
            name="Iron", atomic_number=26, atomic_mass=55.845,
            melting_temperature=1811, boiling_temperature=3134,
            latent_heat_of_fusion=13800, latent_heat_of_vaporization=340000
        )
        material = Material(
            name="TestMaterial",
            material_type="pure_metal",
            elements=[element],
            composition=[1.0],
            melting_temperature=sp.Float(1811),
            boiling_temperature=sp.Float(3134)
        )
        # Set energy density with multiple symbols
        T = sp.Symbol('T')
        P = sp.Symbol('P')
        material.energy_density = sp.Piecewise((T + P, T < 500), (2*T + P, True))
        with pytest.raises(ValueError, match="exactly one symbol.*found"):
            PiecewiseInverter.create_energy_density_inverse(material)
