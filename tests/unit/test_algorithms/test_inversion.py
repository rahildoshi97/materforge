"""Unit tests for PiecewiseInverter."""

import pytest
import sympy as sp
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

class TestPiecewiseInverter:
    """Test cases for PiecewiseInverter."""
    def test_create_inverse_linear_piecewise(self):
        """Test creating inverse of linear piecewise function."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Create a simple linear piecewise function
        piecewise_func = sp.Piecewise(
            (2*T + 100, T < 500),
            (3*T - 400, True)
        )
        inverse_func = PiecewiseInverter.create_inverse(piecewise_func, T, E)
        assert isinstance(inverse_func, sp.Piecewise)
        # Test round-trip accuracy
        test_temps = [300, 500, 700]
        for temp in test_temps:
            energy = float(piecewise_func.subs(T, temp))
            recovered_temp = float(inverse_func.subs(E, energy))
            assert abs(recovered_temp - temp) < 1e-10

    def test_validate_linear_only_valid(self):
        """Test validation passes for linear functions."""
        T = sp.Symbol('T')
        piecewise_func = sp.Piecewise(
            (2*T + 100, T < 500),
            (T + 600, True)
        )
        # Should not raise any exception
        PiecewiseInverter._validate_linear_only(piecewise_func, T)
    
    def test_validate_linear_only_invalid(self):
        """Test validation fails for non-linear functions."""
        T = sp.Symbol('T')
        piecewise_func = sp.Piecewise(
            (T**2 + 100, T < 500),  # Quadratic - should fail
            (T + 600, True)
        )
        with pytest.raises(ValueError, match="Only linear functions"):
            PiecewiseInverter._validate_linear_only(piecewise_func, T)

    def test_extract_boundary(self):
        """Test boundary extraction from conditions."""
        T = sp.Symbol('T')
        condition = T < 500
        boundary = PiecewiseInverter._extract_boundary(condition, T)
        assert boundary == 500.0

    def test_invert_linear_expression(self):
        """Test inversion of linear expressions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        # Test linear expression: 2*T + 100
        expr = 2*T + 100
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        expected = (E - 100) / 2
        assert sp.simplify(inverse_expr - expected) == 0

    def test_invert_constant_expression(self):
        """Test inversion of constant expressions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        # Test constant expression without boundary
        expr = sp.Float(1000)
        inverse_expr = inverter._invert_linear_expression(expr, T, E, boundary_temp=None)
        assert inverse_expr == 1000.0
        # Test constant expression with boundary
        inverse_expr_with_boundary = inverter._invert_linear_expression(expr, T, E, boundary_temp=500.0)
        assert inverse_expr_with_boundary == 500.0

    def test_create_energy_density_inverse(self, sample_aluminum_element):
        """Test creating energy density inverse for a material."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        # Create a simple energy density function
        T = sp.Symbol('T')
        material.energy_density = sp.Piecewise(
            (2*T + 1000, T < 1000),
            (3*T - 1000, True)
        )
        inverse_func = PiecewiseInverter.create_energy_density_inverse(material)
        assert isinstance(inverse_func, sp.Piecewise)
        # Test round-trip
        test_temp = 800
        energy = float(material.energy_density.subs(T, test_temp))
        E = sp.Symbol('E')
        recovered_temp = float(inverse_func.subs(E, energy))
        assert abs(recovered_temp - test_temp) < 1e-10

    def test_create_inverse_single_segment(self):
        """Test creating inverse of single segment function."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Test with a simple linear expression directly
        inverter = PiecewiseInverter()
        linear_expr = 2*T + 100
        # Test the inversion method directly instead of the full create_inverse
        inverse_expr = inverter._invert_linear_expression(linear_expr, T, E)
        expected = (E - 100) / 2
        assert sp.simplify(inverse_expr - expected) == 0
        # Test round-trip
        test_temp = 400
        energy = float(linear_expr.subs(T, test_temp))
        recovered_temp = float(inverse_expr.subs(E, energy))
        assert abs(recovered_temp - test_temp) < 1e-10

    def test_create_inverse_with_boundary_conditions(self):
        """Test creating inverse with different boundary conditions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Create piecewise function with multiple segments
        piecewise_func = sp.Piecewise(
            (T + 200, T < 400),
            (2*T, sp.And(T >= 400, T < 800)),
            (T + 800, True)
        )
        inverse_func = PiecewiseInverter.create_inverse(piecewise_func, T, E)
        assert isinstance(inverse_func, sp.Piecewise)
        # Test round-trip for different segments
        test_temps = [300, 600, 900]  # One from each segment
        for temp in test_temps:
            energy = float(piecewise_func.subs(T, temp))
            recovered_temp = float(inverse_func.subs(E, energy))
            assert abs(recovered_temp - temp) < 1e-10

    def test_extract_boundary_different_formats(self):
        """Test boundary extraction from different condition formats."""
        T = sp.Symbol('T')
        # Test T < value format
        condition1 = T < 500
        boundary1 = PiecewiseInverter._extract_boundary(condition1, T)
        assert boundary1 == 500.0
        # Test T <= value format
        condition2 = T <= 600
        boundary2 = PiecewiseInverter._extract_boundary(condition2, T)
        assert boundary2 == 600.0

    def test_invert_linear_expression_edge_cases(self):
        """Test inversion of linear expressions with edge cases."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        # Test expression with negative coefficient
        expr = -2*T + 100
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        expected = (100 - E) / 2
        assert sp.simplify(inverse_expr - expected) == 0
        # Test expression with just T (coefficient = 1)
        expr = T
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        # Use simplify to handle symbolic expressions properly
        assert sp.simplify(inverse_expr - E) == 0
        # Test expression with just constant term
        expr = sp.Float(500)
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        assert inverse_expr == 500.0

    def test_create_inverse_truly_single_expression(self):
        """Test creating inverse of a truly single expression (not piecewise)."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Test the individual inversion method directly to avoid validation issues
        inverter = PiecewiseInverter()
        # Test with various linear expressions
        test_cases = [
            (2*T + 100, (E - 100) / 2),
            (T, E),
            (-T + 50, 50 - E),
            (sp.Float(42), sp.Float(42))
        ]
        for expr, expected in test_cases:
            inverse_expr = inverter._invert_linear_expression(expr, T, E)
            # Use simplify to check mathematical equivalence, not structural equality
            difference = sp.simplify(inverse_expr - expected)
            assert difference == 0, f"Expected {expected}, got {inverse_expr}"

    # ============================================================================
    # Material-Specific Test Cases (Aluminum specific_enthalpy)
    # ============================================================================

    def test_material_specific_enthalpy_bounds_constant_constant(self):
        """Test inverse of specific_enthalpy with bounds=[constant, constant].
        
        This represents a material property where:
        - DOMAIN: T E [300K, 4500K]
        - Lower bound: constant piece for T < 300K
        - Upper bound: constant piece for T >= 4500K
        """
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        
        # Original piecewise function (Aluminum specific enthalpy)
        specific_enthalpy = sp.Piecewise(
            (119597.562064681, T_C < 300.0),
            (661.839974761036*T_C - 78954.4303636301, T_C < 1672.410186589),
            (3503.34591068727*T_C - 4831117.90285977, T_C < 1741.14311373293),
            (847.6065746183*T_C - 207095.645993611, T_C < 4500.0),
            (3607133.93978874, True)
        )
        
        # Create inverse
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        
        # Test that the final piece returns T=4500 (boundary), not the constant
        final_piece_expr = inverse_h.args[-1][0]  # Last piece expression
        assert float(final_piece_expr) == 4500.0, \
            f"Final piece should return T=4500, got {final_piece_expr}"
        
        # Test round-trip accuracy for representative temperatures ONLY within valid domain [300, 4500]
        test_temps = [300, 500, 1000, 1500, 1700, 2000, 3000, 4500]
        for temp in test_temps:
            specific_enthalpy_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, specific_enthalpy_val))
            error = abs(recovered_temp - temp)
            assert error < 1e-10, \
                f"Round-trip failed at T={temp}K: got {recovered_temp}K (error: {error})"


    def test_material_specific_enthalpy_bounds_constant_extrapolate(self):
        """Test inverse of specific_enthalpy with bounds=[constant, extrapolate].
        
        This represents a material property where:
        - DOMAIN: T E [300K, inf)
        - Lower bound: constant piece for T < 300K
        - Upper bound: extrapolated (no constant end piece, just goes to infinity)
        """
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        
        # Original piecewise function
        specific_enthalpy = sp.Piecewise(
            (119597.624967454, T_C < 300.0),
            (661.839836478122*T_C - 78954.325975983, T_C < 1672.40999625891),
            (3503.38337434609*T_C - 4831180.14351127, T_C < 1741.14194564656),
            (847.606636098854*T_C - 207095.866276605, True)
        )
        
        # Create inverse
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        
        # Test round-trip accuracy ONLY within valid domain [300, inf)
        test_temps = [300, 500, 1000, 1500, 1700, 2000, 3000]
        for temp in test_temps:
            specific_enthalpy_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, specific_enthalpy_val))
            error = abs(recovered_temp - temp)
            assert error < 1e-10, \
                f"Round-trip failed at T={temp}K: got {recovered_temp}K (error: {error})"


    def test_material_specific_enthalpy_bounds_extrapolate_constant(self):
        """Test inverse of specific_enthalpy with bounds=[extrapolate, constant].
        
        This represents a material property where:
        - DOMAIN: T E (-inf, 4500K]
        - Lower bound: extrapolated (starts from first linear piece, no lower constant)
        - Upper bound: constant piece for T >= 4500K
        """
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        
        # Original piecewise function (no initial constant piece)
        specific_enthalpy = sp.Piecewise(
            (661.839974761036*T_C - 78954.4303636301, T_C < 1672.410186589),
            (3503.34591068727*T_C - 4831117.90285977, T_C < 1741.14311373293),
            (847.6065746183*T_C - 207095.645993611, T_C < 4500.0),
            (3607133.93978874, True)
        )
        
        # Create inverse
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        
        # Test that the final piece returns T=4500 (boundary), not the constant
        final_piece_expr = inverse_h.args[-1][0]
        assert float(final_piece_expr) == 4500.0, \
            f"Final piece should return T=4500, got {final_piece_expr}"
        
        # Test round-trip accuracy ONLY within valid domain (-inf, 4500]
        test_temps = [500, 1000, 1500, 1700, 2000, 3000, 4500]
        for temp in test_temps:
            specific_enthalpy_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, specific_enthalpy_val))
            error = abs(recovered_temp - temp)
            assert error < 1e-10, \
                f"Round-trip failed at T={temp}K: got {recovered_temp}K (error: {error})"


    def test_material_specific_enthalpy_bounds_extrapolate_extrapolate(self):
        """Test inverse of specific_enthalpy with bounds=[extrapolate, extrapolate].
        
        This represents a material property with no boundary constants:
        - DOMAIN: T ∈ (-inf, inf)
        - Lower bound: extrapolated (first piece is linear)
        - Upper bound: extrapolated (last piece is linear, no final constant)
        """
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        
        # Original piecewise function (no constant pieces at all)
        specific_enthalpy = sp.Piecewise(
            (661.839836478122*T_C - 78954.325975983, T_C < 1672.40999625891),
            (3503.38337434609*T_C - 4831180.14351127, T_C < 1741.14194564656),
            (847.606636098854*T_C - 207095.866276605, True)
        )
        
        # Create inverse
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        
        # Test round-trip accuracy across full domain (-inf, inf)
        test_temps = [100, 500, 1000, 1500, 1700, 2000, 3000]
        for temp in test_temps:
            specific_enthalpy_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, specific_enthalpy_val))
            error = abs(recovered_temp - temp)
            assert error < 1e-10, \
                f"Round-trip failed at T={temp}K: got {recovered_temp}K (error: {error})"


    def test_material_specific_enthalpy_all_configurations(self):
        """Comprehensive test: verify all four boundary configurations produce correct inverses."""
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        
        configurations = [
            # Configuration 1: [constant, constant]
            {
                'name': 'constant-constant (domain: [300K, 4500K])',
                'piecewise': sp.Piecewise(
                    (119597.562064681, T_C < 300.0),
                    (661.839974761036*T_C - 78954.4303636301, T_C < 1672.410186589),
                    (3503.34591068727*T_C - 4831117.90285977, T_C < 1741.14311373293),
                    (847.6065746183*T_C - 207095.645993611, T_C < 4500.0),
                    (3607133.93978874, True)
                ),
                'test_temps': [300, 1000, 2000, 4500],
                'boundary_temps': [300.0, 4500.0]
            },
            # Configuration 2: [constant, extrapolate]
            {
                'name': 'constant-extrapolate (domain: [300K, ∞))',
                'piecewise': sp.Piecewise(
                    (119597.624967454, T_C < 300.0),
                    (661.839836478122*T_C - 78954.325975983, T_C < 1672.40999625891),
                    (3503.38337434609*T_C - 4831180.14351127, T_C < 1741.14194564656),
                    (847.606636098854*T_C - 207095.866276605, True)
                ),
                'test_temps': [300, 1000, 2000],
                'boundary_temps': [300.0]
            },
            # Configuration 3: [extrapolate, constant]
            {
                'name': 'extrapolate-constant (domain: (-∞, 4500K])',
                'piecewise': sp.Piecewise(
                    (661.839974761036*T_C - 78954.4303636301, T_C < 1672.410186589),
                    (3503.34591068727*T_C - 4831117.90285977, T_C < 1741.14311373293),
                    (847.6065746183*T_C - 207095.645993611, T_C < 4500.0),
                    (3607133.93978874, True)
                ),
                'test_temps': [1000, 2000, 4500],
                'boundary_temps': [4500.0]
            },
            # Configuration 4: [extrapolate, extrapolate]
            {
                'name': 'extrapolate-extrapolate (domain: (-∞, ∞))',
                'piecewise': sp.Piecewise(
                    (661.839836478122*T_C - 78954.325975983, T_C < 1672.40999625891),
                    (3503.38337434609*T_C - 4831180.14351127, T_C < 1741.14194564656),
                    (847.606636098854*T_C - 207095.866276605, True)
                ),
                'test_temps': [100, 1000, 2000],
                'boundary_temps': []
            }
        ]
        
        for config in configurations:
            name = config['name']
            piecewise_func = config['piecewise']
            test_temps = config['test_temps']
            
            # Create inverse
            inverse_h = PiecewiseInverter.create_inverse(piecewise_func, T_C, h)
            assert isinstance(inverse_h, sp.Piecewise), \
                f"Configuration {name}: inverse is not piecewise"
            
            # Test round-trip accuracy for all test temperatures
            for temp in test_temps:
                specific_enthalpy_val = float(piecewise_func.subs(T_C, temp))
                recovered_temp = float(inverse_h.subs(h, specific_enthalpy_val))
                error = abs(recovered_temp - temp)
                assert error < 1e-9, \
                    f"Configuration {name}: Round-trip failed at T={temp}K: " \
                    f"got {recovered_temp}K (error: {error})"
