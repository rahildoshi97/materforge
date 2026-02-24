"""Unit tests for PiecewiseInverter."""
import pytest
import sympy as sp
from materforge.algorithms.piecewise_inverter import PiecewiseInverter
from materforge.core.materials import Material

class TestPiecewiseInverter:
    """Test cases for PiecewiseInverter."""
    def test_create_inverse_linear_piecewise(self):
        """Test creating inverse of linear piecewise function."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        piecewise_func = sp.Piecewise(
            (2 * T + 100, T < 500),
            (3 * T - 400, True),
        )
        inverse_func = PiecewiseInverter.create_inverse(piecewise_func, T, E)
        assert isinstance(inverse_func, sp.Piecewise)
        test_temps = [300, 500, 700]
        for temp in test_temps:
            energy = float(piecewise_func.subs(T, temp))
            recovered_temp = float(inverse_func.subs(E, energy))
            assert abs(recovered_temp - temp) < 1e-10

    def test_validate_linear_only_valid(self):
        """Test validation passes for linear functions."""
        T = sp.Symbol('T')
        piecewise_func = sp.Piecewise(
            (2 * T + 100, T < 500),
            (T + 600, True),
        )
        PiecewiseInverter._validate_linear_only(piecewise_func, T)
    def test_validate_linear_only_invalid(self):
        """Test validation fails for non-linear functions."""
        T = sp.Symbol('T')
        piecewise_func = sp.Piecewise(
            (T**2 + 100, T < 500),
            (T + 600, True),
        )
        with pytest.raises(ValueError, match="Only linear functions"):
            PiecewiseInverter._validate_linear_only(piecewise_func, T)

    def test_extract_boundary(self):
        """Test boundary extraction from conditions."""
        T = sp.Symbol('T')
        boundary = PiecewiseInverter._extract_boundary(T < 500, T)
        assert boundary == 500.0

    def test_invert_linear_expression(self):
        """Test inversion of linear expressions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        expr = 2 * T + 100
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        expected = (E - 100) / 2
        assert sp.simplify(inverse_expr - expected) == 0

    def test_invert_constant_expression(self):
        """Test inversion of constant expressions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        expr = sp.Float(1000)
        inverse_expr = inverter._invert_linear_expression(expr, T, E, boundary_temp=None)
        assert float(inverse_expr) == 1000.0
        inverse_expr_with_boundary = inverter._invert_linear_expression(expr, T, E, boundary_temp=500.0)
        assert float(inverse_expr_with_boundary) == 500.0

    def test_create_inverse_single_segment(self):
        """Test creating inverse of single segment function."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        linear_expr = 2 * T + 100
        inverse_expr = inverter._invert_linear_expression(linear_expr, T, E)
        expected = (E - 100) / 2
        assert sp.simplify(inverse_expr - expected) == 0
        test_temp = 400
        energy = float(linear_expr.subs(T, test_temp))
        recovered_temp = float(inverse_expr.subs(E, energy))
        assert abs(recovered_temp - test_temp) < 1e-10

    def test_create_inverse_with_boundary_conditions(self):
        """Test creating inverse with different boundary conditions."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        piecewise_func = sp.Piecewise(
            (T + 200, T < 400),
            (2 * T, sp.And(T >= 400, T < 800)),
            (T + 800, True),
        )
        inverse_func = PiecewiseInverter.create_inverse(piecewise_func, T, E)
        assert isinstance(inverse_func, sp.Piecewise)
        test_temps = [300, 600, 900]
        for temp in test_temps:
            energy = float(piecewise_func.subs(T, temp))
            recovered_temp = float(inverse_func.subs(E, energy))
            assert abs(recovered_temp - temp) < 1e-10

    def test_extract_boundary_different_formats(self):
        """Test boundary extraction from different condition formats."""
        T = sp.Symbol('T')
        assert PiecewiseInverter._extract_boundary(T < 500, T) == 500.0
        assert PiecewiseInverter._extract_boundary(T <= 600, T) == 600.0

    def test_invert_linear_expression_edge_cases(self):
        """Test inversion of linear expressions with edge cases."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        expr = -2 * T + 100
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        expected = (100 - E) / 2
        assert sp.simplify(inverse_expr - expected) == 0
        expr = T
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        assert sp.simplify(inverse_expr - E) == 0
        expr = sp.Float(500)
        inverse_expr = inverter._invert_linear_expression(expr, T, E)
        assert float(inverse_expr) == 500.0

    def test_create_inverse_truly_single_expression(self):
        """Test creating inverse of a truly single expression (not piecewise)."""
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        inverter = PiecewiseInverter()
        test_cases = [
            (2 * T + 100, (E - 100) / 2),
            (T, E),
            (-T + 50, 50 - E),
            (sp.Float(42), sp.Float(42)),
        ]
        for expr, expected in test_cases:
            inverse_expr = inverter._invert_linear_expression(expr, T, E)
            assert sp.simplify(inverse_expr - expected) == 0

    # ============================================================================
    # Material-specific enthalpy inverse tests
    # ============================================================================

    def test_material_specific_enthalpy_bounds_constant_constant(self):
        """Inverse of specific_enthalpy with bounds=[constant, constant]."""
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        specific_enthalpy = sp.Piecewise(
            (119597.562064681, T_C < 300.0),
            (661.839974761036 * T_C - 78954.4303636301, T_C < 1672.410186589),
            (3503.34591068727 * T_C - 4831117.90285977, T_C < 1741.14311373293),
            (847.6065746183 * T_C - 207095.645993611, T_C < 4500.0),
            (3607133.93978874, True),
        )
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        final_piece_expr = inverse_h.args[-1][0]
        assert float(final_piece_expr) == 4500.0
        test_temps = [300, 500, 1000, 1500, 1700, 2000, 3000, 4500]
        for temp in test_temps:
            h_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, h_val))
            assert abs(recovered_temp - temp) < 1e-10

    def test_material_specific_enthalpy_bounds_constant_extrapolate(self):
        """Inverse of specific_enthalpy with bounds=[constant, extrapolate]."""
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        specific_enthalpy = sp.Piecewise(
            (119597.624967454, T_C < 300.0),
            (661.839836478122 * T_C - 78954.325975983, T_C < 1672.40999625891),
            (3503.38337434609 * T_C - 4831180.14351127, T_C < 1741.14194564656),
            (847.606636098854 * T_C - 207095.866276605, True),
        )
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        test_temps = [300, 500, 1000, 1500, 1700, 2000, 3000]
        for temp in test_temps:
            h_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, h_val))
            assert abs(recovered_temp - temp) < 1e-10

    def test_material_specific_enthalpy_bounds_extrapolate_constant(self):
        """Inverse of specific_enthalpy with bounds=[extrapolate, constant]."""
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        specific_enthalpy = sp.Piecewise(
            (661.839974761036 * T_C - 78954.4303636301, T_C < 1672.410186589),
            (3503.34591068727 * T_C - 4831117.90285977, T_C < 1741.14311373293),
            (847.6065746183 * T_C - 207095.645993611, T_C < 4500.0),
            (3607133.93978874, True),
        )
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        final_piece_expr = inverse_h.args[-1][0]
        assert float(final_piece_expr) == 4500.0
        test_temps = [500, 1000, 1500, 1700, 2000, 3000, 4500]
        for temp in test_temps:
            h_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, h_val))
            assert abs(recovered_temp - temp) < 1e-10

    def test_material_specific_enthalpy_bounds_extrapolate_extrapolate(self):
        """Inverse of specific_enthalpy with bounds=[extrapolate, extrapolate]."""
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        specific_enthalpy = sp.Piecewise(
            (661.839836478122 * T_C - 78954.325975983, T_C < 1672.40999625891),
            (3503.38337434609 * T_C - 4831180.14351127, T_C < 1741.14194564656),
            (847.606636098854 * T_C - 207095.866276605, True),
        )
        inverse_h = PiecewiseInverter.create_inverse(specific_enthalpy, T_C, h)
        assert isinstance(inverse_h, sp.Piecewise)
        test_temps = [100, 500, 1000, 1500, 1700, 2000, 3000]
        for temp in test_temps:
            h_val = float(specific_enthalpy.subs(T_C, temp))
            recovered_temp = float(inverse_h.subs(h, h_val))
            assert abs(recovered_temp - temp) < 1e-10

    def test_material_specific_enthalpy_all_configurations(self):
        """All four boundary configurations produce correct inverses."""
        T_C = sp.Symbol('T_C')
        h = sp.Symbol('h')
        configurations = [
            {
                'name': 'constant-constant',
                'piecewise': sp.Piecewise(
                    (119597.562064681, T_C < 300.0),
                    (661.839974761036 * T_C - 78954.4303636301, T_C < 1672.410186589),
                    (3503.34591068727 * T_C - 4831117.90285977, T_C < 1741.14311373293),
                    (847.6065746183 * T_C - 207095.645993611, T_C < 4500.0),
                    (3607133.93978874, True),
                ),
                'test_temps': [300, 1000, 2000, 4500],
            },
            {
                'name': 'constant-extrapolate',
                'piecewise': sp.Piecewise(
                    (119597.624967454, T_C < 300.0),
                    (661.839836478122 * T_C - 78954.325975983, T_C < 1672.40999625891),
                    (3503.38337434609 * T_C - 4831180.14351127, T_C < 1741.14194564656),
                    (847.606636098854 * T_C - 207095.866276605, True),
                ),
                'test_temps': [300, 1000, 2000],
            },
            {
                'name': 'extrapolate-constant',
                'piecewise': sp.Piecewise(
                    (661.839974761036 * T_C - 78954.4303636301, T_C < 1672.410186589),
                    (3503.34591068727 * T_C - 4831117.90285977, T_C < 1741.14311373293),
                    (847.6065746183 * T_C - 207095.645993611, T_C < 4500.0),
                    (3607133.93978874, True),
                ),
                'test_temps': [1000, 2000, 4500],
            },
            {
                'name': 'extrapolate-extrapolate',
                'piecewise': sp.Piecewise(
                    (661.839836478122 * T_C - 78954.325975983, T_C < 1672.40999625891),
                    (3503.38337434609 * T_C - 4831180.14351127, T_C < 1741.14194564656),
                    (847.606636098854 * T_C - 207095.866276605, True),
                ),
                'test_temps': [100, 1000, 2000],
            },
        ]
        for cfg in configurations:
            inverse_h = PiecewiseInverter.create_inverse(cfg['piecewise'], T_C, h)
            assert isinstance(inverse_h, sp.Piecewise), f"{cfg['name']}: inverse not piecewise"
            for temp in cfg['test_temps']:
                h_val = float(cfg['piecewise'].subs(T_C, temp))
                recovered_temp = float(inverse_h.subs(h, h_val))
                assert abs(recovered_temp - temp) < 1e-9, (f"{cfg['name']}: round-trip failed at T={temp}, got {recovered_temp}")
