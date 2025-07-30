"""Unit tests for PiecewiseBuilder."""

import pytest
import numpy as np
import sympy as sp
from materforge.algorithms.piecewise_builder import PiecewiseBuilder

class TestPiecewiseBuilder:
    """Test cases for PiecewiseBuilder."""
    def test_build_from_data_basic(self, temp_symbol):
        """Test basic piecewise function creation from data."""
        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950, 1000])
        config = {'bounds': ['constant', 'constant']}
        piecewise_func = PiecewiseBuilder.build_from_data(
            temp_array, prop_array, temp_symbol, config, "test_property"
        )
        assert isinstance(piecewise_func, sp.Piecewise)
        # Test evaluation at known points
        assert float(piecewise_func.subs(temp_symbol, 350)) == 925.0  # Interpolated
        assert float(piecewise_func.subs(temp_symbol, 250)) == 900.0  # Constant below
        assert float(piecewise_func.subs(temp_symbol, 600)) == 1000.0  # Constant above

    def test_build_from_data_extrapolation(self, temp_symbol):
        """Test piecewise function with extrapolation boundaries."""
        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950, 1000])
        config = {'bounds': ['extrapolate', 'extrapolate']}
        piecewise_func = PiecewiseBuilder.build_from_data(
            temp_array, prop_array, temp_symbol, config, "test_property"
        )
        # When both boundaries are extrapolate, result might be a simple expression
        assert isinstance(piecewise_func, (sp.Piecewise, sp.Expr))
        # Test extrapolation - the key is that it should work mathematically
        # Below range: slope = (950-900)/(400-300) = 0.5
        # At T=200: 900 + 0.5*(200-300) = 850
        result_200 = float(piecewise_func.subs(temp_symbol, 200))
        assert abs(result_200 - 850.0) < 1e-10  # Use small tolerance for floating point
        # Test within range
        result_350 = float(piecewise_func.subs(temp_symbol, 350))
        assert abs(result_350 - 925.0) < 1e-10
        # Test above range: slope = (1000-950)/(500-400) = 0.5
        # At T=600: 1000 + 0.5*(600-500) = 1050
        result_600 = float(piecewise_func.subs(temp_symbol, 600))
        assert abs(result_600 - 1050.0) < 1e-10

    def test_build_from_formulas(self, temp_symbol):
        """Test piecewise function creation from symbolic formulas."""
        temp_points = np.array([300, 400, 500])
        equations = ["2*T + 300", "T + 550"]
        piecewise_func = PiecewiseBuilder.build_from_formulas(
            temp_points, equations, temp_symbol
        )
        assert isinstance(piecewise_func, sp.Piecewise)
        # Test evaluation
        # At T=350 (first segment): 2*350 + 300 = 1000
        assert float(piecewise_func.subs(temp_symbol, 350)) == 1000.0
        # At T=450 (second segment): 450 + 550 = 1000
        assert float(piecewise_func.subs(temp_symbol, 450)) == 1000.0

    def test_build_from_data_mismatched_arrays(self, temp_symbol):
        """Test error handling for mismatched array lengths."""
        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950])  # Different length
        config = {'bounds': ['constant', 'constant']}
        with pytest.raises(ValueError, match="must have same length"):
            PiecewiseBuilder.build_from_data(
                temp_array, prop_array, temp_symbol, config, "test_property"
            )
    def test_build_from_formulas_mismatched_count(self, temp_symbol):
        """Test error handling for mismatched formula count."""
        temp_points = np.array([300, 400, 500])
        equations = ["2*T + 300"]  # Should have 2 equations for 3 points
        with pytest.raises(ValueError, match="Number of equations"):
            PiecewiseBuilder.build_from_formulas(temp_points, equations, temp_symbol)

    def test_build_from_data_mixed_boundaries(self, temp_symbol):
        """Test piecewise function with mixed boundary conditions."""
        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950, 1000])
        config = {'bounds': ['constant', 'extrapolate']}
        piecewise_func = PiecewiseBuilder.build_from_data(
            temp_array, prop_array, temp_symbol, config, "test_property"
        )
        # This should definitely be a Piecewise function with mixed boundaries
        assert isinstance(piecewise_func, sp.Piecewise)
        # Test constant boundary below
        assert float(piecewise_func.subs(temp_symbol, 250)) == 900.0
        # Test extrapolation above
        result_600 = float(piecewise_func.subs(temp_symbol, 600))
        assert abs(result_600 - 1050.0) < 1e-10  # Should extrapolate
