"""Unit tests for RegressionProcessor."""

import pytest
import numpy as np
import sympy as sp
from materforge.algorithms.regression_processor import RegressionProcessor

class TestRegressionProcessor:
    """Test cases for RegressionProcessor."""
    def test_process_regression_params_valid(self):
        """Test processing valid regression parameters."""
        config = {
            'regression': {
                'simplify': 'pre',
                'degree': 2,
                'segments': 3
            }
        }
        has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
            config, "test_property", 10
        )
        assert has_regression is True
        assert simplify_type == 'pre'
        assert degree == 2
        assert segments == 3

    def test_process_regression_params_no_regression(self):
        """Test processing config without regression."""
        config = {'bounds': ['constant', 'constant']}
        has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
            config, "test_property", 10
        )
        assert has_regression is False
        assert simplify_type is None
        assert degree is None
        assert segments is None

    def test_process_regression_with_pwlf(self):
        """Test regression processing using pwlf."""
        temp_array = np.array([300, 400, 500, 600])
        prop_array = np.array([900, 1000, 1100, 1200])
        T = sp.Symbol('T')
        result = RegressionProcessor.process_regression(
            temp_array, prop_array, T, 'constant', 'constant',
            degree=1, segments=1
        )
        assert isinstance(result, sp.Piecewise)

    def test_get_symbolic_conditions(self):
        """Test symbolic conditions generation."""
        import pwlf
        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950, 1000])
        v_pwlf = pwlf.PiecewiseLinFit(temp_array, prop_array, degree=1)
        v_pwlf.fit(n_segments=1)
        T = sp.Symbol('T')
        conditions = RegressionProcessor.get_symbolic_conditions(
            v_pwlf, T, 'constant', 'constant'
        )
        assert isinstance(conditions, list)
        assert len(conditions) > 0
        for expr, cond in conditions:
            assert isinstance(expr, (sp.Expr, sp.Float, sp.Integer))

    def test_get_symbolic_eqn(self):
        """Test symbolic equation generation."""
        import pwlf

        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950, 1000])
        v_pwlf = pwlf.PiecewiseLinFit(temp_array, prop_array, degree=1)
        v_pwlf.fit(n_segments=1)
        T = sp.Symbol('T')
        eqn = RegressionProcessor.get_symbolic_eqn(v_pwlf, 1, T)
        assert isinstance(eqn, sp.Expr)
        eqn_numeric = RegressionProcessor.get_symbolic_eqn(v_pwlf, 1, 350.0)
        assert isinstance(eqn_numeric, (float, int))

    def test_regression_validation_errors(self):
        """Test validation errors in regression processing."""
        config = {
            'regression': {
                'simplify': 'pre',
                'degree': 2,
                'segments': 15  # Too many segments
            }
        }
        with pytest.raises(ValueError, match="Number of segments"):
            RegressionProcessor.process_regression_params(
                config, "test_property", 10
            )
