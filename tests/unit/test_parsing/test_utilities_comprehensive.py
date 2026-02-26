"""Comprehensive tests for parsing utilities."""

import pytest
import numpy as np
import sympy as sp
from unittest.mock import Mock
from materforge.core.materials import Material
from materforge.parsing.utils.utilities import create_step_visualization_data, handle_numeric_dependency


class TestUtilitiesComprehensive:
    """Edge case tests for utility functions."""

    @pytest.fixture
    def sample_material(self):
        mat = Material(name="TestMaterial")
        mat.melting_temperature = 1811.0
        mat.boiling_temperature = 3134.0
        return mat

    def test_handle_numeric_dependency_symbolic(self, sample_material):
        """Symbolic dependency -> returns False, property not added."""
        processor = Mock()
        processor.processed_properties = set()
        T = sp.Symbol('T')
        piecewise_expr = sp.Piecewise((100, T < 500), (200, True))
        result = handle_numeric_dependency(processor, sample_material, 'test_prop', piecewise_expr, T)
        assert result is False
        assert 'test_prop' not in processor.processed_properties

    def test_handle_numeric_dependency_numeric_success(self, sample_material):
        """Numeric dependency -> returns True, property assigned to material."""
        processor = Mock()
        processor.processed_properties = set()
        T = sp.Symbol('T')
        piecewise_expr = sp.Piecewise((100, T < 500), (200, True))
        result = handle_numeric_dependency(processor, sample_material, 'test_prop', piecewise_expr, 400.0)
        assert result is True
        assert hasattr(sample_material, 'test_prop')
        assert 'test_prop' in processor.processed_properties

    def test_handle_numeric_dependency_evaluation_error(self, sample_material):
        """Expression with an unresolvable symbol raises ValueError."""
        processor = Mock()
        processor.processed_properties = set()
        T = sp.Symbol('T')
        # After substituting T=400, 'undef' remains symbolic -> float() raises TypeError
        bad_expr = T + sp.Symbol('undef')
        with pytest.raises(ValueError, match="Failed to evaluate.*at dependency="):
            handle_numeric_dependency(
                processor, sample_material, 'test_prop', bad_expr, 400.0)

    def test_handle_numeric_dependency_complex_expression(self, sample_material):
        """Complex symbolic expression evaluates to a finite float."""
        processor = Mock()
        processor.processed_properties = set()
        T = sp.Symbol('T')
        expr = sp.sin(T) * sp.exp(T / 1000) + sp.log(T)
        result = handle_numeric_dependency(processor, sample_material, 'complex_prop', expr, 400.0)
        assert result is True
        assert hasattr(sample_material, 'complex_prop')
        assert np.isfinite(float(getattr(sample_material, 'complex_prop')))

    def test_create_step_visualization_data_normal_range(self):
        """Normal range: 5-point output, extends beyond range bounds."""
        x_data, y_data = create_step_visualization_data(
        500.0, [100.0, 200.0], np.array([300, 400, 500, 600, 700]))
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert x_data[0] < 300
        assert x_data[-1] > 700
        assert y_data[0] == 100.0
        assert y_data[-1] == 200.0

    def test_create_step_visualization_data_narrow_range(self):
        """Very narrow range still produces correct 5-point step."""
        x_data, y_data = create_step_visualization_data(500.0, [100.0, 200.0], np.array([499, 501]))
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert y_data[0] == 100.0
        assert y_data[-1] == 200.0

    def test_create_step_visualization_data_edge_values(self):
        """Very small transition temp and extreme value scales produce finite output."""
        x_data, y_data = create_step_visualization_data(
            0.1, [1e-10, 1e10], np.array([0.05, 0.1, 0.15]))
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert all(np.isfinite(x_data))
        assert all(np.isfinite(y_data))
        assert y_data[0] == 1e-10
        assert y_data[-1] == 1e10

    def test_create_step_visualization_data_single_point_range(self):
        """Single-point range handled gracefully."""
        x_data, y_data = create_step_visualization_data(500.0, [100.0, 200.0], np.array([500]))
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert y_data[0] == 100.0
        assert y_data[-1] == 200.0

    def test_create_step_visualization_data_zero_values(self):
        """Both values zero -> all y output is zero."""
        x_data, y_data = create_step_visualization_data(500.0, [0.0, 0.0], np.array([300, 400, 500, 600, 700]))
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert all(y == 0.0 for y in y_data)

    def test_create_step_visualization_data_negative_values(self):
        """Negative property values produce finite, correctly ordered output."""
        x_data, y_data = create_step_visualization_data(500.0, [-100.0, -50.0], np.array([300, 400, 500, 600, 700]))
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert y_data[0] == -100.0
        assert y_data[-1] == -50.0
        assert all(np.isfinite(y_data))
