"""Comprehensive tests for parsing utilities."""

import pytest
import numpy as np
import sympy as sp
from unittest.mock import Mock
from materforge.parsing.utils.utilities import handle_numeric_temperature, create_step_visualization_data
from materforge.core.materials import Material
from materforge.core.elements import ChemicalElement


class TestUtilitiesComprehensive:
    """Test edge cases in utility functions."""

    @pytest.fixture
    def sample_material(self):
        """Create a sample material for testing."""
        element = ChemicalElement(
            name="Iron", atomic_number=26, atomic_mass=55.845,
            melting_temperature=1811, boiling_temperature=3134,
            latent_heat_of_fusion=13800, latent_heat_of_vaporization=340000
        )
        return Material(
            name="TestMaterial",
            material_type="pure_metal",
            elements=[element],
            composition=[1.0],
            melting_temperature=sp.Float(1811),
            boiling_temperature=sp.Float(3134)
        )

    def test_handle_numeric_temperature_symbolic(self, sample_material):
        """Test numeric temperature handler with symbolic temperature."""
        processor = Mock()
        processor.processed_properties = set()
        T = sp.Symbol('T')
        piecewise_expr = sp.Piecewise((100, T < 500), (200, True))
        result = handle_numeric_temperature(processor, sample_material, 'test_prop', piecewise_expr, T)
        assert result is False  # Should return False for symbolic temperature
        assert 'test_prop' not in processor.processed_properties

    def test_handle_numeric_temperature_numeric_success(self, sample_material):
        """Test numeric temperature handler with successful evaluation."""
        processor = Mock()
        processor.processed_properties = set()
        T_val = 400.0
        piecewise_expr = sp.Piecewise((100, sp.Symbol('T') < 500), (200, True))
        result = handle_numeric_temperature(processor, sample_material, 'test_prop', piecewise_expr, T_val)
        assert result is True  # Should return True for numeric temperature
        assert hasattr(sample_material, 'test_prop')
        assert 'test_prop' in processor.processed_properties

    def test_handle_numeric_temperature_evaluation_error(self, sample_material):
        """Test numeric temperature handler with evaluation errors."""
        processor = Mock()
        processor.processed_properties = set()
        T_val = 400.0
        # Create expression that will fail evaluation
        piecewise_expr = sp.Symbol('undefined_symbol') / 0
        with pytest.raises(ValueError, match="Failed to evaluate.*at T="):
            handle_numeric_temperature(processor, sample_material, 'test_prop', piecewise_expr, T_val)

    def test_handle_numeric_temperature_complex_expression(self, sample_material):
        """Test numeric temperature handler with complex expressions."""
        processor = Mock()
        processor.processed_properties = set()
        T_val = 400.0
        T = sp.Symbol('T')
        # Complex expression with multiple operations
        piecewise_expr = sp.sin(T) * sp.exp(T/1000) + sp.log(T)
        result = handle_numeric_temperature(processor, sample_material, 'complex_prop', piecewise_expr, T_val)
        assert result is True
        assert hasattr(sample_material, 'complex_prop')
        # Verify the value is reasonable
        value = float(getattr(sample_material, 'complex_prop'))
        assert np.isfinite(value)

    def test_create_step_visualization_data_normal_range(self):
        """Test step visualization data creation with normal temperature range."""
        transition_temp = 500.0
        val_array = [100.0, 200.0]
        temp_range = np.array([300, 400, 500, 600, 700])
        x_data, y_data = create_step_visualization_data(transition_temp, val_array, temp_range)
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert x_data[0] < temp_range[0]  # Should extend below range
        assert x_data[-1] > temp_range[-1]  # Should extend above range
        assert y_data[0] == val_array[0]  # Before transition
        assert y_data[-1] == val_array[1]  # After transition

    def test_create_step_visualization_data_narrow_range(self):
        """Test step visualization with very narrow temperature range."""
        transition_temp = 500.0
        val_array = [100.0, 200.0]
        temp_range = np.array([499, 501])  # Very narrow range
        x_data, y_data = create_step_visualization_data(transition_temp, val_array, temp_range)
        assert len(x_data) == 5
        assert len(y_data) == 5
        # Should still create proper step visualization
        assert y_data[0] == val_array[0]
        assert y_data[-1] == val_array[1]

    def test_create_step_visualization_data_edge_values(self):
        """Test step visualization with edge case values."""
        transition_temp = 0.1  # Very low temperature
        val_array = [1e-10, 1e10]  # Very different scales
        temp_range = np.array([0.05, 0.1, 0.15])
        x_data, y_data = create_step_visualization_data(transition_temp, val_array, temp_range)
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert all(np.isfinite(x_data))
        assert all(np.isfinite(y_data))
        assert y_data[0] == val_array[0]
        assert y_data[-1] == val_array[1]

    def test_create_step_visualization_data_single_point_range(self):
        """Test step visualization with single point temperature range."""
        transition_temp = 500.0
        val_array = [100.0, 200.0]
        temp_range = np.array([500])  # Single point
        x_data, y_data = create_step_visualization_data(transition_temp, val_array, temp_range)
        assert len(x_data) == 5
        assert len(y_data) == 5
        # Should handle single point gracefully
        assert y_data[0] == val_array[0]
        assert y_data[-1] == val_array[1]

    def test_create_step_visualization_data_zero_values(self):
        """Test step visualization with zero values."""
        transition_temp = 500.0
        val_array = [0.0, 0.0]  # Both values zero
        temp_range = np.array([300, 400, 500, 600, 700])
        x_data, y_data = create_step_visualization_data(transition_temp, val_array, temp_range)
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert all(y == 0.0 for y in y_data)

    def test_create_step_visualization_data_negative_values(self):
        """Test step visualization with negative values."""
        transition_temp = 500.0
        val_array = [-100.0, -50.0]  # Negative values
        temp_range = np.array([300, 400, 500, 600, 700])
        x_data, y_data = create_step_visualization_data(transition_temp, val_array, temp_range)
        assert len(x_data) == 5
        assert len(y_data) == 5
        assert y_data[0] == val_array[0]
        assert y_data[-1] == val_array[1]
        assert all(np.isfinite(y_data))
