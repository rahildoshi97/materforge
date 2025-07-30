"""Unit tests for interpolation algorithms."""

import pytest
import numpy as np
from materforge.algorithms.interpolation import interpolate_value, ensure_ascending_order

class TestInterpolation:
    """Test cases for interpolation functions."""
    def test_interpolate_value_within_range(self):
        """Test interpolation within data range."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([900, 950, 1000])
        # Test interpolation at midpoint
        result = interpolate_value(350, x_array, y_array, 'constant', 'constant')
        expected = 925.0  # Linear interpolation between 900 and 950
        assert np.isclose(result, expected)

    def test_interpolate_value_below_range_constant(self):
        """Test interpolation below range with constant boundary."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([900, 950, 1000])
        result = interpolate_value(250, x_array, y_array, 'constant', 'constant')
        assert result == 900.0  # Should return first value

    def test_interpolate_value_above_range_constant(self):
        """Test interpolation above range with constant boundary."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([900, 950, 1000])
        result = interpolate_value(600, x_array, y_array, 'constant', 'constant')
        assert result == 1000.0  # Should return last value

    def test_interpolate_value_below_range_extrapolate(self):
        """Test interpolation below range with extrapolation."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([900, 950, 1000])
        result = interpolate_value(200, x_array, y_array, 'extrapolate', 'constant')
        # Slope between first two points: (950-900)/(400-300) = 0.5
        # Extrapolated value: 900 + 0.5 * (200-300) = 900 - 50 = 850
        expected = 850.0
        assert np.isclose(result, expected)

    def test_ensure_ascending_order_already_ascending(self):
        """Test ensure_ascending_order with already ascending array."""
        temp_array = np.array([300, 400, 500])
        prop_array = np.array([900, 950, 1000])
        result_temp, result_prop = ensure_ascending_order(temp_array, prop_array)
        np.testing.assert_array_equal(result_temp, temp_array)
        np.testing.assert_array_equal(result_prop, prop_array)

    def test_ensure_ascending_order_descending(self):
        """Test ensure_ascending_order with descending array."""
        temp_array = np.array([500, 400, 300])
        prop_array = np.array([1000, 950, 900])
        result_temp, result_prop = ensure_ascending_order(temp_array, prop_array)
        expected_temp = np.array([300, 400, 500])
        expected_prop = np.array([900, 950, 1000])
        np.testing.assert_array_equal(result_temp, expected_temp)
        np.testing.assert_array_equal(result_prop, expected_prop)

    def test_ensure_ascending_order_non_monotonic(self):
        """Test ensure_ascending_order with non-monotonic array."""
        temp_array = np.array([300, 500, 400])
        prop_array = np.array([900, 1000, 950])
        with pytest.raises(ValueError, match="not strictly ascending or strictly descending"):
            ensure_ascending_order(temp_array, prop_array)

    def test_interpolate_value_exact_match(self):
        """Test interpolation at exact data points."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([900, 950, 1000])
        # Test exact matches
        assert interpolate_value(300, x_array, y_array, 'constant', 'constant') == 900.0
        assert interpolate_value(400, x_array, y_array, 'constant', 'constant') == 950.0
        assert interpolate_value(500, x_array, y_array, 'constant', 'constant') == 1000.0

    def test_interpolate_value_edge_cases(self):
        """Test interpolation edge cases."""
        x_array = np.array([300, 400])
        y_array = np.array([900, 950])
        # Test with minimum data points
        result = interpolate_value(350, x_array, y_array, 'constant', 'constant')
        expected = 925.0
        assert np.isclose(result, expected)
