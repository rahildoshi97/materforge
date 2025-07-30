"""Comprehensive edge case tests for interpolation module."""

import pytest
import numpy as np
from materforge.algorithms.interpolation import interpolate_value, ensure_ascending_order


class TestInterpolationComprehensive:
    """Test comprehensive edge cases in interpolation functionality."""

    def test_interpolate_value_equal_boundary_points(self):
        """Test interpolation when consecutive points are equal."""
        x_array = np.array([300, 300, 400])  # Duplicate first point
        y_array = np.array([100, 100, 150])
        with pytest.raises(ValueError, match="Cannot extrapolate.*equal"):
            interpolate_value(250, x_array, y_array, 'extrapolate', 'constant')

    def test_interpolate_value_equal_last_points(self):
        """Test interpolation when last two points are equal."""
        x_array = np.array([300, 400, 400])  # Duplicate last point
        y_array = np.array([100, 150, 150])
        with pytest.raises(ValueError, match="Cannot extrapolate.*equal"):
            interpolate_value(450, x_array, y_array, 'constant', 'extrapolate')

    def test_interpolate_value_at_exact_boundary(self):
        """Test interpolation exactly at boundary points."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        # Exactly at lower boundary
        result = interpolate_value(300, x_array, y_array, 'constant', 'constant')
        assert result == 100.0
        # Exactly at upper boundary
        result = interpolate_value(500, x_array, y_array, 'constant', 'constant')
        assert result == 200.0

    def test_interpolate_value_very_close_to_boundary(self):
        """Test interpolation very close to boundaries."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        # Very close to lower boundary (should extrapolate)
        result = interpolate_value(299.999, x_array, y_array, 'extrapolate', 'constant')
        expected = 100 + (150 - 100) * (299.999 - 300) / (400 - 300)
        assert abs(result - expected) < 1e-10
        # Very close to upper boundary (should extrapolate)
        result = interpolate_value(500.001, x_array, y_array, 'constant', 'extrapolate')
        expected = 200 + (200 - 150) * (500.001 - 500) / (500 - 400)
        assert abs(result - expected) < 1e-10

    def test_interpolate_value_extreme_values(self):
        """Test interpolation with extreme numerical values."""
        x_array = np.array([1e-10, 1e-5, 1e10])
        y_array = np.array([1e-20, 1e20, 1e-20])
        # Test interpolation in the middle
        result = interpolate_value(1e-6, x_array, y_array, 'constant', 'constant')
        assert np.isfinite(result)

    def test_interpolate_value_negative_values(self):
        """Test interpolation with negative values."""
        x_array = np.array([-500, -300, -100])
        y_array = np.array([-200, -150, -100])
        # Test within range
        result = interpolate_value(-400, x_array, y_array, 'constant', 'constant')
        expected = np.interp(-400, x_array, y_array)
        assert abs(result - expected) < 1e-10
        # Test extrapolation
        result = interpolate_value(-600, x_array, y_array, 'extrapolate', 'constant')
        slope = (-150 - (-200)) / (-300 - (-500))
        expected = -200 + slope * (-600 - (-500))
        assert abs(result - expected) < 1e-10

    def test_interpolate_value_within_range_standard(self):
        """Test standard interpolation within range."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        # Test interpolation at midpoint
        result = interpolate_value(350, x_array, y_array, 'constant', 'constant')
        expected = 100 + (150 - 100) * (350 - 300) / (400 - 300)
        assert abs(result - expected) < 1e-10

    def test_interpolate_value_below_range_constant(self):
        """Test interpolation below range with constant boundary."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        result = interpolate_value(250, x_array, y_array, 'constant', 'constant')
        assert result == 100.0  # Should return first value

    def test_interpolate_value_above_range_constant(self):
        """Test interpolation above range with constant boundary."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        result = interpolate_value(550, x_array, y_array, 'constant', 'constant')
        assert result == 200.0  # Should return last value

    def test_interpolate_value_below_range_extrapolate(self):
        """Test interpolation below range with extrapolation."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        result = interpolate_value(250, x_array, y_array, 'extrapolate', 'constant')
        # Should extrapolate using first two points
        slope = (150 - 100) / (400 - 300)
        expected = 100 + slope * (250 - 300)
        assert abs(result - expected) < 1e-10

    def test_interpolate_value_above_range_extrapolate(self):
        """Test interpolation above range with extrapolation."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        result = interpolate_value(550, x_array, y_array, 'constant', 'extrapolate')
        # Should extrapolate using last two points
        slope = (200 - 150) / (500 - 400)
        expected = 200 + slope * (550 - 500)
        assert abs(result - expected) < 1e-10

    def test_ensure_ascending_order_already_ascending(self):
        """Test ensure_ascending_order with already ascending arrays."""
        temp_array = np.array([100, 200, 300])
        value_array1 = np.array([10, 20, 30])
        value_array2 = np.array([1, 2, 3])
        result = ensure_ascending_order(temp_array, value_array1, value_array2)
        # Should return unchanged
        np.testing.assert_array_equal(result[0], temp_array)
        np.testing.assert_array_equal(result[1], value_array1)
        np.testing.assert_array_equal(result[2], value_array2)

    def test_ensure_ascending_order_descending(self):
        """Test ensure_ascending_order with descending arrays."""
        temp_array = np.array([300, 200, 100])
        value_array1 = np.array([30, 20, 10])
        value_array2 = np.array([3, 2, 1])
        result = ensure_ascending_order(temp_array, value_array1, value_array2)
        # Should be flipped
        np.testing.assert_array_equal(result[0], np.array([100, 200, 300]))
        np.testing.assert_array_equal(result[1], np.array([10, 20, 30]))
        np.testing.assert_array_equal(result[2], np.array([1, 2, 3]))

    def test_ensure_ascending_order_mixed_order(self):
        """Test ensure_ascending_order with mixed order (should raise error)."""
        temp_array = np.array([100, 300, 200])  # Neither ascending nor descending
        value_array = np.array([10, 30, 20])
        with pytest.raises(ValueError, match="not strictly ascending or strictly descending"):
            ensure_ascending_order(temp_array, value_array)

    def test_ensure_ascending_order_single_element(self):
        """Test ensure_ascending_order with single element."""
        temp_array = np.array([100])
        value_array = np.array([10])
        result = ensure_ascending_order(temp_array, value_array)
        # Should return unchanged
        np.testing.assert_array_equal(result[0], temp_array)
        np.testing.assert_array_equal(result[1], value_array)

    def test_ensure_ascending_order_empty_array(self):
        """Test ensure_ascending_order with empty arrays."""
        temp_array = np.array([])
        value_array = np.array([])
        result = ensure_ascending_order(temp_array, value_array)
        # Should return unchanged
        np.testing.assert_array_equal(result[0], temp_array)
        np.testing.assert_array_equal(result[1], value_array)

    def test_ensure_ascending_order_constant_array(self):
        """Test ensure_ascending_order with constant values."""
        temp_array = np.array([100, 100, 100])  # All same values
        value_array = np.array([10, 10, 10])
        with pytest.raises(ValueError, match="not strictly ascending or strictly descending"):
            ensure_ascending_order(temp_array, value_array)

    def test_ensure_ascending_order_floating_point_precision(self):
        """Test ensure_ascending_order with floating point precision issues."""
        # Values that are very close but not exactly equal
        temp_array = np.array([100.0, 100.0 + 1e-12, 100.0 + 2e-12])  # Use larger differences
        value_array = np.array([10, 20, 30])
        # Should be treated as ascending (differences > tolerance)
        result = ensure_ascending_order(temp_array, value_array)
        np.testing.assert_array_equal(result[0], temp_array)
        np.testing.assert_array_equal(result[1], value_array)

    def test_interpolate_value_two_point_array(self):
        """Test interpolation with minimum two-point array."""
        x_array = np.array([300, 500])
        y_array = np.array([100, 200])
        # Test interpolation at midpoint
        result = interpolate_value(400, x_array, y_array, 'constant', 'constant')
        expected = 100 + (200 - 100) * (400 - 300) / (500 - 300)
        assert abs(result - expected) < 1e-10
        # Test extrapolation below
        result = interpolate_value(200, x_array, y_array, 'extrapolate', 'constant')
        slope = (200 - 100) / (500 - 300)
        expected = 100 + slope * (200 - 300)
        assert abs(result - expected) < 1e-10

    def test_interpolate_value_nan_input(self):
        """Test interpolation with NaN input value."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        # NaN input should raise error
        with pytest.raises(ValueError, match="Temperature T must be finite, got nan"):
            interpolate_value(np.nan, x_array, y_array, 'constant', 'constant')

    def test_interpolate_value_inf_input(self):
        """Test interpolation with infinite input value."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        # Infinite input should raise validation error
        with pytest.raises(ValueError, match="Temperature T must be finite"):
            interpolate_value(np.inf, x_array, y_array, 'constant', 'extrapolate')
        # Also test negative infinity
        with pytest.raises(ValueError, match="Temperature T must be finite"):
            interpolate_value(-np.inf, x_array, y_array, 'constant', 'extrapolate')

    def test_interpolate_value_mismatched_array_lengths(self):
        """Test interpolation with mismatched array lengths."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150])  # Different length
        with pytest.raises(ValueError):
            interpolate_value(350, x_array, y_array, 'constant', 'constant')

    def test_interpolate_value_empty_arrays(self):
        """Test interpolation with empty arrays."""
        x_array = np.array([])
        y_array = np.array([])
        with pytest.raises(ValueError):
            interpolate_value(350, x_array, y_array, 'constant', 'constant')

    def test_interpolate_value_single_point_array(self):
        """Test interpolation with single point array."""
        x_array = np.array([300])
        y_array = np.array([100])
        # Single point should return that value regardless of input
        result = interpolate_value(350, x_array, y_array, 'constant', 'constant')
        assert result == 100.0
        result = interpolate_value(250, x_array, y_array, 'extrapolate', 'extrapolate')
        assert result == 100.0

    def test_ensure_ascending_order_two_elements(self):
        """Test ensure_ascending_order with two elements."""
        # Ascending case
        temp_array = np.array([100, 200])
        value_array = np.array([10, 20])
        result = ensure_ascending_order(temp_array, value_array)
        np.testing.assert_array_equal(result[0], temp_array)
        np.testing.assert_array_equal(result[1], value_array)
        # Descending case
        temp_array = np.array([200, 100])
        value_array = np.array([20, 10])
        result = ensure_ascending_order(temp_array, value_array)
        np.testing.assert_array_equal(result[0], np.array([100, 200]))
        np.testing.assert_array_equal(result[1], np.array([10, 20]))

    def test_interpolate_value_boundary_conditions_comprehensive(self):
        """Comprehensive test of all boundary condition combinations."""
        x_array = np.array([300, 400, 500])
        y_array = np.array([100, 150, 200])
        # Test all boundary combinations
        boundary_combinations = [
            ('constant', 'constant'),
            ('constant', 'extrapolate'),
            ('extrapolate', 'constant'),
            ('extrapolate', 'extrapolate')
        ]
        for lower_bound, upper_bound in boundary_combinations:
            # Test below range
            result_below = interpolate_value(250, x_array, y_array, lower_bound, upper_bound)
            assert np.isfinite(result_below), f"Non-finite result for bounds ({lower_bound}, {upper_bound})"
            # Test above range
            result_above = interpolate_value(550, x_array, y_array, lower_bound, upper_bound)
            assert np.isfinite(result_above), f"Non-finite result for bounds ({lower_bound}, {upper_bound})"
            # Test within range (should be same regardless of bounds)
            result_within = interpolate_value(350, x_array, y_array, lower_bound, upper_bound)
            expected_within = 100 + (150 - 100) * (350 - 300) / (400 - 300)
            assert abs(result_within - expected_within) < 1e-10
