"""Unit tests for array validation functions."""

import pytest
import numpy as np
from materforge.parsing.validation.array_validator import is_monotonic

class TestArrayValidator:
    """Test cases for array validation functions."""
    def test_is_monotonic_strictly_increasing(self):
        """Test strictly increasing array validation."""
        array = np.array([1, 2, 3, 4, 5])
        assert is_monotonic(array, mode="strictly_increasing") is True
        array_float = np.array([1.1, 2.2, 3.3, 4.4])
        assert is_monotonic(array_float, mode="strictly_increasing") is True

    def test_is_monotonic_not_strictly_increasing(self):
        """Test arrays that are not strictly increasing."""
        array_equal = np.array([1, 2, 2, 3])
        assert is_monotonic(array_equal, mode="strictly_increasing", raise_error=False) is False
        array_decreasing = np.array([1, 3, 2, 4])
        assert is_monotonic(array_decreasing, mode="strictly_increasing", raise_error=False) is False

    def test_is_monotonic_with_error_raising(self):
        """Test that errors are raised when expected."""
        array = np.array([1, 2, 2, 3])
        with pytest.raises(ValueError, match="not strictly increasing"):
            is_monotonic(array, mode="strictly_increasing", raise_error=True)

    def test_is_monotonic_single_element(self):
        """Test monotonicity with single element array."""
        array = np.array([42])
        assert is_monotonic(array, mode="strictly_increasing") is True

    def test_is_monotonic_empty_array(self):
        """Test monotonicity with empty array."""
        array = np.array([])
        assert is_monotonic(array, mode="strictly_increasing") is True

    def test_is_monotonic_temperature_array(self):
        """Test monotonicity with realistic temperature arrays."""
        temp_array = np.array([273.15, 373.15, 473.15, 573.15, 673.15])
        assert is_monotonic(temp_array, mode="strictly_increasing") is True
        temp_array_equal = np.array([273.15, 373.15, 373.15, 473.15])
        assert is_monotonic(temp_array_equal, mode="strictly_increasing", raise_error=False) is False
