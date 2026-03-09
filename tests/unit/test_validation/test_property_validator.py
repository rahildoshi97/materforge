"""Unit tests for property validation functions."""

import pytest
import numpy as np
from materforge.parsing.validation.property_validator import validate_monotonic_property


class TestValidateMonotonicProperty:
    """Tests for the pure validate_monotonic_property validator."""

    # --- Strictly increasing (default mode) ---

    def test_strictly_increasing_valid(self):
        dep = np.array([300, 400, 500, 600])
        prop = np.array([1000, 1500, 2000, 2500])
        validate_monotonic_property("energy_density", dep, prop)

    def test_strictly_increasing_invalid(self):
        dep = np.array([300, 400, 500, 600])
        prop = np.array([1000, 1500, 1200, 2500])
        with pytest.raises(ValueError, match="violates strictly increasing constraint"):
            validate_monotonic_property("energy_density", dep, prop)

    def test_constant_values_fail_strictly_increasing(self):
        dep = np.array([300, 400, 500])
        prop = np.array([100, 100, 100])
        with pytest.raises(ValueError, match="violates strictly increasing constraint"):
            validate_monotonic_property("some_prop", dep, prop)

    def test_constant_values_pass_non_decreasing(self):
        dep = np.array([300, 400, 500])
        prop = np.array([100, 100, 100])
        validate_monotonic_property("some_prop", dep, prop, mode="non_decreasing")

    # --- Strictly decreasing ---

    def test_strictly_decreasing_valid(self):
        dep = np.array([300, 400, 500])
        prop = np.array([200, 150, 100])
        validate_monotonic_property("density", dep, prop, mode="strictly_decreasing")

    def test_strictly_decreasing_invalid(self):
        dep = np.array([300, 400, 500])
        prop = np.array([200, 150, 160])
        with pytest.raises(ValueError, match="violates strictly decreasing constraint"):
            validate_monotonic_property("density", dep, prop, mode="strictly_decreasing")

    # --- Non-increasing ---

    def test_non_increasing_valid(self):
        dep = np.array([300, 400, 500])
        prop = np.array([200, 200, 100])
        validate_monotonic_property("some_prop", dep, prop, mode="non_increasing")

    # --- Tolerance ---

    def test_tolerance_respected(self):
        dep = np.array([300, 400, 500])
        # Tiny dip well within tolerance - should pass
        prop = np.array([100.0, 200.0, 200.0 - 1e-13])
        validate_monotonic_property("some_prop", dep, prop,
                                    mode="non_decreasing", tolerance=1e-10)

    # --- Edge cases ---

    def test_single_element(self):
        dep = np.array([300.0])
        prop = np.array([100.0])
        validate_monotonic_property("single_point", dep, prop)

    def test_two_elements_increasing(self):
        dep = np.array([300, 400])
        prop = np.array([100, 150])
        validate_monotonic_property("two_points", dep, prop)

    def test_two_elements_decreasing_fails_default(self):
        dep = np.array([300, 400])
        prop = np.array([150, 100])
        with pytest.raises(ValueError, match="violates strictly increasing constraint"):
            validate_monotonic_property("two_points", dep, prop)

    def test_error_message_contains_ranges(self):
        dep = np.array([300.0, 400.0, 500.0])
        prop = np.array([100.0, 80.0, 120.0])
        with pytest.raises(ValueError) as exc_info:
            validate_monotonic_property("heat_capacity", dep, prop)
        msg = str(exc_info.value)
        assert "Dependency range" in msg
        assert "Property range" in msg
        assert "Validation details" in msg

    # --- Realistic material data ---

    def test_energy_density_realistic(self):
        dep = np.array([273.15, 373.15, 473.15, 573.15, 673.15])
        energy_density = np.array([0.0, 90.0, 190.0, 300.0, 420.0])
        validate_monotonic_property("energy_density", dep, energy_density)

    def test_heat_capacity_non_monotone_passes_by_default(self):
        """heat_capacity is not monotone - the validator is pure and raises.
        The calling code is responsible for deciding whether to call it at all.
        """
        dep = np.array([273.15, 373.15, 473.15, 573.15, 673.15])
        heat_capacity = np.array([900, 950, 1000, 920, 1100])  # dip at index 3
        with pytest.raises(ValueError, match="violates strictly increasing constraint"):
            validate_monotonic_property("heat_capacity", dep, heat_capacity)

    def test_density_decreasing_fails_default_mode(self):
        dep = np.array([273.15, 373.15, 473.15, 573.15, 673.15])
        density = np.array([2700, 2690, 2680, 2670, 2660])
        with pytest.raises(ValueError, match="violates strictly increasing constraint"):
            validate_monotonic_property("density", dep, density)

    def test_density_decreasing_passes_correct_mode(self):
        dep = np.array([273.15, 373.15, 473.15, 573.15, 673.15])
        density = np.array([2700, 2690, 2680, 2670, 2660])
        validate_monotonic_property("density", dep, density, mode="strictly_decreasing")
