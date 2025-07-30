"""Unit tests for TemperatureResolver."""

import pytest
import numpy as np
import sympy as sp
from materforge.parsing.processors.temperature_resolver import TemperatureResolver

class TestTemperatureResolver:
    """Test cases for TemperatureResolver."""
    def test_resolve_numeric_temperature(self):
        """Test resolution of numeric temperature."""
        result = TemperatureResolver.resolve_temperature_definition(500.0)
        expected = np.array([500.0])
        np.testing.assert_array_equal(result, expected)

    def test_resolve_temperature_list(self):
        """Test resolution of temperature list."""
        temp_list = [300, 400, 500]
        result = TemperatureResolver.resolve_temperature_definition(temp_list)
        expected = np.array([300.0, 400.0, 500.0])
        np.testing.assert_array_equal(result, expected)

    def test_resolve_equidistant_format(self):
        """Test resolution of equidistant format."""
        result = TemperatureResolver.resolve_temperature_definition("(300, 50)", n_values=5)
        expected = np.array([300.0, 350.0, 400.0, 450.0, 500.0])
        np.testing.assert_array_equal(result, expected)

    def test_resolve_range_format(self):
        """Test resolution of range format."""
        result = TemperatureResolver.resolve_temperature_definition("(300, 500, 50)")
        expected = np.linspace(300, 500, 50)
        np.testing.assert_array_almost_equal(result, expected)
        assert result.shape == (50,)
        assert result[0] == 300.0
        assert result[-1] == 500.0

    def test_resolve_temperature_reference(self, sample_aluminum_element):
        """Test resolution of temperature reference."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        result = TemperatureResolver.resolve_temperature_reference("melting_temperature", material)
        assert result == 933.47

    def test_resolve_temperature_arithmetic(self, sample_aluminum_element):
        """Test resolution of temperature arithmetic expressions."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        result = TemperatureResolver.resolve_temperature_reference("melting_temperature + 50", material)
        assert result == 983.47
        result = TemperatureResolver.resolve_temperature_reference("melting_temperature - 100", material)
        assert result == 833.47

    def test_invalid_temperature_below_absolute_zero(self):
        """Test that temperatures below absolute zero raise error."""
        with pytest.raises(ValueError, match="Temperature must be above absolute zero"):
            TemperatureResolver.resolve_temperature_definition(-10.0)

    def test_invalid_temperature_reference(self, sample_aluminum_element):
        """Test invalid temperature reference."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        with pytest.raises(ValueError, match="Unknown temperature reference"):
            TemperatureResolver.resolve_temperature_reference("invalid_temperature", material)
