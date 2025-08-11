"""Unit tests for TemperatureResolver."""

import pytest
import numpy as np
import sympy as sp
from materforge.parsing.processors.dependency_resolver import DependencyResolver

class TestDependencyResolver:
    """Test cases for DependencyResolver."""
    def test_resolve_numeric_temperature(self):
        """Test resolution of numeric temperature."""
        result = DependencyResolver.resolve_dependency_definition(500.0)
        expected = np.array([500.0])
        np.testing.assert_array_equal(result, expected)

    def test_resolve_temperature_list(self):
        """Test resolution of dependency list."""
        temp_list = [300, 400, 500]
        result = DependencyResolver.resolve_dependency_definition(temp_list)
        expected = np.array([300.0, 400.0, 500.0])
        np.testing.assert_array_equal(result, expected)

    def test_resolve_equidistant_format(self):
        """Test resolution of equidistant format."""
        result = DependencyResolver.resolve_dependency_definition("(300, 50)", n_values=5)
        expected = np.array([300.0, 350.0, 400.0, 450.0, 500.0])
        np.testing.assert_array_equal(result, expected)

    def test_resolve_range_format(self):
        """Test resolution of range format."""
        result = DependencyResolver.resolve_dependency_definition("(300, 500, 50)")
        expected = np.linspace(300, 500, 50)
        np.testing.assert_array_almost_equal(result, expected)
        assert result.shape == (50,)
        assert result[0] == 300.0
        assert result[-1] == 500.0

    def test_resolve_dependency_reference(self, sample_aluminum_element):
        """Test resolution of dependency reference."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        result = DependencyResolver.resolve_dependency_reference("melting_temperature", material)
        assert result == 933.47

    def test_resolve_dependency_arithmetic(self, sample_aluminum_element):
        """Test resolution of dependency arithmetic expressions."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        result = DependencyResolver.resolve_dependency_reference("melting_temperature + 50", material)
        assert result == 983.47
        result = DependencyResolver.resolve_dependency_reference("melting_temperature - 100", material)
        assert result == 833.47

    def test_invalid_dependency_below_absolute_zero(self):
        """Test that temperatures below absolute zero raise error."""
        with pytest.raises(ValueError, match="Dependency must be above absolute zero"):
            DependencyResolver.resolve_dependency_definition(-10.0)

    def test_invalid_dependency_reference(self, sample_aluminum_element):
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
        with pytest.raises(ValueError, match="Unknown dependency reference"):
            DependencyResolver.resolve_dependency_reference("invalid_temperature", material)
