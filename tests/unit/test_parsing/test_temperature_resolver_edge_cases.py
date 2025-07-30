"""Additional edge case tests for temperature resolver."""

import pytest
import numpy as np
import sympy as sp
from materforge.parsing.processors.temperature_resolver import TemperatureResolver
from materforge.core.materials import Material
from materforge.core.elements import ChemicalElement


class TestTemperatureResolverEdgeCases:
    """Test edge cases and error conditions."""

    @pytest.fixture
    def sample_material(self):
        """Create a sample material for testing."""
        element = ChemicalElement(
            name="Iron", atomic_number=26, atomic_mass=55.845,
            melting_temperature=1811, boiling_temperature=3134,
            latent_heat_of_fusion=13800, latent_heat_of_vaporization=340000
        )
        return Material(
            name="TestSteel",
            material_type="pure_metal",
            elements=[element],
            composition=[1.0],
            melting_temperature=sp.Float(1811),
            boiling_temperature=sp.Float(3134)
        )

    def test_resolve_negative_temperature(self):
        """Test error handling for negative temperatures."""
        with pytest.raises(ValueError, match="above absolute zero"):
            TemperatureResolver.resolve_temperature_definition(-100.0)

    def test_resolve_zero_kelvin(self):
        """Test error handling for absolute zero."""
        with pytest.raises(ValueError, match="above absolute zero"):
            TemperatureResolver.resolve_temperature_definition(0.0)

    def test_resolve_very_small_positive_temperature(self):
        """Test handling of very small positive temperatures."""
        result = TemperatureResolver.resolve_temperature_definition(0.1)
        assert result[0] == 0.1

    def test_resolve_equidistant_zero_increment(self):
        """Test error handling for zero increment in equidistant format."""
        with pytest.raises(ValueError, match="increment.*cannot be zero"):
            TemperatureResolver.resolve_temperature_definition("(300, 0)", n_values=5)

    def test_resolve_equidistant_insufficient_points(self):
        """Test error handling for insufficient points."""
        with pytest.raises(ValueError, match="Number of values must be at least.*got"):
            TemperatureResolver.resolve_temperature_definition("(300, 50)", n_values=1)

    def test_resolve_range_invalid_step(self):
        """Test error handling for invalid step in range format."""
        with pytest.raises(ValueError):
            TemperatureResolver.resolve_temperature_definition("(300, 500, 0)")

    def test_resolve_invalid_string_format(self):
        """Test error handling for invalid string formats."""
        with pytest.raises(ValueError):
            TemperatureResolver.resolve_temperature_definition("invalid_format")

    def test_resolve_temperature_reference_missing_material(self):
        """Test error handling when material is required but not provided."""
        with pytest.raises(ValueError, match="require material"):
            TemperatureResolver.resolve_temperature_definition("melting_temperature")

    def test_resolve_invalid_temperature_reference(self, sample_material):
        """Test error handling for invalid temperature references."""
        with pytest.raises(ValueError, match="Unknown temperature reference"):
            TemperatureResolver.resolve_temperature_reference("invalid_temp_ref", sample_material)

    def test_resolve_arithmetic_invalid_base_reference(self, sample_material):
        """Test error handling for invalid base reference in arithmetic."""
        with pytest.raises(ValueError, match="Unknown temperature reference"):
            TemperatureResolver.resolve_temperature_reference("invalid_ref + 50", sample_material)

    def test_resolve_list_with_invalid_reference(self, sample_material):
        """Test error handling for invalid references in lists."""
        with pytest.raises(ValueError):
            TemperatureResolver.resolve_temperature_definition(
                [300, "invalid_reference", 500], material=sample_material
            )

    def test_validate_temperature_array_empty(self):
        """Test validation of empty temperature arrays."""
        with pytest.raises(ValueError, match="empty"):
            TemperatureResolver.validate_temperature_array(np.array([]), "test")

    def test_validate_temperature_array_insufficient_points(self):
        """Test validation of arrays with insufficient points."""
        with pytest.raises(ValueError, match="at least.*points"):
            TemperatureResolver.validate_temperature_array(np.array([300]), "test")

    def test_validate_temperature_array_below_absolute_zero(self):
        """Test validation of arrays with sub-zero temperatures."""
        with pytest.raises(ValueError, match="above absolute zero"):
            TemperatureResolver.validate_temperature_array(np.array([300, -100, 500]), "test")

    def test_validate_temperature_array_non_finite(self):
        """Test validation of arrays with non-finite values."""
        with pytest.raises(ValueError, match="non-finite"):
            TemperatureResolver.validate_temperature_array(
                np.array([300, np.inf, 500]), "test"
            )
