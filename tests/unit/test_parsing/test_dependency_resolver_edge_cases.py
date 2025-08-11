"""Additional edge case tests for temperature resolver."""

import pytest
import numpy as np
import sympy as sp
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.core.materials import Material
from materforge.core.elements import ChemicalElement


class TestDependencyResolverEdgeCases:
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

    def test_resolve_negative_dependency(self):
        """Test error handling for negative dependencies."""
        with pytest.raises(ValueError, match="above absolute zero"):
            DependencyResolver.resolve_dependency_definition(-100.0)

    def test_resolve_zero_kelvin(self):
        """Test error handling for absolute zero."""
        with pytest.raises(ValueError, match="above absolute zero"):
            DependencyResolver.resolve_dependency_definition(0.0)

    def test_resolve_very_small_positive_dependency(self):
        """Test handling of very small positive dependencies."""
        result = DependencyResolver.resolve_dependency_definition(0.1)
        assert result[0] == 0.1

    def test_resolve_equidistant_zero_increment(self):
        """Test error handling for zero increment in equidistant format."""
        with pytest.raises(ValueError, match="increment.*cannot be zero"):
            DependencyResolver.resolve_dependency_definition("(300, 0)", n_values=5)

    def test_resolve_equidistant_insufficient_points(self):
        """Test error handling for insufficient points."""
        with pytest.raises(ValueError, match="Number of values must be at least.*got"):
            DependencyResolver.resolve_dependency_definition("(300, 50)", n_values=1)

    def test_resolve_range_invalid_step(self):
        """Test error handling for invalid step in range format."""
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition("(300, 500, 0)")

    def test_resolve_invalid_string_format(self):
        """Test error handling for invalid string formats."""
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition("invalid_format")

    def test_resolve_dependency_reference_missing_material(self):
        """Test error handling when material is required but not provided."""
        with pytest.raises(ValueError, match="require material"):
            DependencyResolver.resolve_dependency_definition("melting_temperature")

    def test_resolve_invalid_dependency_reference(self, sample_material):
        """Test error handling for invalid dependency references."""
        with pytest.raises(ValueError, match="Unknown dependency reference"):
            DependencyResolver.resolve_dependency_reference("invalid_temp_ref", sample_material)

    def test_resolve_arithmetic_invalid_base_reference(self, sample_material):
        """Test error handling for invalid base reference in arithmetic."""
        with pytest.raises(ValueError, match="Unknown dependency reference"):
            DependencyResolver.resolve_dependency_reference("invalid_ref + 50", sample_material)

    def test_resolve_list_with_invalid_reference(self, sample_material):
        """Test error handling for invalid references in lists."""
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition(
                [300, "invalid_reference", 500], material=sample_material
            )

    def test_validate_dependency_array_empty(self):
        """Test validation of empty dependency arrays."""
        with pytest.raises(ValueError, match="empty"):
            DependencyResolver.validate_dependency_array(np.array([]), "test")

    def test_validate_dependency_array_insufficient_points(self):
        """Test validation of arrays with insufficient points."""
        with pytest.raises(ValueError, match="at least.*points"):
            DependencyResolver.validate_dependency_array(np.array([300]), "test")

    def test_validate_dependency_array_below_absolute_zero(self):
        """Test validation of arrays with sub-zero temperatures."""
        with pytest.raises(ValueError, match="above absolute zero"):
            DependencyResolver.validate_dependency_array(np.array([300, -100, 500]), "test")

    def test_validate_temperature_array_non_finite(self):
        """Test validation of arrays with non-finite values."""
        with pytest.raises(ValueError, match="non-finite"):
            DependencyResolver.validate_dependency_array(
                np.array([300, np.inf, 500]), "test"
            )
