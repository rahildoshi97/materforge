"""Additional edge case tests for dependency resolver."""
import pytest
import numpy as np
from materforge.core.materials import Material
from materforge.parsing.processors.dependency_resolver import DependencyResolver


class TestDependencyResolverEdgeCases:
    """Edge cases and error conditions for DependencyResolver."""

    @pytest.fixture
    def sample_material(self):
        mat = Material(name="TestSteel")
        mat.melting_temperature = 1811.0
        mat.boiling_temperature = 3134.0
        return mat

    def test_resolve_negative_dependency(self):
        with pytest.raises(ValueError, match="above absolute zero"):
            DependencyResolver.resolve_dependency_definition(-100.0)

    def test_resolve_zero_kelvin(self):
        with pytest.raises(ValueError, match="above absolute zero"):
            DependencyResolver.resolve_dependency_definition(0.0)

    def test_resolve_very_small_positive_dependency(self):
        result = DependencyResolver.resolve_dependency_definition(0.1)
        assert result[0] == 0.1

    def test_resolve_equidistant_zero_increment(self):
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition("(300, 0)", n_values=5)

    def test_resolve_equidistant_insufficient_points(self):
        with pytest.raises(ValueError, match="Number of values must be at least.*got"):
            DependencyResolver.resolve_dependency_definition("(300, 50)", n_values=1)

    def test_resolve_range_invalid_step(self):
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition("(300, 500, 0)")

    def test_resolve_invalid_string_format(self):
        """Bare non-reference string without material raises ValueError."""
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition("invalid_format")

    def test_resolve_dependency_reference_missing_material(self):
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition("melting_temperature")

    def test_resolve_valid_scalar_property(self, sample_material):
        """Scalar property defined on material resolves correctly."""
        result = DependencyResolver.resolve_dependency_reference(
            "melting_temperature", sample_material)
        assert result == 1811.0

    def test_resolve_arithmetic_on_scalar_property(self, sample_material):
        """Arithmetic expression on a valid scalar property resolves correctly."""
        result = DependencyResolver.resolve_dependency_reference(
            "melting_temperature + 50", sample_material)
        assert result == 1861.0

    def test_resolve_invalid_dependency_reference(self, sample_material):
        """Unknown property name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown dependency reference"):
            DependencyResolver.resolve_dependency_reference("invalid_temp_ref", sample_material)

    def test_resolve_arithmetic_invalid_base_reference(self, sample_material):
        """Arithmetic expression with unknown base name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown dependency reference"):
            DependencyResolver.resolve_dependency_reference("invalid_ref + 50", sample_material)

    def test_resolve_list_with_invalid_reference(self, sample_material):
        """Invalid string reference inside a list raises ValueError."""
        with pytest.raises(ValueError):
            DependencyResolver.resolve_dependency_definition(
                [300, "invalid_reference", 500], material=sample_material)

    def test_resolve_list_with_valid_reference(self, sample_material):
        """Valid scalar reference inside a list resolves correctly."""
        result = DependencyResolver.resolve_dependency_definition(
            [300, "melting_temperature", 3000], material=sample_material)
        assert result[1] == 1811.0

    def test_validate_dependency_array_empty(self):
        with pytest.raises(ValueError, match="empty"):
            DependencyResolver.validate_dependency_array(np.array([]), "test")

    def test_validate_dependency_array_insufficient_points(self):
        with pytest.raises(ValueError, match="at least.*points"):
            DependencyResolver.validate_dependency_array(np.array([300]), "test")

    def test_validate_dependency_array_below_absolute_zero(self):
        with pytest.raises(ValueError, match="above absolute zero"):
            DependencyResolver.validate_dependency_array(np.array([300, -100, 500]), "test")

    def test_validate_dependency_array_non_finite(self):
        with pytest.raises(ValueError, match="non-finite"):
            DependencyResolver.validate_dependency_array(
                np.array([300, np.inf, 500]), "test")
