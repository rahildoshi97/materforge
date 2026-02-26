"""Unit tests for DependencyResolver."""

import pytest
import numpy as np
from materforge.core.materials import Material
from materforge.parsing.processors.dependency_resolver import DependencyResolver

def _make_material(**props) -> Material:
    mat = Material(name="TestMaterial")
    for k, v in props.items():
        setattr(mat, k, v)
    return mat


class TestDependencyResolver:
    """Tests for DependencyResolver."""

    def test_resolve_numeric_temperature(self):
        result = DependencyResolver.resolve_dependency_definition(500.0)
        np.testing.assert_array_equal(result, np.array([500.0]))

    def test_resolve_temperature_list(self):
        result = DependencyResolver.resolve_dependency_definition([300, 400, 500])
        np.testing.assert_array_equal(result, np.array([300.0, 400.0, 500.0]))

    def test_resolve_equidistant_format(self):
        result = DependencyResolver.resolve_dependency_definition("(300, 50)", n_values=5)
        np.testing.assert_array_equal(
            result, np.array([300.0, 350.0, 400.0, 450.0, 500.0]))

    def test_resolve_range_format(self):
        result = DependencyResolver.resolve_dependency_definition("(300, 500, 50)")
        expected = np.linspace(300, 500, 50)
        np.testing.assert_array_almost_equal(result, expected)
        assert result[0] == 300.0
        assert result[-1] == 500.0

    def test_resolve_dependency_reference(self):
        """Named scalar property on material resolves to its float value."""
        mat = _make_material(melting_temperature=933.47, boiling_temperature=2792.0)
        result = DependencyResolver.resolve_dependency_reference("melting_temperature", mat)
        assert result == pytest.approx(933.47)

    def test_resolve_dependency_arithmetic(self):
        """Arithmetic expressions on scalar properties resolve correctly."""
        mat = _make_material(melting_temperature=933.47, boiling_temperature=2792.0)
        assert DependencyResolver.resolve_dependency_reference("melting_temperature + 50", mat) == pytest.approx(983.47)
        assert DependencyResolver.resolve_dependency_reference("melting_temperature - 100", mat) == pytest.approx(833.47)

    def test_invalid_dependency_reference(self):
        """Unknown property name raises ValueError."""
        mat = _make_material(melting_temperature=933.47, boiling_temperature=2792.0)
        with pytest.raises(ValueError, match="Unknown dependency reference"):
            DependencyResolver.resolve_dependency_reference("invalid_temperature", mat)
