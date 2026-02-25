"""Unit tests for Material class."""
import pytest
import sympy as sp

from materforge.core.materials import Material

class TestMaterial:
    """Tests for the dynamic-property Material model."""

    def test_construction_with_name(self):
        mat = Material(name="Aluminum")
        assert mat.name == "Aluminum"

    def test_construction_empty_properties(self):
        mat = Material(name="Empty")
        assert len(mat.property_names()) == 0  # type-agnostic

    def test_name_stored_correctly(self):
        mat = Material(name="Steel 1.4301")
        assert mat.name == "Steel 1.4301"

    def test_assign_scalar_float(self):
        mat = Material(name="Test")
        mat.density = 7850.0
        assert mat.density == 7850.0

    def test_assign_sympy_float(self):
        mat = Material(name="Test")
        mat.melting_temperature = sp.Float(933.47)
        assert float(mat.melting_temperature) == pytest.approx(933.47)

    def test_assign_sympy_piecewise(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.heat_capacity = sp.Piecewise((450 + 0.1*T, T < 1000), (550.0, True))
        assert isinstance(mat.heat_capacity, sp.Piecewise)

    def test_assign_sympy_expression(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.viscosity = 2.5*T + 100
        assert mat.viscosity.free_symbols == {T}

    def test_overwrite_property(self):
        mat = Material(name="Test")
        mat.density = 7000.0
        mat.density = 7850.0
        assert mat.density == 7850.0

    def test_multiple_properties_independent(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.density = 7850.0
        mat.heat_capacity = 450 + 0.1*T
        mat.melting_temperature = 1811.0
        assert mat.density == 7850.0
        assert mat.melting_temperature == 1811.0
        assert mat.heat_capacity.free_symbols == {T}

    def test_property_names_empty(self):
        mat = Material(name="Test")
        assert len(mat.property_names()) == 0  # type-agnostic

    def test_property_names_single(self):
        mat = Material(name="Test")
        mat.density = 7850.0
        assert 'density' in mat.property_names()

    def test_property_names_multiple(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.density = 7850.0
        mat.heat_capacity = 450 + 0.1*T
        mat.melting_temperature = 933.47
        names = mat.property_names()
        for expected in ('density', 'heat_capacity', 'melting_temperature'):
            assert expected in names

    def test_property_names_does_not_include_name(self):
        """'name' is a constructor field, not a dynamic property."""
        mat = Material(name="Test")
        assert 'name' not in mat.property_names()

    def test_property_names_returns_copy(self):
        """Mutating the returned collection must not affect the material."""
        mat = Material(name="Test")
        mat.density = 7850.0
        names = mat.property_names()
        if isinstance(names, set):
            names.add('fake_property')
        else:
            names.append('fake_property')
        assert 'fake_property' not in mat.property_names()

    def test_evaluate_symbolic_property(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.heat_capacity = 450 + 0.1*T
        result = mat.evaluate_properties_at_temperature(T, 500.0)
        assert result['heat_capacity'] == pytest.approx(500.0)

    def test_evaluate_constant_property(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.density = 7850.0
        result = mat.evaluate_properties_at_temperature(T, 500.0)
        assert result['density'] == pytest.approx(7850.0)

    def test_evaluate_piecewise_below_boundary(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.heat_capacity = sp.Piecewise((400.0, T < 1000), (600.0, True))
        result = mat.evaluate_properties_at_temperature(T, 500.0)
        assert result['heat_capacity'] == pytest.approx(400.0)

    def test_evaluate_piecewise_above_boundary(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.heat_capacity = sp.Piecewise((400.0, T < 1000), (600.0, True))
        result = mat.evaluate_properties_at_temperature(T, 1500.0)
        assert result['heat_capacity'] == pytest.approx(600.0)

    def test_evaluate_no_properties_returns_empty(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        result = mat.evaluate_properties_at_temperature(T, 500.0)
        assert result == {}

    def test_evaluate_returns_float_values(self):
        T = sp.Symbol('T')
        mat = Material(name="Test")
        mat.heat_capacity = 450 + 0.1*T
        result = mat.evaluate_properties_at_temperature(T, 500.0)
        assert isinstance(result['heat_capacity'], float)

    def test_sample_valid_material_fixture(self, sample_valid_material):
        assert sample_valid_material.name == "Test Aluminum"
        assert sample_valid_material.melting_temperature == pytest.approx(933.47)
        assert sample_valid_material.boiling_temperature == pytest.approx(2792.0)

    def test_sample_valid_alloy_fixture(self, sample_valid_alloy):
        assert sample_valid_alloy.name == "Test Steel"
        assert sample_valid_alloy.solidus_temperature == pytest.approx(1400.0)
        assert sample_valid_alloy.liquidus_temperature == pytest.approx(1450.0)
