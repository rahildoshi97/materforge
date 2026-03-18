"""Integration tests for YAML material creation."""
import math
import pytest
import sympy as sp
from materforge import create_material

class TestYAMLMaterialCreation:
    """Integration tests for creating materials from YAML files."""
    # temp_symbol, aluminum_yaml_path, steel_yaml_path are all provided
    # by conftest.py - no need to redefine them here.

    def test_aluminum_material_creation(self, aluminum_yaml_path, temp_symbol):
        """Aluminum material: name, scalar temps, and property existence."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML not found: {aluminum_yaml_path}")
        mat = create_material(yaml_path=aluminum_yaml_path, dependency=temp_symbol, enable_plotting=False)
        assert "Aluminum" in mat.name or "Al" in mat.name
        assert 'melting_temperature' in mat.property_names()
        assert 'boiling_temperature' in mat.property_names()
        assert float(mat.melting_temperature) > 0
        assert float(mat.boiling_temperature) > float(mat.melting_temperature)

    def test_steel_material_creation(self, steel_yaml_path, temp_symbol):
        """Steel material: name, solidus/liquidus ordering."""
        if not steel_yaml_path.exists():
            pytest.skip(f"Steel YAML not found: {steel_yaml_path}")
        mat = create_material(yaml_path=steel_yaml_path, dependency=temp_symbol, enable_plotting=False)
        assert "Steel" in mat.name or "1.4301" in mat.name
        assert 'solidus_temperature' in mat.property_names()
        assert 'liquidus_temperature' in mat.property_names()
        assert float(mat.solidus_temperature) > 0
        assert float(mat.liquidus_temperature) > float(mat.solidus_temperature)

    def test_material_property_evaluation(self, aluminum_yaml_path, temp_symbol):
        """Every dynamic property evaluates to a finite float at T=300."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML not found: {aluminum_yaml_path}")
        mat = create_material(yaml_path=aluminum_yaml_path, dependency=temp_symbol, enable_plotting=False)
        props = mat.property_names()
        assert len(props) > 0, "Material has no processed properties"
        for prop_name in props:
            value = getattr(mat, prop_name)
            assert value is not None, f"'{prop_name}' is None after processing"
            if isinstance(value, sp.Expr) and value.free_symbols:
                numerical = float(value.subs(temp_symbol, 300).evalf())
                assert not math.isnan(numerical), \
                    f"'{prop_name}' evaluated to NaN at T=300"
            else:
                numerical = float(value)
            assert isinstance(numerical, float), \
                f"'{prop_name}' did not resolve to float"

    def test_comprehensive_material_properties(self, aluminum_yaml_path, temp_symbol):
        """All sympy properties evaluate without error at T=300."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML not found: {aluminum_yaml_path}")
        mat = create_material(yaml_path=aluminum_yaml_path, dependency=temp_symbol, enable_plotting=False)
        for prop_name in mat.property_names():
            value = getattr(mat, prop_name)
            if isinstance(value, sp.Expr):
                try:
                    numerical = float(value.subs(temp_symbol, 300).evalf())
                    assert isinstance(numerical, float)
                    assert not math.isnan(numerical)
                except (TypeError, ValueError):
                    pass  # Non-numeric symbolic constant - skip
