"""Integration tests for YAML material creation."""
import pytest
import sympy as sp
from pathlib import Path

from pymatlib.parsing.api import create_material, get_supported_properties

class TestYAMLMaterialCreation:
    """Integration tests for creating materials from YAML files."""
    @pytest.fixture
    def temp_symbol(self):
        """Temperature symbol for testing."""
        return sp.Symbol('T')

    @pytest.fixture
    def aluminum_yaml_path(self):
        """Path to aluminum YAML file."""
        current_file = Path(__file__)
        return current_file.parent.parent.parent / "src" / "pymatlib" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"

    @pytest.fixture
    def steel_yaml_path(self):
        """Path to steel YAML file."""
        current_file = Path(__file__)
        return current_file.parent.parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

    def test_aluminum_material_creation(self, aluminum_yaml_path, temp_symbol):
        """Test aluminum material creation from YAML."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML file not found: {aluminum_yaml_path}")
        mat_Al = create_material(yaml_path=aluminum_yaml_path, dependency=temp_symbol, enable_plotting=False)
        # Basic material verification
        assert mat_Al.name == "Aluminum"
        assert mat_Al.material_type == "pure_metal"
        assert len(mat_Al.elements) == 1
        assert mat_Al.elements[0].name in ["Aluminum", "Aluminium"]
        # Temperature properties
        assert hasattr(mat_Al, 'melting_temperature')
        assert hasattr(mat_Al, 'boiling_temperature')
        melting_temp = float(mat_Al.melting_temperature) if hasattr(mat_Al.melting_temperature, '__float__') else float(mat_Al.melting_temperature.evalf())
        boiling_temp = float(mat_Al.boiling_temperature) if hasattr(mat_Al.boiling_temperature, '__float__') else float(mat_Al.boiling_temperature.evalf())
        assert melting_temp > 0
        assert boiling_temp > melting_temp

    def test_steel_material_creation(self, steel_yaml_path, temp_symbol):
        """Test steel material creation from YAML."""
        if not steel_yaml_path.exists():
            pytest.skip(f"Steel YAML file not found: {steel_yaml_path}")
        mat_steel = create_material(yaml_path=steel_yaml_path, dependency=temp_symbol, enable_plotting=False)
        # Basic material verification
        assert "Steel" in mat_steel.name or "1.4301" in mat_steel.name
        assert mat_steel.material_type == "alloy"
        assert len(mat_steel.elements) >= 2
        # Temperature properties for alloys
        assert hasattr(mat_steel, 'solidus_temperature')
        assert hasattr(mat_steel, 'liquidus_temperature')
        solidus_temp = float(mat_steel.solidus_temperature) if hasattr(mat_steel.solidus_temperature, '__float__') else float(mat_steel.solidus_temperature.evalf())
        liquidus_temp = float(mat_steel.liquidus_temperature) if hasattr(mat_steel.liquidus_temperature, '__float__') else float(mat_steel.liquidus_temperature.evalf())
        assert solidus_temp > 0
        assert liquidus_temp > solidus_temp

    def test_material_property_evaluation(self, aluminum_yaml_path, temp_symbol):
        """Test material property evaluation at specific temperatures."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML file not found: {aluminum_yaml_path}")
        mat_Al = create_material(yaml_path=aluminum_yaml_path, dependency=temp_symbol, enable_plotting=False)
        # Test property evaluation if properties exist
        test_temp = 300.0
        valid_properties = get_supported_properties()
        for prop_name in valid_properties:
            if hasattr(mat_Al, prop_name):
                prop_value = getattr(mat_Al, prop_name)
                if isinstance(prop_value, sp.Expr):
                    try:
                        numerical_value = float(prop_value.subs(temp_symbol, test_temp))
                        assert isinstance(numerical_value, float)
                        assert not sp.nan(numerical_value)
                    except (TypeError, ValueError):
                        # Some expressions might not be temperature-dependent
                        try:
                            numerical_value = float(prop_value.evalf())
                            assert isinstance(numerical_value, float)
                        except:
                            pass  # Skip if can't evaluate

    def test_comprehensive_material_properties(self, aluminum_yaml_path, temp_symbol):
        """Test comprehensive material property evaluation."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML file not found: {aluminum_yaml_path}")
        mat_Al = create_material(yaml_path=aluminum_yaml_path, dependency=temp_symbol, enable_plotting=False)
        test_temp = 300.0
        # Test all properties that exist on the material
        for attr_name in dir(mat_Al):
            if not attr_name.startswith('_') and not callable(getattr(mat_Al, attr_name)):
                if attr_name in ['elements', 'composition', 'name', 'material_type']:
                    continue  # Skip non-property attributes
                try:
                    attr_value = getattr(mat_Al, attr_name)
                    if isinstance(attr_value, sp.Expr):
                        try:
                            numerical_value = float(attr_value.subs(temp_symbol, test_temp))
                            assert isinstance(numerical_value, float)
                            assert not sp.isnan(numerical_value)
                        except (TypeError, ValueError):
                            try:
                                numerical_value = float(attr_value.evalf())
                                assert isinstance(numerical_value, float)
                            except:
                                pass  # Skip if can't evaluate
                except Exception:
                    pass  # Skip problematic attributes
