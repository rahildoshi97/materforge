"""End-to-end integration tests."""

import pytest
import tempfile
from ruamel.yaml import YAML
from pathlib import Path

from pymatlib.parsing.api import create_material

class TestEndToEnd:
    """End-to-end integration tests."""
    def test_create_aluminum_material_from_existing_yaml(self, temp_symbol):
        """Test material creation from existing aluminum YAML file."""
        # Try multiple possible paths for the aluminum YAML file
        possible_paths = [
            Path("src/pymatlib/data/materials/pure_metals/Al/Al.yaml"),
            Path("../src/pymatlib/data/materials/pure_metals/Al/Al.yaml"),
            Path("pymatlib/data/materials/pure_metals/Al/Al.yaml"),
        ]
        aluminum_yaml_path = None
        for path in possible_paths:
            if path.exists():
                aluminum_yaml_path = path
                break
        # Skip test if file doesn't exist
        if aluminum_yaml_path is None or not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML file not found. Tried paths: {[str(p) for p in possible_paths]}")
        # Create material
        material = create_material(aluminum_yaml_path, temp_symbol, enable_plotting=False)
        # Verify material properties
        assert material.name == "Aluminum"
        assert material.material_type == "pure_metal"
        assert len(material.elements) == 1
        assert material.elements[0].name in ["Aluminum", "Aluminium"]  # Handle both spellings
        # Verify processed properties
        assert hasattr(material, 'density')
        # Handle SymPy expressions properly
        if hasattr(material.density, 'evalf'):
            try:
                density_value = float(material.density.evalf())
            except (TypeError, ValueError):
                # If evalf() fails, try substituting a temperature value
                density_value = float(material.density.subs(temp_symbol, 300))
        elif hasattr(material.density, 'subs'):
            # If it's a symbolic expression, substitute a temperature value
            density_value = float(material.density.subs(temp_symbol, 300))
        else:
            density_value = float(material.density)
        assert density_value > 0

    def test_create_steel_material_from_existing_yaml(self, temp_symbol):
        """Test material creation from existing steel YAML file."""
        # Try multiple possible paths for the steel YAML file
        possible_paths = [
            Path("src/pymatlib/data/materials/alloys/1.4301/1.4301.yaml"),
            Path("../src/pymatlib/data/materials/alloys/1.4301/1.4301.yaml"),
            Path("pymatlib/data/materials/alloys/1.4301/1.4301.yaml"),
        ]
        steel_yaml_path = None
        for path in possible_paths:
            if path.exists():
                steel_yaml_path = path
                break
        # Skip test if file doesn't exist
        if steel_yaml_path is None or not steel_yaml_path.exists():
            pytest.skip(f"Steel YAML file not found. Tried paths: {[str(p) for p in possible_paths]}")
        # Create material
        material = create_material(steel_yaml_path, temp_symbol, enable_plotting=False)
        # Verify material properties
        assert "Steel" in material.name or "1.4301" in material.name
        assert material.material_type == "alloy"
        assert len(material.elements) >= 2
        # Verify temperature properties
        assert hasattr(material, 'solidus_temperature')
        assert hasattr(material, 'liquidus_temperature')
        # Handle SymPy expressions for temperature properties
        if hasattr(material.solidus_temperature, 'evalf'):
            solidus_temp = float(material.solidus_temperature.evalf())
        else:
            solidus_temp = float(material.solidus_temperature)
        if hasattr(material.liquidus_temperature, 'evalf'):
            liquidus_temp = float(material.liquidus_temperature.evalf())
        else:
            liquidus_temp = float(material.liquidus_temperature)
        assert solidus_temp > 0
        assert liquidus_temp > solidus_temp

    def test_create_aluminum_material_from_temp_yaml(self, temp_symbol):
        """Test complete material creation workflow with temporary YAML."""
        sample_config = {
            'name': 'Test Aluminum',
            'material_type': 'pure_metal',
            'composition': {'Al': 1.0},
            'melting_temperature': 933.47,
            'boiling_temperature': 2792.0,
            'properties': {
                'density': 2700.0,
                'heat_capacity': {
                    'temperature': [300, 400, 500],
                    'value': [900, 950, 1000],
                    'bounds': ['constant', 'constant']
                }
            }
        }
        # Create temporary YAML file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(sample_config, f)
            yaml_path = Path(f.name)
        try:
            # Create material
            material = create_material(yaml_path, temp_symbol, enable_plotting=False)
            # Verify material properties
            assert material.name == "Test Aluminum"
            assert material.material_type == "pure_metal"
            assert len(material.elements) == 1
            # Verify processed properties
            assert hasattr(material, 'density')
            assert hasattr(material, 'heat_capacity')
            # Test property evaluation - handle SymPy expressions
            if hasattr(material.density, 'evalf'):
                density_value = float(material.density.evalf())
            else:
                density_value = float(material.density)
            assert density_value == 2700.0
            # Test piecewise property evaluation
            heat_cap_at_400K = float(material.heat_capacity.subs(temp_symbol, 400))
            assert 900 < heat_cap_at_400K < 1000  # Should be interpolated
        finally:
            yaml_path.unlink()  # Clean up

    def test_create_steel_alloy_from_temp_yaml(self, temp_symbol):
        """Test alloy creation workflow with temporary YAML."""
        sample_config = {
            'name': 'Test Stainless Steel 304L',
            'material_type': 'alloy',
            'composition': {'Fe': 0.68, 'Cr': 0.20, 'Ni': 0.10, 'C': 0.02},
            'solidus_temperature': 1400.0,
            'liquidus_temperature': 1450.0,
            'initial_boiling_temperature': 2800.0,
            'final_boiling_temperature': 2900.0,
            'properties': {
                'density': 7850.0,
                'heat_capacity': {
                    'temperature': [300, 600, 900, 1200],
                    'value': [500, 550, 600, 650],
                    'bounds': ['constant', 'constant']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(sample_config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, temp_symbol, enable_plotting=False)

            assert material.name == "Test Stainless Steel 304L"
            assert material.material_type == "alloy"
            assert len(material.elements) >= 2
            # Verify temperature properties - handle SymPy expressions
            if hasattr(material.solidus_temperature, 'evalf'):
                solidus_temp = float(material.solidus_temperature.evalf())
            else:
                solidus_temp = float(material.solidus_temperature)
            if hasattr(material.liquidus_temperature, 'evalf'):
                liquidus_temp = float(material.liquidus_temperature.evalf())
            else:
                liquidus_temp = float(material.liquidus_temperature)
            assert solidus_temp == 1400.0
            assert liquidus_temp == 1450.0
        finally:
            yaml_path.unlink()

    def test_error_handling_invalid_yaml(self, temp_symbol):
        """Test error handling for invalid YAML."""
        invalid_config = {
            'name': 'Invalid Material',
            'material_type': 'pure_metal',
            'composition': {'Al': 0.8},  # Sum is not 1.0
            'melting_temperature': 933.47,
            'boiling_temperature': 2792.0,
            'properties': {}
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(invalid_config, f)
            yaml_path = Path(f.name)
        try:
            # Update the regex pattern to match the actual error message
            with pytest.raises(ValueError, match=r"Composition fractions must sum to 1\.0"):
                create_material(yaml_path, temp_symbol, enable_plotting=False)
        finally:
            yaml_path.unlink()
