"""Unit tests for YAML parser components."""

import pytest
import tempfile
from pathlib import Path
from ruamel.yaml import YAML
import sympy as sp
from pymatlib.parsing.config.material_yaml_parser import MaterialYAMLParser

class TestMaterialYAMLParser:
    """Test cases for MaterialYAMLParser."""
    def test_yaml_parser_initialization_valid_file(self):
        """Test parser initialization with valid YAML file."""
        config = {
            'name': 'Test Material',
            'material_type': 'pure_metal',
            'composition': {'Al': 1.0},
            'melting_temperature': 933.47,
            'boiling_temperature': 2792.0,
            'properties': {'density': 2700.0}
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            parser = MaterialYAMLParser(yaml_path)
            assert parser.config_path == yaml_path
            assert parser.config['name'] == 'Test Material'
            assert parser.config['material_type'] == 'pure_metal'
        finally:
            yaml_path.unlink()

    def test_yaml_parser_invalid_file_path(self):
        """Test parser initialization with invalid file path."""
        invalid_path = Path("nonexistent_file.yaml")
        with pytest.raises(FileNotFoundError, match="YAML file not found"):
            MaterialYAMLParser(invalid_path)

    def test_yaml_parser_invalid_yaml_syntax(self):
        """Test parser with invalid YAML syntax."""
        invalid_yaml = """
        name: Test Material
        material_type: pure_metal
        composition:
          Al: 1.0
        melting_temperature: 933.47
        boiling_temperature: 2792.0
        properties:
          density: 2700.0
        invalid_syntax: [unclosed list
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(invalid_yaml)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(Exception):  # YAML syntax errors can be various types
                MaterialYAMLParser(yaml_path)
        finally:
            yaml_path.unlink()

    def test_yaml_parser_missing_required_fields(self):
        """Test parser with missing required fields."""
        incomplete_config = {
            'name': 'Incomplete Material',
            # Missing material_type, composition, etc.
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(incomplete_config, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="Missing required field"):
                MaterialYAMLParser(yaml_path)
        finally:
            yaml_path.unlink()

    def test_create_material_from_parser(self):
        """Test complete material creation through parser."""
        config = {
            'name': 'Parser Test Material',
            'material_type': 'pure_metal',
            'composition': {'Cu': 1.0},
            'melting_temperature': 1357.77,
            'boiling_temperature': 2835.0,
            'properties': {
                'density': 8960.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900],
                    'value': [385, 420, 455],
                    'bounds': ['constant', 'constant']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            T = sp.Symbol('T')
            parser = MaterialYAMLParser(yaml_path)
            material = parser.create_material(dependency=T, enable_plotting=False)
            assert material.name == "Parser Test Material"
            assert material.material_type == "pure_metal"
            assert len(material.elements) == 1
            assert material.elements[0].atomic_number == 29.0  # Copper
            # Test property evaluation
            if hasattr(material.density, 'evalf'):
                density_val = float(material.density.evalf())
            else:
                density_val = float(material.density)
            assert density_val == 8960.0
        finally:
            yaml_path.unlink()

    def test_yaml_parser_composition_validation(self):
        """Test composition validation in parser."""
        invalid_config = {
            'name': 'Invalid Composition',
            'material_type': 'pure_metal',
            'composition': {'Al': 0.8},  # Sum != 1.0
            'melting_temperature': 933.47,
            'boiling_temperature': 2792.0,
            'properties': {}
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(invalid_config, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="Composition fractions must sum to 1.0"):
                MaterialYAMLParser(yaml_path)
        finally:
            yaml_path.unlink()

    def test_yaml_parser_property_validation(self):
        """Test property name validation."""
        config_with_invalid_property = {
            'name': 'Invalid Property Material',
            'material_type': 'pure_metal',
            'composition': {'Al': 1.0},
            'melting_temperature': 933.47,
            'boiling_temperature': 2792.0,
            'properties': {
                'invalid_property_name': 2700.0  # Not in VALID_YAML_PROPERTIES
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config_with_invalid_property, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="Invalid properties found"):
                MaterialYAMLParser(yaml_path)
        finally:
            yaml_path.unlink()
