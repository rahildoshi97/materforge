"""Unit tests for YAML parser components."""

import pytest
import tempfile
from pathlib import Path
from ruamel.yaml import YAML
import sympy as sp
from materforge.parsing.config.material_yaml_parser import MaterialYAMLParser

class TestMaterialYAMLParser:
    """Test cases for MaterialYAMLParser."""
    def test_yaml_parser_initialization_valid_file(self):
        """Parser reads name and properties from a valid YAML file."""
        config = {
            'name': 'Test Material',
            'properties': {
                'melting_temperature': 933.47,
                'boiling_temperature': 2792.0,
                'density': 2700.0,
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            parser = MaterialYAMLParser(yaml_path)
            assert parser.config_path == yaml_path
            assert parser.config['name'] == 'Test Material'
        finally:
            yaml_path.unlink()

    def test_yaml_parser_invalid_file_path(self):
        """Non-existent path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="YAML file not found"):
            MaterialYAMLParser(Path("nonexistent_file.yaml"))

    def test_yaml_parser_invalid_yaml_syntax(self):
        """Malformed YAML syntax raises a parsing exception."""
        invalid_yaml = (
            "name: Test Material\n"
            "properties:\n"
            "  density: 2700.0\n"
            "  invalid_syntax: [unclosed list\n"
        )
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(invalid_yaml)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(Exception):
                MaterialYAMLParser(yaml_path)
        finally:
            yaml_path.unlink()

    def test_yaml_parser_missing_required_fields(self):
        """Config without a 'properties' key raises ValueError."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump({'name': 'Incomplete Material'}, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="Missing required field"):
                MaterialYAMLParser(yaml_path)
        finally:
            yaml_path.unlink()

    def test_create_material_from_parser(self):
        """Parser produces a Material with correct name and dynamic properties."""
        config = {
            'name': 'Parser Test Material',
            'properties': {
                'melting_temperature': 1357.77,
                'boiling_temperature': 2835.0,
                'density': 8960.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900],
                    'value': [385, 420, 455],
                    'bounds': ['constant', 'constant'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            T = sp.Symbol('T')
            material = MaterialYAMLParser(yaml_path).create_material(dependency=T, enable_plotting=False)
            assert material.name == "Parser Test Material"
            assert float(material.density) == pytest.approx(8960.0)
            assert float(material.melting_temperature) == pytest.approx(1357.77)
            assert 'heat_capacity' in material.property_names()
        finally:
            yaml_path.unlink()

    def test_yaml_parser_empty_properties_validation(self):
        """Empty properties block raises ValueError."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump({'name': 'Empty Props', 'properties': {}}, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="cannot be empty"):
                MaterialYAMLParser(yaml_path).create_material(
                    dependency=sp.Symbol('T'), enable_plotting=False)
        finally:
            yaml_path.unlink()

    def test_yaml_parser_missing_bounds_validation(self):
        """Tabular property without 'bounds' key raises ValueError."""
        config = {
            'name': 'Missing Bounds',
            'properties': {
                'heat_capacity': {
                    'dependency': [300, 400, 500],
                    'value': [900, 950, 1000],
                    # bounds intentionally omitted
                },
            },
        }
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.yaml', delete=False
        ) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="missing required keys"):
                MaterialYAMLParser(yaml_path).create_material(dependency=sp.Symbol('T'), enable_plotting=False)
        finally:
            yaml_path.unlink()
