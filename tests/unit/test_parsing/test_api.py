"""Tests for the main API module."""

import pytest
import tempfile
from pathlib import Path
import sympy as sp
from materforge.parsing.api import create_material, validate_yaml_file, get_supported_properties


class TestCreateMaterial:
    """Test the main create_material function."""

    def test_create_material_with_symbolic_temperature(self):
        """Test material creation with symbolic temperature."""
        # Create a minimal test YAML
        yaml_content = """
name: TestMaterial
material_type: pure_metal
composition:
  Fe: 1.0
melting_temperature: 1811.0
boiling_temperature: 3134.0
properties:
  density: 7874.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            yaml_path = Path(f.name)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "TestMaterial"
            assert material.material_type == "pure_metal"
            assert hasattr(material, 'density')
        finally:
            yaml_path.unlink()

    def test_create_material_file_not_found(self):
        """Test error handling for missing files."""
        with pytest.raises(FileNotFoundError):
            create_material("nonexistent.yaml", sp.Symbol('T'), enable_plotting=False)

    def test_create_material_invalid_yaml(self):
        """Test error handling for invalid YAML."""
        yaml_content = """
name: TestMaterial
invalid_structure: [
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(Exception):  # Could be various YAML errors
                create_material(yaml_path, 500.0)
        finally:
            yaml_path.unlink()


class TestValidateYamlFile:
    """Test YAML validation functionality."""

    def test_validate_valid_yaml(self):
        """Test validation of valid YAML file."""
        yaml_content = """
name: TestMaterial
material_type: pure_metal
composition:
  Fe: 1.0
melting_temperature: 1811.0
boiling_temperature: 3134.0
properties:
  density: 7874.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            yaml_path = Path(f.name)
        try:
            result = validate_yaml_file(yaml_path)
            assert result is True
        finally:
            yaml_path.unlink()

    def test_validate_missing_file(self):
        """Test validation error for missing file."""
        with pytest.raises(FileNotFoundError):
            validate_yaml_file("nonexistent.yaml")

    def test_validate_invalid_yaml_content(self):
        """Test validation error for invalid content."""
        yaml_content = """
name: TestMaterial
# Missing required fields
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError):
                validate_yaml_file(yaml_path)
        finally:
            yaml_path.unlink()


class TestGetSupportedProperties:
    """Test supported properties listing."""

    def test_get_supported_properties_returns_list(self):
        """Test that supported properties returns a list."""
        props = get_supported_properties()
        assert isinstance(props, list)
        assert len(props) > 0

    def test_get_supported_properties_contains_expected(self):
        """Test that expected properties are in the list."""
        props = get_supported_properties()
        expected_props = ['density', 'heat_capacity', 'thermal_diffusivity']
        for prop in expected_props:
            assert prop in props
