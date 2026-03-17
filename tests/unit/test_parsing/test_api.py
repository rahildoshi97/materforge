"""Tests for the main API module."""

import pytest
import tempfile
from pathlib import Path
import sympy as sp
from materforge import create_material, validate_yaml_file, get_material_property_names
from materforge.parsing.validation.errors import MaterialConfigError


VALID_PURE_METAL_YAML = """
name: TestMaterial
properties:
  melting_temperature: 1811.0
  boiling_temperature: 3134.0
  density: 7874.0
"""

INVALID_STRUCTURE_YAML = """
name: TestMaterial
invalid_structure: [
"""

INCOMPLETE_YAML = """
name: TestMaterial
# Missing properties block entirely
"""


def _write_temp_yaml(content: str) -> Path:
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(content)
        return Path(f.name)


class TestCreateMaterial:
    """Tests for the main create_material function."""

    def test_create_material_with_symbolic_temperature(self):
        """Material creation with a symbolic dependency succeeds."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "TestMaterial"
            assert 'density' in material.property_names()
        finally:
            yaml_path.unlink()

    def test_create_material_properties_are_not_none(self):
        """All names returned by property_names() must have non-None values."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            for name in material.property_names():
                assert getattr(material, name) is not None, (
                    f"Property '{name}' is tracked but its value is None")
        finally:
            yaml_path.unlink()

    def test_create_material_file_not_found(self):
        """Missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            create_material("nonexistent.yaml", sp.Symbol('T'), enable_plotting=False)

    def test_create_material_non_symbol_dependency_raises_type_error(self):
        """Passing a float as dependency raises TypeError."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            with pytest.raises(TypeError):
                create_material(yaml_path, 500.0, enable_plotting=False)  # type: ignore
        finally:
            yaml_path.unlink()

    def test_create_material_invalid_yaml(self):
        """Malformed YAML raises a scanner error during parsing."""
        yaml_path = _write_temp_yaml(INVALID_STRUCTURE_YAML)
        try:
            with pytest.raises(Exception):   # scanner.ScannerError - keeping Exception is fine here
                create_material(yaml_path, sp.Symbol('T'), enable_plotting=False)
        finally:
            yaml_path.unlink()


class TestValidateYamlFile:
    """Test YAML validation functionality."""

    def test_validate_valid_yaml(self):
        """Valid YAML file returns True."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            assert validate_yaml_file(yaml_path) is True
        finally:
            yaml_path.unlink()

    def test_validate_missing_file(self):
        """Missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            validate_yaml_file("nonexistent.yaml")

    def test_validate_incomplete_yaml_raises_material_config_error(self):
        """YAML missing the 'properties' block raises MaterialConfigError."""
        yaml_path = _write_temp_yaml(INCOMPLETE_YAML)
        try:
            with pytest.raises(MaterialConfigError):
                validate_yaml_file(yaml_path)
        finally:
            yaml_path.unlink()


class TestDynamicPropertyTracker:
    """
    Property names are no longer a fixed global list - any name defined in
    a YAML file is valid. These tests verify the dynamic tracker.
    """

    def test_property_names_returns_set(self):
        """property_names() returns a non-empty set after processing."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            names = material.property_names()
            assert isinstance(names, set)
            assert len(names) > 0, "Material has no processed properties"
        finally:
            yaml_path.unlink()

    def test_get_material_property_names_matches_instance(self):
        """API get_material_property_names() agrees with material.property_names()."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            assert set(get_material_property_names(material)) == material.property_names()
        finally:
            yaml_path.unlink()

    def test_name_field_not_in_property_names(self):
        """'name' is a constructor field and must never appear in property_names()."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            assert 'name' not in material.property_names()
        finally:
            yaml_path.unlink()

    def test_scalar_temps_are_in_property_names(self):
        """Scalar temperatures defined in YAML are tracked as dynamic properties."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            assert 'melting_temperature' in material.property_names()
            assert 'boiling_temperature' in material.property_names()
        finally:
            yaml_path.unlink()

    def test_get_material_property_names_rejects_non_material(self):
        """get_material_property_names() raises TypeError for wrong input type."""
        with pytest.raises(TypeError):
            get_material_property_names("not_a_material")  # type: ignore
