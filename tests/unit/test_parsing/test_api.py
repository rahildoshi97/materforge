"""Tests for the main API module."""

import pytest
import tempfile
from pathlib import Path
import sympy as sp
from materforge.parsing.api import create_material, validate_yaml_file, get_material_property_names


# ------------------------------------------------------------------
# Shared YAML fixtures
# ------------------------------------------------------------------

VALID_PURE_METAL_YAML = """
name: TestMaterial
material_type: pure_metal
composition:
  Fe: 1.0
melting_temperature: 1811.0
boiling_temperature: 3134.0
properties:
  density: 7874.0
"""

INVALID_STRUCTURE_YAML = """
name: TestMaterial
invalid_structure: [
"""

INCOMPLETE_YAML = """
name: TestMaterial
# Missing all required fields
"""


def _write_temp_yaml(content: str) -> Path:
    """Write content to a named temp file and return its Path."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(content)
        return Path(f.name)


# ------------------------------------------------------------------
# TestCreateMaterial
# ------------------------------------------------------------------

class TestCreateMaterial:
    """Test the main create_material function."""

    def test_create_material_with_symbolic_temperature(self):
        """Material creation with a symbolic dependency succeeds."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)

            assert material.name == "TestMaterial"
            assert material.material_type == "pure_metal"

            # Use property_names() — hasattr always returns True on the dataclass
            # even for fields that were never set by the property processor
            assert 'density' in material.property_names(), (
                "Expected 'density' in dynamic property tracker after processing")
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
                    f"Property '{name}' is in property_names() but its value is None")
        finally:
            yaml_path.unlink()

    def test_create_material_file_not_found(self):
        """Missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            create_material("nonexistent.yaml", sp.Symbol('T'), enable_plotting=False)

    def test_create_material_non_symbol_dependency_raises_type_error(self):
        """Passing a float as dependency raises TypeError, not a generic Exception."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            with pytest.raises(TypeError):
                create_material(yaml_path, 500.0, enable_plotting=False)  # type: ignore
        finally:
            yaml_path.unlink()

    def test_create_material_invalid_yaml(self):
        """Malformed YAML raises an exception during parsing."""
        yaml_path = _write_temp_yaml(INVALID_STRUCTURE_YAML)
        try:
            with pytest.raises(Exception):
                create_material(yaml_path, sp.Symbol('T'), enable_plotting=False)
        finally:
            yaml_path.unlink()


# ------------------------------------------------------------------
# TestValidateYamlFile
# ------------------------------------------------------------------

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

    def test_validate_incomplete_yaml_raises_value_error(self):
        """YAML missing required fields raises ValueError."""
        yaml_path = _write_temp_yaml(INCOMPLETE_YAML)
        try:
            with pytest.raises(ValueError):
                validate_yaml_file(yaml_path)
        finally:
            yaml_path.unlink()


# ------------------------------------------------------------------
# TestDynamicPropertyTracker
# ------------------------------------------------------------------

class TestDynamicPropertyTracker:
    """
    Replaces TestGetSupportedProperties.

    Property names are no longer a fixed global list — any name defined in
    a YAML file is valid. These tests verify the dynamic tracker that replaced
    the old predefined field approach.
    """

    def test_property_names_returns_list(self):
        """property_names() returns a non-empty list after processing."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            names = material.property_names()
            assert isinstance(names, list)
            assert len(names) > 0, "Material has no processed properties"
        finally:
            yaml_path.unlink()

    def test_get_material_property_names_matches_instance(self):
        """API function get_material_property_names() agrees with material.property_names()."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            assert set(get_material_property_names(material)) == set(material.property_names())
        finally:
            yaml_path.unlink()

    def test_structural_fields_not_in_property_names(self):
        """Temperature and identity fields must never appear in property_names()."""
        yaml_path = _write_temp_yaml(VALID_PURE_METAL_YAML)
        try:
            T = sp.Symbol('T')
            material = create_material(yaml_path, T, enable_plotting=False)
            structural = {
                'name', 'material_type', 'elements', 'composition',
                'melting_temperature', 'boiling_temperature',
                'solidus_temperature', 'liquidus_temperature',
                'initial_boiling_temperature', 'final_boiling_temperature',
                'atomic_number', 'atomic_mass',
            }
            overlap = structural & set(material.property_names())
            assert not overlap, (
                f"Structural fields must not be tracked as properties: {overlap}")
        finally:
            yaml_path.unlink()

    def test_get_material_property_names_rejects_non_material(self):
        """get_material_property_names() raises ValueError for wrong input type."""
        with pytest.raises(ValueError):
            get_material_property_names("not_a_material")  # type: ignore
