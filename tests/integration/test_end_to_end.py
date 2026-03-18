"""End-to-end integration tests."""
import tempfile
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from materforge import create_material

class TestEndToEnd:
    """End-to-end integration tests."""
    def test_create_aluminum_material_from_existing_yaml(self, aluminum_yaml_path, temp_symbol):
        """Material creation from existing aluminum YAML file."""
        if not aluminum_yaml_path.exists():
            pytest.skip(f"Aluminum YAML not found: {aluminum_yaml_path}")
        material = create_material(aluminum_yaml_path, temp_symbol, enable_plotting=False)
        assert "Aluminum" in material.name or "Al" in material.name
        assert 'density' in material.property_names()
        density_val = float(material.density) if isinstance(
            material.density, (int, float)) else float(
            material.density.subs(temp_symbol, 300))
        assert density_val > 0

    def test_create_steel_material_from_existing_yaml(
        self, steel_yaml_path, temp_symbol
    ):
        """Material creation from existing steel YAML file."""
        if not steel_yaml_path.exists():
            pytest.skip(f"Steel YAML not found: {steel_yaml_path}")
        material = create_material(steel_yaml_path, temp_symbol, enable_plotting=False)
        assert "Steel" in material.name or "1.4301" in material.name
        # Temperature scalars are dynamic properties now
        assert 'solidus_temperature' in material.property_names()
        assert 'liquidus_temperature' in material.property_names()
        assert float(material.solidus_temperature) > 0
        assert float(material.liquidus_temperature) > float(
            material.solidus_temperature)
    def test_create_material_from_temp_yaml(self, temp_symbol):
        """Complete material creation workflow with a temporary YAML."""
        config = {
            'name': 'Test Aluminum',
            'properties': {
                'melting_temperature': 933.47,
                'boiling_temperature': 2792.0,
                'density': 2700.0,
                'heat_capacity': {
                    'dependency': [300, 400, 500],
                    'value': [900, 950, 1000],
                    'bounds': ['constant', 'constant'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, temp_symbol, enable_plotting=False)
            assert material.name == "Test Aluminum"
            assert 'density' in material.property_names()
            assert 'heat_capacity' in material.property_names()
            assert float(material.density) == 2700.0
            assert float(material.melting_temperature) == pytest.approx(933.47)
            heat_cap_at_400 = float(material.heat_capacity.subs(temp_symbol, 400))
            assert 900 < heat_cap_at_400 < 1000
        finally:
            yaml_path.unlink()

    def test_create_alloy_from_temp_yaml(self, temp_symbol):
        """Alloy creation workflow with a temporary YAML."""
        config = {
            'name': 'Test Stainless Steel 304L',
            'properties': {
            'solidus_temperature': 1400.0,
            'liquidus_temperature': 1450.0,
            'initial_boiling_temperature': 2800.0,
            'final_boiling_temperature': 2900.0,
                'density': 7850.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900, 1200],
                    'value': [500, 550, 600, 650],
                    'bounds': ['constant', 'constant'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, temp_symbol, enable_plotting=False)

            assert material.name == "Test Stainless Steel 304L"
            assert float(material.solidus_temperature) == pytest.approx(1400.0)
            assert float(material.liquidus_temperature) == pytest.approx(1450.0)
            assert float(material.liquidus_temperature) > float(
                material.solidus_temperature)
        finally:
            yaml_path.unlink()

    def test_error_handling_empty_properties(self, temp_symbol):
        """Empty properties block raises ValueError."""
        config = {'name': 'Invalid Material', 'properties': {}}
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.yaml', delete=False
        ) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="cannot be empty"):
                create_material(yaml_path, temp_symbol, enable_plotting=False)
        finally:
            yaml_path.unlink()

    def test_error_handling_missing_bounds(self, temp_symbol):
        """Tabular property missing required 'bounds' key raises ValueError."""
        config = {
            'name': 'Invalid Material',
            'properties': {
                'heat_capacity': {
                    'dependency': [300, 400, 500],
                    'value': [900, 950, 1000],
                    # bounds intentionally omitted
                },
            },
        }
        with tempfile.NamedTemporaryFile( mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="missing required keys"):
                create_material(yaml_path, temp_symbol, enable_plotting=False)
        finally:
            yaml_path.unlink()

    def test_error_handling_nonexistent_file(self, temp_symbol):
        """Non-existent YAML path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            create_material(
                Path("nonexistent.yaml"), temp_symbol, enable_plotting=False)
