"""Integration tests for material creation workflows."""
import tempfile
from pathlib import Path
import pytest
import sympy as sp
from ruamel.yaml import YAML
from materforge.parsing.api import create_material

class TestMaterialCreation:
    """Integration tests for complete material creation workflows."""
    def test_create_pure_metal_complete_workflow(self):
        """Complete pure-metal creation: scalar constants + tabular properties."""
        T = sp.Symbol('T')
        config = {
            'name': 'Test Titanium',
            'properties': {
                'melting_temperature': 1941.0,
                'boiling_temperature': 3560.0,
                'density': 4506.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900, 1200],
                    'value': [523, 565, 590, 615],
                    'bounds': ['constant', 'constant'],
                },
                'heat_conductivity': {
                    'dependency': [300, 600, 900],
                    'value': [21.9, 24.5, 27.1],
                    'bounds': ['extrapolate', 'extrapolate'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "Test Titanium"
            assert float(material.melting_temperature) == pytest.approx(1941.0)
            assert float(material.boiling_temperature) == pytest.approx(3560.0)
            assert 'density' in material.property_names()
            assert 'heat_capacity' in material.property_names()
            assert 'heat_conductivity' in material.property_names()
            assert float(material.density) == pytest.approx(4506.0)
            heat_cap_500 = float(material.heat_capacity.subs(T, 500))
            assert 523 < heat_cap_500 < 615
            cond_500 = float(material.heat_conductivity.subs(T, 500))
            assert cond_500 > 0
        finally:
            yaml_path.unlink()

    def test_create_complex_alloy_workflow(self):
        """Complex alloy with scalar temperatures + piecewise equation property."""
        T = sp.Symbol('T')
        config = {
            'name': 'Test Inconel 718',
            'properties': {
            'solidus_temperature': 1533.0,
            'liquidus_temperature': 1609.0,
            'initial_boiling_temperature': 3000.0,
            'final_boiling_temperature': 3200.0,
                'density': 8220.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900, 1200, 1500],
                    'value': [435, 485, 520, 555, 590],
                    'bounds': ['constant', 'constant'],
                },
                'heat_conductivity': {
                    'dependency': [300, 600, 900],
                    'equation': ['211000 - 45*T', '180000 - 20*T'],
                    'bounds': ['constant', 'extrapolate'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "Test Inconel 718"
            assert float(material.solidus_temperature) == pytest.approx(1533.0)
            assert float(material.liquidus_temperature) == pytest.approx(1609.0)
            assert float(material.liquidus_temperature) > float(
                material.solidus_temperature)
            assert 'density' in material.property_names()
            assert 'heat_capacity' in material.property_names()
            assert 'heat_conductivity' in material.property_names()
            assert float(material.density) == pytest.approx(8220.0)
            # Piecewise equation: first piece active at T=500 (300 ≤ T < 600)
            cond_500 = float(material.heat_conductivity.subs(T, 500))
            assert cond_500 == pytest.approx(211000 - 45 * 500)
        finally:
            yaml_path.unlink()

    def test_create_material_with_step_function(self):
        """Material with a step function referencing a scalar constant."""
        T = sp.Symbol('T')
        config = {
            'name': 'Test Step Material',
            'properties': {
                'solidus_temperature': 1400.0,
                'latent_heat': {
                    'dependency': 'solidus_temperature',
                    'value': [0.0, 250000.0],
                    'bounds': ['constant', 'constant'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.yaml', delete=False
        ) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "Test Step Material"
            assert 'latent_heat' in material.property_names()
            # Below transition: value should be 0.0
            assert float(material.latent_heat.subs(T, 1000)) == pytest.approx(0.0)
            # Above transition: value should be 250000.0
            assert float(material.latent_heat.subs(T, 1800)) == pytest.approx(250000.0)
        finally:
            yaml_path.unlink()
