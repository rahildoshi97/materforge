"""Integration tests for material creation workflows."""

import pytest
import tempfile
from pathlib import Path
import sympy as sp
from ruamel.yaml import YAML
from materforge.parsing.api import create_material

class TestMaterialCreation:
    """Integration tests for complete material creation workflows."""
    def test_create_pure_metal_complete_workflow(self):
        """Test complete pure metal creation workflow."""
        T = sp.Symbol('T')
        config = {
            'name': 'Test Titanium',
            'material_type': 'pure_metal',
            'composition': {'Ti': 1.0},
            'melting_temperature': 1941.0,
            'boiling_temperature': 3560.0,
            'properties': {
                'density': 4506.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900, 1200],
                    'value': [523, 565, 590, 615],
                    'bounds': ['constant', 'constant']
                },
                'heat_conductivity': {
                    'dependency': [300, 600, 900],
                    'value': [21.9, 24.5, 27.1],
                    'bounds': ['extrapolate', 'extrapolate']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            # Verify basic properties
            assert material.name == "Test Titanium"
            assert material.material_type == "pure_metal"
            assert len(material.elements) == 1
            # Verify temperature properties
            assert float(material.melting_temperature) == 1941.0
            assert float(material.boiling_temperature) == 3560.0
            # Verify processed properties
            assert hasattr(material, 'density')
            assert hasattr(material, 'heat_capacity')
            assert hasattr(material, 'heat_conductivity')
            # Test property evaluation
            density_val = float(material.density.evalf()) if hasattr(material.density, 'evalf') else float(material.density)
            assert density_val == 4506.0
            # Test temperature-dependent properties
            heat_cap_500K = float(material.heat_capacity.subs(T, 500))
            assert 523 < heat_cap_500K < 615
            thermal_cond_500K = float(material.heat_conductivity.subs(T, 500))
            assert thermal_cond_500K > 0
        finally:
            yaml_path.unlink()

    def test_create_complex_alloy_workflow(self):
        """Test complex alloy creation with multiple properties."""
        T = sp.Symbol('T')
        config = {
            'name': 'Test Inconel 718',
            'material_type': 'alloy',
            'composition': {
                'Ni': 0.525,
                'Cr': 0.19,
                'Fe': 0.185,
                'Mn': 0.053,
                'Mo': 0.031,
                'Ti': 0.009,
                'Al': 0.005,
                'Cu': 0.002
            },
            'solidus_temperature': 1533.0,
            'liquidus_temperature': 1609.0,
            'initial_boiling_temperature': 3000.0,
            'final_boiling_temperature': 3200.0,
            'properties': {
                'density': 8220.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900, 1200, 1500],
                    'value': [435, 485, 520, 555, 590],
                    'bounds': ['constant', 'constant']
                },
                'heat_conductivity': {
                    'dependency': [300, 600, 900],
                    'equation': ['211000 - 45*T', '180000 - 20*T'],
                    'bounds': ['constant', 'extrapolate']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            # Verify alloy properties
            assert material.name == "Test Inconel 718"
            assert material.material_type == "alloy"
            assert len(material.elements) == 8
            # Verify composition sums to 1
            assert abs(sum(material.composition) - 1.0) < 1e-10
            # Verify temperature properties
            assert float(material.solidus_temperature) == 1533.0
            assert float(material.liquidus_temperature) == 1609.0
            # Test piecewise equation property
            elastic_mod_500K = float(material.heat_conductivity.subs(T, 500))
            expected_500K = 211000 - 45*500  # First equation
            assert abs(elastic_mod_500K - expected_500K) < 1.0
        finally:
            yaml_path.unlink()
