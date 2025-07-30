"""Regression tests to ensure bugs don't reappear."""

import sympy as sp
import tempfile
from pathlib import Path
from ruamel.yaml import YAML
from materforge.parsing.api import create_material

class TestRegressionSuite:
    """Regression tests for known issues and bug fixes."""
    def test_composition_sum_validation_regression(self):
        """Regression test for composition sum validation bug."""
        # This test ensures that composition validation works correctly
        T = sp.Symbol('T')
        # Configuration that should be valid but might fail due to floating point precision
        config = {
            'name': 'Precision Test Alloy',
            'material_type': 'alloy',
            'composition': {
                'Fe': 0.333333333,
                'Cr': 0.333333333,
                'Ni': 0.333333334  # Sum = 1.0 with floating point precision
            },
            'solidus_temperature': 1400.0,
            'liquidus_temperature': 1500.0,
            'initial_boiling_temperature': 2800.0,
            'final_boiling_temperature': 2900.0,
            'properties': {'density': 7850.0}
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            # This should not raise a composition error
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "Precision Test Alloy"
            assert abs(sum(material.composition) - 1.0) < 1e-10
        finally:
            yaml_path.unlink()

    def test_temperature_boundary_regression(self):
        """Regression test for temperature boundary handling."""
        T = sp.Symbol('T')
        config = {
            'name': 'Boundary Test Material',
            'material_type': 'pure_metal',
            'composition': {'Al': 1.0},
            'melting_temperature': 933.47,
            'boiling_temperature': 2792.0,
            'properties': {
                'heat_capacity': {
                    'dependency': [300, 600, 900],
                    'value': [900, 950, 1000],
                    'bounds': ['constant', 'constant']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            # Test evaluation exactly at boundaries
            heat_cap_300 = float(material.heat_capacity.subs(T, 300))
            heat_cap_600 = float(material.heat_capacity.subs(T, 600))
            heat_cap_900 = float(material.heat_capacity.subs(T, 900))
            assert heat_cap_300 == 900.0
            assert heat_cap_600 == 950.0
            assert heat_cap_900 == 1000.0
        finally:
            yaml_path.unlink()

    def test_empty_properties_regression(self):
        """Regression test for handling materials with minimal properties."""
        T = sp.Symbol('T')
        config = {
            'name': 'Minimal Material',
            'material_type': 'pure_metal',
            'composition': {'Cu': 1.0},
            'melting_temperature': 1357.77,
            'boiling_temperature': 2835.0,
            'properties': {}  # No additional properties
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            # This should work without errors
            material = create_material(yaml_path, T, enable_plotting=False)
            assert material.name == "Minimal Material"
            assert len(material.elements) == 1
            assert material.elements[0].name == "Copper"

        finally:
            yaml_path.unlink()

    def test_symbolic_temperature_consistency_regression(self):
        """Regression test for symbolic temperature consistency."""
        # Bug: Previously, different temperature symbols could cause issues
        T1 = sp.Symbol('T')
        T2 = sp.Symbol('Temperature')
        config = {
            'name': 'Symbol Test Material',
            'material_type': 'pure_metal',
            'composition': {'Fe': 1.0},
            'melting_temperature': 1811.0,
            'boiling_temperature': 3134.0,
            'properties': {
                'thermal_expansion_coefficient': {
                    'dependency': [300, 600, 900],
                    'value': [1.2e-5, 1.4e-5, 1.6e-5],
                    'bounds': ['constant', 'constant']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            # Create with first symbol
            material1 = create_material(yaml_path, T1, enable_plotting=False)
            # Create with second symbol
            material2 = create_material(yaml_path, T2, enable_plotting=False)
            # Both should work and give consistent results when evaluated
            result1 = float(material1.thermal_expansion_coefficient.subs(T1, 500))
            result2 = float(material2.thermal_expansion_coefficient.subs(T2, 500))
            assert abs(result1 - result2) < 1e-15
        finally:
            yaml_path.unlink()
