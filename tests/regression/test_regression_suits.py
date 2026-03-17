"""Regression tests to ensure bugs don't reappear."""
import tempfile
from pathlib import Path
import pytest
import sympy as sp
from ruamel.yaml import YAML
from materforge import create_material

class TestRegressionSuite:
    """Regression tests for known issues and bug fixes."""
    def test_step_function_scalar_reference_regression(self):
        """Regression: step function referencing a scalar constant resolves correctly.
        Previously, scalar constants were not guaranteed to be processed before
        dependent step-function properties, causing KeyError in the resolver.
        """
        T = sp.Symbol('T')
        config = {
            'name': 'Step Reference Test',
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
            assert 'latent_heat' in material.property_names()
            assert float(material.latent_heat.subs(T, 1000)) == pytest.approx(0.0)
            assert float(material.latent_heat.subs(T, 1800)) == pytest.approx(250000.0)
        finally:
            yaml_path.unlink()

    def test_temperature_boundary_regression(self):
        """Regression: piecewise property evaluates exactly at tabular boundaries."""
        T = sp.Symbol('T')
        config = {
            'name': 'Boundary Test Material',
            'properties': {
                'melting_temperature': 933.47,
                'boiling_temperature': 2792.0,
                'heat_capacity': {
                    'dependency': [300, 600, 900],
                    'value': [900, 950, 1000],
                    'bounds': ['constant', 'constant'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            assert float(material.heat_capacity.subs(T, 300)) == pytest.approx(900.0)
            assert float(material.heat_capacity.subs(T, 600)) == pytest.approx(950.0)
            assert float(material.heat_capacity.subs(T, 900)) == pytest.approx(1000.0)
        finally:
            yaml_path.unlink()

    def test_empty_properties_regression(self):
        """Regression: empty properties block raises a clear ValueError."""
        T = sp.Symbol('T')
        config = {'name': 'Minimal Material', 'properties': {}}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            with pytest.raises(ValueError, match="cannot be empty"):
                create_material(yaml_path, T, enable_plotting=False)
        finally:
            yaml_path.unlink()

    def test_symbolic_temperature_consistency_regression(self):
        """Regression: different symbol names produce identical numerical results."""
        T1 = sp.Symbol('T')
        T2 = sp.Symbol('Temperature')
        config = {
            'name': 'Symbol Test Material',
            'properties': {
                'melting_temperature': 1811.0,
                'boiling_temperature': 3134.0,
                'thermal_expansion_coefficient': {
                    'dependency': [300, 600, 900],
                    'value': [1.2e-5, 1.4e-5, 1.6e-5],
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
            mat1 = create_material(yaml_path, T1, enable_plotting=False)
            mat2 = create_material(yaml_path, T2, enable_plotting=False)

            result1 = float(mat1.thermal_expansion_coefficient.subs(T1, 500))
            result2 = float(mat2.thermal_expansion_coefficient.subs(T2, 500))
            assert result1 == pytest.approx(result2)
        finally:
            yaml_path.unlink()

    def test_multiple_scalar_constants_regression(self):
        """Regression: multiple scalar constants all land on material before any
        dependent property is processed, regardless of YAML declaration order."""
        T = sp.Symbol('T')
        config = {
            'name': 'Multi-Scalar Test',
            'properties': {
                'solidus_temperature': 1400.0,
                'liquidus_temperature': 1450.0,
                'latent_heat_solidus': {
                    'dependency': 'solidus_temperature',
                    'value': [0.0, 100000.0],
                    'bounds': ['constant', 'constant'],
                },
                'latent_heat_liquidus': {
                    'dependency': 'liquidus_temperature',
                    'value': [100000.0, 250000.0],
                    'bounds': ['constant', 'constant'],
                },
            },
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            YAML().dump(config, f)
            yaml_path = Path(f.name)
        try:
            material = create_material(yaml_path, T, enable_plotting=False)
            assert float(material.latent_heat_solidus.subs(T, 1000)) == pytest.approx(0.0)
            assert float(material.latent_heat_solidus.subs(T, 1420)) == pytest.approx(100000.0)
            assert float(material.latent_heat_liquidus.subs(T, 1420)) == pytest.approx(100000.0)
            assert float(material.latent_heat_liquidus.subs(T, 1500)) == pytest.approx(250000.0)
        finally:
            yaml_path.unlink()
