"""Unit tests for ComputedPropertyProcessor."""
import pytest
import sympy as sp
from materforge.core.materials import Material
from materforge.parsing.processors.computed_property_processor import ComputedPropertyProcessor
from materforge.parsing.validation.errors import CircularDependencyError

def _make_material(**scalar_props) -> Material:
    mat = Material(name="TestMaterial")
    for k, v in scalar_props.items():
        setattr(mat, k, v)
    return mat


class TestComputedPropertyProcessor:
    """Tests for ComputedPropertyProcessor."""

    def test_computed_property_processor_initialization(self):
        properties = {
            'density': 2700.0,
            'heat_capacity': {'temperature': [300, 400], 'value': [900, 950]},
        }
        processed = set()
        processor = ComputedPropertyProcessor(properties, processed)
        assert processor.properties == properties
        assert processor.processed_properties == processed

    def test_extract_equation_dependencies_simple(self):
        deps = ComputedPropertyProcessor._extract_equation_dependencies(
            "density * heat_capacity")
        assert set(deps) == {'density', 'heat_capacity'}

    def test_extract_equation_dependencies_complex(self):
        deps = ComputedPropertyProcessor._extract_equation_dependencies(
            "density * heat_capacity + thermal_conductivity / viscosity")
        assert set(deps) == {
            'density', 'heat_capacity', 'thermal_conductivity', 'viscosity'}

    def test_extract_equation_dependencies_excludes_temperature_symbol(self):
        """Temperature symbol T must not appear in extracted dependencies."""
        deps = ComputedPropertyProcessor._extract_equation_dependencies(
            "density * T + heat_capacity")
        assert set(deps) == {'density', 'heat_capacity'}
        assert 'T' not in deps

    def test_validate_circular_dependencies_no_cycle(self):
        properties = {
            'density': 2700.0,
            'volume': {'equation': 'mass / density'},
            'mass': 1000.0,
        }
        processor = ComputedPropertyProcessor(properties, set())
        processor._validate_circular_dependencies('volume', ['mass', 'density'], set())

    def test_validate_circular_dependencies_with_cycle(self):
        """Cycle a -> b -> c -> a raises CircularDependencyError."""
        properties = {
            'prop_a': {'equation': 'prop_b + 100'},
            'prop_b': {'equation': 'prop_c * 2'},
            'prop_c': {'equation': 'prop_a / 3'},
        }
        processor = ComputedPropertyProcessor(properties, set())
        with pytest.raises(CircularDependencyError):
            processor._validate_circular_dependencies('prop_a', ['prop_b'], set())

    def test_validate_circular_dependencies_self_reference(self):
        properties = {'recursive_prop': {'equation': 'recursive_prop + 1'}}
        processor = ComputedPropertyProcessor(properties, set())
        with pytest.raises(CircularDependencyError):
            processor._validate_circular_dependencies('recursive_prop', ['recursive_prop'], set())

    def test_process_computed_property_simple(self):
        """Product of two scalar constants is computed and tracked."""
        mat = _make_material(density=sp.Float(2700.0), volume=sp.Float(0.001))
        properties = {
            'density': 2700.0,
            'volume': 0.001,
            'mass': {'equation': 'density * volume', 'dependency': [300, 400, 500]},
        }
        processed = {'density', 'volume'}
        processor = ComputedPropertyProcessor(properties, processed)
        processor.process_computed_property(mat, 'mass', sp.Symbol('T'))
        assert hasattr(mat, 'mass')
        assert 'mass' in processed

    def test_process_computed_property_missing_dependency(self):
        """Referencing undefined symbols raises ValueError."""
        mat = Material(name="TestMaterial")
        properties = {
            'mass': {
                'equation': 'density * volume',  # neither assigned to mat
                'dependency': [300, 400, 500],
            },
        }
        processor = ComputedPropertyProcessor(properties, set())
        with pytest.raises(ValueError, match="Missing dependencies in expression"):
            processor.process_computed_property(mat, 'mass', sp.Symbol('T'))

    def test_process_computed_property_transitive_dependencies(self):
        """Processor auto-resolves transitive computed dependencies."""
        mat = _make_material(density=sp.Float(2700.0))
        properties = {
            'density': 2700.0,
            'specific_volume': {
                'equation': '1 / density',
                'dependency': [300, 400, 500],
            },
            'normalized_volume': {
                'equation': 'specific_volume * 1000',
                'dependency': [300, 400, 500],
            },
        }
        processed = {'density'}
        processor = ComputedPropertyProcessor(properties, processed)
        processor.process_computed_property(mat, 'normalized_volume', sp.Symbol('T'))
        assert hasattr(mat, 'specific_volume')
        assert hasattr(mat, 'normalized_volume')
        assert 'specific_volume' in processed
        assert 'normalized_volume' in processed

    def test_parse_and_process_expression_returns_sympy_expr(self):
        mat = _make_material(
            density=sp.Float(2700.0),
            heat_capacity=sp.Float(900.0),
        )
        properties = {
            'density': 2700.0,
            'heat_capacity': 900.0,
            'thermal_mass': {
                'equation': 'density * heat_capacity',
                'dependency': [300, 400, 500],
            },
        }
        processor = ComputedPropertyProcessor(
            properties, {'density', 'heat_capacity'})
        result = processor._parse_and_process_expression("density * heat_capacity", mat, sp.Symbol('T'), 'thermal_mass')
        assert isinstance(result, sp.Expr)

    def test_equation_referencing_scalar_material_property(self):
        """Equation referencing a scalar already on the material resolves correctly."""
        mat = _make_material(melting_temperature=933.47, boiling_temperature=2792.0)
        properties = {
            'test_property': {
                'equation': 'melting_temperature * 2',
                'dependency': [300, 400, 500],
            },
        }
        processed = set()
        processor = ComputedPropertyProcessor(properties, processed)
        processor.process_computed_property(mat, 'test_property', sp.Symbol('T'))
        assert hasattr(mat, 'test_property')
        assert 'test_property' in processed
