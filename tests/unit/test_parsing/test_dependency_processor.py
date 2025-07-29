"""Unit tests for DependencyProcessor."""

import pytest
import sympy as sp
from pymatlib.parsing.processors.dependency_processor import DependencyProcessor
from pymatlib.parsing.validation.errors import DependencyError, CircularDependencyError

class TestDependencyProcessor:
    """Test cases for DependencyProcessor."""
    def test_dependency_processor_initialization(self):
        """Test dependency processor initialization."""
        properties = {
            'density': 2700.0,
            'heat_capacity': {'temperature': [300, 400], 'value': [900, 950]}
        }
        processed_properties = set()
        processor = DependencyProcessor(properties, processed_properties)
        assert processor.properties == properties
        assert processor.processed_properties == processed_properties

    def test_extract_equation_dependencies_simple(self):
        """Test extracting dependencies from simple equations."""
        equation = "density * heat_capacity"
        dependencies = DependencyProcessor._extract_equation_dependencies(equation)
        expected = ['density', 'heat_capacity']
        assert set(dependencies) == set(expected)

    def test_extract_equation_dependencies_complex(self):
        """Test extracting dependencies from complex equations."""
        equation = "density * heat_capacity + thermal_conductivity / viscosity"
        dependencies = DependencyProcessor._extract_equation_dependencies(equation)
        expected = ['density', 'heat_capacity', 'thermal_conductivity', 'viscosity']
        assert set(dependencies) == set(expected)

    def test_extract_equation_dependencies_with_temperature(self):
        """Test that temperature symbol is excluded from dependencies."""
        equation = "density * T + heat_capacity"
        dependencies = DependencyProcessor._extract_equation_dependencies(equation)
        # T should be excluded, only other symbols included
        expected = ['density', 'heat_capacity']
        assert set(dependencies) == set(expected)

    def test_validate_circular_dependencies_no_cycle(self):
        """Test validation with no circular dependencies."""
        properties = {
            'density': 2700.0,
            'volume': {'equation': 'mass / density'},
            'mass': 1000.0
        }
        processed_properties = set()
        processor = DependencyProcessor(properties, processed_properties)
        # Should not raise any exception
        processor._validate_circular_dependencies('volume', ['mass', 'density'], set())

    def test_validate_circular_dependencies_with_cycle(self):
        """Test validation with circular dependencies."""
        properties = {
            'prop_a': {'equation': 'prop_b + 100'},
            'prop_b': {'equation': 'prop_c * 2'},
            'prop_c': {'equation': 'prop_a / 3'}  # Creates cycle: a -> b -> c -> a
        }
        processed_properties = set()
        processor = DependencyProcessor(properties, processed_properties)
        with pytest.raises(CircularDependencyError):
            processor._validate_circular_dependencies('prop_a', ['prop_b'], set())

    def test_validate_circular_dependencies_self_reference(self):
        """Test validation with self-referencing property."""
        properties = {
            'recursive_prop': {'equation': 'recursive_prop + 1'}
        }
        processed_properties = set()
        processor = DependencyProcessor(properties, processed_properties)
        with pytest.raises(CircularDependencyError):
            processor._validate_circular_dependencies('recursive_prop', ['recursive_prop'], set())

    def test_process_computed_property_simple(self, sample_aluminum_element):
        """Test processing simple computed property."""
        from pymatlib.core.materials import Material

        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        # Set base properties
        material.density = sp.Float(2700.0)
        material.volume = sp.Float(0.001)  # 1 liter
        properties = {
            'density': 2700.0,
            'volume': 0.001,
            'mass': {
                'equation': 'density * volume',
                'dependency': [300, 400, 500]
            }
        }
        processed_properties = {'density', 'volume'}
        processor = DependencyProcessor(properties, processed_properties)
        T = sp.Symbol('T')
        processor.process_computed_property(material, 'mass', T)
        assert hasattr(material, 'mass')
        assert 'mass' in processed_properties

    def test_process_computed_property_missing_dependency(self, sample_aluminum_element):
        """Test processing computed property with missing dependency."""
        from pymatlib.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        properties = {
            'mass': {
                'equation': 'density * volume',  # density and volume not defined
                'dependency': [300, 400, 500]
            }
        }
        processed_properties = set()
        processor = DependencyProcessor(properties, processed_properties)
        T = sp.Symbol('T')
        with pytest.raises(ValueError, match="Missing dependencies in expression"):
            processor.process_computed_property(material, 'mass', T)

    def test_process_computed_property_with_dependencies(self, sample_aluminum_element):
        """Test processing computed property that depends on other computed properties."""
        from pymatlib.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        # Set base property
        material.density = sp.Float(2700.0)
        properties = {
            'density': 2700.0,
            'specific_volume': {
                'equation': '1 / density',
                'dependency': [300, 400, 500]
            },
            'normalized_volume': {
                'equation': 'specific_volume * 1000',
                'dependency': [300, 400, 500]
            }
        }
        processed_properties = {'density'}
        processor = DependencyProcessor(properties, processed_properties)
        T = sp.Symbol('T')
        # Process the dependent property
        processor.process_computed_property(material, 'normalized_volume', T)
        assert hasattr(material, 'specific_volume')
        assert hasattr(material, 'normalized_volume')
        assert 'specific_volume' in processed_properties
        assert 'normalized_volume' in processed_properties

    def test_parse_and_process_expression_basic(self, sample_aluminum_element):
        """Test basic expression parsing and processing."""
        from pymatlib.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        # Set required properties
        material.density = sp.Float(2700.0)
        material.heat_capacity = sp.Float(900.0)
        properties = {
            'density': 2700.0,
            'heat_capacity': 900.0,
            'thermal_mass': {
                'equation': 'density * heat_capacity',
                'dependency': [300, 400, 500]
            }
        }
        processed_properties = {'density', 'heat_capacity'}
        processor = DependencyProcessor(properties, processed_properties)
        T = sp.Symbol('T')
        # Test the expression parsing method directly
        expression = "density * heat_capacity"
        result = processor._parse_and_process_expression(expression, material, T, 'thermal_mass')
        assert isinstance(result, sp.Expr)

    def test_get_temperature_value_from_material(self, sample_aluminum_element):
        """Test getting temperature values from material references."""
        from pymatlib.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        properties = {
            'test_property': {
                'equation': 'melting_temperature * 2',
                'dependency': [300, 400, 500]
            }
        }
        processed_properties = set()
        processor = DependencyProcessor(properties, processed_properties)
        T = sp.Symbol('T')
        # This should work since melting_temperature is available on the material
        processor.process_computed_property(material, 'test_property', T)
        assert hasattr(material, 'test_property')
        assert 'test_property' in processed_properties
