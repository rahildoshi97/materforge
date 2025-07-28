"""Shared pytest fixtures for PyMatLib tests."""
import pytest
import numpy as np
import sympy as sp
from pathlib import Path

from pymatlib.core.materials import Material
from pymatlib.core.elements import ChemicalElement
from pymatlib.data.elements.element_data import element_map

@pytest.fixture
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent.parent / "src" / "pymatlib" / "data" / "materials"

@pytest.fixture
def aluminum_yaml_path():
    """Path to aluminum YAML file."""
    current_file = Path(__file__)
    return current_file.parent.parent / "src" / "pymatlib" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"

@pytest.fixture
def steel_yaml_path():
    """Path to steel YAML file."""
    current_file = Path(__file__)
    return current_file.parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

@pytest.fixture
def sample_aluminum_element():
    """Sample aluminum element for testing."""
    return element_map['Al']

@pytest.fixture
def sample_steel_elements():
    """Sample steel alloy elements."""
    return [element_map['Fe'], element_map['C'], element_map['Cr'], element_map['Ni']]

@pytest.fixture
def temp_symbol():
    """Temperature symbol for testing."""
    return sp.Symbol('T')

@pytest.fixture
def sample_temperature_array():
    """Sample temperature array."""
    return np.linspace(300, 1000, 50)

@pytest.fixture
def sample_property_array():
    """Sample property array."""
    return np.linspace(900, 1200, 50)

@pytest.fixture
def mock_visualizer():
    """Mock visualizer for testing."""
    class MockVisualizer:
        def __init__(self):
            self.visualization_enabled = True
            self.visualized_properties = []

        def is_visualization_enabled(self):
            return self.visualization_enabled

        def visualize_property(self, **kwargs):
            self.visualized_properties.append(kwargs)

        def initialize_plots(self):
            pass

        def save_property_plots(self):
            pass

        def reset_visualization_tracking(self):
            self.visualized_properties = []

    return MockVisualizer()

@pytest.fixture
def sample_valid_material(sample_aluminum_element):
    """Sample valid material for testing."""
    return Material(
        name="Test Aluminum",
        material_type="pure_metal",
        elements=[sample_aluminum_element],
        composition=[1.0],
        melting_temperature=sp.Float(933.47),
        boiling_temperature=sp.Float(2792.0)
    )

@pytest.fixture
def sample_valid_alloy(sample_steel_elements):
    """Sample valid alloy for testing."""
    return Material(
        name="Test Steel",
        material_type="alloy",
        elements=sample_steel_elements,
        composition=[0.68, 0.02, 0.20, 0.10],
        solidus_temperature=sp.Float(1400.0),
        liquidus_temperature=sp.Float(1450.0),
        initial_boiling_temperature=sp.Float(2800.0),
        final_boiling_temperature=sp.Float(2900.0)
    )
