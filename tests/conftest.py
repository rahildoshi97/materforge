"""Shared pytest fixtures for MaterForge tests."""
import matplotlib
matplotlib.use('Agg')

import pytest
import numpy as np
import sympy as sp
from pathlib import Path

from materforge.core.materials import Material

@pytest.fixture
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials"

@pytest.fixture
def aluminum_yaml_path():
    """Path to aluminum YAML file."""
    current_file = Path(__file__)
    return current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "Al.yaml"

@pytest.fixture
def steel_yaml_path():
    """Path to steel YAML file."""
    return (Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials" / "1.4301.yaml")

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
def sample_valid_material():
    """Minimal pure-metal material with scalar temperature constants."""
    mat = Material(name="Test Aluminum")
    mat.melting_temperature = 933.47
    mat.boiling_temperature = 2792.0
    return mat

@pytest.fixture
def sample_valid_alloy():
    """Minimal alloy material with scalar temperature constants."""
    mat = Material(name="Test Steel")
    mat.solidus_temperature = 1400.0
    mat.liquidus_temperature = 1450.0
    mat.initial_boiling_temperature = 2800.0
    mat.final_boiling_temperature = 2900.0
    return mat
