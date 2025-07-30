"""Error handling tests for visualization module."""

import pytest
import tempfile
import sympy as sp
from pathlib import Path
from unittest.mock import Mock, patch
from materforge.visualization.plotters import PropertyVisualizer
from materforge.core.materials import Material
from materforge.core.elements import ChemicalElement


class TestPropertyVisualizerErrorHandling:
    """Test error handling in property visualization."""

    @pytest.fixture
    def mock_parser(self):
        """Create a mock parser for testing."""
        parser = Mock()
        parser.base_dir = Path(tempfile.gettempdir())
        parser.config_path = "test.yaml"
        parser.categorized_properties = {
            'CONSTANT': [('density', 7850.0)],
            'KEY_VAL': [],
            'FILE': [],
            'STEP_FUNCTION': [],
            'PIECEWISE_EQUATION': [],
            'COMPUTE': []
        }
        parser.config = {
            'name': 'TestMaterial',
            'material_type': 'pure_metal'
        }
        return parser

    @pytest.fixture
    def sample_material(self):
        """Create a sample material for testing."""
        element = ChemicalElement(
            name="Iron", atomic_number=26, atomic_mass=55.845,
            melting_temperature=1811, boiling_temperature=3134,
            latent_heat_of_fusion=13800, latent_heat_of_vaporization=340000
        )

        material = Material(
            name="TestMaterial",
            material_type="pure_metal",
            elements=[element],
            composition=[1.0],
            melting_temperature=sp.Float(1811),
            boiling_temperature=sp.Float(3134)
        )
        material.density = sp.Float(7850)
        return material

    def test_visualizer_no_figure_available(self, mock_parser, sample_material):
        """Test visualization when no figure is available."""
        visualizer = PropertyVisualizer(mock_parser)
        # Don't initialize plots
        T = sp.Symbol('T')
        # Should not raise error, just log warning
        visualizer.visualize_property(
            material=sample_material,
            prop_name='density',
            T=T,
            prop_type='CONSTANT'
        )
        # Property should not be added to visualized set
        assert 'density' not in visualizer.visualized_properties

    def test_visualizer_numeric_temperature_skip(self, mock_parser, sample_material):
        """Test that visualization is skipped for numeric temperature."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        # Should skip visualization for numeric temperature
        visualizer.visualize_property(
            material=sample_material,
            prop_name='density',
            T=500.0,  # Numeric temperature
            prop_type='CONSTANT'
        )
        assert 'density' not in visualizer.visualized_properties

    def test_visualizer_disabled(self, mock_parser, sample_material):
        """Test visualization when disabled."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.is_enabled = False
        # Should not initialize plots
        visualizer.initialize_plots()
        assert visualizer.fig is None

    def test_visualizer_duplicate_property(self, mock_parser, sample_material):
        """Test handling of duplicate property visualization."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        visualizer.visualized_properties.add('density')
        T = sp.Symbol('T')
        # Should skip already visualized property
        visualizer.visualize_property(
            material=sample_material,
            prop_name='density',
            T=T,
            prop_type='CONSTANT'
        )
        # Should still be in set (not added twice)
        assert 'density' in visualizer.visualized_properties

    @patch('matplotlib.pyplot.tight_layout')
    def test_save_plots_layout_error(self, mock_tight_layout, mock_parser):
        """Test handling of layout errors during plot saving."""
        mock_tight_layout.side_effect = Exception("Layout error")
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        # Should handle layout error gracefully
        visualizer.save_property_plots()

    def test_visualizer_invalid_property_evaluation(self, mock_parser, sample_material):
        """Test handling of property evaluation errors."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        # Create a material with an invalid property
        material = sample_material
        material.invalid_prop = sp.Symbol('undefined_symbol') / 0  # Division by zero
        T = sp.Symbol('T')
        # Should handle evaluation error gracefully
        visualizer.visualize_property(
            material=material,
            prop_name='invalid_prop',
            T=T,
            prop_type='COMPUTE'
        )

    def test_visualizer_missing_property(self, mock_parser, sample_material):
        """Test handling of missing properties."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        T = sp.Symbol('T')
        # Should handle missing property gracefully with descriptive error
        with pytest.raises(ValueError, match="Unexpected error in property.*has no attribute"):
            visualizer.visualize_property(
                material=sample_material,
                prop_name='nonexistent_property',
                T=T,
                prop_type='CONSTANT'
            )

    def test_reset_visualization_tracking(self, mock_parser):
        """Test resetting visualization tracking."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.visualized_properties.add('test_prop')
        assert len(visualizer.visualized_properties) == 1
        visualizer.reset_visualization_tracking()
        assert len(visualizer.visualized_properties) == 0
