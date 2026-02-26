"""Error handling tests for visualization module."""

import pytest
import tempfile
import sympy as sp
from pathlib import Path
from unittest.mock import Mock, patch
from materforge.core.materials import Material
from materforge.visualization.plotters import PropertyVisualizer


class TestPropertyVisualizerErrorHandling:
    """Error handling tests for PropertyVisualizer."""

    @pytest.fixture
    def mock_parser(self):
        parser = Mock()
        parser.base_dir = Path(tempfile.gettempdir())
        parser.config_path = "test.yaml"
        parser.categorized_properties = {
            'CONSTANT': [('density', 7850.0)],
            'KEY_VAL': [],
            'FILE': [],
            'STEP_FUNCTION': [],
            'PIECEWISE_EQUATION': [],
            'COMPUTE': [],
        }
        parser.config = {'name': 'TestMaterial'}
        return parser

    @pytest.fixture
    def sample_material(self):
        mat = Material(name="TestMaterial")
        mat.melting_temperature = 1811.0
        mat.boiling_temperature = 3134.0
        mat.density = sp.Float(7850)
        return mat

    def test_visualizer_no_figure_available(self, mock_parser, sample_material):
        """Property skipped silently when plots are not initialised."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.visualize_property(
            material=sample_material,
            prop_name='density',
            dependency=sp.Symbol('T'),
            prop_type='CONSTANT',
        )
        assert 'density' not in visualizer.visualized_properties

    def test_visualizer_numeric_temperature_skip(self, mock_parser, sample_material):
        """Numeric temperature causes visualizer to skip the property."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        visualizer.visualize_property(
            material=sample_material,
            prop_name='density',
            dependency=500.0,
            prop_type='CONSTANT',
        )
        assert 'density' not in visualizer.visualized_properties

    def test_visualizer_disabled(self, mock_parser):
        """Disabled visualizer does not initialise figure."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.is_enabled = False
        visualizer.initialize_plots()
        assert visualizer.fig is None

    def test_visualizer_duplicate_property(self, mock_parser, sample_material):
        """Already-visualised property is not added a second time."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        visualizer.visualized_properties.add('density')
        visualizer.visualize_property(
            material=sample_material,
            prop_name='density',
            dependency=sp.Symbol('T'),
            prop_type='CONSTANT',
        )
        assert 'density' in visualizer.visualized_properties

    @patch('matplotlib.pyplot.tight_layout')
    def test_save_plots_layout_error(self, mock_tight_layout, mock_parser):
        """Layout error during save is handled gracefully."""
        mock_tight_layout.side_effect = Exception("Layout error")
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        visualizer.save_property_plots()  # must not raise

    def test_visualizer_invalid_property_evaluation(self, mock_parser, sample_material):
        """Expression that cannot be evaluated is handled gracefully."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        sample_material.invalid_prop = sp.Symbol('undefined_symbol') / 0
        visualizer.visualize_property(
            material=sample_material,
            prop_name='invalid_prop',
            dependency=sp.Symbol('T'),
            prop_type='COMPUTE',
        )

    def test_visualizer_missing_property(self, mock_parser, sample_material):
        """Non-existent property name raises ValueError."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.initialize_plots()
        with pytest.raises(ValueError, match="Unexpected error visualizing.*has no attribute"):
            visualizer.visualize_property(
                material=sample_material,
                prop_name='nonexistent_property',
                dependency=sp.Symbol('T'),
                prop_type='CONSTANT',
            )

    def test_reset_visualization_tracking(self, mock_parser):
        """Reset clears all tracked properties."""
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.visualized_properties.add('test_prop')
        visualizer.reset_visualization_tracking()
        assert len(visualizer.visualized_properties) == 0
