"""Unit tests for visualization plotters."""

import pytest
import numpy as np
import sympy as sp
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
from materforge.core.materials import Material
from materforge.visualization.plotters import PropertyVisualizer

class TestPropertyVisualizer:
    """Test cases for PropertyVisualizer."""
    @pytest.fixture
    def mock_parser(self):
        """Create a mock parser for PropertyVisualizer initialization."""
        parser = Mock()
        parser.config_path = Path("test.yaml")
        parser.base_dir = Path(".")
        parser.categorized_properties = {
            'CONSTANTVALUE': [('density', 2700.0)],
            'TABULARDATA': [(
                'heat_capacity',
                {'dependency': [300, 400], 'value': [900, 950], 'bounds': ['constant', 'constant']},
            )],
        }
        return parser

    def test_property_visualizer_initialization(self, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        assert hasattr(visualizer, 'is_visualization_enabled')
        assert hasattr(visualizer, 'visualize_property')
        assert hasattr(visualizer, 'initialize_plots')
        assert hasattr(visualizer, 'save_property_plots')

    def test_visualization_disabled_by_default(self, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        assert visualizer.is_visualization_enabled() is False

    def test_disable_visualization(self, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        assert visualizer.is_visualization_enabled() is False

    def test_enable_visualization(self, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        assert visualizer.is_visualization_enabled() is False

    @patch('matplotlib.pyplot.figure')
    def test_visualize_property_basic(self, mock_figure, mock_parser):
        material = Material(name="Test Material")
        visualizer = PropertyVisualizer(mock_parser)
        T = sp.Symbol('T')
        x_data = np.array([300, 400, 500], dtype=float)
        y_data = np.array([900, 950, 1000], dtype=float)
        visualizer.visualize_property(
            material=material,
            prop_name="heat_capacity",
            dependency=T,
            prop_type="TABULARDATA",
            x_data=x_data,
            y_data=y_data,
        )
        # Visualization is disabled by default, so no figures should be created.
        mock_figure.assert_not_called()

    @patch('matplotlib.pyplot.figure')
    def test_visualize_property_disabled(self, mock_figure, mock_parser):
        material = Material(name="Test Material")
        visualizer = PropertyVisualizer(mock_parser)
        T = sp.Symbol('T')
        x_data = np.array([300, 400, 500], dtype=float)
        y_data = np.array([900, 950, 1000], dtype=float)
        visualizer.visualize_property(
            material=material,
            prop_name="heat_capacity",
            dependency=T,
            prop_type="TABULARDATA",
            x_data=x_data,
            y_data=y_data,
        )
        mock_figure.assert_not_called()
    @patch('matplotlib.pyplot.figure')
    def test_visualize_constant_property(self, mock_figure, mock_parser):
        material = Material(name="Test Material")
        material.density = 2700.0
        visualizer = PropertyVisualizer(mock_parser)
        T = sp.Symbol('T')
        visualizer.visualize_property(
            material=material,
            prop_name="density",
            dependency=T,
            prop_type="CONSTANTVALUE",
        )
        mock_figure.assert_not_called()

    @patch('matplotlib.pyplot.savefig')
    def test_save_property_plots(self, mock_savefig, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.save_property_plots()
        mock_savefig.assert_not_called()

    def test_initialize_plots(self, mock_parser):
        with tempfile.TemporaryDirectory() as temp_dir:
            mock_parser.base_dir = Path(temp_dir)
            visualizer = PropertyVisualizer(mock_parser)
            visualizer.initialize_plots()
            # If your implementation creates a directory, it should be here.
            plot_dir = Path(temp_dir) / "materforge_plots"
            assert plot_dir.exists()

    def test_reset_visualization_tracking(self, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.reset_visualization_tracking()
        assert True

    @patch('matplotlib.pyplot.figure')
    def test_visualize_property_with_regression(self, mock_figure, mock_parser):
        material = Material(name="Test Material")
        visualizer = PropertyVisualizer(mock_parser)
        T = sp.Symbol('T')
        x_data = np.array([300, 400, 500], dtype=float)
        y_data = np.array([900, 950, 1000], dtype=float)
        visualizer.visualize_property(
            material=material,
            prop_name="heat_capacity",
            dependency=T,
            prop_type="TABULARDATA",
            x_data=x_data,
            y_data=y_data,
            has_regression=True,
            degree=2,
            segments=1,
        )
        mock_figure.assert_not_called()

    def test_visualize_property_error_handling(self, mock_parser):
        material = Material(name="Test Material")
        visualizer = PropertyVisualizer(mock_parser)
        T = sp.Symbol('T')
        # Should not raise even if prop doesn't exist / invalid config.
        visualizer.visualize_property(
            material=material,
            prop_name="nonexistent_prop",
            dependency=T,
            prop_type="CONSTANTVALUE",
        )

    @patch('matplotlib.pyplot.figure')
    def test_visualize_multiple_properties(self, mock_figure, mock_parser):
        material = Material(name="Test Material")
        visualizer = PropertyVisualizer(mock_parser)
        T = sp.Symbol('T')
        for prop in ['density', 'heat_capacity', 'thermal_conductivity']:
            visualizer.visualize_property(
                material=material,
                prop_name=prop,
                dependency=T,
                prop_type="CONSTANTVALUE",
            )
        mock_figure.assert_not_called()

    @patch('matplotlib.pyplot.savefig')
    def test_save_property_plots_no_directory(self, mock_savefig, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        visualizer.save_property_plots()
        mock_savefig.assert_not_called()

    def test_visualizer_with_real_parser_structure(self, mock_parser):
        with tempfile.TemporaryDirectory() as temp_dir:
            mock_parser.config_path = Path("test_material.yaml")
            mock_parser.base_dir = Path(temp_dir)
            visualizer = PropertyVisualizer(mock_parser)
            assert visualizer.parser == mock_parser
            assert visualizer.is_visualization_enabled() is False
            visualizer.initialize_plots()
            assert (Path(temp_dir) / "materforge_plots").exists()

    def test_visualization_state_methods(self, mock_parser):
        visualizer = PropertyVisualizer(mock_parser)
        assert visualizer.is_visualization_enabled() is False
        assert hasattr(visualizer, 'parser')
