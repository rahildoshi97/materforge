"""Comprehensive tests for property processor base functionality."""
import pytest
import numpy as np
import sympy as sp
from pathlib import Path
from unittest.mock import Mock
from materforge.parsing.processors.property_processor_base import PropertyProcessorBase
from materforge.core.materials import Material


class TestPropertyProcessorBaseComprehensive:
    """Comprehensive tests for property processor base functionality."""

    @pytest.fixture
    def processor(self):
        return PropertyProcessorBase()

    @pytest.fixture
    def sample_material(self):
        mat = Material(name="TestMaterial")
        mat.melting_temperature = 1811.0
        mat.boiling_temperature = 3134.0
        return mat

    def test_initialization(self, processor):
        assert processor.processed_properties == set()
        assert processor.base_dir is None
        assert processor.visualizer is None

    def test_set_processing_context(self, processor):
        base_dir = Path("/test/path")
        visualizer = Mock()
        processed_props = {'prop1', 'prop2'}

        processor.set_processing_context(base_dir, visualizer, processed_props)

        assert processor.base_dir == base_dir
        assert processor.visualizer == visualizer
        assert processor.processed_properties == processed_props

    def test_set_visualizer(self, processor):
        visualizer = Mock()
        processor.set_visualizer(visualizer)
        assert processor.visualizer == visualizer

    def test_finalize_with_data_arrays_symbolic(self, processor, sample_material):
        """Symbolic dependency -> returns False, assigns piecewise to material."""
        T = sp.Symbol('T')
        config = {'bounds': ['constant', 'constant']}

        result = processor.finalize_with_data_arrays(
            material=sample_material,
            prop_name='test_prop',
            dep_array=np.array([300, 400, 500]),
            prop_array=np.array([100, 150, 200]),
            dependency=T,
            config=config,
            prop_type='KEY_VAL',
        )

        assert result is False
        assert hasattr(sample_material, 'test_prop')
        assert 'test_prop' in processor.processed_properties

    def test_finalize_with_data_arrays_numeric_temperature(self, processor, sample_material):
        """Numeric dependency -> returns True, assigns interpolated sp.Float."""
        config = {'bounds': ['constant', 'constant']}

        result = processor.finalize_with_data_arrays(
            material=sample_material,
            prop_name='test_prop',
            dep_array=np.array([300, 400, 500]),
            prop_array=np.array([100, 150, 200]),
            dependency=400.0,
            config=config,
            prop_type='KEY_VAL',
        )

        assert result is True
        assert hasattr(sample_material, 'test_prop')
        assert isinstance(getattr(sample_material, 'test_prop'), sp.Float)

    def test_finalize_with_data_arrays_none_input(self, processor, sample_material):
        """None array raises ValueError."""
        T = sp.Symbol('T')
        config = {'bounds': ['constant', 'constant']}

        with pytest.raises(ValueError, match="cannot be None"):
            processor.finalize_with_data_arrays(
                material=sample_material,
                prop_name='test_prop',
                dep_array=None,
                prop_array=np.array([100, 150, 200]),
                dependency=T,
                config=config,
                prop_type='KEY_VAL',
            )

    def test_finalize_with_data_arrays_mismatched_lengths(self, processor, sample_material):
        """Mismatched array lengths raise ValueError."""
        T = sp.Symbol('T')
        config = {'bounds': ['constant', 'constant']}

        with pytest.raises(ValueError, match="equal length"):
            processor.finalize_with_data_arrays(
                material=sample_material,
                prop_name='test_prop',
                dep_array=np.array([300, 400]),
                prop_array=np.array([100, 150, 200]),
                dependency=T,
                config=config,
                prop_type='KEY_VAL',
            )

    def test_finalize_with_piecewise_function(self, processor, sample_material):
        """Piecewise function assigned verbatim for symbolic dependency."""
        T = sp.Symbol('T')
        piecewise_func = sp.Piecewise((100, T < 400), (200, True))
        config = {'bounds': ['constant', 'constant']}

        result = processor.finalize_with_piecewise_function(
            material=sample_material,
            prop_name='test_prop',
            piecewise_func=piecewise_func,
            dependency=T,
            config=config,
            prop_type='PIECEWISE_EQUATION',
        )

        assert result is False
        assert hasattr(sample_material, 'test_prop')
        assert getattr(sample_material, 'test_prop') == piecewise_func

    def test_reset_processing_state(self, processor):
        processor.processed_properties.add('test_prop')
        processor.reset_processing_state()
        assert len(processor.processed_properties) == 0

    def test_get_processed_properties_returns_copy(self, processor):
        processor.processed_properties.update({'prop1', 'prop2'})
        props = processor.get_processed_properties()
        assert props == {'prop1', 'prop2'}
        props.add('prop3')
        assert 'prop3' not in processor.processed_properties

    def test_is_property_processed(self, processor):
        assert not processor.is_property_processed('test_prop')
        processor.processed_properties.add('test_prop')
        assert processor.is_property_processed('test_prop')


    def test_visualize_if_enabled_no_visualizer(self, processor, sample_material):
        """No error when visualizer is None."""
        processor._visualize_if_enabled(
            material=sample_material,
            prop_name='test_prop',
            dependency=sp.Symbol('T'),
            prop_type='CONSTANT',
        )

    def test_visualize_if_enabled_numeric_temperature(self, processor, sample_material):
        """Visualizer is not called for numeric temperature."""
        processor.visualizer = Mock()
        processor._visualize_if_enabled(
            material=sample_material,
            prop_name='test_prop',
            dependency=400.0,
            prop_type='CONSTANT',
        )
        assert not processor.visualizer.visualize_property.called

    def test_visualize_if_enabled_disabled_visualizer(self, processor, sample_material):
        """Visualizer is not called when disabled."""
        visualizer = Mock()
        visualizer.is_visualization_enabled.return_value = False
        processor.visualizer = visualizer

        processor._visualize_if_enabled(
            material=sample_material,
            prop_name='test_prop',
            dependency=sp.Symbol('T'),
            prop_type='CONSTANT',
        )

        visualizer.is_visualization_enabled.assert_called_once()
        assert not visualizer.visualize_property.called

    def test_visualize_if_enabled_with_regression_config(self, processor, sample_material):
        """Regression config is forwarded correctly to the visualizer."""
        visualizer = Mock()
        visualizer.is_visualization_enabled.return_value = True
        processor.visualizer = visualizer

        config = {
            'bounds': ['constant', 'extrapolate'],
            'regression': {'simplify': 'pre', 'degree': 2, 'segments': 3},
        }

        processor._visualize_if_enabled(
            material=sample_material,
            prop_name='test_prop',
            dependency=sp.Symbol('T'),
            prop_type='KEY_VAL',
            config=config,
            bounds=(300, 500),
        )

        visualizer.visualize_property.assert_called_once()
        call_kwargs = visualizer.visualize_property.call_args[1]
        assert call_kwargs['has_regression'] is True
        assert call_kwargs['simplify_type'] == 'pre'
        assert call_kwargs['degree'] == 2
        assert call_kwargs['segments'] == 3
        assert call_kwargs['lower_bound_type'] == 'constant'
        assert call_kwargs['upper_bound_type'] == 'extrapolate'
