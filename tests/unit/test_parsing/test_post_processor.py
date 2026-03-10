"""Comprehensive tests for post processor functionality with complete isolation."""

import pytest
import numpy as np
import sympy as sp
import sys
import importlib
import logging
from unittest.mock import patch
from materforge.core.materials import Material
from materforge.parsing.processors.post_processor import PropertyPostProcessor
from materforge.parsing.validation.property_type_detector import PropertyType


class TestPropertyPostProcessor:
    """Post-processing functionality tests with complete isolation."""

    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        """Ensure complete test isolation by reloading processor modules after each test."""
        yield
        for module_name in (
            'materforge.parsing.processors.post_processor',
            'materforge.parsing.processors.temperature_resolver',
        ):
            if module_name in sys.modules:
                importlib.reload(sys.modules[module_name])

    @pytest.fixture
    def sample_material(self):
        T = sp.Symbol('T')
        mat = Material(name="TestMaterial")
        mat.melting_temperature = 1811.0
        mat.boiling_temperature = 3134.0
        mat.heat_capacity = 450 + 0.1 * T
        mat.density = sp.Float(7850)
        return mat

    @pytest.fixture
    def post_processor(self):
        return PropertyPostProcessor()

    def test_post_process_properties_numeric_temperature(self, post_processor, sample_material):
        """Numeric temperature causes the post-processor to skip silently."""
        post_processor.post_process_properties(
            sample_material, 500.0, {}, {PropertyType.CONSTANT_VALUE: []}, set()
        )

    def test_post_process_properties_no_regression(self, post_processor, sample_material):
        """Properties without a regression block are skipped without error."""
        T = sp.Symbol('T')
        config = {
            'dependency': [300, 400, 500],
            'value': [450, 460, 470],
            'bounds': ['constant', 'constant'],
        }
        post_processor.post_process_properties(
            sample_material, T,
            {'heat_capacity': config},
            {PropertyType.TABULAR_DATA: [('heat_capacity', config)]},
            {'heat_capacity'},
        )

    def test_post_process_properties_with_post_regression(self, post_processor, sample_material):
        """Post-regression is applied and the property remains on the material."""
        T = sp.Symbol('T')
        config = {
            'dependency': [300, 400, 500, 600],
            'value': [450, 460, 470, 480],
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'post', 'degree': 1, 'segments': 2},
        }
        sample_material.heat_capacity = 450 + 0.1 * T
        post_processor.post_process_properties(
            sample_material, T,
            {'heat_capacity': config},
            {PropertyType.TABULAR_DATA: [('heat_capacity', config)]},
            {'heat_capacity'},
        )
        assert hasattr(sample_material, 'heat_capacity')

    def test_post_process_properties_missing_property(self, post_processor, sample_material):
        """Property absent from processed_properties set is skipped gracefully."""
        T = sp.Symbol('T')
        config = {
            'dependency': [300, 400, 500],
            'value': [100, 110, 120],
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'post', 'degree': 1, 'segments': 2},
        }
        post_processor.post_process_properties(
            sample_material, T,
            {'missing_property': config},
            {PropertyType.TABULAR_DATA: [('missing_property', config)]},
            set(),  # not in processed set -> skipped
        )

    def test_post_process_properties_integral_property(self, post_processor, sample_material):
        """Integral properties are skipped during post-processing."""
        T = sp.Symbol('T')
        sample_material.integral_prop = sp.Integral(T, T)
        config = {
            'dependency': [300, 400, 500],
            'value': [100, 110, 120],
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'post', 'degree': 1, 'segments': 2},
        }
        post_processor.post_process_properties(
            sample_material, T,
            {'integral_prop': config},
            {PropertyType.COMPUTED_PROPERTY: [('integral_prop', config)]},
            {'integral_prop'},
        )

    def test_post_process_unprocessed_properties_warning(self, post_processor, sample_material, caplog):
        """Warning logged when some properties were never processed."""
        T = sp.Symbol('T')
        props = {
            'prop1': {'dependency': [300, 400], 'value': [100, 110]},
            'prop2': {'dependency': [300, 400], 'value': [200, 210]},
        }
        with caplog.at_level(logging.WARNING):
            post_processor.post_process_properties(
                sample_material, T,
                props,
                {PropertyType.TABULAR_DATA: [('prop1', props['prop1']), ('prop2', props['prop2'])]},
                {'prop1'},  # prop2 intentionally unprocessed
            )
        assert "Some properties were not processed" in caplog.text
        assert "prop2" in caplog.text

    def test_post_process_error_handling(self, post_processor, sample_material):
        """ValueError raised with summary when post-processing a broken property."""
        T = sp.Symbol('T')
        sample_material.error_prop = sp.Symbol('undefined_symbol')
        config = {
            'dependency': [300, 400, 500],
            'value': [100, 110, 120],
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'post', 'degree': 1, 'segments': 2},
        }
        with pytest.raises(ValueError, match="Post-processing errors occurred"):
            post_processor.post_process_properties(
                sample_material, T,
                {'error_prop': config},
                {PropertyType.COMPUTED_PROPERTY: [('error_prop', config)]},
                {'error_prop'},
            )

    def test_apply_post_regression_invalid_dependency_string(self, post_processor, sample_material):
        """Unresolvable string dependency raises ValueError."""
        T = sp.Symbol('T')
        sample_material.test_prop = 450 + 0.1 * T
        config = {
            'dependency': "invalid_string",
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'pre', 'degree': 1, 'segments': 2},
        }
        with pytest.raises(ValueError, match="Failed to extract dependency array."):
            post_processor._apply_post_regression(sample_material, 'test_prop', config, T)

    def test_apply_post_regression_non_numeric_list(self, post_processor, sample_material):
        """Non-numeric list values raise ValueError."""
        T = sp.Symbol('T')
        sample_material.test_prop = 450 + 0.1 * T
        config = {
            'dependency': ['not_a_number', 'also_not_a_number'],
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'pre', 'degree': 1, 'segments': 2},
        }
        with pytest.raises(ValueError):
            post_processor._apply_post_regression(sample_material, 'test_prop', config, T)

    def test_apply_post_regression_invalid_config_type(self, post_processor, sample_material):
        """Dict-typed dependency raises ValueError."""
        T = sp.Symbol('T')
        sample_material.test_prop = 450 + 0.1 * T
        config = {
            'dependency': {'invalid': 'config'},
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'pre', 'degree': 1, 'segments': 2},
        }
        with pytest.raises(ValueError):
            post_processor._apply_post_regression(sample_material, 'test_prop', config, T)

    def test_apply_post_regression_string_dtype_array(self, post_processor, sample_material):
        """String-dtype arrays are coerced to float without error."""
        T = sp.Symbol('T')
        sample_material.test_prop = 450 + 0.1 * T
        config = {
            'dependency': [300, 400, 500],
            'bounds': ['constant', 'constant'],
            'regression': {'simplify': 'pre', 'degree': 1, 'segments': 2},
        }
        with patch('materforge.parsing.processors.post_processor.DependencyResolver.extract_from_config') as mock_extract:
            mock_extract.return_value = np.array(['300', '400', '500'], dtype='U10')
            post_processor._apply_post_regression(
                sample_material, 'test_prop', config, T)
