"""Unit tests for constants modules."""

import scipy.constants as sc
from materforge.data.constants.processing_constants import ProcessingConstants


class TestProcessingConstants:
    """Test cases for ProcessingConstants."""
    def test_processing_constants_exist(self):
        """Test that processing constants are defined."""
        assert hasattr(ProcessingConstants, 'STEP_FUNCTION_OFFSET')
        assert hasattr(ProcessingConstants, 'DEFAULT_TOLERANCE')

    def test_step_function_offset_value(self):
        """Test step function offset is positive."""
        assert ProcessingConstants.STEP_FUNCTION_OFFSET > 0
        assert isinstance(ProcessingConstants.STEP_FUNCTION_OFFSET, (int, float))

    def test_default_tolerance_value(self):
        """Test default tolerance is a small positive fraction."""
        assert 0 < ProcessingConstants.DEFAULT_TOLERANCE < 1.0
        assert isinstance(ProcessingConstants.DEFAULT_TOLERANCE, (int, float))

    def test_temperature_epsilon(self):
        """Test temperature epsilon is a small positive fraction."""
        if hasattr(ProcessingConstants, 'DEPENDENCY_EPSILON'):
            assert 0 < ProcessingConstants.DEPENDENCY_EPSILON < 1.0

    def test_min_data_points(self):
        """Test minimum data points is a positive integer."""
        if hasattr(ProcessingConstants, 'MIN_DATA_POINTS'):
            assert ProcessingConstants.MIN_DATA_POINTS > 0
            assert isinstance(ProcessingConstants.MIN_DATA_POINTS, int)

    def test_regression_constants(self):
        """Test regression-related constants have correct types."""
        if hasattr(ProcessingConstants, 'DEFAULT_REGRESSION_SEED'):
            assert isinstance(ProcessingConstants.DEFAULT_REGRESSION_SEED, int)
        if hasattr(ProcessingConstants, 'MAX_REGRESSION_SEGMENTS'):
            assert ProcessingConstants.MAX_REGRESSION_SEGMENTS > 0
            assert isinstance(ProcessingConstants.MAX_REGRESSION_SEGMENTS, int)
