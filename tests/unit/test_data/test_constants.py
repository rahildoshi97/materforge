"""Unit tests for constants modules."""

import pytest
from materforge.data.constants.processing_constants import ProcessingConstants
from materforge.data.constants.physical_constants import PhysicalConstants

class TestProcessingConstants:
    """Test cases for ProcessingConstants."""
    def test_processing_constants_exist(self):
        """Test that processing constants are defined."""
        assert hasattr(ProcessingConstants, 'STEP_FUNCTION_OFFSET')
        assert hasattr(ProcessingConstants, 'DEFAULT_TOLERANCE')

    def test_step_function_offset_value(self):
        """Test step function offset value."""
        assert ProcessingConstants.STEP_FUNCTION_OFFSET > 0
        assert isinstance(ProcessingConstants.STEP_FUNCTION_OFFSET, (int, float))

    def test_default_tolerance_value(self):
        """Test default tolerance value."""
        assert ProcessingConstants.DEFAULT_TOLERANCE > 0
        assert ProcessingConstants.DEFAULT_TOLERANCE < 1.0
        assert isinstance(ProcessingConstants.DEFAULT_TOLERANCE, (int, float))

    def test_temperature_epsilon(self):
        """Test temperature epsilon constant."""
        if hasattr(ProcessingConstants, 'TEMPERATURE_EPSILON'):
            assert ProcessingConstants.TEMPERATURE_EPSILON > 0
            assert ProcessingConstants.TEMPERATURE_EPSILON < 1.0

    def test_min_data_points(self):
        """Test minimum data points constant."""
        if hasattr(ProcessingConstants, 'MIN_DATA_POINTS'):
            assert ProcessingConstants.MIN_DATA_POINTS > 0
            assert isinstance(ProcessingConstants.MIN_DATA_POINTS, int)

    def test_regression_constants(self):
        """Test regression-related constants."""
        if hasattr(ProcessingConstants, 'DEFAULT_REGRESSION_SEED'):
            assert isinstance(ProcessingConstants.DEFAULT_REGRESSION_SEED, int)
        if hasattr(ProcessingConstants, 'MAX_REGRESSION_SEGMENTS'):
            assert ProcessingConstants.MAX_REGRESSION_SEGMENTS > 0
            assert isinstance(ProcessingConstants.MAX_REGRESSION_SEGMENTS, int)

class TestPhysicalConstants:
    """Test cases for PhysicalConstants."""

    def test_physical_constants_exist(self):
        """Test that physical constants are defined."""
        assert hasattr(PhysicalConstants, 'ABSOLUTE_ZERO')
        assert hasattr(PhysicalConstants, 'BOLTZMANN_CONSTANT')
        assert hasattr(PhysicalConstants, 'AVOGADRO_NUMBER')
        assert hasattr(PhysicalConstants, 'GAS_CONSTANT')

    def test_absolute_zero_value(self):
        """Test absolute zero constant value."""
        assert PhysicalConstants.ABSOLUTE_ZERO == 0.0
        assert isinstance(PhysicalConstants.ABSOLUTE_ZERO, (int, float))

    def test_boltzmann_constant_value(self):
        """Test Boltzmann constant value."""
        # Boltzmann constant in J/K
        expected_kb = 1.380649e-23
        assert abs(PhysicalConstants.BOLTZMANN_CONSTANT - expected_kb) < 1e-28

    def test_avogadro_number_value(self):
        """Test Avogadro number value."""
        # Avogadro number in mol^-1
        expected_na = 6.02214076e23
        assert abs(PhysicalConstants.AVOGADRO_NUMBER - expected_na) < 1e18

    def test_gas_constant_value(self):
        """Test gas constant value."""
        # Gas constant in J/(mol·K)
        expected_r = 8.314462618
        assert abs(PhysicalConstants.GAS_CONSTANT - expected_r) < 1e-9

    def test_constants_are_numeric(self):
        """Test that all constants are numeric values."""
        assert isinstance(PhysicalConstants.ABSOLUTE_ZERO, (int, float))
        assert isinstance(PhysicalConstants.BOLTZMANN_CONSTANT, (int, float))
        assert isinstance(PhysicalConstants.AVOGADRO_NUMBER, (int, float))
        assert isinstance(PhysicalConstants.GAS_CONSTANT, (int, float))
        assert isinstance(ProcessingConstants.STEP_FUNCTION_OFFSET, (int, float))
        assert isinstance(ProcessingConstants.DEFAULT_TOLERANCE, (int, float))

    def test_physical_constants_relationships(self):
        """Test relationships between physical constants."""
        # Gas constant should equal Boltzmann constant times Avogadro number
        calculated_r = PhysicalConstants.BOLTZMANN_CONSTANT * PhysicalConstants.AVOGADRO_NUMBER
        assert abs(calculated_r - PhysicalConstants.GAS_CONSTANT) < 1e-9
