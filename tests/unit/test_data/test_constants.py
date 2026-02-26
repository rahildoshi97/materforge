"""Unit tests for constants modules."""

import scipy.constants as sc
from materforge.data.constants.processing_constants import ProcessingConstants
from materforge.data.constants.physical_constants import PhysicalConstants

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

class TestPhysicalConstants:
    """Test cases for PhysicalConstants sourced from scipy.constants (CODATA 2022)."""

    def test_physical_constants_exist(self):
        """Test all expected constants are present."""
        expected = [
            'ABSOLUTE_ZERO', 'ROOM_TEMPERATURE', 'TRIPLE_POINT_WATER',
            'AMU', 'N_A', 'AVOGADRO_NUMBER',
            'BOLTZMANN_CONSTANT', 'GAS_CONSTANT',
            'STEFAN_BOLTZMANN_CONSTANT', 'GRAVITY',
        ]
        for attr in expected:
            assert hasattr(PhysicalConstants, attr), f"Missing constant: {attr}"

    def test_constants_are_numeric(self):
        """Test all constants are numeric."""
        attrs = [
            'ABSOLUTE_ZERO', 'ROOM_TEMPERATURE', 'TRIPLE_POINT_WATER',
            'AMU', 'N_A', 'AVOGADRO_NUMBER',
            'BOLTZMANN_CONSTANT', 'GAS_CONSTANT',
            'STEFAN_BOLTZMANN_CONSTANT', 'GRAVITY',
        ]
        for attr in attrs:
            assert isinstance(getattr(PhysicalConstants, attr), (int, float)), \
                f"{attr} is not numeric"

    # -------------------------------------------------------------------------
    # Temperature references (hardcoded conventional values)
    # -------------------------------------------------------------------------

    def test_absolute_zero(self):
        """Absolute zero is exactly 0.0 K."""
        assert PhysicalConstants.ABSOLUTE_ZERO == 0.0

    def test_room_temperature(self):
        """Standard room temperature is 298.15 K."""
        assert PhysicalConstants.ROOM_TEMPERATURE == 298.15

    def test_triple_point_water(self):
        """Triple point of water is 273.16 K."""
        assert PhysicalConstants.TRIPLE_POINT_WATER == 273.16

    # -------------------------------------------------------------------------
    # scipy.constants delegation
    # -------------------------------------------------------------------------

    def test_amu_matches_scipy(self):
        """AMU delegates to scipy.constants.atomic_mass (within numerical tolerance)."""
        assert abs(PhysicalConstants.AMU - sc.atomic_mass) < 1e-30

    def test_avogadro_matches_scipy(self):
        """N_A and AVOGADRO_NUMBER delegate to scipy.constants.Avogadro."""
        assert PhysicalConstants.N_A == sc.Avogadro
        assert PhysicalConstants.AVOGADRO_NUMBER == sc.Avogadro

    def test_boltzmann_matches_scipy(self):
        """BOLTZMANN_CONSTANT delegates to scipy.constants.Boltzmann."""
        assert PhysicalConstants.BOLTZMANN_CONSTANT == sc.Boltzmann

    def test_gas_constant_matches_scipy(self):
        """GAS_CONSTANT delegates to scipy.constants.R (within numerical tolerance)."""
        assert abs(PhysicalConstants.GAS_CONSTANT - sc.R) < 1e-9

    def test_stefan_boltzmann_matches_scipy(self):
        """STEFAN_BOLTZMANN_CONSTANT delegates to scipy.constants.Stefan_Boltzmann."""
        assert abs(PhysicalConstants.STEFAN_BOLTZMANN_CONSTANT - sc.Stefan_Boltzmann) < 1e-15

    def test_GRAVITY_matches_scipy(self):
        """GRAVITY delegates to scipy.constants.g."""
        assert PhysicalConstants.GRAVITY == sc.g

    # -------------------------------------------------------------------------
    # Physical relationships - internal consistency checks
    # -------------------------------------------------------------------------

    def test_gas_constant_is_kb_times_na(self):
        """R = k_B * N_A must hold within floating point precision."""
        calculated_r = PhysicalConstants.BOLTZMANN_CONSTANT * PhysicalConstants.AVOGADRO_NUMBER
        assert abs(calculated_r - PhysicalConstants.GAS_CONSTANT) < 1e-9

    def test_avogadro_alias_consistent(self):
        """N_A and AVOGADRO_NUMBER must be identical."""
        assert PhysicalConstants.N_A == PhysicalConstants.AVOGADRO_NUMBER
