"""Unit tests for PropertyTypeDetector."""

import pytest
from pymatlib.parsing.validation.property_type_detector import PropertyType, PropertyTypeDetector

class TestPropertyTypeDetector:
    """Test cases for PropertyTypeDetector."""
    @pytest.mark.parametrize("config,expected_type", [
        (5.0, PropertyType.CONSTANT_VALUE),
        ("2.5", PropertyType.CONSTANT_VALUE),
        ({"file_path": "data.csv", "dependency_header": "T", "value_header": "rho", "bounds": ["constant", "constant"]}, PropertyType.FILE_IMPORT),
        ({"dependency": "melting_temperature", "value": [900, 1000]}, PropertyType.STEP_FUNCTION),
        ({"dependency": [300, 400, 500], "value": [900, 950, 1000]}, PropertyType.TABULAR_DATA),
        ({"dependency": [300, 400, 500], "equation": ["2*T + 100", "3*T - 50"], "bounds": ["constant", "constant"]}, PropertyType.PIECEWISE_EQUATION),
        ({"dependency": [300, 400, 500], "equation": "density * heat_capacity"}, PropertyType.COMPUTED_PROPERTY),
    ])
    def test_determine_property_type(self, config, expected_type):
        """Test property type detection for various configurations."""
        result = PropertyTypeDetector.determine_property_type("test_prop", config)
        assert result == expected_type

    def test_determine_property_type_invalid_integer(self):
        """Test that integer constants raise appropriate error."""
        with pytest.raises(ValueError, match="must be defined as a float"):
            PropertyTypeDetector.determine_property_type("test_prop", 5)

    def test_determine_property_type_unknown_pattern(self):
        """Test unknown configuration pattern."""
        config = {"unknown_key": "unknown_value"}
        with pytest.raises(ValueError, match="doesn't match any known configuration pattern"):
            PropertyTypeDetector.determine_property_type("test_prop", config)

    def test_validate_constant_property_valid(self):
        """Test validation of valid constant property."""
        # Should not raise any exception
        PropertyTypeDetector.validate_property_config("test_prop", 5.0, PropertyType.CONSTANT_VALUE)

    def test_validate_constant_property_invalid(self):
        """Test validation of invalid constant property."""
        with pytest.raises(ValueError, match="could not be converted to a float"):
            PropertyTypeDetector.validate_property_config("test_prop", "invalid", PropertyType.CONSTANT_VALUE)

    def test_validate_file_property_missing_keys(self):
        """Test validation of file property with missing keys."""
        config = {"file_path": "data.csv"}  # Missing required keys
        with pytest.raises(ValueError, match="Invalid configuration"):
            PropertyTypeDetector.validate_property_config("test_prop", config, PropertyType.FILE_IMPORT)

    def test_validate_step_function_invalid_values(self):
        """Test validation of step function with invalid values."""
        config = {
            "dependency": "melting_temperature",
            "value": [900]  # Should have exactly 2 values
        }
        with pytest.raises(ValueError, match="must be a list of exactly two numbers"):
            PropertyTypeDetector.validate_property_config("test_prop", config, PropertyType.STEP_FUNCTION)
