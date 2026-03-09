"""Unit tests for PropertyTypeDetector."""
import pytest
from materforge.parsing.validation.property_type_detector import PropertyType, PropertyTypeDetector

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
        """Property type detection covers all known configuration shapes."""
        result = PropertyTypeDetector.determine_property_type("test_prop", config)
        assert result == expected_type

    def test_determine_property_type_invalid_integer(self):
        """Bare integer constant raises ValueError (must be float)."""
        with pytest.raises(ValueError, match="must be defined as a float"):
            PropertyTypeDetector.determine_property_type("test_prop", 5)

    def test_determine_property_type_unknown_pattern(self):
        """Unrecognised config dict raises ValueError."""
        with pytest.raises(ValueError, match="does not match any known configuration pattern"):
            PropertyTypeDetector.determine_property_type("test_prop", {"unknown_key": "unknown_value"})

    def test_validate_constant_property_valid(self):
        """Valid float constant passes without raising."""
        PropertyTypeDetector.validate_property_config("test_prop", 5.0, PropertyType.CONSTANT_VALUE)

    def test_validate_constant_property_invalid(self):
        """Non-numeric string constant raises ValueError."""
        with pytest.raises(ValueError):
            PropertyTypeDetector.validate_property_config("test_prop", "invalid", PropertyType.CONSTANT_VALUE)

    def test_validate_file_property_missing_keys(self):
        """File-import config missing required keys raises ValueError."""
        with pytest.raises(ValueError, match="Invalid configuration"):
            PropertyTypeDetector.validate_property_config("test_prop", {"file_path": "data.csv"}, PropertyType.FILE_IMPORT)

    def test_validate_step_function_invalid_values(self):
        """Step function with fewer than two values raises ValueError."""
        config = {"dependency": "melting_temperature", "value": [900]}
        with pytest.raises(ValueError):
            PropertyTypeDetector.validate_property_config("test_prop", config, PropertyType.STEP_FUNCTION)
