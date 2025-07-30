from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class ProcessingConstants:
    """Processing constants used throughout the YAML parser."""
    # Tolerance and precision
    DEFAULT_TOLERANCE: Final[float] = 1e-8
    TEMPERATURE_EPSILON: Final[float] = 1e-10
    FLOATING_POINT_TOLERANCE: Final[float] = 1e-12
    MONOTONICITY_THRESHOLD: Final[float] = 1e-8
    # Composition
    COMPOSITION_THRESHOLD: Final[float] = 1e-10
    # Temperature limits
    MIN_TEMPERATURE: Final[float] = 0.0
    # Data validation
    MIN_DATA_POINTS: Final[int] = 2
    MIN_TEMPERATURE_POINTS: Final[int] = 2
    # Regression limits
    MAX_REGRESSION_SEGMENTS: Final[int] = 8
    WARNING_REGRESSION_SEGMENTS: Final[int] = 6
    DEFAULT_REGRESSION_SEED: Final[int] = 13579
    # Visualization
    DEFAULT_VISUALIZATION_POINTS: Final[int] = 1000
    TEMPERATURE_PADDING_FACTOR: Final[float] = 0.1
    STEP_FUNCTION_OFFSET: Final[float] = 100.0
    # File processing
    MAX_MISSING_VALUE_PERCENTAGE: Final[float] = 50.0
    # Default ranges
    DEFAULT_TEMP_LOWER: Final[float] = 273.15
    DEFAULT_TEMP_UPPER: Final[float] = 3000.0
    # Shared regex for temperature arithmetic
    TEMP_ARITHMETIC_REGEX: Final[str] = r'^(\w+)\s*([+-])\s*(\d+(?:\.\d+)?)$'


@dataclass(frozen=True)
class ErrorMessages:
    """Standardized error message templates."""
    TEMPERATURE_BELOW_ZERO: Final[str] = "Temperature must be above absolute zero ({min_temp}K), got {temp}K"
    INSUFFICIENT_DATA_POINTS: Final[str] = "Insufficient data points ({count}), minimum required: {min_points}"
    INVALID_PROPERTY_TYPE: Final[str] = "Invalid property type for '{prop_name}': {prop_type}"
    MISSING_DEPENDENCY: Final[str] = "Missing dependency '{dep}' for property '{prop_name}'"
    CIRCULAR_DEPENDENCY: Final[str] = "Circular dependency detected: {cycle_path}"


@dataclass(frozen=True)
class FileConstants:
    """File processing related constants."""
    SUPPORTED_EXTENSIONS: Final[tuple] = ('.csv', '.xlsx', '.txt')
    MAX_FILE_SIZE_MB: Final[int] = 100
    DEFAULT_ENCODING: Final[str] = 'utf-8'
    # Missing value indicators
    NA_VALUES: Final[tuple] = ('', ' ', '  ', '   ', 'nan', 'NaN', 'NULL', 'null', 'N/A', 'n/a', 'NA')
