"""Constants used for YAML parsing and property processing."""

# ── TOP-LEVEL KEYS ───────────────────────────────────────────────

INDEPENDENT_VARIABLES_KEY = "independent_variables"  # Maps physical dependencies to symbols

# Property structure keys (Multi-dependency format)
DEPENDENCIES_KEY = "dependencies"  # List of dependencies for a property
RANGES_KEY = "ranges"  # Dictionary mapping dependency names to ranges
BOUNDS_KEY = "bounds"  # Dictionary mapping dependency names to boundary conditions
COLUMNS_KEY = "columns"  # Dictionary mapping dependency names to column names

MAX_STEP_FUNCTION_DEPENDENCIES = 1  # Maximum number of dependencies for step functions

# File property keys
FILE_PATH_KEY = "file_path"
DEPENDENCY_COLUMN_KEY = "dependency_column"  # Legacy for single dependency
PROPERTY_COLUMN_KEY = "property_column"

# Legacy keys (for backward compatibility - will be removed)
# TEMPERATURE_COLUMN_KEY = "temperature_column"
# TEMPERATURE_KEY = "temperature"
# DEPENDENCY_KEY = "dependency"  # Legacy single dependency key

# Value and equation keys
VALUE_KEY = "value"
EQUATION_KEY = "equation"

# Boundary condition keys
CONSTANT_KEY = "constant"
EXTRAPOLATE_KEY = "extrapolate"

# Regression keys
REGRESSION_KEY = "regression"
SIMPLIFY_KEY = "simplify"
DEGREE_KEY = "degree"
SEGMENTS_KEY = "segments"
PRE_KEY = "pre"
POST_KEY = "post"

# Material type keys
MATERIAL_TYPE_KEY = "material_type"
PURE_METAL_KEY = "pure_metal"
ALLOY_KEY = "alloy"

# Composition key
COMPOSITION_KEY = "composition"

# Temperature reference keys
MELTING_TEMPERATURE_KEY = "melting_temperature"
BOILING_TEMPERATURE_KEY = "boiling_temperature"
SOLIDUS_TEMPERATURE_KEY = "solidus_temperature"
LIQUIDUS_TEMPERATURE_KEY = "liquidus_temperature"
INITIAL_BOILING_TEMPERATURE_KEY = "initial_boiling_temperature"
FINAL_BOILING_TEMPERATURE_KEY = "final_boiling_temperature"

# Properties and material name
PROPERTIES_KEY = "properties"
NAME_KEY = "name"

# Supported material types
SUPPORTED_MATERIAL_TYPES = {
    PURE_METAL_KEY, ALLOY_KEY
}

# Supported dependency names
SUPPORTED_DEPENDENCY_NAMES = {
    "temperature", "strain_rate", "concentration", "pressure", "time"
}

# Default dependency name for backward compatibility
# DEFAULT_DEPENDENCY_NAME = "temperature"

MAX_DEPENDENCIES_KEY = 1

# Automatically export all constants
__all__ = [name for name in globals() if not name.startswith('_') and name.isupper()]
