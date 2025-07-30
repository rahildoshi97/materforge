"""Constants used for YAML parsing and property processing."""

# ── TOP-LEVEL KEYS ───────────────────────────────────────────────
INDEPENDENT_VARIABLES_KEY = "independent_variables"  # Note: plural for multi-dependency

# Property structure keys
DEPENDENCIES_KEY = "dependencies"
RANGES_KEY = "ranges"
COLUMNS_KEY = "columns"

# File property keys
FILE_PATH_KEY = "file_path"
DEPENDENCY_COLUMN_KEY = "dependency_column"  # Generic column key
PROPERTY_COLUMN_KEY = "property_column"

# Legacy keys (for backward compatibility)
TEMPERATURE_COLUMN_KEY = "temperature_column"
TEMPERATURE_KEY = "temperature"

# Generic dependency keys
DEPENDENCY_KEY = "dependency"  # Generic dependency key for legacy format

# Value and equation keys
VALUE_KEY = "value"
EQUATION_KEY = "equation"

# Boundary condition keys
BOUNDS_KEY = "bounds"
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

# Supported dependency names
SUPPORTED_DEPENDENCY_NAMES = {
    "temperature", "strain_rate", "concentration", "pressure", "time"
}

# Automatically export all constants
__all__ = [name for name in globals() if not name.startswith('_') and name.isupper()]
