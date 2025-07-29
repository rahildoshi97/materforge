"""Constants used for YAML parsing and property processing."""

# File property keys
FILE_PATH_KEY = "file_path"
X_DATA_COLUMN_KEY = "dependency_column"
PROPERTY_COLUMN_KEY = "property_column"

# Temperature and value keys
X_DATA_KEY = "dependency"
VALUE_KEY = "value"

# Equation keys
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
# Extend to include other material types as needed

# Composition key
COMPOSITION_KEY = "composition"

# Pure metal temperature points
MELTING_TEMPERATURE_KEY = "melting_temperature"
BOILING_TEMPERATURE_KEY = "boiling_temperature"

# Alloy temperature points
SOLIDUS_TEMPERATURE_KEY = "solidus_temperature"
LIQUIDUS_TEMPERATURE_KEY = "liquidus_temperature"
INITIAL_BOILING_TEMPERATURE_KEY = "initial_boiling_temperature"
FINAL_BOILING_TEMPERATURE_KEY = "final_boiling_temperature"

# Properties key
PROPERTIES_KEY = "properties"

# Material name key
NAME_KEY = "name"

# Automatically export all constants (those that don't start with underscore)
__all__ = [name for name in globals() if not name.startswith('_') and name.isupper()]
