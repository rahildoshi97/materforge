# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import sympy as sp

# Property structure keys
DEPENDENCIES_KEY = "dependencies"
RANGES_KEY = "ranges"
COLUMNS_KEY = "columns"

# File property keys
FILE_PATH_KEY = "file_path"
DEPENDENCY_COLUMN_KEY = "dependency_column"  # Generic column key
PROPERTY_COLUMN_KEY = "property_column"

# Generic dependency keys
DEPENDENCY_KEY = "dependency"  # Generic dependency key for legacy format

# Value and equation keys
VALUE_KEY = "value"
EQUATION_KEY = "equation"

# Boundary condition keys
BOUNDS_KEY = "bounds"
CONSTANT_KEY = "constant"
LINEAR_KEY = "linear"

# Regression keys
REGRESSION_KEY = "regression"
SIMPLIFY_KEY = "simplify"
DEGREE_KEY = "degree"
SEGMENTS_KEY = "segments"
PRE_KEY = "pre"
POST_KEY = "post"

# Properties and material name
PROPERTIES_KEY = "properties"
NAME_KEY = "name"

# The placeholder symbol used in YAML equation strings.
# All expressions parsed from YAML are built in terms of this symbol;
# it is substituted with the caller-supplied dependency symbol at runtime.
YAML_PLACEHOLDER = sp.Symbol('T')

# Automatically export all constants
__all__ = [name for name in globals() if not name.startswith('_') and name.isupper()]
