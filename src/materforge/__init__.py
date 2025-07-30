"""
MaterForge - A Python library for temperature-dependent material property modeling.

This library provides tools for defining, processing, and evaluating material properties
as functions of temperature, with support for various property definition formats
including constants, piecewise functions, file-based data, and computed properties.

Key Features:
- Temperature-dependent material property modeling
- Multiple property definition formats (YAML-based)
- Symbolic mathematics integration with SymPy
- Piecewise function creation and evaluation
- Material property visualization
- Integration with numerical simulation frameworks

Main Components:
- Core: Material definitions and fundamental data structures
- Parsing: YAML configuration parsing and property processing
- Algorithms: Mathematical operations and property computations
- Visualization: Property plotting and analysis tools
- Data: Material databases and physical constants
"""

# Enhanced version handling with multiple fallbacks
try:
    from ._version import version as __version__
except ImportError:
    try:
        from importlib.metadata import version
        __version__ = version("materforge")
    except ImportError:
        try:
            from importlib_metadata import version
            __version__ = version("materforge")
        except ImportError:
            __version__ = "0.3.0+unknown"  # Fallback version

# Core material definitions
from .core.materials import Material
from .core.elements import ChemicalElement
from .core.symbol_registry import SymbolRegistry

# Main API functions
from .parsing.api import (
    create_material,
    get_supported_properties,
    validate_yaml_file
)

# Property processing
from .parsing.processors.property_processor import PropertyProcessor
from .parsing.validation.property_type_detector import PropertyType

# Algorithms
from .algorithms.interpolation import interpolate_value, ensure_ascending_order
from .algorithms.piecewise_builder import PiecewiseBuilder
from .algorithms.piecewise_inverter import PiecewiseInverter

# Visualization
from .visualization.plotters import PropertyVisualizer

__all__ = [
    # Version
    '__version__',

    # Core classes
    'Material',
    'ChemicalElement',
    'SymbolRegistry',

    # Main API
    'create_material',
    'get_supported_properties',
    'validate_yaml_file',

    # Processing
    'PropertyProcessor',
    'PropertyType',

    # Algorithms
    'interpolate_value',
    'ensure_ascending_order',
    'PiecewiseBuilder',
    'PiecewiseInverter',

    # Visualization
    'PropertyVisualizer'
]

# Package metadata
__author__ = "Rahil Doshi"
__email__ = "rahil.doshi@fau.de"
__description__ = "Temperature-dependent material property modeling library"
