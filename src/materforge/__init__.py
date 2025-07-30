"""
MaterForge - Materials Formulation Engine with Python

A high-performance Python library for material property modeling
and analysis. MaterForge provides comprehensive tools for defining, processing, and
evaluating material properties as functions of temperature and other dependencies.

Key Features:
- Temperature-dependent material property modeling
- Multiple property definition formats (YAML-based)
- Symbolic mathematics integration with SymPy
- Piecewise function creation and evaluation
- Material property visualization
- Regression and data analysis capabilities
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
            __version__ = "0.5.5+unknown"  # Updated fallback version

# Core material definitions
from .core.materials import Material
from .core.elements import ChemicalElement
from .core.symbol_registry import SymbolRegistry

# Main API functions
from .parsing.api import (
    create_material,
    get_supported_properties,
    validate_yaml_file,
    get_material_info
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
    'get_material_info',

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
__description__ = "Materials Formulation Engine with Python"
