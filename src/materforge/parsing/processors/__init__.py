"""
Property processing modules for PyMatLib.

This package contains the core property processing functionality,
including specialized handlers for different property types,
dependency resolution, and post-processing.
"""

from .property_processor import PropertyProcessor
from .property_handlers import (
    BasePropertyHandler,
    ConstantValuePropertyHandler,
    StepFunctionPropertyHandler,
    FileImportPropertyHandler,
    TabularDataPropertyHandler,
    PiecewiseEquationPropertyHandler,
    ComputedPropertyHandler
)
from .dependency_processor import DependencyProcessor
from .post_processor import PropertyPostProcessor
from .temperature_resolver import TemperatureResolver

__all__ = [
    'PropertyProcessor',
    'BasePropertyHandler',
    'ConstantValuePropertyHandler',
    'StepFunctionPropertyHandler',
    'FileImportPropertyHandler',
    'TabularDataPropertyHandler',
    'PiecewiseEquationPropertyHandler',
    'ComputedPropertyHandler',
    'DependencyProcessor',
    'PropertyPostProcessor',
    'TemperatureResolver'
]
