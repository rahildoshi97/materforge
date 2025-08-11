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
from .computed_property_processor import ComputedPropertyProcessor
from .post_processor import PropertyPostProcessor
from .dependency_resolver import DependencyResolver

__all__ = [
    'PropertyProcessor',
    'BasePropertyHandler',
    'ConstantValuePropertyHandler',
    'StepFunctionPropertyHandler',
    'FileImportPropertyHandler',
    'TabularDataPropertyHandler',
    'PiecewiseEquationPropertyHandler',
    'ComputedPropertyHandler',
    'ComputedPropertyProcessor',
    'PropertyPostProcessor',
    'DependencyResolver'
]
