"""Validation utilities for PyMatLib."""

from .array_validator import is_monotonic
from .errors import PropertyError, DependencyError, CircularDependencyError
from .property_type_detector import PropertyType, PropertyTypeDetector
from .property_validator import validate_monotonic_energy_density, validate_monotonic_property

__all__ = [
    "PropertyError",
    "DependencyError",
    "CircularDependencyError",
    "PropertyType",
    "PropertyTypeDetector",
    "validate_monotonic_energy_density",
    "validate_monotonic_property",
    "is_monotonic"
]
