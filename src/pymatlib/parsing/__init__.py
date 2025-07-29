"""
Parsing and configuration modules for PyMatLib.

This package handles YAML parsing, property type detection, validation,
and material creation from configuration files.
"""

from .api import create_material, get_supported_properties, validate_yaml_file, get_material_info, evaluate_material_properties, get_material_property_names
from .config.material_yaml_parser import MaterialYAMLParser
from .validation.property_type_detector import PropertyType, PropertyTypeDetector
from .processors.property_processor import PropertyProcessor

__all__ = [
    'create_material',
    'get_supported_properties',
    'validate_yaml_file',
    'MaterialYAMLParser',
    'PropertyType',
    'PropertyTypeDetector',
    'PropertyProcessor',
    'get_material_info',
    'evaluate_material_properties',
    'get_material_property_names',
]
