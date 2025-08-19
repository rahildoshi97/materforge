"""
Main API module for MaterForge material property library.

This module provides the primary interface for creating and working with material objects
from YAML configuration files. It includes functions for material creation, validation,
property evaluation, and information extraction.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.config.material_yaml_parser import MaterialYAMLParser
from materforge.parsing.config.yaml_keys import (
    NAME_KEY, MATERIAL_TYPE_KEY, COMPOSITION_KEY, PROPERTIES_KEY,
    PURE_METAL_KEY, MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
    SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
    INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY,
    ALLOY_KEY
)

logger = logging.getLogger(__name__)


# ====================================================================
# CORE MATERIAL CREATION AND VALIDATION
# ====================================================================

def create_material(yaml_path: Union[str, Path], dependency: sp.Symbol,
                    enable_plotting: bool = True) -> Material:
    """Create material instance from YAML configuration file.

    This function serves as the main entry point for creating material
    objects from YAML configuration files. It handles the parsing of the configuration
    and creation of the material with the specified temperature.

    Parameters
    ----------
    yaml_path : Union[str, Path]
        Path to the YAML configuration file
    dependency : sp.Symbol
        Sympy symbol for property evaluation. Use a symbolic variable
        (e.g., sp.Symbol('T') or sp.Symbol('u_C')) for symbolic temperature expressions
    enable_plotting : bool, optional
        Whether to generate visualization plots (default: True)

    Returns
    -------
    Material
        The material instance with all properties initialized

    Raises
    ------
    FileNotFoundError
        If the YAML file doesn't exist
    ValueError
        If the YAML content is invalid or material creation fails
    TypeError
        If temperature parameter has invalid type

    Notes
    -----
    In YAML files, always use 'T' as the temperature variable in equations.
    The system will automatically substitute this with your provided symbol.

    Examples
    --------
    Create a material with symbolic temperature expressions:

    >>> import sympy as sp
    >>> T = sp.Symbol('T')
    >>> material = create_material('steel.yaml', T)
    >>> print(material.name)
    Steel

    Create a material with a custom temperature symbol:

    >>> u_C = sp.Symbol('u_C')
    >>> material_copper = create_material('copper.yaml', u_C)
    """
    logger.info("Creating material from: %s with dependency=%s, plotting=%s", yaml_path, dependency, enable_plotting)
    try:
        # Accept symbolic temperatures only
        if not isinstance(dependency, sp.Symbol):
            raise TypeError(f"Dependency '{dependency}' must be a sympy Symbol, got {type(dependency)}")
        parser = MaterialYAMLParser(yaml_path=yaml_path)
        material = parser.create_material(dependency=dependency, enable_plotting=enable_plotting)
        logger.info("Successfully created material: %s with %d properties",
                    material.name, len([attr for attr in dir(material)
                                        if not attr.startswith('_') and hasattr(material, attr)]))
        return material
    except Exception as e:
        logger.error("Failed to create material from %s: %s", yaml_path, e, exc_info=True)
        raise


def validate_yaml_file(yaml_path: Union[str, Path]) -> bool:
    """
    Validate a YAML file without creating the material.
    Args:
        yaml_path: Path to the YAML configuration file to validate
    Returns:
        bool: True if the file is valid
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the YAML content is invalid
    Example:
        try:
            is_valid = validate_yaml_file('steel.yaml')
            print(f"YAML file is valid: {is_valid}")
        except ValueError as e:
            print(f"Validation failed: {e}")
    """
    logger.info("Validating YAML file: %s", yaml_path)
    try:
        _ = MaterialYAMLParser(yaml_path)
        logger.info("YAML validation successful for: %s", yaml_path)
        return True
    except FileNotFoundError as e:
        logger.error("YAML file not found: %s", yaml_path)
        raise FileNotFoundError(f"YAML file not found: {yaml_path}") from e
    except ValueError as e:
        logger.error("YAML validation failed for %s: %s", yaml_path, e)
        raise ValueError(f"YAML validation failed: {str(e)}") from e
    except Exception as e:
        logger.error("Unexpected error validating YAML %s: %s", yaml_path, e, exc_info=True)
        raise ValueError(f"Unexpected error validating YAML: {str(e)}") from e


# ====================================================================
# MATERIAL INFORMATION AND PROPERTIES
# ====================================================================

def get_material_info(yaml_path: Union[str, Path]) -> Dict:
    """
    Get basic information about a material configuration without full processing.
    Args:
        yaml_path: Path to the YAML configuration file
    Returns:
        Dict: Dictionary containing material information including:
            - name: Material name
            - material_type: Type of material (pure_metal or alloy)
            - composition: Element composition dictionary
            - properties: List of available property names
            - total_properties: Number of properties defined
            - property_types: Count of each property type
            - Temperature properties based on material type
    Raises:
        FileNotFoundError: If the YAML file doesn't exist
        ValueError: If the YAML content is invalid
    Example:
        info = get_material_info('steel.yaml')
        print(f"Material: {info['name']}")
        print(f"Properties: {info['total_properties']}")
        print(f"Type: {info['material_type']}")
    """
    logger.info("Extracting material info from: %s", yaml_path)
    try:
        parser = MaterialYAMLParser(yaml_path=yaml_path)
        config = parser.config
        # Base information
        info = {
            'name': config.get(NAME_KEY, 'Unknown'),
            'material_type': config.get(MATERIAL_TYPE_KEY, 'Unknown'),
            'composition': config.get(COMPOSITION_KEY, {}),
        }
        # Add properties information
        properties = config.get(PROPERTIES_KEY, {})
        info['properties'] = list(properties.keys())
        info['total_properties'] = len(properties)
        # Add temperature-specific properties based on material type
        if info['material_type'] == PURE_METAL_KEY:
            info['melting_temperature'] = config.get(MELTING_TEMPERATURE_KEY, 'Undefined')
            info['boiling_temperature'] = config.get(BOILING_TEMPERATURE_KEY, 'Undefined')
        elif info['material_type'] == ALLOY_KEY:
            info['solidus_temperature'] = config.get(SOLIDUS_TEMPERATURE_KEY, 'Undefined')
            info['liquidus_temperature'] = config.get(LIQUIDUS_TEMPERATURE_KEY, 'Undefined')
            info['initial_boiling_temperature'] = config.get(INITIAL_BOILING_TEMPERATURE_KEY, 'Undefined')
            info['final_boiling_temperature'] = config.get(FINAL_BOILING_TEMPERATURE_KEY, 'Undefined')
        # Add property categorization info if available
        if hasattr(parser, 'categorized_properties') and parser.categorized_properties:
            info['property_types'] = {
                prop_type.name: len(props)
                for prop_type, props in parser.categorized_properties.items()
                if len(props) > 0
            }
        logger.info("Successfully extracted info for material: %s", info['name'])
        return info
    except FileNotFoundError as e:
        logger.error("YAML file not found: %s", yaml_path)
        raise FileNotFoundError(f"YAML file not found: {yaml_path}") from e
    except KeyError as e:
        logger.error("Missing required field in YAML %s: %s", yaml_path, e)
        raise ValueError(f"Missing required field in YAML: {str(e)}") from e
    except Exception as e:
        logger.error("Failed to extract material info from %s: %s", yaml_path, e, exc_info=True)
        raise ValueError(f"Failed to extract material info: {str(e)}") from e


def get_supported_properties() -> List[str]:
    """Get a list of all supported material properties.

    Returns
    -------
    List[str]
        List of strings representing valid property names that can be
        defined in YAML files
    Examples
    --------
    >>> props = get_supported_properties()
    >>> print(f"Supported properties: {len(props)}")
    Supported properties: 12
    >>> for prop in props[:3]:
    ...     print(f"  - {prop}")
      - density
      - heat_capacity
      - thermal_conductivity
    """
    return sorted(list(MaterialYAMLParser.VALID_YAML_PROPERTIES))


def get_material_property_names(material: Material) -> List[str]:
    """
    Get list of all available property names for a material instance.
    Args:
        material: Material instance
    Returns:
        List[str]: List of property names that exist (are not None) on the material
    Raises:
        ValueError: If material is not a Material instance
    Example:
        material = create_material('steel.yaml', dependency=sp.Symbol('T'))
        available = get_material_property_names(material)
        print(f"Available properties: {available}")
    """
    if not isinstance(material, Material):
        raise ValueError(f"Expected Material instance, got {type(material).__name__}")
    # Define all possible property names
    property_names = [
        'density', 'dynamic_viscosity', 'energy_density', 'heat_capacity',
        'heat_conductivity', 'kinematic_viscosity', 'latent_heat_of_fusion',
        'latent_heat_of_vaporization', 'specific_enthalpy', 'surface_tension',
        'thermal_diffusivity', 'thermal_expansion_coefficient'
    ]
    # Return only properties that exist and are not None
    return [name for name in property_names if getattr(material, name, None) is not None]


# ====================================================================
# PROPERTY EVALUATION
# ====================================================================

def evaluate_material_properties(material: Material, temperature: Union[float, int],
                                 properties: Optional[List[str]] = None,
                                 include_constants: bool = True) -> Dict[str, float]:
    """
    Convenience function to evaluate material properties at a specific temperature.

    This is a wrapper around Material.evaluate_properties_at_temperature() for
    functional-style usage.
    Args:
        material: Material instance
        temperature: Temperature value in Kelvin
        properties: List of specific property names to evaluate. If None, evaluates all.
        include_constants: Whether to include constant properties in the result
    Returns:
        Dict[str, float]: Dictionary mapping property names to their evaluated values
    Raises:
        ValueError: If material is not a Material instance or temperature is invalid
    Examples:
        # Evaluate all properties
        values = evaluate_material_properties(material, 500.0)
        # Evaluate specific properties
        values = evaluate_material_properties(material, 500.0, ['density', 'heat_capacity'])
        # Get only temperature-dependent properties
        values = evaluate_material_properties(material, 500.0, include_constants=False)
    """
    logger.info("Evaluating material properties via API function")
    if not isinstance(material, Material):
        raise ValueError(f"Expected Material instance, got {type(material).__name__}")
    return material.evaluate_properties_at_temperature(
        temperature=temperature,
        properties=properties,
        include_constants=include_constants
    )


# ====================================================================
# INTERNAL/TESTING FUNCTIONS
# ====================================================================

def _test_api():
    """
    Internal test function for API validation.

    This function is used for internal testing and should not be called
    by end users.
    """
    try:
        # Test basic validation
        test_path = Path("example.yaml")
        if test_path.exists():
            assert validate_yaml_file(test_path) is True
            logger.info("API test passed")
        else:
            logger.warning("Test file not found, skipping API test")
    except (FileNotFoundError, ValueError, AssertionError) as e:
        logger.error(f"API test failed: {e}")


# ====================================================================
# MODULE EXPORTS
# ====================================================================

__all__ = [
    'create_material',
    'validate_yaml_file',
    'get_material_info',
    'get_supported_properties',
    'get_material_property_names',
    'evaluate_material_properties'
]
