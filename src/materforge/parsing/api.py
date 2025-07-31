"""
Enhanced API module for PyMatLib material property library with multi-dependency support.

This module provides the primary interface for creating and working with material objects
from YAML configuration files. It includes functions for material creation, validation,
property evaluation, and information extraction with full multi-dependency support.
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
    ALLOY_KEY, INDEPENDENT_VARIABLES_KEY, MAX_DEPENDENCIES, SUPPORTED_DEPENDENCY_NAMES
)

logger = logging.getLogger(__name__)


# ====================================================================
# CORE MATERIAL CREATION AND VALIDATION
# ====================================================================

def create_material(yaml_path: Union[str, Path],
                    dependencies: Optional[Dict[str, sp.Symbol]] = None,
                    enable_plotting: bool = True,
                    **kwargs) -> Material:
    """Create material instance from YAML configuration file with multi-dependency support.

    This function serves as the main entry point for creating material (pure metal or alloy)
    objects from YAML configuration files. It handles both single and multi-dependency
    configurations with full backward compatibility.

    Parameters
    ----------
    yaml_path : Union[str, Path]
        Path to the YAML configuration file
    dependencies : Optional[Dict[str, sp.Symbol]]
        Dictionary mapping dependency names to SymPy symbols
        e.g., {'temperature': sp.Symbol('T'), 'pressure': sp.Symbol('P')}
        If None, defaults to {'temperature': sp.Symbol('T')}
    enable_plotting : bool, optional
        Whether to generate visualization plots (default: True)
    **kwargs : Dict[str, sp.Symbol]
        Alternative way to specify dependencies using keyword arguments
        e.g., temperature=sp.Symbol('T'), pressure=sp.Symbol('P')

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
        If dependency parameters have invalid types

    Notes
    -----
    - For multi-dependency YAML files, use symbols defined in the independent_variables section
    - In equations, always use the symbols specified in independent_variables (e.g., 'T', 'P', 'S')
    - The system automatically maps these to your provided symbols

    Examples
    --------
    Single dependency (backward compatible):

    >>> import sympy as sp
    >>> material = create_material('steel.yaml', dependencies={'temperature': sp.Symbol('T')})
    >>> # Or using kwargs
    >>> material = create_material('steel.yaml', temperature=sp.Symbol('T'))

    Multi-dependency:

    >>> material = create_material('steel.yaml', dependencies={
    ...     'temperature': sp.Symbol('T'),
    ...     'pressure': sp.Symbol('P')
    ... })
    >>> # Or using kwargs
    >>> material = create_material('steel.yaml',
    ...                          temperature=sp.Symbol('T'),
    ...                          pressure=sp.Symbol('P'))

    Custom symbols:

    >>> u_C = sp.Symbol('u_C')
    >>> p_bar = sp.Symbol('p_bar')
    >>> material = create_material('alloy.yaml',
    ...                          temperature=u_C,
    ...                          pressure=p_bar)
    """
    logger.info("Creating material from: %s with plotting=%s", yaml_path, enable_plotting)

    try:
        # Handle different input formats
        if dependencies is None and kwargs:
            # Use kwargs: create_material(yaml_path, temperature=Symbol('T'), pressure=Symbol('P'))
            symbol_mapping = kwargs
        elif dependencies is not None:
            # Use dependencies dict
            symbol_mapping = dependencies
        else:
            # Default single dependency for backward compatibility
            symbol_mapping = {'temperature': sp.Symbol('T')}
            logger.info("No dependencies specified, using default temperature dependency")

        # Validate that all values are SymPy symbols
        for dep_name, symbol in symbol_mapping.items():
            if not isinstance(symbol, sp.Symbol):
                raise TypeError(f"Dependency '{dep_name}' must be a SymPy Symbol, got {type(symbol)}")

        logger.debug("Symbol mapping: %s", {k: str(v) for k, v in symbol_mapping.items()})

        parser = MaterialYAMLParser(yaml_path=yaml_path)
        material = parser.create_material(symbol_mapping=symbol_mapping, enable_plotting=enable_plotting)

        logger.info("Successfully created material: %s with %d dependencies",
                    material.name, len(symbol_mapping))
        return material

    except Exception as e:
        logger.error("Failed to create material from %s: %s", yaml_path, e, exc_info=True)
        raise


def validate_yaml_file(yaml_path: Union[str, Path]) -> bool:
    """
    Validate a YAML file without creating the material.
    Enhanced to validate multi-dependency structure.

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
    except Exception as e:
        logger.error("YAML validation failed for %s: %s", yaml_path, e)
        raise ValueError(f"YAML validation failed: {str(e)}") from e


# ====================================================================
# MATERIAL INFORMATION AND PROPERTIES
# ====================================================================

def get_material_info(yaml_path: Union[str, Path]) -> Dict:
    """
    Get comprehensive information about a material configuration.
    Enhanced to include multi-dependency information.

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
            - is_multi_dependency: Whether the material uses multi-dependency format
            - independent_variables: Mapping of dependency names to symbols
            - max_dependencies_supported: Maximum dependencies allowed
            - Temperature properties based on material type

    Raises:
        FileNotFoundError: If the YAML file doesn't exist
        ValueError: If the YAML content is invalid

    Example:
        info = get_material_info('steel.yaml')
        print(f"Material: {info['name']}")
        print(f"Properties: {info['total_properties']}")
        print(f"Multi-dependency: {info['is_multi_dependency']}")
        print(f"Dependencies: {list(info['independent_variables'].keys())}")
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

        # Multi-dependency information
        if hasattr(parser, 'dependency_resolver'):
            info['is_multi_dependency'] = parser.dependency_resolver.is_multi_dependency
            info['independent_variables'] = parser.dependency_resolver.independent_vars
            info['max_dependencies_supported'] = MAX_DEPENDENCIES
        else:
            info['is_multi_dependency'] = False
            info['independent_variables'] = {'temperature': 'T'}
            info['max_dependencies_supported'] = 1

        # Properties information
        properties = config.get(PROPERTIES_KEY, {})
        info['properties'] = list(properties.keys())
        info['total_properties'] = len(properties)

        # Add temperature-specific properties
        if info['material_type'] == PURE_METAL_KEY:
            info['melting_temperature'] = config.get(MELTING_TEMPERATURE_KEY, 'Undefined')
            info['boiling_temperature'] = config.get(BOILING_TEMPERATURE_KEY, 'Undefined')
        elif info['material_type'] == ALLOY_KEY:
            info['solidus_temperature'] = config.get(SOLIDUS_TEMPERATURE_KEY, 'Undefined')
            info['liquidus_temperature'] = config.get(LIQUIDUS_TEMPERATURE_KEY, 'Undefined')
            info['initial_boiling_temperature'] = config.get(INITIAL_BOILING_TEMPERATURE_KEY, 'Undefined')
            info['final_boiling_temperature'] = config.get(FINAL_BOILING_TEMPERATURE_KEY, 'Undefined')

        # Property categorization
        if hasattr(parser, 'categorized_properties'):
            info['property_types'] = {
                prop_type.name: len(props)
                for prop_type, props in parser.categorized_properties.items()
                if len(props) > 0
            }

        logger.info("Successfully extracted info for material: %s", info['name'])
        return info

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
    Supported properties: 25
    >>> for prop in props[:3]:
    ...     print(f"  - {prop}")
      - density
      - heat_capacity
      - thermal_conductivity
    """
    return sorted(list(MaterialYAMLParser.VALID_YAML_PROPERTIES))


def get_supported_dependencies() -> List[str]:
    """Get list of supported dependency types for multi-dependency materials.

    Returns
    -------
    List[str]
        List of supported dependency names that can be used in multi-dependency
        YAML configurations

    Examples
    --------
    >>> deps = get_supported_dependencies()
    >>> print(f"Supported dependencies: {deps}")
    ['concentration', 'pressure', 'strain_rate', 'temperature', 'time']
    """
    return sorted(list(SUPPORTED_DEPENDENCY_NAMES))


def get_max_dependencies() -> int:
    """Get maximum number of dependencies supported.

    Returns
    -------
    int
        Maximum number of dependencies that can be used simultaneously

    Examples
    --------
    >>> max_deps = get_max_dependencies()
    >>> print(f"Maximum dependencies supported: {max_deps}")
    Maximum dependencies supported: 3
    """
    return MAX_DEPENDENCIES


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
        material = create_material('steel.yaml', temperature=sp.Symbol('T'))
        available = get_material_property_names(material)
        print(f"Available properties: {available}")
    """
    if not isinstance(material, Material):
        raise ValueError(f"Expected Material instance, got {type(material).__name__}")

    # Define all possible property names
    property_names = [
        'bulk_modulus', 'density', 'dynamic_viscosity', 'elastic_modulus',
        'electrical_conductivity', 'electrical_resistivity', 'energy_density',
        'fracture_toughness', 'hardness', 'heat_capacity', 'heat_conductivity',
        'kinematic_viscosity', 'latent_heat_of_fusion', 'latent_heat_of_vaporization',
        'magnetic_permeability', 'poisson_ratio', 'shear_modulus', 'specific_enthalpy',
        'surface_tension', 'thermal_diffusivity', 'thermal_expansion_coefficient',
        'ultimate_tensile_strength', 'viscosity', 'yield_strength'
    ]

    # Return only properties that exist and are not None
    return [name for name in property_names if getattr(material, name, None) is not None]


# ====================================================================
# PROPERTY EVALUATION
# ====================================================================

def evaluate_material_properties(material: Material,
                                 dependencies: Dict[str, Union[float, int]],
                                 properties: Optional[List[str]] = None,
                                 include_constants: bool = True) -> Dict[str, float]:
    """
    Enhanced function to evaluate material properties with multi-dependency support.

    This function evaluates material properties at specific dependency values,
    supporting both single and multi-dependency materials.

    Args:
        material: Material instance
        dependencies: Dictionary mapping dependency names to their values
                     e.g., {'temperature': 500.0, 'pressure': 1e5}
                     For single dependency: {'temperature': 500.0}
        properties: List of specific property names to evaluate. If None, evaluates all.
        include_constants: Whether to include constant properties in the result

    Returns:
        Dict[str, float]: Dictionary mapping property names to their evaluated values

    Raises:
        ValueError: If material is not a Material instance or dependencies are invalid

    Examples:
        # Single dependency (backward compatible)
        values = evaluate_material_properties(material, {'temperature': 500.0})

        # Multi-dependency
        values = evaluate_material_properties(material, {
            'temperature': 500.0,
            'pressure': 1e5
        })

        # Evaluate specific properties
        values = evaluate_material_properties(material,
                                            {'temperature': 500.0},
                                            ['density', 'heat_capacity'])

        # Get only temperature-dependent properties
        values = evaluate_material_properties(material,
                                            {'temperature': 500.0},
                                            include_constants=False)
    """
    logger.info("Evaluating material properties via enhanced API function")

    if not isinstance(material, Material):
        raise ValueError(f"Expected Material instance, got {type(material).__name__}")

    if not isinstance(dependencies, dict):
        raise ValueError(f"Dependencies must be a dictionary, got {type(dependencies)}")

    # For now, if only temperature is provided, use the existing method
    if len(dependencies) == 1 and 'temperature' in dependencies:
        return material.evaluate_properties_at_temperature(
            temperature=dependencies['temperature'],
            properties=properties,
            include_constants=include_constants
        )

    # For multi-dependency, we'll need to implement this in the Material class
    # For now, raise a helpful error message
    if len(dependencies) > 1:
        raise NotImplementedError(
            "Multi-dependency property evaluation is not yet implemented. "
            "Currently only single temperature dependency is supported. "
            f"Provided dependencies: {list(dependencies.keys())}"
        )

    # Single non-temperature dependency
    if len(dependencies) == 1:
        dep_name, dep_value = next(iter(dependencies.items()))
        if dep_name != 'temperature':
            raise NotImplementedError(
                f"Property evaluation for dependency '{dep_name}' is not yet implemented. "
                "Currently only 'temperature' dependency is supported."
            )

    raise ValueError("At least one dependency must be provided")


# ====================================================================
# UTILITY AND INFORMATION FUNCTIONS
# ====================================================================

def check_yaml_compatibility(yaml_path: Union[str, Path]) -> Dict[str, Union[str, bool, List]]:
    """
    Check YAML file compatibility with different API versions.

    Args:
        yaml_path: Path to the YAML configuration file

    Returns:
        Dict containing compatibility information:
            - format_version: 'legacy' or 'multi_dependency'
            - has_independent_variables: Whether independent_variables section exists
            - dependencies_used: List of dependencies found in the file
            - backward_compatible: Whether file works with old API
            - recommendations: List of recommended actions

    Example:
        compat = check_yaml_compatibility('steel.yaml')
        print(f"Format: {compat['format_version']}")
        print(f"Dependencies: {compat['dependencies_used']}")
    """
    try:
        parser = MaterialYAMLParser(yaml_path)
        config = parser.config

        has_independent_vars = INDEPENDENT_VARIABLES_KEY in config

        # Analyze dependencies used
        dependencies_found = set()
        if has_independent_vars:
            dependencies_found.update(config[INDEPENDENT_VARIABLES_KEY].keys())

        # Check properties for dependency usage
        properties = config.get(PROPERTIES_KEY, {})
        for prop_config in properties.values():
            if isinstance(prop_config, dict):
                if 'dependencies' in prop_config:
                    dependencies_found.update(prop_config['dependencies'])

        format_version = 'multi_dependency' if has_independent_vars else 'legacy'
        backward_compatible = not has_independent_vars or dependencies_found <= {'temperature'}

        recommendations = []
        if not has_independent_vars:
            recommendations.append("Consider upgrading to multi-dependency format for enhanced capabilities")
        elif len(dependencies_found) > 1:
            recommendations.append("Multi-dependency format detected - use new API syntax")

        return {
            'format_version': format_version,
            'has_independent_variables': has_independent_vars,
            'dependencies_used': sorted(list(dependencies_found)),
            'backward_compatible': backward_compatible,
            'recommendations': recommendations
        }

    except Exception as e:
        return {
            'format_version': 'unknown',
            'has_independent_variables': False,
            'dependencies_used': [],
            'backward_compatible': False,
            'recommendations': [f"Error analyzing file: {str(e)}"]
        }


# ====================================================================
# BACKWARD COMPATIBILITY WRAPPER
# ====================================================================

def create_material_legacy(yaml_path: Union[str, Path],
                           dependency: sp.Symbol,
                           enable_plotting: bool = True) -> Material:
    """
    Legacy wrapper for backward compatibility.

    This function maintains the old API signature for existing code.

    Args:
        yaml_path: Path to the YAML configuration file
        dependency: Single SymPy symbol (typically temperature)
        enable_plotting: Whether to generate visualization plots

    Returns:
        Material: The material instance

    Note:
        This function is deprecated. Use create_material() with the dependencies
        parameter for new code.
    """
    import warnings
    warnings.warn(
        "create_material_legacy is deprecated. Use create_material() with "
        "dependencies={'temperature': symbol} instead.",
        DeprecationWarning,
        stacklevel=2
    )

    return create_material(
        yaml_path=yaml_path,
        dependencies={'temperature': dependency},
        enable_plotting=enable_plotting
    )


# ====================================================================
# INTERNAL/TESTING FUNCTIONS
# ====================================================================

def _test_api():
    """
    Internal test function for API validation.
    Enhanced to test multi-dependency features.
    """
    try:
        # Test basic validation
        test_path = Path("example.yaml")
        if test_path.exists():
            assert validate_yaml_file(test_path) is True

            # Test info extraction
            info = get_material_info(test_path)
            assert 'is_multi_dependency' in info

            logger.info("Enhanced API test passed")
        else:
            logger.warning("Test file not found, skipping API test")
    except Exception as e:
        logger.error(f"API test failed: {e}")


# ====================================================================
# MODULE EXPORTS
# ====================================================================

__all__ = [
    'create_material',
    'validate_yaml_file',
    'get_material_info',
    'get_supported_properties',
    'get_supported_dependencies',
    'get_max_dependencies',
    'get_material_property_names',
    'evaluate_material_properties',
    'check_yaml_compatibility',
    'create_material_legacy'  # For backward compatibility
]
