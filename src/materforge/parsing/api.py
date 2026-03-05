# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Main API module for MaterForge material property library."""

import logging
from pathlib import Path
from typing import Dict, Set, Union
import sympy as sp
from materforge.core.materials import Material
from materforge.parsing.config.material_yaml_parser import MaterialYAMLParser
from materforge.parsing.config.yaml_keys import NAME_KEY, PROPERTIES_KEY

logger = logging.getLogger(__name__)


# ====================================================================
# CORE MATERIAL CREATION AND VALIDATION
# ====================================================================

def create_material(yaml_path: Union[str, Path], dependency: sp.Symbol,
                    enable_plotting: bool = True) -> Material:
    """Creates a Material from a YAML configuration file.

    Args:
        yaml_path:       Path to the YAML configuration file.
        dependency:      SymPy symbol used as the independent variable.
                         YAML equations always use the placeholder 'T';
                         it is substituted with this symbol at runtime.
        enable_plotting: Generate visualisation plots (default: True).
    Returns:
        Fully initialised Material instance.
    Raises:
        FileNotFoundError: YAML file does not exist.
        TypeError:         dependency is not a sp.Symbol.
        ValueError:        YAML content is invalid or material creation fails.
    Example:
        >>> material = create_material('steel.yaml', sp.Symbol('T'))
        >>> material = create_material('copper.yaml', sp.Symbol('u_C'), enable_plotting=False)
    """
    logger.info("Creating material from: %s (dependency=%s, plotting=%s)",
                yaml_path, dependency, enable_plotting)
    if not isinstance(dependency, sp.Symbol):
        raise TypeError(
            f"dependency '{dependency}' must be a sympy Symbol, got {type(dependency).__name__}")
    try:
        parser = MaterialYAMLParser(yaml_path=yaml_path)
        material = parser.create_material(dependency=dependency, enable_plotting=enable_plotting)
        logger.info("Successfully created material '%s' with %d properties",
                    material.name, len(material.property_names()))
        return material
    except Exception as e:
        # Top of the materforge stack - log once here, then re-raise as-is
        logger.error("Failed to create material from %s: %s", yaml_path, e)
        raise


def validate_yaml_file(yaml_path: Union[str, Path]) -> bool:
    """Validates a YAML file without creating the material.

    Args:
        yaml_path: Path to the YAML configuration file to validate.
    Returns:
        True if the file is structurally valid.
    Raises:
        FileNotFoundError: File does not exist.
        ValueError:        Content is invalid.
    Example:
        >>> is_valid = validate_yaml_file('steel.yaml')
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
    """Extracts material metadata from a YAML file without full processing.

    Args:
        yaml_path: Path to the YAML configuration file.
    Returns:
        Dict with keys: name, properties, total_properties, property_types,
        and any additional top-level YAML fields.
    Raises:
        FileNotFoundError: File does not exist.
        ValueError:        Content is invalid or a required field is missing.
    Example:
        >>> info = get_material_info('steel.yaml')
        >>> print(info['name'], info['total_properties'])
    """
    logger.info("Extracting material info from: %s", yaml_path)
    try:
        parser = MaterialYAMLParser(yaml_path=yaml_path)
        config = parser.config
        info: Dict = {'name': config.get(NAME_KEY, 'Unknown')}
        properties = config.get(PROPERTIES_KEY, {})
        info['properties'] = list(properties.keys())
        info['total_properties'] = len(properties)
        reserved_keys = {NAME_KEY, PROPERTIES_KEY}
        for key, value in config.items():
            if key not in reserved_keys:
                info[key] = value
        if hasattr(parser, "categorized_properties") and parser.categorized_properties:
            info['property_types'] = {
                pt.name: len(props)
                for pt, props in parser.categorized_properties.items()
                if props
            }
        logger.info("Successfully extracted info for material: '%s'", info['name'])
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


def get_material_property_names(material: Material) -> Set[str]:
    """Returns all property names dynamically assigned to a material instance.

    Args:
        material: A fully processed Material instance.
    Returns:
        Property names assigned during processing.
    Raises:
        ValueError: Argument is not a Material instance.
    Example:
        >>> material = create_material('steel.yaml', sp.Symbol('T'))
        >>> print(get_material_property_names(material))
    """
    if not isinstance(material, Material):
        raise ValueError(f"Expected Material instance, got {type(material).__name__}")
    return material.property_names()


# ====================================================================
# PROPERTY EVALUATION
# ====================================================================

def evaluate_material_properties(material: Material, symbol: sp.Symbol, value) -> Material:
    """Evaluates all symbolic properties at a given numeric value.

    Args:
        material: A fully processed Material instance.
        symbol:   SymPy symbol to substitute.
        value:    Numeric value to substitute.
    Returns:
        New Material with all properties evaluated to numeric values.
    Raises:
        ValueError: If material is not a Material instance.
    Example:
        >>> evaluate_material_properties(material, T, 500.0)
    """
    if not isinstance(material, Material):
        raise ValueError(f"Expected Material instance, got {type(material).__name__}")
    return material.evaluate(symbol, value)


# ====================================================================
# INTERNAL/TESTING FUNCTIONS
# ====================================================================

def _test_api() -> None:
    """Internal test function. Not intended for end-user use."""
    try:
        test_path = Path("example.yaml")
        if test_path.exists():
            assert validate_yaml_file(test_path) is True
            logger.info("API test passed")
        else:
            logger.warning("Test file not found, skipping API test")
    except (FileNotFoundError, ValueError, AssertionError) as e:
        logger.error("API test failed: %s", e)


# ====================================================================
# MODULE EXPORTS
# ====================================================================

__all__ = [
    'create_material',
    'validate_yaml_file',
    'get_material_info',
    'get_material_property_names',
    'evaluate_material_properties',
]
