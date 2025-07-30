import logging
import sympy as sp
from dataclasses import dataclass
from typing import List, Union

logger = logging.getLogger(__name__)


@dataclass
class ChemicalElement:
    name: str
    atomic_number: float
    atomic_mass: float
    melting_temperature: float
    boiling_temperature: float
    latent_heat_of_fusion: float
    latent_heat_of_vaporization: float


# Utility functions for element property interpolation
def interpolate(values: List[float], composition: List[float]) -> Union[float, sp.Expr]:
    """
    Interpolates a property based on its values and composition.
    Args:
        values (list): List of property values.
        composition (list): List of composition percentages.
    Returns:
        float: Interpolated property value.
    """
    if len(values) != len(composition):
        logger.error("Length mismatch: values=%d, composition=%d", len(values), len(composition))
        raise ValueError(f"Values and composition arrays must have same length: {len(values)} vs {len(composition)}")
    logger.debug("Interpolating property with %d components", len(values))
    result = 0.
    for i, (v, c) in enumerate(zip(values, composition)):
        if c < 0 or c > 1:
            logger.warning("Composition value %d out of range [0,1]: %f", i, c)
        result += v * c
        logger.debug("Component %d: value=%.3f, composition=%.3f, contribution=%.3f", i, v, c, v*c)
    logger.debug("Interpolation result: %.6f", result)
    return result


def interpolate_atomic_number(elements: List[ChemicalElement], composition: List[float]) -> float:
    """
    Interpolates the atomic number based on the elements and their composition.
    Args:
        elements (list[ChemicalElement]): List of elements.
        composition (list[float]): List of composition percentages.
    Returns:
        float: Interpolated atomic number.
    """
    logger.debug("Interpolating atomic number for %d elements", len(elements))
    values = [element.atomic_number for element in elements]
    result = interpolate(values, composition)
    logger.info("Interpolated atomic number: %.3f", result)
    return result


def interpolate_atomic_mass(elements: List[ChemicalElement], composition: List[float]) -> float:
    """
    Interpolates the atomic mass based on the elements and their composition.
    Args:
        elements (list[ChemicalElement]): List of elements.
        composition (list[float]): List of composition percentages.
    Returns:
        float: Interpolated atomic mass.
    """
    logger.debug("Interpolating atomic mass for %d elements", len(elements))
    values = [element.atomic_mass for element in elements]
    result = interpolate(values, composition)
    logger.info("Interpolated atomic mass: %.3f g/mol", result)
    return result


def interpolate_melting_temperature(elements: List[ChemicalElement], composition: List[float]) -> float:
    """
    Interpolates the melting temperature based on the elements and their composition.
    """
    logger.debug("Interpolating melting temperature for %d elements", len(elements))
    values = [element.melting_temperature for element in elements]
    result = interpolate(values, composition)
    logger.info("Interpolated melting temperature: %.1f K", result)
    return result


def interpolate_boiling_temperature(elements: List[ChemicalElement], composition: List[float]) -> float:
    """
    Interpolates the boiling temperature based on the elements and their composition.
    Args:
        elements (list[ChemicalElement]): List of elements.
        composition (list[float]): List of composition percentages.
    Returns:
        float: Interpolated boiling temperature.
    """
    logger.debug("Interpolating boiling temperature for %d elements", len(elements))
    values = [element.boiling_temperature for element in elements]
    result = interpolate(values, composition)
    logger.info("Interpolated boiling temperature: %.1f K", result)
    return result
