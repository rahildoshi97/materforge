# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
import sympy as sp
from dataclasses import dataclass
from typing import List, Union, Optional

logger = logging.getLogger(__name__)


@dataclass
class ChemicalElement:
    name: str
    atomic_number: float
    atomic_mass: float
    melting_temperature: Optional[float] = None
    boiling_temperature: Optional[float] = None
    latent_heat_of_fusion: Optional[float] = None
    latent_heat_of_vaporization: Optional[float] = None


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


def interpolate_melting_temperature(elements: List[ChemicalElement], composition: List[float]) -> Optional[float]:
    logger.debug("Interpolating melting temperature for %d elements", len(elements))
    if any(e.melting_temperature is None for e in elements):
        logger.warning("One or more elements have no melting temperature data; cannot interpolate")
        return None
    values = [element.melting_temperature for element in elements]
    result = interpolate(values, composition)
    logger.info("Interpolated melting temperature: %.1f K", result)
    return result


def interpolate_boiling_temperature(elements: List[ChemicalElement], composition: List[float]) -> Optional[float]:
    logger.debug("Interpolating boiling temperature for %d elements", len(elements))
    if any(e.boiling_temperature is None for e in elements):
        logger.warning("One or more elements have no boiling temperature data; cannot interpolate")
        return None
    values = [element.boiling_temperature for element in elements]
    result = interpolate(values, composition)
    logger.info("Interpolated boiling temperature: %.1f K", result)
    return result

def validate_element_completeness(element: ChemicalElement, context: str = "pure_metal") -> None:
    """
    Validate that a ChemicalElement has all properties required for a given context.

    Args:
        element: The ChemicalElement to validate.
        context: Usage context — 'pure_metal' requires full phase data;
                 'alloy' only requires atomic properties.
    Raises:
        ValueError: If required properties are missing for the given context.
    """
    missing = []
    if context == "pure_metal":
        if element.melting_temperature is None:
            missing.append("melting_temperature")
        if element.boiling_temperature is None:
            missing.append("boiling_temperature")
        if element.latent_heat_of_fusion is None:
            missing.append("latent_heat_of_fusion")
        if element.latent_heat_of_vaporization is None:
            missing.append("latent_heat_of_vaporization")
    if missing:
        raise ValueError(
            f"Element '{element.name}' ({element.atomic_number}) is missing required properties "
            f"for context '{context}': {missing}. "
            f"This element cannot be used as a pure metal in MaterForge. "
            f"Check mendeleev data: https://mendeleev.readthedocs.io/en/stable/data.html"
        )
