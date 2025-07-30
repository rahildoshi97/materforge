"""
Material data, constants, and element definitions.

This package provides access to physical constants, processing constants,
chemical element data, and material property databases used throughout PyMatLib.
"""

from .constants.physical_constants import PhysicalConstants
from .constants.processing_constants import ProcessingConstants, ErrorMessages, FileConstants
from .elements.element_data import (
    element_map,
    CARBON, NITROGEN, ALUMINIUM, SILICON, PHOSPHORUS, SULFUR,
    TITANIUM, VANADIUM, CHROMIUM, MANGANESE, IRON, NICKEL, MOLYBDENUM
)

__all__ = [
    "PhysicalConstants",
    "ProcessingConstants",
    "ErrorMessages",
    "FileConstants",
    "element_map",
    # Individual elements for direct access
    "CARBON", "NITROGEN", "ALUMINIUM", "SILICON", "PHOSPHORUS", "SULFUR",
    "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE", "IRON", "NICKEL", "MOLYBDENUM"
]
