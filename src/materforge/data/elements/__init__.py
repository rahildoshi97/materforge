"""Chemical element data and definitions."""

from .element_data import (
    element_map, get_element,
    CARBON, NITROGEN, ALUMINIUM, SILICON, PHOSPHORUS, SULFUR,
    TITANIUM, VANADIUM, CHROMIUM, MANGANESE, IRON, NICKEL, MOLYBDENUM
)

__all__ = [
    "element_map", "get_element",
    # Individual elements for direct access
    "CARBON", "NITROGEN", "ALUMINIUM", "SILICON", "PHOSPHORUS", "SULFUR",
    "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE", "IRON", "NICKEL", "MOLYBDENUM"
]
