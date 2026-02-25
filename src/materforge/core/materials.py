# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Set, Union

import sympy as sp

logger = logging.getLogger(__name__)

_CORE_FIELDS = frozenset({'name', '_dynamic_properties'})


@dataclass
class Material:
    """Generic material container with fully dynamic property tracking.
    All thermophysical properties are assigned dynamically
    via setattr and tracked automatically. No material_type, composition, or
    elements required.

    Attributes:
        name: Human-readable material identifier.
    """
    name: str
    _dynamic_properties: Set[str] = field(default_factory=set, init=False, repr=False)
    # --- Dynamic property tracking ---
    def __setattr__(self, name: str, value) -> None:
        object.__setattr__(self, name, value)
        # Auto-track anything that is not a structural field
        if name not in _CORE_FIELDS and not name.startswith('_'):
            try:
                object.__getattribute__(self, '_dynamic_properties').add(name)
            except AttributeError:
                # Fires during dataclass __init__ before _dynamic_properties exists
                pass

    def __getattr__(self, name: str) -> Any:
        # Only fires when normal lookup fails - dynamic properties won't hit this
        # since they are in __dict__ via __setattr__
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def property_names(self) -> set:
        """Returns all dynamically assigned property names."""
        return set(self._dynamic_properties)

    # --- Property evaluation ---
    def evaluate_properties_at_temperature(self, symbol: sp.Symbol, value: Union[float, int]) -> Dict[str, float]:
        """Evaluates all properties by substituting symbol=value.

        Args:
            symbol: SymPy symbol to substitute (e.g. sp.Symbol('T')).
            value:  Value to substitute.
        Returns:
            All property names mapped to evaluated float values.
            Properties that fail evaluation are silently excluded.
        Raises:
            ValueError: If symbol is not sp.Symbol, or value is non-numeric.
        """
        if not isinstance(symbol, sp.Symbol):
            raise ValueError(f"symbol must be sp.Symbol, got {type(symbol).__name__}")
        if not isinstance(value, (int, float)):
            raise ValueError(f"value must be numeric, got {type(value).__name__}")
        logger.info("Evaluating '%s' at %s=%.2f", self.name, symbol, value)
        results: Dict[str, float] = {}
        for prop_name in self.property_names():
            prop_value = getattr(self, prop_name)
            if prop_value is None:
                continue
            try:
                if hasattr(prop_value, 'free_symbols') and prop_value.free_symbols:
                    substituted = prop_value.subs(symbol, value).doit()
                    if hasattr(substituted, 'free_symbols') and substituted.free_symbols:
                        raise ValueError(f"symbol '{symbol}' not found in expression; "
                            f"expression requires: {sorted(str(s) for s in prop_value.free_symbols)}")
                    results[prop_name] = float(substituted.evalf())
                elif isinstance(prop_value, (sp.Float, sp.Integer, int, float)):
                    results[prop_name] = float(prop_value)
                else:
                    logger.warning("Unrecognised type for '%s': %s - skipping",
                                prop_name, type(prop_value).__name__)
            except Exception as e:
                logger.error("Failed to evaluate '%s': %s", prop_name, e)
        logger.info("Evaluated %d/%d properties at %s=%.2f",
                    len(results), len(self.property_names()), symbol, value)
        return results

    # --- Repr ---
    def __str__(self) -> str:
        return f"Material: {self.name} ({len(self._dynamic_properties)} properties)"

    def __repr__(self) -> str:
        return (f"Material(name='{self.name}', "
                f"properties={sorted(self.property_names())})")
