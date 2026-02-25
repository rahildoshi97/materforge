# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Union

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
    def evaluate_properties_at_temperature(self, temperature: Union[float, int],
                                           properties: Optional[List[str]] = None,
                                           include_constants: bool = True) -> Dict[str, float]:
        """Evaluates all (or selected) properties at a given temperature.
        Args:
            temperature: Temperature value in Kelvin. Must be positive.
            properties: Subset of property names to evaluate. Evaluates all if None.
            include_constants: Whether to include scalar/constant properties.
        Returns:
            Property names mapped to their evaluated float values.
            Properties that fail evaluation are silently excluded.
        Raises:
            ValueError: If temperature is non-numeric or non-positive, or if any
                requested property name does not exist on this material.
        """
        logger.info("Evaluating properties at T=%.2f K for '%s'", temperature, self.name)
        if not isinstance(temperature, (int, float)):
            raise ValueError(f"Temperature must be numeric, got {type(temperature).__name__}")
        if temperature <= 0:
            raise ValueError(f"Temperature must be positive, got {temperature}")
        existing = {n: getattr(self, n) for n in self.property_names() if getattr(self, n) is not None}
        if properties is not None:
            if not isinstance(properties, list):
                raise ValueError("'properties' must be a list of strings")
            invalid = set(properties) - set(existing)
            if invalid:
                raise ValueError(f"Unknown properties: {invalid}. Available: {sorted(existing)}")
            existing = {k: v for k, v in existing.items() if k in properties}
        logger.debug("Evaluating %d properties: %s", len(existing), sorted(existing))
        temp_symbol = None
        for v in existing.values():
            if hasattr(v, 'free_symbols') and v.free_symbols:
                temp_symbol = next(iter(v.free_symbols))
                break
        results: Dict[str, float] = {}
        for prop_name, prop_value in existing.items():
            try:
                if hasattr(prop_value, 'free_symbols') and prop_value.free_symbols:
                    evaluated = (prop_value.subs(temp_symbol, temperature) if temp_symbol else prop_value)
                    results[prop_name] = float(evaluated.evalf() if hasattr(evaluated, 'evalf') else evaluated)
                elif isinstance(prop_value, (sp.Float, sp.Integer, int, float)):
                    if include_constants:
                        results[prop_name] = float(prop_value)
                else:
                    logger.warning("Unrecognised type for property '%s': %s - skipping", prop_name, type(prop_value).__name__)
            except Exception as e:
                logger.error("Failed to evaluate property '%s': %s", prop_name, e)
        logger.info("Evaluated %d/%d properties at T=%.2f K", len(results), len(existing), temperature)
        return results

    # --- Repr ---
    def __str__(self) -> str:
        return f"Material: {self.name} ({len(self._dynamic_properties)} properties)"

    def __repr__(self) -> str:
        return (f"Material(name='{self.name}', "
                f"properties={sorted(self.property_names())})")
