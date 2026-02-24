# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import sympy as sp

from materforge.core.elements import ChemicalElement, interpolate_atomic_mass, interpolate_atomic_number
from materforge.core.exceptions import MaterialCompositionError, MaterialTemperatureError

logger = logging.getLogger(__name__)

# Structural identity fields (never tracked as material properties)
_CORE_FIELDS = frozenset({'name', 'material_type', 'elements', 'composition',
    'melting_temperature', 'boiling_temperature',
    'solidus_temperature', 'liquidus_temperature',
    'initial_boiling_temperature', 'final_boiling_temperature',
    'atomic_number', 'atomic_mass', '_dynamic_properties',
})


@dataclass
class Material:
    """
    Represents a material with its composition and thermophysical properties.

    Supports pure metals, alloys, and future material types (ceramics, polymers, etc.).
    Thermophysical properties are assigned dynamically via setattr and tracked
    automatically - no predefined property fields required.
    """
    # Basic material information
    name: str
    material_type: str
    elements: List[ChemicalElement]
    composition: Union[np.ndarray, List[float], Tuple]
    # Temperature bounds - pure metals
    melting_temperature: Optional[sp.Float] = None
    boiling_temperature: Optional[sp.Float] = None
    # Temperature bounds - alloys
    solidus_temperature: Optional[sp.Float] = None
    liquidus_temperature: Optional[sp.Float] = None
    initial_boiling_temperature: Optional[sp.Float] = None
    final_boiling_temperature: Optional[sp.Float] = None
    # Derived — computed in __post_init__, not user-supplied
    atomic_number: Optional[float] = field(default=None, init=False)
    atomic_mass: Optional[float] = field(default=None, init=False)
    # Dynamic property tracker - automatically populated by __setattr__
    _dynamic_properties: Set[str] = field(default_factory=set, init=False, repr=False)
    # Validation bounds - pure metals (K)
    MIN_MELTING_TEMP = 302.0    # Cs
    MAX_MELTING_TEMP = 3695.0   # W
    MIN_BOILING_TEMP = 630.0    # Hg
    MAX_BOILING_TEMP = 6203.0   # W
    # Validation bounds - alloys (K)
    MIN_SOLIDUS_TEMP = 250.0
    MAX_SOLIDUS_TEMP = 2000.0
    MIN_LIQUIDUS_TEMP = 300.0
    MAX_LIQUIDUS_TEMP = 2200.0

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

    def property_names(self) -> List[str]:
        """Return all dynamically assigned thermophysical property names."""
        return list(self._dynamic_properties)

    @classmethod
    def property_field_names(cls) -> List[str]:
        """Deprecated. Use instance method property_names() instead."""
        import warnings
        warnings.warn(
            "property_field_names() is deprecated. Use instance method property_names().",
            DeprecationWarning,
            stacklevel=2,
        )
        return []

    # --- Initialisation ---
    def __post_init__(self) -> None:
        logger.info("Initializing material: %s (type: %s)", self.name, self.material_type)
        self._validate_composition()
        self._validate_temperatures()
        self._calculate_properties()

    # --- Validation ---
    def _validate_composition(self) -> None:
        logger.debug("Validating composition for material: %s", self.name)
        if not self.elements:
            raise ValueError("Elements list cannot be empty")
        if len(self.elements) != len(self.composition):
            raise ValueError(
                f"Number of elements ({len(self.elements)}) must match composition length ({len(self.composition)})")
        if not np.isclose(sum(self.composition), 1.0, atol=1e-10):
            raise MaterialCompositionError(f"The sum of the composition array must be 1.0, got {sum(self.composition)}")
        if self.material_type == 'pure_metal' and len(self.elements) != 1:
            raise MaterialCompositionError(f"Pure metals must have exactly 1 element, got {len(self.elements)}")
        if self.material_type == 'alloy' and len(self.elements) < 2:
            raise MaterialCompositionError(f"Alloys must have at least 2 elements, got {len(self.elements)}")

    def _validate_temperatures(self) -> None:
        logger.debug(f"Starting temperature validation for {self.material_type}: {self.name}")
        try:
            if self.material_type == 'pure_metal':
                self._validate_pure_metal_temperatures()
            elif self.material_type == 'alloy':
                self._validate_alloy_temperatures()
            else:
                # Future material types: no structural temperature validation yet
                logger.warning(
                    "No temperature validation implemented for material_type '%s'. Skipping.",
                    self.material_type)
            logger.debug(f"Temperature validation passed for {self.material_type}: {self.name}")
        except MaterialTemperatureError as e:
            logger.error(f"Temperature validation failed for {self.name}: {e}")
            raise

    def _validate_pure_metal_temperatures(self) -> None:
        if self.melting_temperature is None:
            raise MaterialTemperatureError("Pure metals must specify melting_temperature")
        if self.boiling_temperature is None:
            raise MaterialTemperatureError("Pure metals must specify boiling_temperature")
        if not isinstance(self.melting_temperature, sp.Float):
            raise MaterialTemperatureError(
                f"melting_temperature must be a SymPy Float, got {type(self.melting_temperature).__name__}"
            )
        if not isinstance(self.boiling_temperature, sp.Float):
            raise MaterialTemperatureError(
                f"boiling_temperature must be a SymPy Float, got {type(self.boiling_temperature).__name__}"
            )
        # Validate temperature ranges
        self._validate_temperature_value(
            self.melting_temperature, "melting_temperature",
            self.MIN_MELTING_TEMP, self.MAX_MELTING_TEMP)
        self._validate_temperature_value(
            self.boiling_temperature, "boiling_temperature",
            self.MIN_BOILING_TEMP, self.MAX_BOILING_TEMP)
        if float(self.melting_temperature) >= float(self.boiling_temperature):
            raise MaterialTemperatureError(
                f"melting_temperature ({float(self.melting_temperature)}K) must be less than "
                f"boiling_temperature ({float(self.boiling_temperature)}K)")

    def _validate_alloy_temperatures(self) -> None:
        required_temps = [
            ('solidus_temperature', self.solidus_temperature),
            ('liquidus_temperature', self.liquidus_temperature),
            ('initial_boiling_temperature', self.initial_boiling_temperature),
            ('final_boiling_temperature', self.final_boiling_temperature),
        ]
        missing_temps = [name for name, temp in required_temps if temp is None]
        if missing_temps:
            raise MaterialTemperatureError(
                f"Alloys must specify all temperature properties. "
                f"Missing: {', '.join(missing_temps)}")
        for temp_name, temp_val in required_temps:
            if temp_val is not None and not isinstance(temp_val, sp.Float):
                raise MaterialTemperatureError(
                    f"{temp_name} must be a SymPy Float, got {type(temp_val).__name__}")
        self._validate_temperature_value(
            self.solidus_temperature, "solidus_temperature",
            self.MIN_SOLIDUS_TEMP, self.MAX_SOLIDUS_TEMP)
        self._validate_temperature_value(
            self.liquidus_temperature, "liquidus_temperature",
            self.MIN_LIQUIDUS_TEMP, self.MAX_LIQUIDUS_TEMP)
        self._validate_temperature_value(
            self.initial_boiling_temperature, "initial_boiling_temperature",
            self.MIN_BOILING_TEMP, self.MAX_BOILING_TEMP)
        self._validate_temperature_value(
            self.final_boiling_temperature, "final_boiling_temperature",
            self.MIN_BOILING_TEMP, self.MAX_BOILING_TEMP)
        if float(self.solidus_temperature) > float(self.liquidus_temperature):
            raise MaterialTemperatureError(
                f"solidus_temperature ({float(self.solidus_temperature)}K) must be <= "
                f"liquidus_temperature ({float(self.liquidus_temperature)}K)")
        if float(self.initial_boiling_temperature) > float(self.final_boiling_temperature):
            raise MaterialTemperatureError(
                f"initial_boiling_temperature ({float(self.initial_boiling_temperature)}K) must be <= "
                f"final_boiling_temperature ({float(self.final_boiling_temperature)}K)")
        if float(self.liquidus_temperature) >= float(self.initial_boiling_temperature):
            raise MaterialTemperatureError(
                f"liquidus_temperature ({float(self.liquidus_temperature)}K) must be < "
                f"initial_boiling_temperature ({float(self.initial_boiling_temperature)}K)")

    @staticmethod
    def _validate_temperature_value(temperature: Union[float, sp.Float], temp_name: str,
                                    min_temp: Optional[float] = None, max_temp: Optional[float] = None) -> None:
        if temperature is None:
            raise MaterialTemperatureError(f"{temp_name} cannot be None")
        temp_val = float(temperature)
        from materforge.data.constants import PhysicalConstants
        if temp_val <= PhysicalConstants.ABSOLUTE_ZERO:
            raise MaterialTemperatureError(f"{temp_name} must be above absolute zero, got {temp_val}K")
        if min_temp is not None and temp_val < min_temp:
            raise MaterialTemperatureError(f"{temp_name} {temp_val}K is below minimum allowed value ({min_temp}K)")
        if max_temp is not None and temp_val > max_temp:
            raise MaterialTemperatureError(f"{temp_name} {temp_val}K is above maximum allowed value ({max_temp}K)")

    # --- Derived property calculation ---
    def _calculate_properties(self) -> None:
        logger.debug("Calculating interpolated properties for %s", self.name)
        if self.material_type == 'pure_metal':
            element = self.elements[0]
            self.atomic_number = float(element.atomic_number)
            self.atomic_mass = float(element.atomic_mass)
            logger.debug("Pure metal - atomic_number: %.1f, atomic_mass: %.3f", self.atomic_number, self.atomic_mass)
        elif self.material_type == 'alloy':
            self.atomic_number = interpolate_atomic_number(self.elements, self.composition)
            self.atomic_mass = interpolate_atomic_mass(self.elements, self.composition)
            logger.debug("Alloy - atomic_number: %.3f, atomic_mass: %.3f", self.atomic_number, self.atomic_mass)
        else:
            # Future material types: atomic properties unavailable without element_map support
            logger.warning("No atomic property calculation for material_type '%s'. "
                "atomic_number and atomic_mass will be None.", self.material_type)

    # --- Public interface ---
    def solidification_interval(self) -> Tuple[sp.Float, sp.Float]:
        if self.material_type != 'alloy':
            raise ValueError("Solidification interval is only applicable to alloys")
        return self.solidus_temperature, self.liquidus_temperature

    def melting_point(self) -> sp.Float:
        if self.material_type == 'pure_metal':
            return self.melting_temperature
        return self.solidus_temperature

    def boiling_point(self) -> sp.Float:
        if self.material_type == 'pure_metal':
            return self.boiling_temperature
        return self.initial_boiling_temperature

    def evaluate_properties_at_temperature(self, temperature: Union[float, int],
                                           properties: Optional[List[str]] = None,
                                           include_constants: bool = True) -> Dict[str, float]:
        logger.info("Evaluating properties at T=%.1f K for material: %s", temperature, self.name)
        if not isinstance(temperature, (int, float)):
            raise ValueError(f"Temperature must be numeric, got {type(temperature).__name__}")
        if temperature <= 0:
            raise ValueError(f"Temperature must be positive, got {temperature}")
        # Use dynamic tracker - works for any property name including user-defined ones
        all_properties = {name: getattr(self, name) for name in self.property_names()}
        existing_properties = {k: v for k, v in all_properties.items() if v is not None}
        if properties is not None:
            if not isinstance(properties, list):
                raise ValueError("Properties must be a list of strings")
            invalid_props = set(properties) - set(existing_properties.keys())
            if invalid_props:
                raise ValueError(f"Invalid properties: {invalid_props}. Available: {list(existing_properties.keys())}")
            existing_properties = {k: v for k, v in existing_properties.items() if k in properties}
        logger.debug("Evaluating %d properties: %s", len(existing_properties), list(existing_properties.keys()))
        # Find the temperature symbol used across all symbolic properties
        temp_symbol = None
        for prop_value in existing_properties.values():
            if hasattr(prop_value, 'free_symbols') and prop_value.free_symbols:
                temp_symbol = list(prop_value.free_symbols)[0]
                break
        results = {}
        for prop_name, prop_value in existing_properties.items():
            try:
                if hasattr(prop_value, 'free_symbols') and prop_value.free_symbols:
                    evaluated = prop_value.subs(temp_symbol, temperature) if temp_symbol else prop_value
                    results[prop_name] = float(evaluated.evalf() if hasattr(evaluated, 'evalf') else evaluated)
                elif isinstance(prop_value, (sp.Float, sp.Integer)):
                    if include_constants:
                        results[prop_name] = float(prop_value)
                elif isinstance(prop_value, (int, float)):
                    if include_constants:
                        results[prop_name] = float(prop_value)
                else:
                    logger.warning("Unknown property type for '%s': %s", prop_name, type(prop_value))
            except Exception as e:
                logger.error("Failed to evaluate property '%s': %s", prop_name, e)
                results[prop_name] = None
        results = {k: v for k, v in results.items() if v is not None}
        logger.info("Successfully evaluated %d properties at T=%.1f K", len(results), temperature)
        return results

    # --- Repr ---
    def __str__(self) -> str:
        element_names = [e.name for e in self.elements]
        if self.material_type == 'pure_metal':
            return f"Pure Metal: {self.name} ({element_names[0]})"
        composition_str = ", ".join(f"{elem}: {comp:.3f}" for elem, comp in zip(element_names, self.composition))
        return f"Alloy: {self.name} ({composition_str})"

    def __repr__(self) -> str:
        return (f"Material(name='{self.name}', material_type='{self.material_type}', "
                f"elements={len(self.elements)}, composition={self.composition})")
