# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
import warnings
from typing import Optional
from mendeleev import element as mendeleev_element
from materforge.data.constants import PhysicalConstants
from materforge.core.elements import ChemicalElement

logger = logging.getLogger(__name__)

# Silence mendeleev's allotrope UserWarnings globally (P, S, Se, Sn, etc.)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    module=r"mendeleev\.models",
)

def _build_element(symbol: str) -> Optional[ChemicalElement]:
    """
    Build a ChemicalElement from mendeleev data.

    Unit conversions applied:
        - atomic_weight: Da (= g/mol) → kg via PhysicalConstants.AMU
        - melting_point / boiling_point: K (used as-is)
        - fusion_heat / evaporation_heat: kJ/mol → J/kg
          Formula: J/kg = kJ/mol × 1e6 / atomic_weight[g/mol]
    Args:
        symbol: IUPAC chemical symbol (e.g. 'Fe', 'Al')
    Returns:
        ChemicalElement instance, or None if mendeleev data is incomplete.
    Sources:
        - Atomic weights: IUPAC 2021 via mendeleev
        - Thermophysical data: mendeleev database
          (see https://mendeleev.readthedocs.io/en/stable/data.html)
    """
    try:
        e = mendeleev_element(symbol)
        if e.atomic_weight is None:
            logger.warning("Element %s: missing atomic_weight, skipping", symbol)
            return None
        atomic_weight_g_per_mol = float(e.atomic_weight)  # g/mol = Da numerically
        atomic_mass_kg = atomic_weight_g_per_mol * PhysicalConstants.AMU
        melting_temperature = e.melting_point   # K
        boiling_temperature = e.boiling_point   # K
        def _to_j_per_kg(heat_kj_per_mol: Optional[float]) -> Optional[float]:
            """Convert latent heat from kJ/mol to J/kg."""
            if heat_kj_per_mol is None:
                return None
            return heat_kj_per_mol * 1e6 / atomic_weight_g_per_mol
        return ChemicalElement(
            name=e.name,
            atomic_number=e.atomic_number,
            atomic_mass=atomic_mass_kg,
            melting_temperature=melting_temperature,
            boiling_temperature=boiling_temperature,
            latent_heat_of_fusion=_to_j_per_kg(e.fusion_heat),
            latent_heat_of_vaporization=_to_j_per_kg(e.evaporation_heat),
        )
    except Exception as ex:
        logger.error("Failed to build element '%s' from mendeleev: %s", symbol, ex)
        return None

# ---------------------------------------------------------------------------
# Cache: populated lazily on first access via get_element()
# ---------------------------------------------------------------------------
element_map: dict[str, ChemicalElement] = {}

def get_element(symbol: str) -> ChemicalElement:
    """
    Get a ChemicalElement by its IUPAC symbol.

    Elements are built from mendeleev on first access and cached in
    element_map for all subsequent lookups. This keeps import time
    near-zero regardless of how many elements exist in the periodic table.
    Args:
        symbol: IUPAC chemical symbol (e.g. 'Fe', 'Al', 'Cr')
    Returns:
        ChemicalElement with properties sourced from mendeleev.
    Raises:
        KeyError: If the symbol is invalid or mendeleev data is incomplete.
    Examples:
        >>> iron = get_element('Fe')
        >>> iron.atomic_number
        26
    """
    symbol = symbol.strip()
    logger.debug("Looking up element: %s", symbol)
    if symbol in element_map:
        el = element_map[symbol]
        logger.debug("Cache hit: %s (Z=%d)", el.name, el.atomic_number)
        return el
    el = _build_element(symbol)
    if el is None:
        raise KeyError(
            f"Element '{symbol}' not found or has incomplete mendeleev data. "
            f"Currently cached elements: {list(element_map.keys())}"
        )
    element_map[symbol] = el
    logger.info(
        "Loaded element %s (Z=%d) from mendeleev. Cache size: %d",
        el.name, el.atomic_number, len(element_map),
    )
    return el
