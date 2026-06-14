# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Access to the example material files bundled with MaterForge.

MaterForge is built around YAML configs that users write themselves. The handful
of materials shipped with the package are reference examples, not a curated
database. These helpers load them by name so you don't have to dig the file out
of the installed package - handy for trying the API, demos, and tests.
"""

from __future__ import annotations

import logging
from importlib.resources import files
from pathlib import Path
from typing import Dict, List

import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.api import create_material

logger = logging.getLogger(__name__)

_MATERIALS_PACKAGE = "materforge.data.materials"
_YAML_SUFFIXES = (".yaml", ".yml")


def list_materials() -> List[str]:
    """Returns the names of the bundled example materials, sorted alphabetically."""
    return sorted(_discover_materials())


def get_material_path(name: str) -> Path:
    """Returns the path to a bundled material's YAML file.

    Matching is case-insensitive. Useful when you want to read or copy the YAML
    rather than load it directly.

    Args:
        name: Name of a bundled material (see :func:`list_materials`).
    Returns:
        Filesystem path to the material's YAML file.
    Raises:
        ValueError: If name does not match a bundled material; the message lists
            the available names.
    """
    catalog = _discover_materials()
    target = name.strip().casefold()
    for material_name, path in catalog.items():
        if material_name.casefold() == target:
            return path
    available = ", ".join(sorted(catalog)) or "(none found)"
    raise ValueError(f"Unknown material '{name}'. Available materials: {available}")


def load_material(name: str, dependency: sp.Symbol, enable_plotting: bool = False) -> Material:
    """Loads a bundled example material by name.

    A thin wrapper around :func:`materforge.create_material` for the YAML files
    shipped with the package. Plotting is off by default so loading an example
    never writes plot files next to the installed package.

    Args:
        name:            Name of a bundled material (see :func:`list_materials`).
        dependency:      SymPy symbol used as the independent variable, exactly
                         as in :func:`materforge.create_material`.
        enable_plotting: Generate visualisation plots (default: False).
    Returns:
        Fully initialised Material instance.
    Raises:
        ValueError: If name does not match a bundled material.
        TypeError:  If dependency is not a sympy Symbol.
    Example:
        >>> import sympy as sp
        >>> from materforge import load_material
        >>> steel = load_material('1.4301', sp.Symbol('T'))
    """
    path = get_material_path(name)
    logger.info("Loading bundled material '%s' from %s", name, path)
    return create_material(path, dependency, enable_plotting=enable_plotting)


def _discover_materials() -> Dict[str, Path]:
    """Maps each bundled material name to its YAML path.

    Names are the YAML file stems (e.g. 'Al', '1.4301'). Companion files such as
    the Excel tables referenced by FILE_IMPORT properties live alongside the YAML
    and are resolved relative to it at load time.
    """
    catalog: Dict[str, Path] = {}
    for entry in files(_MATERIALS_PACKAGE).iterdir():
        path = Path(str(entry))
        if path.suffix.lower() in _YAML_SUFFIXES:
            catalog[path.stem] = path
    return catalog
