# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""
Core data structures and material definitions.

This module contains the fundamental classes and interfaces that define
materials, and the core abstractions used throughout the MaterForge library.
"""

from .materials import Material
from .symbol_registry import SymbolRegistry

__all__ = [
    "Material",
    "SymbolRegistry",
]
