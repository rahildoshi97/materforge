# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""
Material data, constants, and element definitions.

This package provides access to physical constants, processing constants,
chemical element data, and material property databases used throughout MaterForge.
"""

from .constants.processing_constants import ProcessingConstants, ErrorMessages, FileConstants

__all__ = [
    "ProcessingConstants",
    "ErrorMessages",
    "FileConstants",
]
