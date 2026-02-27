# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Property-specific validation functions."""

import logging
import numpy as np
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


def validate_monotonic_property(prop_name: str, dep_array: np.ndarray,
                                prop_array: np.ndarray, mode: str = "strictly_increasing",
                                tolerance: float = ProcessingConstants.DEFAULT_TOLERANCE) -> None:
    """Validates that a property array satisfies a monotonicity constraint.

    Args:
        prop_name:  Name of the property (used in error messages only).
        dep_array:  Dependency variable values.
        prop_array: Corresponding property values to validate.
        mode:       One of 'strictly_increasing', 'non_decreasing',
                    'strictly_decreasing', 'non_increasing'.
        tolerance:  Numerical tolerance for comparisons.
    Raises:
        ValueError: If the property array violates the specified monotonicity constraint.
    """
    from materforge.parsing.validation.array_validator import is_monotonic
    try:
        is_monotonic(prop_array, f"Property '{prop_name}'", mode, tolerance, raise_error=True)
    except ValueError as e:
        raise ValueError(
            f"Property '{prop_name}' violates {mode.replace('_', ' ')} constraint.\n"
            f"Dependency range: {np.min(dep_array):.6e} - {np.max(dep_array):.6e}\n"
            f"Property range:   {np.min(prop_array):.6e} - {np.max(prop_array):.6e}\n"
            f"Validation details: {str(e)}"
        ) from e
