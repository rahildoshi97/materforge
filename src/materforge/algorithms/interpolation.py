# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
import numpy as np
from typing import Tuple

from materforge.parsing.config.yaml_keys import LINEAR_KEY
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


def interpolate_value(dep_value, dep_array: np.ndarray, prop_array: np.ndarray,
                      lower_bound_type: str, upper_bound_type: str):
    """Interpolates a property value at dep_value using the provided data arrays.

    Args:
        dep_value:        Dependency value at which to interpolate.
        dep_array:        Dependency axis data (must be sorted ascending).
        prop_array:       Corresponding property values.
        lower_bound_type: Behaviour below dep_array[0]: CONSTANT_KEY or LINEAR_KEY.
        upper_bound_type: Behaviour above dep_array[-1]: CONSTANT_KEY or LINEAR_KEY.
    Returns:
        Interpolated value.
    Raises:
        ValueError: If dep_value is non-finite, arrays are empty or mismatched,
                    or extrapolation is requested on a degenerate interval.
    """
    if not np.isfinite(dep_value):
        raise ValueError(f"dep_value must be finite, got {dep_value}")
    if len(dep_array) == 0 or len(prop_array) == 0:
        raise ValueError("Input arrays cannot be empty")
    if len(dep_array) != len(prop_array):
        raise ValueError(f"Array length mismatch: dep_array({len(dep_array)}) != prop_array({len(prop_array)})")

    logger.debug("Interpolating at dep_value=%.4g (lower=%s, upper=%s)",
                 dep_value, lower_bound_type, upper_bound_type)
    logger.debug("Data range: dep=[%.4g, %.4g], prop=[%.3e, %.3e]",
                 dep_array[0], dep_array[-1], np.min(prop_array), np.max(prop_array))

    try:
        if dep_value < dep_array[0] and lower_bound_type == LINEAR_KEY:
            logger.debug("dep_value below data range - applying lower bound: %s", lower_bound_type)
            if len(dep_array) < 2:
                return prop_array[0]
            denom = dep_array[1] - dep_array[0]
            if denom == 0:
                raise ValueError("Cannot extrapolate: first two dependency values are equal")
            slope = (prop_array[1] - prop_array[0]) / denom
            return prop_array[0] + slope * (dep_value - dep_array[0])
        if dep_value >= dep_array[-1] and upper_bound_type == LINEAR_KEY:
            logger.debug("dep_value above data range - applying upper bound: %s", upper_bound_type)
            if len(dep_array) < 2:
                return prop_array[-1]
            denom = dep_array[-1] - dep_array[-2]
            if denom == 0:
                raise ValueError("Cannot extrapolate: last two dependency values are equal")
            slope = (prop_array[-1] - prop_array[-2]) / denom
            return float(prop_array[-1] + slope * (dep_value - dep_array[-1]))
        return float(np.interp(dep_value, dep_array, prop_array))
    except Exception as e:
        raise ValueError(f"Interpolation failed at dep_value={dep_value}: {str(e)}") from e


def ensure_ascending_order(dep_array: np.ndarray, *value_arrays: np.ndarray) -> Tuple[np.ndarray, ...]:
    """Ensures dep_array is in ascending order, flipping all arrays together if needed.

    Args:
        dep_array:     Dependency axis data.
        *value_arrays: Any number of paired value arrays to flip alongside dep_array.
    Returns:
        Tuple of (dep_array, *value_arrays), all in ascending dep order.
    Raises:
        ValueError: If the array is neither strictly ascending nor strictly descending.
    """
    if len(dep_array) < 2:
        return (dep_array,) + value_arrays
    try:
        diffs = np.diff(dep_array)
        tol = ProcessingConstants.FLOATING_POINT_TOLERANCE
        n = len(diffs)
        ascending_count  = int(np.sum(diffs >  tol))
        descending_count = int(np.sum(diffs < -tol))
        constant_count   = int(np.sum(np.abs(diffs) <= tol))
        logger.debug("Array diffs: %d ascending, %d descending, %d constant",
                     ascending_count, descending_count, constant_count)
        is_ascending  = ascending_count  == n or (ascending_count  > 0 and constant_count == n - ascending_count)
        is_descending = descending_count == n or (descending_count > 0 and constant_count == n - descending_count)
        if is_ascending:
            return (dep_array,) + value_arrays
        if is_descending:
            logger.debug("Array is descending - flipping all arrays")
            return (np.flip(dep_array),) + tuple(np.flip(arr) for arr in value_arrays)
        summary = (str(dep_array.tolist()) if len(dep_array) <= 20
                   else f"[{dep_array[0]}, ..., {dep_array[-1]}] (length={len(dep_array)})")
        raise ValueError(f"Array is neither strictly ascending nor descending: {summary}")
    except Exception as e:
        raise ValueError(f"Failed to ensure ascending order: {str(e)}") from e
