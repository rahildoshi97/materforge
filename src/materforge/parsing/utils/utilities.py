# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import List, Tuple
import numpy as np

from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


# --- Core Utility Functions ---
def create_step_visualization_data(transition_point: float, val_array: List[float],
                                   dep_range: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Creates step function visualization data around a transition point.

    Args:
        transition_point: The dependency value at which the step occurs.
        val_array:        Two-element list [value_before, value_after].
        dep_range:        Full dependency range array for margin calculation.
    Returns:
        Tuple of (x_data, y_data) arrays suitable for plotting.
    """
    margin = (np.max(dep_range) - np.min(dep_range)) * ProcessingConstants.DEPENDENCY_PADDING_FACTOR
    epsilon = ProcessingConstants.DEPENDENCY_EPSILON
    x_data = np.array([
        np.min(dep_range) - margin,
        transition_point - epsilon,
        transition_point,
        transition_point + epsilon,
        np.max(dep_range) + margin
    ])
    y_data = np.array([
        val_array[0],
        val_array[0],
        val_array[0],
        val_array[1],
        val_array[1]
    ])
    return x_data, y_data


def ensure_sympy_compatible(value):
    """Converts a value to a SymPy-compatible Python scalar or list.

    Handles NumPy scalars and arrays that would otherwise cause type errors
    in SymPy expression construction.
    """
    if hasattr(value, 'item'):  # NumPy scalar
        return float(value.item())
    elif isinstance(value, (np.float64, np.int64, np.float32, np.int32, np.number)):
        return float(value)
    elif isinstance(value, (list, np.ndarray)):
        return [float(x) for x in value]
    else:
        return float(value)
