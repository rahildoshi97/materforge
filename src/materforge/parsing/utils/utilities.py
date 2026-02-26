# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import List, Tuple, Union
import numpy as np
import sympy as sp
from materforge.core.materials import Material
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)

# --- Core Utility Functions ---
def handle_numeric_dependency(processor_instance, material: Material,
                              prop_name: str, piecewise_expr: sp.Expr,
                              dependency: Union[float, sp.Symbol]) -> bool:
    """Evaluates a piecewise expression at a numeric dependency value and assigns
    the result to the material.

    Args:
        processor_instance: Processor instance owning processed_properties tracking.
        material:           Material to assign the evaluated value to.
        prop_name:          Name of the property being evaluated.
        piecewise_expr:     SymPy expression to evaluate.
        dependency:         Dependency value or symbol. If sp.Symbol, returns False
                            and no evaluation is performed.
    Returns:
        True if numeric evaluation was performed, False if dependency is symbolic.
    Raises:
        ValueError: If numeric evaluation fails.
    """
    if isinstance(dependency, sp.Symbol):
        return False
    try:
        placeholder = sp.Symbol('T')
        value = float(piecewise_expr.subs(placeholder, dependency).evalf())
        setattr(material, prop_name, sp.Float(value))
        processor_instance.processed_properties.add(prop_name)
        logger.debug("Numeric evaluation completed for '%s': %s", prop_name, value)
        return True
    except Exception as e:
        raise ValueError(f"Failed to evaluate '{prop_name}' at dependency={dependency}: {str(e)}") from e

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
