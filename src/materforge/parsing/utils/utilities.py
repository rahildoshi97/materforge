import logging
from typing import List, Tuple, Union

import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


# --- Core Utility Functions ---
def handle_numeric_temperature(processor_instance, material: Material,
                               prop_name: str, piecewise_expr: sp.Expr,
                               T: Union[float, sp.Symbol]) -> bool:
    """
    Handle numeric temperature evaluation with consistent error handling.
    Returns:
        bool: True if numeric evaluation was performed, False if symbolic
    """
    if isinstance(T, sp.Symbol):
        return False
    try:
        T_standard = sp.Symbol('T')
        value = float(piecewise_expr.subs(T_standard, T).evalf())
        setattr(material, prop_name, sp.Float(value))
        processor_instance.processed_properties.add(prop_name)
        logger.debug(f"Numeric evaluation completed for '{prop_name}': {value}")
        return True
    except Exception as e:
        raise ValueError(f"Failed to evaluate {prop_name} at T={T}: {str(e)}")


def create_step_visualization_data(transition_temp: float, val_array: List[float],
                                   temp_range: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Create proper step function visualization data."""
    margin = (np.max(temp_range) - np.min(temp_range)) * ProcessingConstants.TEMPERATURE_PADDING_FACTOR
    epsilon = ProcessingConstants.TEMPERATURE_EPSILON
    x_data = np.array([
        np.min(temp_range) - margin,
        transition_temp - epsilon,
        transition_temp,
        transition_temp + epsilon,
        np.max(temp_range) + margin
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
    """Ensure value is compatible with SymPy operations."""
    if hasattr(value, 'item'):  # NumPy scalar
        return float(value.item())
    elif isinstance(value, (np.float64, np.int64, np.float32, np.int32, np.number)):
        return float(value)
    elif isinstance(value, (list, np.ndarray)):
        return [float(x) for x in value]
    else:
        return float(value)
