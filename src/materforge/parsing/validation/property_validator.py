"""Property-specific validation functions."""

import logging
import numpy as np
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


def validate_monotonic_energy_density(prop_name: str, temp_array: np.ndarray,
                                      prop_array: np.ndarray,
                                      tolerance: float = ProcessingConstants.DEFAULT_TOLERANCE) -> None:
    """Validate that energy_density is monotonically non-decreasing with temperature."""
    if prop_name != 'energy_density':
        return
    validate_monotonic_property(prop_name, temp_array, prop_array,
                                mode="strictly_increasing", tolerance=tolerance)


def validate_monotonic_property(prop_name: str, temp_array: np.ndarray,
                                prop_array: np.ndarray, mode: str = "strictly_increasing",
                                tolerance: float = ProcessingConstants.DEFAULT_TOLERANCE) -> None:
    """
    Generalized property monotonicity validation with enhanced error reporting.
    Args:
        prop_name: Name of the property being validated
        temp_array: Temperature data points
        prop_array: Property values corresponding to temperatures
        mode: Monotonicity mode ('strictly_increasing', 'non_decreasing',
              'strictly_decreasing', 'non_increasing')
        tolerance: Numerical tolerance for comparisons
    """
    from materforge.parsing.validation.array_validator import is_monotonic
    try:
        is_monotonic(prop_array, f"Property '{prop_name}'", mode, tolerance, raise_error=True)
    except ValueError as e:
        # Property-specific context
        enhanced_msg = (
            f"Property '{prop_name}' violates {mode.replace('_', ' ')} constraint with temperature.\n"
            f"Temperature range: {np.min(temp_array):.2f}K - {np.max(temp_array):.2f}K\n"
            f"Property range: {np.min(prop_array):.6e} - {np.max(prop_array):.6e}\n"
            f"Validation details: {str(e)}"
        )
        raise ValueError(enhanced_msg) from e
