import logging
import numpy as np
from typing import Tuple

from materforge.parsing.config.yaml_keys import CONSTANT_KEY
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


def interpolate_value(T: float, x_array: np.ndarray, y_array: np.ndarray,
                      lower_bound_type: str, upper_bound_type: str) -> float:
    """Interpolate a value at temperature T using the provided data arrays."""
    # Input validation
    if not np.isfinite(T):
        raise ValueError(f"Temperature T must be finite, got {T}")
    if len(x_array) == 0 or len(y_array) == 0:
        raise ValueError("Input arrays cannot be empty")
    if len(x_array) != len(y_array):
        raise ValueError(f"Array length mismatch: x_array({len(x_array)}) != y_array({len(y_array)})")
    logger.debug("Interpolating value at T=%.1f with bounds: lower=%s, upper=%s",
                 T, lower_bound_type, upper_bound_type)
    logger.debug("Data range: T∈[%.1f, %.1f], y∈[%.3e, %.3e]",
                 x_array[0], x_array[-1], np.min(y_array), np.max(y_array))
    try:
        # Handle single-point arrays
        if len(x_array) == 1:
            logger.debug("Single-point array: returning constant value %.6f", y_array[0])
            return float(y_array[0])
        if T < x_array[0]:
            logger.debug("T below data range, applying lower bound: %s", lower_bound_type)
            if lower_bound_type == CONSTANT_KEY:
                result = float(y_array[0])
                logger.debug("Lower constant extrapolation: %.6f", result)
                return result
            else:  # 'extrapolate'
                denominator = x_array[1] - x_array[0]
                if denominator == 0:
                    logger.error("Cannot extrapolate: first two temperature values are equal")
                    raise ValueError("Cannot extrapolate: first two temperature values are equal.")
                slope = (y_array[1] - y_array[0]) / denominator
                result = float(y_array[0] + slope * (T - x_array[0]))
                logger.debug("Lower linear extrapolation: slope=%.6f, result=%.6f", slope, result)
                return result
        elif T >= x_array[-1]:
            logger.debug("T above data range, applying upper bound: %s", upper_bound_type)
            if upper_bound_type == CONSTANT_KEY:
                result = float(y_array[-1])
                logger.debug("Upper constant extrapolation: %.6f", result)
                return result
            else:  # 'extrapolate'
                denominator = x_array[-1] - x_array[-2]
                if denominator == 0:
                    logger.error("Cannot extrapolate: last two temperature values are equal")
                    raise ValueError("Cannot extrapolate: last two temperature values are equal.")
                slope = (y_array[-1] - y_array[-2]) / denominator
                result = float(y_array[-1] + slope * (T - x_array[-1]))
                logger.debug("Upper linear extrapolation: slope=%.6f, result=%.6f", slope, result)
                return result
        else:
            logger.debug("T within data range, using linear interpolation")
            result = float(np.interp(T, x_array, y_array))
            logger.debug("Linear interpolation result: %.6f", result)
            return result
    except Exception as e:
        logger.error("Interpolation failed at T=%.1f: %s", T, e, exc_info=True)
        raise ValueError(f"Interpolation failed at T={T}: {str(e)}") from e


def ensure_ascending_order(temp_array: np.ndarray, *value_arrays: np.ndarray) -> Tuple[np.ndarray, ...]:
    """Ensure temperature array is in ascending order, flipping all provided arrays if needed."""
    logger.debug("Checking array order for %d arrays (temp + %d value arrays)",
                 len(value_arrays) + 1, len(value_arrays))
    if len(temp_array) < 2:
        logger.debug("Array too short for order check: length=%d", len(temp_array))
        return (temp_array,) + value_arrays
    try:
        diffs = np.diff(temp_array)
        # Use tolerance for floating-point comparisons
        tolerance = ProcessingConstants.FLOATING_POINT_TOLERANCE
        ascending_count = np.sum(diffs > tolerance)
        descending_count = np.sum(diffs < -tolerance)
        constant_count = np.sum(np.abs(diffs) <= tolerance)
        logger.debug("Array differences: %d ascending, %d descending, %d constant",
                     ascending_count, descending_count, constant_count)
        # Check for strictly ascending (with tolerance)
        if ascending_count == len(diffs) or (ascending_count > 0 and constant_count == len(diffs) - ascending_count):
            logger.debug("Array is ascending (within tolerance)")
            return (temp_array,) + value_arrays
        # Check for strictly descending (with tolerance)
        elif descending_count == len(diffs) or (descending_count > 0 and constant_count == len(diffs) - descending_count):
            logger.debug("Array is descending, flipping all arrays")
            flipped_temp = np.flip(temp_array)
            flipped_values = tuple(np.flip(arr) for arr in value_arrays)
            logger.debug("Flipped arrays: temp range [%.1f, %.1f] -> [%.1f, %.1f]",
                         temp_array[0], temp_array[-1], flipped_temp[0], flipped_temp[-1])
            return (flipped_temp,) + flipped_values
        else:
            logger.error("Array is neither strictly ascending nor descending")
            logger.error("Temperature array: %s", temp_array.tolist() if len(temp_array) <= 20
                         else f"[{temp_array[0]}, ..., {temp_array[-1]}] (length={len(temp_array)})")
            raise ValueError(f"Array is not strictly ascending or strictly descending: {temp_array}")
    except Exception as e:
        logger.error("Error checking array order: %s", e, exc_info=True)
        raise ValueError(f"Failed to ensure ascending order: {str(e)}") from e
