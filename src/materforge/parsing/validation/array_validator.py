"""General data validation utilities."""

import logging
import numpy as np
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


def is_monotonic(arr: np.ndarray, name: str = "Array",
                 mode: str = "strictly_increasing",
                 threshold: float = ProcessingConstants.MONOTONICITY_THRESHOLD,
                 raise_error: bool = True) -> bool:
    """Universal monotonicity checker supporting multiple modes."""
    for i in range(1, len(arr)):
        diff = arr[i] - arr[i - 1]
        violation = False
        if mode == "strictly_increasing" and diff <= threshold:
            violation = True
        elif mode == "non_decreasing" and diff < -threshold:
            violation = True
        elif mode == "strictly_decreasing" and diff >= -threshold:
            violation = True
        elif mode == "non_increasing" and diff > threshold:
            violation = True
        if violation:
            start_idx = max(0, i - 2)
            end_idx = min(len(arr), i + 3)
            context = "\nSurrounding values:\n"
            for j in range(start_idx, end_idx):
                context += f"Index {j}: {arr[j]:.10e}\n"
            error_msg = (
                f"{name} is not {mode.replace('_', ' ')} at index {i}:\n"
                f"Previous value ({i - 1}): {arr[i - 1]:.10e}\n"
                f"Current value ({i}): {arr[i]:.10e}\n"
                f"Difference: {diff:.10e}\n"
                f"{context}"
            )
            if raise_error:
                raise ValueError(error_msg)
            else:
                logger.warning(f"Warning: {error_msg}")
                return False
    logger.info(f"{name} is {mode.replace('_', ' ')}")
    return True
