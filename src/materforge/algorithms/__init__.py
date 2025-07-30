"""
Core computational algorithms for materials property processing.

This module provides mathematical algorithms for processing material properties,
including interpolation, regression, piecewise function construction, and
property inversion operations.
"""

from .interpolation import interpolate_value, ensure_ascending_order
from .regression_processor import RegressionProcessor
from .piecewise_builder import PiecewiseBuilder
from .piecewise_inverter import PiecewiseInverter

__all__ = [
    "interpolate_value",
    "ensure_ascending_order",
    "RegressionProcessor",
    "PiecewiseBuilder",
    "PiecewiseInverter"
]
