"""
Base processor class for material property processing with common finalization logic.

This module provides the PropertyProcessorBase class that serves as a foundation
for all property processors in the MaterForge library, ensuring consistent
processing patterns and reducing code duplication.
"""

import logging
import numpy as np
import sympy as sp
from typing import Dict, Union, Tuple, Optional
from pathlib import Path

from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.core.materials import Material
from materforge.parsing.utils.utilities import handle_numeric_temperature

logger = logging.getLogger(__name__)


class PropertyProcessorBase:
    """
    Base class for material property processors with common finalization logic.

    This class provides shared functionality for processing material properties,
    including piecewise function creation, visualization, and property assignment.
    It serves as the foundation for all property handlers in the MaterForge library.
    Attributes:
        processed_properties (set): Set of property names that have been processed
        base_dir (Path): Base directory for file operations
        visualizer: Optional visualizer instance for property plotting
    """

    def __init__(self):
        """Initialize the base processor with empty state."""
        self.processed_properties = set()
        self.base_dir: Path = None # type: ignore
        self.visualizer = None
        logger.debug("PropertyProcessorBase initialized")

    def set_processing_context(self, base_dir: Path, visualizer, processed_properties: set):
        """Set processing context shared across handlers."""
        self.base_dir = base_dir
        self.visualizer = visualizer
        self.processed_properties = processed_properties

    def finalize_with_piecewise_function(self, material: Material, prop_name: str,
                                         piecewise_func: sp.Piecewise, T: Union[float, sp.Symbol],
                                         config: Dict, prop_type: str,
                                         x_data: Optional[np.ndarray] = None,
                                         y_data: Optional[np.ndarray] = None) -> bool:
        """
        Finalize property processing when you already have a piecewise function.
        This method handles:
        1. Numeric temperature evaluation
        2. Property assignment
        3. Visualization
        4. Processing tracking
        Returns True if processing is complete (numeric case), False if symbolic.
        """
        logger.debug(f"Finalizing property '{prop_name}' with existing piecewise function")
        # Handle numeric temperature case using utility function
        if handle_numeric_temperature(self, material, prop_name, piecewise_func, T):
            return True
        # Assign symbolic piecewise function
        setattr(material, prop_name, piecewise_func)
        # Generate visualization
        bounds = None
        if x_data is not None and len(x_data) > 0:
            bounds = (np.min(x_data), np.max(x_data))
        self._visualize_if_enabled(material=material, prop_name=prop_name, T=T, prop_type=prop_type,
                                   x_data=x_data, y_data=y_data, config=config, bounds=bounds)
        # Track processed property
        self.processed_properties.add(prop_name)
        logger.debug(f"Successfully finalized property '{prop_name}'")
        return False

    def finalize_with_data_arrays(self, material: Material, prop_name: str,
                                  temp_array: np.ndarray, prop_array: np.ndarray,
                                  T: Union[float, sp.Symbol], config: Dict, prop_type: str) -> bool:
        """
        Finalize property processing when you have data arrays.
        This method handles:
        1. Numeric temperature evaluation (interpolation)
        2. Piecewise function creation
        3. Property assignment
        4. Visualization
        5. Processing tracking
        Returns True if processing is complete (numeric case), False if symbolic.
        """
        logger.debug(f"Finalizing property '{prop_name}' with data arrays")
        # Validate input arrays
        if temp_array is None or prop_array is None:
            raise ValueError(f"Temperature and property arrays cannot be None for '{prop_name}'")
        if len(temp_array) != len(prop_array):
            raise ValueError(f"Temperature and property arrays must have same length for '{prop_name}'")
        # Extract boundary configuration
        lower_bound_type, upper_bound_type = config.get('bounds', ['constant', 'constant'])
        # Handle numeric temperature case (interpolation)
        if not isinstance(T, sp.Symbol):
            from materforge.algorithms.interpolation import interpolate_value
            try:
                value = interpolate_value(T, temp_array, prop_array, lower_bound_type, upper_bound_type)
                setattr(material, prop_name, sp.Float(value))
                self.processed_properties.add(prop_name)
                logger.debug(f"Numeric interpolation completed for '{prop_name}': {value}")
                return True
            except Exception as e:
                raise ValueError(f"Failed to interpolate {prop_name} at T={T}: {str(e)}")
        # Create symbolic piecewise function
        try:
            piecewise_func = PiecewiseBuilder.build_from_data(temp_array, prop_array, T, config, prop_name)
            setattr(material, prop_name, piecewise_func)
            logger.debug(f"Created piecewise function for property '{prop_name}'")
        except Exception as e:
            logger.error(f"Failed to create piecewise function for '{prop_name}': {e}")
            raise ValueError(f"Failed to finalize property '{prop_name}': {str(e)}") from e
        # Generate visualization
        self._visualize_if_enabled(material=material, prop_name=prop_name, T=T, prop_type=prop_type,
                                   x_data=temp_array, y_data=prop_array, config=config,
                                   bounds=(np.min(temp_array), np.max(temp_array)))
        # Track processed property
        self.processed_properties.add(prop_name)
        logger.debug(f"Successfully finalized property '{prop_name}'")
        return False

    def _visualize_if_enabled(self, material: Material, prop_name: str,
                              T: Union[float, sp.Symbol], prop_type: str,
                              x_data: Optional[np.ndarray] = None,
                              y_data: Optional[np.ndarray] = None,
                              config: Optional[Dict] = None,
                              bounds: Optional[Tuple[float, float]] = None) -> None:
        """
        Generate visualization if visualizer is available and enabled.

        This method handles the visualization of processed properties by extracting
        visualization parameters from the configuration and calling the appropriate
        visualizer methods. It gracefully handles cases where visualization is
        disabled or unavailable.
        Args:
            material: Material object containing the processed property
            prop_name: Name of the property to visualize
            T: Temperature symbol or numeric value
            prop_type: Type of property for visualization context
            x_data: Temperature data points (optional)
            y_data: Property values corresponding to temperatures (optional)
            config: Configuration dictionary containing visualization settings (optional)
            bounds: Temperature bounds tuple (min_temp, max_temp) (optional)
        Note:
            This method logs warnings for visualization failures but does not raise
            exceptions, ensuring that visualization issues don't interrupt property
            processing.
        """
        # Check if visualizer is available
        if self.visualizer is None:
            logger.debug(f"No visualizer available for property '{prop_name}'")
            return
        # Skip visualization for numeric temperature
        if not isinstance(T, sp.Symbol):
            logger.debug(f"Skipping visualization for property '{prop_name}' - numeric temperature")
            return
        # Check if visualization is enabled
        if not hasattr(self.visualizer, 'is_visualization_enabled') or \
                not self.visualizer.is_visualization_enabled():
            logger.debug(f"Visualization disabled for property '{prop_name}'")
            return
        try:
            # Initialize default visualization parameters
            has_regression = False
            simplify_type = None
            degree = 1
            segments = 1
            lower_bound_type = 'constant'
            upper_bound_type = 'constant'
            # Extract visualization parameters from config
            if config:
                bounds_config = config.get('bounds', ['constant', 'constant'])
                lower_bound_type, upper_bound_type = bounds_config
                # Check for regression configuration
                if 'regression' in config:
                    has_regression = True
                    regression_config = config['regression']
                    simplify_type = regression_config.get('simplify', 'pre')
                    degree = regression_config.get('degree', 1)
                    segments = regression_config.get('segments', 1)
            # Set temperature bounds if available
            lower_bound = bounds[0] if bounds else None
            upper_bound = bounds[1] if bounds else None
            # Call visualizer with extracted parameters
            self.visualizer.visualize_property(
                material=material,
                prop_name=prop_name,
                T=T,
                prop_type=prop_type,
                x_data=x_data,
                y_data=y_data,
                has_regression=has_regression,
                simplify_type=simplify_type,
                degree=degree,
                segments=segments,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                lower_bound_type=lower_bound_type,
                upper_bound_type=upper_bound_type
            )
            logger.debug(f"Generated visualization for property '{prop_name}'")
        except Exception as e:
            logger.warning(f"Failed to generate visualization for '{prop_name}': {e}")
            # Don't raise exception - visualization failure shouldn't stop processing

    def set_visualizer(self, visualizer) -> None:
        """
        Set the visualizer instance for property plotting.
        Args:
            visualizer: Visualizer instance that implements visualization methods
        """
        self.visualizer = visualizer
        logger.debug("Visualizer set for PropertyProcessorBase")

    def reset_processing_state(self) -> None:
        """
        Reset the processing state for reuse of the processor instance.

        This method clears the set of processed properties, allowing the same
        processor instance to be used for processing multiple materials.
        """
        self.processed_properties.clear()
        logger.debug("Processing state reset")

    def get_processed_properties(self) -> set:
        """
        Get the set of properties that have been processed.
        Returns:
            set: Set of property names that have been successfully processed
        """
        return self.processed_properties.copy()

    def is_property_processed(self, prop_name: str) -> bool:
        """
        Check if a specific property has been processed.
        Args:
            prop_name: Name of the property to check
        Returns:
            bool: True if the property has been processed, False otherwise
        """
        return prop_name in self.processed_properties
