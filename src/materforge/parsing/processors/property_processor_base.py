# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""
Base processor class for material property processing with common finalization logic.

Provides PropertyProcessorBase, which serves as the foundation for all property
handlers in MaterForge, ensuring consistent processing patterns and reducing
code duplication.
"""

import logging
import numpy as np
import sympy as sp
from typing import Dict, Union, Tuple, Optional
from pathlib import Path
from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.core.materials import Material


logger = logging.getLogger(__name__)


class PropertyProcessorBase:
    """Base class for material property processors with common finalization logic.

    Provides shared functionality for processing material properties, including
    piecewise function creation, visualization, and property assignment. All
    property handlers in MaterForge inherit from this class.

    Attributes:
        processed_properties: Set of property names that have been processed.
        base_dir:             Base directory for resolving relative file paths.
        visualizer:           Optional visualizer instance for property plotting.
    """

    def __init__(self):
        self.processed_properties = set()
        self.base_dir: Path = None  # type: ignore
        self.visualizer = None
        logger.debug("PropertyProcessorBase initialized")

    def set_processing_context(self, base_dir: Path, visualizer, processed_properties: set):
        """Sets processing context shared across all handlers."""
        self.base_dir = base_dir
        self.visualizer = visualizer
        self.processed_properties = processed_properties

    def finalize_with_piecewise_function(self, material: Material, prop_name: str,
                                         piecewise_func: sp.Piecewise,
                                         dependency: sp.Symbol,
                                         config: Dict, prop_type: str,
                                         x_data: Optional[np.ndarray] = None,
                                         y_data: Optional[np.ndarray] = None) -> bool:
        """Finalises property processing when a piecewise function is already available.

        Handles numeric evaluation, property assignment, visualization, and tracking.

        Args:
            material:       Material instance to assign the property to.
            prop_name:      Name of the property.
            piecewise_func: Pre-built SymPy Piecewise expression.
            dependency:     SymPy symbol (symbolic mode).
            config:         Property configuration dict.
            prop_type:      Property type string for visualization context.
            x_data:         Optional dependency values for visualization.
            y_data:         Optional property values for visualization.
        Returns:
            True if numeric evaluation was performed, False if symbolic.
        """
        logger.debug("Finalizing property '%s' with existing piecewise function", prop_name)
        setattr(material, prop_name, piecewise_func)
        bounds = (np.min(x_data), np.max(x_data)) if x_data is not None and len(x_data) > 0 else None
        self._visualize_if_enabled(
            material=material, prop_name=prop_name, dependency=dependency,
            prop_type=prop_type, x_data=x_data, y_data=y_data, config=config, bounds=bounds
        )
        self.processed_properties.add(prop_name)
        logger.debug("Successfully finalized property '%s'", prop_name)
        return False

    def finalize_with_data_arrays(self, material: Material, prop_name: str,
                                  dep_array: np.ndarray, prop_array: np.ndarray,
                                  dependency: sp.Symbol,
                                  config: Dict, prop_type: str) -> bool:
        """Finalises property processing from raw dependency-value arrays.

        In numeric mode, interpolates at the given value and assigns a float.
        In symbolic mode, builds a SymPy Piecewise and assigns it to the material.

        Args:
            material:    Material instance to assign the property to.
            prop_name:   Name of the property.
            dep_array:   Dependency variable values.
            prop_array:  Corresponding property values.
            dependency:  SymPy symbol (symbolic mode).
            config:      Property configuration dict.
            prop_type:   Property type string for visualization context.
        Returns:
            True if numeric evaluation was performed, False if symbolic.
        Raises:
            ValueError: If arrays are None, mismatched in length, interpolation fails,
                or piecewise construction fails.
        """
        logger.debug("Finalizing property '%s' with data arrays", prop_name)
        if dep_array is None or prop_array is None:
            raise ValueError(f"Dependency and property arrays cannot be None for '{prop_name}'")
        if len(dep_array) != len(prop_array):
            raise ValueError(f"Dependency and property arrays must have equal length for '{prop_name}' "
                f"(got {len(dep_array)} and {len(prop_array)})")

        # Symbolic mode: build piecewise function
        try:
            piecewise_func = PiecewiseBuilder.build_from_data(dep_array, prop_array, dependency, config, prop_name)
            setattr(material, prop_name, piecewise_func)
            logger.debug("Created piecewise function for property '%s'", prop_name)
        except Exception as e:
            raise ValueError(f"Failed to finalize property '{prop_name}': {str(e)}") from e
        self._visualize_if_enabled(material=material, prop_name=prop_name, dependency=dependency,
            prop_type=prop_type, x_data=dep_array, y_data=prop_array, config=config,
            bounds=(np.min(dep_array), np.max(dep_array)))
        self.processed_properties.add(prop_name)
        logger.debug("Successfully finalized property '%s'", prop_name)
        return False

    def _visualize_if_enabled(self, material: Material, prop_name: str,
                              dependency: Union[float, sp.Symbol], prop_type: str,
                              x_data: Optional[np.ndarray] = None,
                              y_data: Optional[np.ndarray] = None,
                              config: Optional[Dict] = None,
                              bounds: Optional[Tuple[float, float]] = None) -> None:
        """Generates a visualization if the visualizer is available and enabled.

        Failures are logged as warnings and never propagate - visualization issues
        must not interrupt property processing.

        Args:
            material:    Material instance.
            prop_name:   Name of the property to visualize.
            dependency:  SymPy symbol. Visualization is skipped for numeric values.
            prop_type:   Property type string for visualization context.
            x_data:      Dependency axis data (optional).
            y_data:      Property axis data (optional).
            config:      Property configuration dict (optional).
            bounds:      (min, max) of the dependency range (optional).
        """
        if self.visualizer is None:
            logger.debug("No visualizer available for property '%s'", prop_name)
            return
        if not isinstance(dependency, sp.Symbol):
            logger.debug("Skipping visualization for '%s' - numeric dependency", prop_name)
            return
        if not hasattr(self.visualizer, 'is_visualization_enabled') or                 not self.visualizer.is_visualization_enabled():
            logger.debug("Visualization disabled for property '%s'", prop_name)
            return
        try:
            has_regression = False
            simplify_type = None
            degree = 1
            segments = 1
            lower_bound_type = 'constant'
            upper_bound_type = 'constant'
            if config:
                lower_bound_type, upper_bound_type = config.get('bounds', ['constant', 'constant'])
                if 'regression' in config:
                    has_regression = True
                    reg = config['regression']
                    simplify_type = reg.get('simplify', 'pre')
                    degree = reg.get('degree', 1)
                    segments = reg.get('segments', 1)
            self.visualizer.visualize_property(
                material=material,
                prop_name=prop_name,
                dependency=dependency,
                prop_type=prop_type,
                x_data=x_data,
                y_data=y_data,
                has_regression=has_regression,
                simplify_type=simplify_type,
                degree=degree,
                segments=segments,
                lower_bound=bounds[0] if bounds else None,
                upper_bound=bounds[1] if bounds else None,
                lower_bound_type=lower_bound_type,
                upper_bound_type=upper_bound_type
            )
            logger.debug("Generated visualization for property '%s'", prop_name)
        except Exception as e:
            logger.warning("Failed to generate visualization for '%s': %s", prop_name, e)

    def set_visualizer(self, visualizer) -> None:
        """Sets the visualizer instance for property plotting."""
        self.visualizer = visualizer
        logger.debug("Visualizer set for PropertyProcessorBase")

    def reset_processing_state(self) -> None:
        """Clears processed properties to allow reuse of the processor instance."""
        self.processed_properties.clear()
        logger.debug("Processing state reset")

    def get_processed_properties(self) -> set:
        """Returns a copy of the set of successfully processed property names."""
        return self.processed_properties.copy()

    def is_property_processed(self, prop_name: str) -> bool:
        """Returns True if the given property has already been processed."""
        return prop_name in self.processed_properties
