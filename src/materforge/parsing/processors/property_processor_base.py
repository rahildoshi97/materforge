"""
Enhanced base processor class for material property processing with multi-dependency support.
"""

import logging
import numpy as np
import sympy as sp
from typing import Dict, Union, Tuple, Optional, List
from pathlib import Path

from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.core.materials import Material
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.utils.utilities import handle_numeric_temperature

logger = logging.getLogger(__name__)


class PropertyProcessorBase:
    """
    Enhanced base class for material property processors with multi-dependency support.
    """

    def __init__(self):
        """Initialize the base processor with empty state."""
        self.processed_properties = set()
        self.base_dir: Path = None
        self.visualizer = None
        self.dependency_resolver: Optional[DependencyResolver] = None
        self.symbol_mapping: Optional[Dict[str, sp.Symbol]] = None
        logger.debug("PropertyProcessorBase initialized")

    def set_multi_dependency_context(self, base_dir: Path, visualizer, processed_properties: set,
                                     dependency_resolver, symbol_mapping: Dict[str, sp.Symbol]):
        """Set multi-dependency processing context shared across handlers."""
        self.base_dir = base_dir
        self.visualizer = visualizer
        self.processed_properties = processed_properties
        self.dependency_resolver = dependency_resolver
        self.symbol_mapping = symbol_mapping

        logger.debug("Multi-dependency context set with %d dependencies: %s",
                     len(symbol_mapping), list(symbol_mapping.keys()))

    def finalize_multi_dependency_property(self, material: Material, prop_name: str,
                                           resolved_ranges: Dict[str, np.ndarray],
                                           prop_values: np.ndarray,
                                           dependencies: List[str],
                                           config: Dict) -> bool:
        """Finalize multi-dependency property processing."""
        logger.debug(f"Finalizing multi-dependency property '{prop_name}' with {len(dependencies)} dependencies")

        # For single dependency, use existing logic
        if len(dependencies) == 1:
            dep_name = dependencies[0]
            temp_array = resolved_ranges[dep_name]

            # Get the symbol for this dependency
            if hasattr(self, 'symbol_mapping') and self.symbol_mapping:
                primary_symbol = list(self.symbol_mapping.values())[0]
            else:
                primary_symbol = sp.Symbol('T')

            return self.finalize_with_data_arrays(
                material=material,
                prop_name=prop_name,
                temp_array=temp_array,
                prop_array=prop_values,
                T=primary_symbol,
                config=config,
                prop_type=f'MULTI_DEPENDENCY_{len(dependencies)}D'
            )
        else:
            # Multi-dependency case - set symbolic expression for now
            # This can be enhanced later for true multi-dimensional support
            setattr(material, prop_name, prop_values)  # Assuming prop_values is symbolic
            self.processed_properties.add(prop_name)
            return False

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Dict) -> None:
        """
        Process a property with multi-dependency support.
        This method should be implemented by each specific handler.
        """
        raise NotImplementedError("Subclasses must implement process_multi_dependency_property")

    def resolve_property_dependencies(self, config: Dict, material: Material) -> Dict[str, np.ndarray]:
        """
        Resolve all dependencies for a property configuration.

        Returns:
            Dict mapping dependency names to their resolved arrays
        """
        if not self.dependency_resolver:
            raise ValueError("Dependency resolver not available")

        return self.dependency_resolver.resolve_dependency_ranges(config, material)

    def get_dependency_symbols(self, dependencies: List[str]) -> Dict[str, sp.Symbol]:
        """Get symbol mapping for specific dependencies."""
        if not self.dependency_resolver or not self.symbol_mapping:
            raise ValueError("Dependency resolver or symbol mapping not available")

        return self.dependency_resolver.get_dependency_symbols(dependencies, self.symbol_mapping)

    def finalize_multi_dependency_property1(self, material: Material, prop_name: str,
                                           resolved_ranges: Dict[str, np.ndarray],
                                           prop_values: np.ndarray,
                                           dependencies: List[str],
                                           config: Dict) -> bool:
        """
        Finalize multi-dependency property processing.

        Returns True if processing is complete (numeric case), False if symbolic.
        """
        logger.debug(f"Finalizing multi-dependency property '{prop_name}' with {len(dependencies)} dependencies")

        # Handle single dependency case (similar to existing logic)
        if len(dependencies) == 1:
            dep_name = dependencies[0]
            temp_array = resolved_ranges[dep_name]

            # Get the symbol for this dependency
            dep_symbols = self.get_dependency_symbols(dependencies)
            yaml_symbol = self.dependency_resolver.independent_vars[dep_name]
            actual_symbol = dep_symbols[yaml_symbol]

            return self.finalize_with_data_arrays(
                material=material,
                prop_name=prop_name,
                temp_array=temp_array,
                prop_array=prop_values,
                T=actual_symbol,
                config=config,
                prop_type=f'MULTI_DEPENDENCY_{len(dependencies)}D'
            )

        # Handle multi-dependency case (2+ dependencies)
        else:
            # For now, create a simple symbolic function
            # This can be enhanced later for true multidimensional support
            primary_dep = dependencies[0]
            primary_array = resolved_ranges[primary_dep]

            dep_symbols = self.get_dependency_symbols(dependencies)
            yaml_symbol = self.dependency_resolver.independent_vars[primary_dep]
            primary_symbol = dep_symbols[yaml_symbol]

            # Create piecewise function for primary dependency
            piecewise_func = PiecewiseBuilder.build_from_data(
                primary_array, prop_values, primary_symbol, config, prop_name
            )

            # Set the property
            setattr(material, prop_name, piecewise_func)

            # Generate visualization for primary dependency
            self._visualize_if_enabled(
                material=material,
                prop_name=prop_name,
                T=primary_symbol,
                prop_type=f'MULTI_DEPENDENCY_{len(dependencies)}D',
                x_data=primary_array,
                y_data=prop_values,
                config=config,
                bounds=(np.min(primary_array), np.max(primary_array))
            )

            # Track processed property
            self.processed_properties.add(prop_name)
            logger.debug(f"Successfully finalized multi-dependency property '{prop_name}'")

            return False  # Symbolic processing

    def finalize_with_data_arrays(self, material: Material, prop_name: str,
                                  temp_array: np.ndarray, prop_array: np.ndarray,
                                  T: Union[float, sp.Symbol], config: Dict, prop_type: str) -> bool:
        """
        Finalize property processing when you have data arrays.
        Enhanced for multi-dependency support.
        """
        logger.debug(f"Finalizing property '{prop_name}' with data arrays")

        # Validate input arrays
        if temp_array is None or prop_array is None:
            raise ValueError(f"Temperature and property arrays cannot be None for '{prop_name}'")

        if len(temp_array) != len(prop_array):
            raise ValueError(f"Temperature and property arrays must have same length for '{prop_name}'")

        # Extract boundary configuration
        bounds_config = config.get('bounds', ['constant', 'constant'])

        # Handle multi-dependency bounds format
        if isinstance(bounds_config, dict) and len(bounds_config) == 1:
            # Single dependency with dict bounds format
            bounds_list = list(bounds_config.values())[0]
            lower_bound_type, upper_bound_type = bounds_list
        elif isinstance(bounds_config, list):
            # Legacy list format
            lower_bound_type, upper_bound_type = bounds_config
        else:
            # Default fallback
            lower_bound_type, upper_bound_type = ['constant', 'constant']

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
        self._visualize_if_enabled(
            material=material,
            prop_name=prop_name,
            T=T,
            prop_type=prop_type,
            x_data=temp_array,
            y_data=prop_array,
            config=config,
            bounds=(np.min(temp_array), np.max(temp_array))
        )

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
        """Generate visualization if visualizer is available and enabled."""
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

                # Handle multi-dependency bounds format
                if isinstance(bounds_config, dict):
                    # Use first dependency's bounds for visualization
                    first_bounds = list(bounds_config.values())[0]
                    lower_bound_type, upper_bound_type = first_bounds
                elif isinstance(bounds_config, list):
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

    def reset_processing_state(self) -> None:
        """Reset the processing state for reuse of the processor instance."""
        self.processed_properties.clear()
        logger.debug("Processing state reset")

    def get_processed_properties(self) -> set:
        """Get the set of properties that have been processed."""
        return self.processed_properties.copy()

    def is_property_processed(self, prop_name: str) -> bool:
        """Check if a specific property has been processed."""
        return prop_name in self.processed_properties
