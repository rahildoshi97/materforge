import logging
from typing import Dict, List, Set, Tuple, Union, Any

import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.validation.property_validator import validate_monotonic_energy_density
from materforge.algorithms.interpolation import ensure_ascending_order
from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.algorithms.regression_processor import RegressionProcessor
from materforge.parsing.validation.property_type_detector import PropertyType
from materforge.parsing.config.yaml_keys import (
    REGRESSION_KEY, POST_KEY, DEPENDENCIES_KEY, RANGES_KEY, SIMPLIFY_KEY, PRE_KEY
)

logger = logging.getLogger(__name__)


class PropertyPostProcessor:
    """Handles post-processing of properties after initial processing with multi-dependency support."""

    def post_process_multi_dependency_properties(self, material: Material,
                                                 symbol_mapping: Dict[str, sp.Symbol],
                                                 dependency_resolver: DependencyResolver,
                                                 properties: Dict[str, Any],
                                                 categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                                 processed_properties: Set[str]) -> None:
        """
        Perform post-processing regression on properties with multi-dependency support.

        Args:
            material: Material instance
            symbol_mapping: Dictionary mapping dependency names to SymPy symbols
            dependency_resolver: Resolver for handling multi-dependency configurations
            properties: Raw property configurations
            categorized_properties: Properties organized by type
            processed_properties: Set of already processed properties
        """
        logger.debug("Starting multi-dependency post-processing of properties")

        # Skip post-processing if no symbolic dependencies
        if not any(isinstance(symbol, sp.Symbol) for symbol in symbol_mapping.values()):
            logger.debug("Skipping post-processing - no symbolic dependencies")
            return

        errors = []
        processed_count = len(processed_properties)
        total_count = sum(len(prop_list) for prop_list in categorized_properties.values())

        logger.info(f"Multi-dependency post-processing: {processed_count}/{total_count} properties processed")

        if processed_count < total_count:
            unprocessed = []
            for prop_list in categorized_properties.values():
                for prop_name, _ in prop_list:
                    if prop_name not in processed_properties:
                        unprocessed.append(prop_name)
            logger.warning(f"Some properties were not processed: {unprocessed}")

        # Apply post-processing regression with multi-dependency support
        for prop_name, prop_config in properties.items():
            try:
                if not isinstance(prop_config, dict) or REGRESSION_KEY not in prop_config:
                    continue

                # Check if this property has multi-dependency format
                if DEPENDENCIES_KEY in prop_config:
                    self._apply_multi_dependency_post_regression(
                        material, prop_name, prop_config, symbol_mapping, dependency_resolver
                    )
                else:
                    # Legacy single dependency format
                    self._apply_legacy_post_regression(
                        material, prop_name, prop_config, symbol_mapping, dependency_resolver
                    )

                logger.debug(f"Successfully post-processed property: {prop_name}")

            except Exception as e:
                error_msg = f"Failed to post-process {prop_name}: {str(e)}"
                logger.error(error_msg, exc_info=True)
                errors.append(error_msg)

        if errors:
            error_summary = "\n".join(errors)
            logger.error(f"Post-processing errors occurred: {error_summary}")
            raise ValueError(f"Post-processing errors occurred:\n{error_summary}")

        logger.debug("Multi-dependency post-processing completed successfully")

    def _apply_multi_dependency_post_regression(self, material: Material, prop_name: str,
                                                prop_config: Dict, symbol_mapping: Dict[str, sp.Symbol],
                                                dependency_resolver: DependencyResolver) -> None:
        """Apply post-processing regression for multi-dependency properties."""
        logger.debug(f"Applying multi-dependency post-regression to property: {prop_name}")

        dependencies = prop_config[DEPENDENCIES_KEY]

        # Check for regression configuration
        has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
            prop_config, prop_name, 100  # Default length for validation
        )

        if not has_regression or simplify_type != POST_KEY:
            logger.debug(f"Skipping post-regression for {prop_name} - not configured for post-processing")
            return

        if not hasattr(material, prop_name):
            logger.warning(f"Property '{prop_name}' not found on material during post-processing")
            return

        prop_value = getattr(material, prop_name)

        if isinstance(prop_value, sp.Integral):
            logger.warning(f"Property '{prop_name}' is an integral and cannot be post-processed")
            return

        if not isinstance(prop_value, sp.Expr):
            logger.debug(f"Skipping non-symbolic property: {prop_name}")
            return

        # Resolve dependency ranges
        try:
            resolved_ranges = dependency_resolver.resolve_dependency_ranges(prop_config, material)
        except Exception as e:
            logger.error(f"Failed to resolve dependency ranges for {prop_name}: {e}", exc_info=True)
            raise ValueError(f"Failed to resolve dependency ranges for {prop_name}: {str(e)}") from e

        # For multi-dependency, we focus on the primary dependency for post-regression
        primary_dependency = dependencies[0]
        primary_array = resolved_ranges[primary_dependency]

        # Get the primary symbol
        dep_symbols = dependency_resolver.get_dependency_symbols(dependencies, symbol_mapping)
        yaml_symbol = dependency_resolver.independent_vars[primary_dependency]
        primary_symbol = dep_symbols[yaml_symbol]

        # Apply post-regression to the primary dependency
        self._apply_post_regression_to_dependency(
            material, prop_name, prop_config, prop_value,
            primary_array, primary_symbol, degree, segments
        )

    def _apply_legacy_post_regression(self, material: Material, prop_name: str,
                                      prop_config: Dict, symbol_mapping: Dict[str, sp.Symbol],
                                      dependency_resolver: DependencyResolver) -> None:
        """Apply post-processing regression for legacy single dependency properties."""
        logger.debug(f"Applying legacy post-regression to property: {prop_name}")

        # Check for regression configuration
        try:
            temp_array = dependency_resolver.extract_from_config(prop_config, material, "temperature")
        except Exception as e:
            logger.error(f"Failed to extract dependency array for {prop_name}: {e}", exc_info=True)
            raise ValueError(f"Failed to extract dependency array for {prop_name}: {str(e)}") from e

        has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
            prop_config, prop_name, len(temp_array)
        )

        if not has_regression or simplify_type != POST_KEY:
            return

        if not hasattr(material, prop_name):
            logger.warning(f"Property '{prop_name}' not found on material during post-processing")
            return

        prop_value = getattr(material, prop_name)

        if isinstance(prop_value, sp.Integral):
            logger.warning(f"Property '{prop_name}' is an integral and cannot be post-processed")
            return

        if not isinstance(prop_value, sp.Expr):
            logger.debug(f"Skipping non-symbolic property: {prop_name}")
            return

        # Use first available symbol (for legacy compatibility)
        T = list(symbol_mapping.values())[0]

        # Apply post-regression
        self._apply_post_regression_to_dependency(
            material, prop_name, prop_config, prop_value,
            temp_array, T, degree, segments
        )

    def _apply_post_regression_to_dependency(self, material: Material, prop_name: str,
                                             prop_config: Dict, prop_value: sp.Expr,
                                             dep_array: np.ndarray, symbol: sp.Symbol,
                                             degree: int, segments: int) -> None:
        """Apply post-processing regression to a specific dependency."""
        logger.debug(f"Applying post-regression to dependency for property: {prop_name}")

        # Validate and convert dependency array
        if isinstance(dep_array, str):
            raise ValueError(f"Dependency array for {prop_name} is a string: '{dep_array}'. Expected numpy array.")

        if not isinstance(dep_array, np.ndarray):
            try:
                dep_array = np.array(dep_array, dtype=np.float64)
            except Exception as e:
                raise ValueError(f"Cannot convert dependency data to numpy array for {prop_name}: {str(e)}") from e

        # Validate dtype
        if dep_array.dtype.kind not in ['f', 'i']:
            logger.warning(
                f"Dependency array for {prop_name} has non-numeric dtype {dep_array.dtype}. Converting to float64.")
            try:
                dep_array = np.asarray(dep_array, dtype=np.float64)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert dependency array to numeric format for {prop_name}: {str(e)}") from e

        # Evaluate property at dependency points
        f_prop = sp.lambdify(symbol, prop_value, 'numpy')
        prop_array = f_prop(dep_array)

        # Validate prop_array
        if hasattr(prop_array, 'dtype') and prop_array.dtype.kind not in ['f', 'i']:
            try:
                prop_array = np.asarray(prop_array, dtype=np.float64)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert property array to numeric format for {prop_name}: {str(e)}") from e

        dep_array, prop_array = ensure_ascending_order(dep_array, prop_array)
        validate_monotonic_energy_density(prop_name, dep_array, prop_array)

        # Apply regression
        # Temporarily modify config to force PRE_KEY regression
        original_simplify_type = prop_config[REGRESSION_KEY][SIMPLIFY_KEY]
        prop_config[REGRESSION_KEY][SIMPLIFY_KEY] = PRE_KEY

        try:
            piecewise_func = PiecewiseBuilder.build_from_data(
                dep_array, prop_array, symbol, prop_config, prop_name
            )
            logger.debug(f"Successfully created piecewise function for {prop_name} with post-regression")
        finally:
            # Restore original config
            prop_config[REGRESSION_KEY][SIMPLIFY_KEY] = original_simplify_type

        setattr(material, prop_name, piecewise_func)
        logger.debug(f"Successfully applied post-regression to property: {prop_name}")

    # Legacy method for backward compatibility
    def post_process_properties(self, material: Material, T: Union[float, sp.Symbol],
                                properties: Dict[str, Any],
                                categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                processed_properties: Set[str]) -> None:
        """
        Legacy post-processing method for backward compatibility.
        Delegates to multi-dependency version with default symbol mapping.
        """
        logger.debug("Using legacy post-processing method - converting to multi-dependency format")

        if not isinstance(T, sp.Symbol):
            logger.debug("Skipping post-processing for numeric temperature")
            return

        # Create default symbol mapping for legacy compatibility
        symbol_mapping = {'temperature': T}

        # Create a minimal dependency resolver for legacy support
        config = {'independent_variables': {'temperature': 'T'}}
        dependency_resolver = DependencyResolver(config)

        # Call the new multi-dependency method
        self.post_process_multi_dependency_properties(
            material, symbol_mapping, dependency_resolver, properties,
            categorized_properties, processed_properties
        )
