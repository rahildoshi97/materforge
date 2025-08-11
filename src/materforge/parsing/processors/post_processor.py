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
from materforge.parsing.config.yaml_keys import REGRESSION_KEY, POST_KEY

logger = logging.getLogger(__name__)


class PropertyPostProcessor:
    """Handles post-processing of properties after initial processing."""

    def post_process_properties(self, material: Material, T: Union[float, sp.Symbol],
                                properties: Dict[str, Any],
                                categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                processed_properties: Set[str]) -> None:
        """Perform post-processing regression on properties after all have been initially processed."""
        logger.debug("Starting post-processing of properties")
        if not isinstance(T, sp.Symbol):
            logger.debug("Skipping post-processing for numeric temperature")
            return
        errors = []
        processed_count = len(processed_properties)
        total_count = sum(len(prop_list) for prop_list in categorized_properties.values())
        logger.info(f"Post-processing: {processed_count}/{total_count} properties processed")
        if processed_count < total_count:
            unprocessed = []
            for prop_list in categorized_properties.values():
                for prop_name, _ in prop_list:
                    if prop_name not in processed_properties:
                        unprocessed.append(prop_name)
            logger.warning(f"Some properties were not processed: {unprocessed}")
        # Apply post-processing regression
        for prop_name, prop_config in properties.items():
            try:
                if not isinstance(prop_config, dict) or REGRESSION_KEY not in prop_config:
                    continue
                temp_array = DependencyResolver.extract_from_config(prop_config, material)
                has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
                    prop_config, prop_name, len(temp_array)
                )
                if not has_regression or simplify_type != POST_KEY:
                    continue
                if not hasattr(material, prop_name):
                    logger.warning(f"Property '{prop_name}' not found on material during post-processing")
                    continue
                prop_value = getattr(material, prop_name)
                if isinstance(prop_value, sp.Integral):
                    logger.warning(f"Property '{prop_name}' is an integral and cannot be post-processed")
                    continue
                if not isinstance(prop_value, sp.Expr):
                    logger.debug(f"Skipping non-symbolic property: {prop_name}")
                    continue
                logger.debug(f"Applying post-processing regression to property: {prop_name}")
                self._apply_post_regression(material, prop_name, prop_config, T)
                logger.debug(f"Successfully post-processed property: {prop_name}")
            except Exception as e:
                error_msg = f"Failed to post-process {prop_name}: {str(e)}"
                logger.error(error_msg, exc_info=True)
                errors.append(error_msg)
        if errors:
            error_summary = "\n".join(errors)
            logger.error(f"Post-processing errors occurred: {error_summary}")
            raise ValueError(f"Post-processing errors occurred:\n{error_summary}")
        logger.debug("Post-processing completed successfully")

    @staticmethod
    def _apply_post_regression(material: Material, prop_name: str,
                               prop_config: Dict, T: sp.Symbol) -> None:
        """Apply post-processing regression with proper type validation."""
        logger.debug(f"Applying post-regression to property: {prop_name}")
        prop_value = getattr(material, prop_name)
        try:
            temp_array = DependencyResolver.extract_from_config(prop_config, material)
        except Exception as e:
            logger.error(f"Failed to extract dependency array for {prop_name}: {e}", exc_info=True)
            raise ValueError(f"Failed to extract dependency array for {prop_name}: {str(e)}") from e
        # Validate and convert temp_array
        if isinstance(temp_array, str):
            raise ValueError(f"Temperature array for {prop_name} is a string: '{temp_array}'. Expected numpy array.")
        if not isinstance(temp_array, np.ndarray):
            try:
                temp_array = np.array(temp_array, dtype=np.float64)
            except Exception as e:
                raise ValueError(f"Cannot convert temperature data to numpy array for {prop_name}: {str(e)}") from e
        # Validate dtype
        if temp_array.dtype.kind not in ['f', 'i']:
            logger.warning(
                f"Temperature array for {prop_name} has non-numeric dtype {temp_array.dtype}. Converting to float64.")
            try:
                temp_array = np.asarray(temp_array, dtype=np.float64)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert temperature array to numeric format for {prop_name}: {str(e)}") from e
        # Evaluate property at temperature points
        f_prop = sp.lambdify(T, prop_value, 'numpy')
        prop_array = f_prop(temp_array)
        # Validate prop_array
        if hasattr(prop_array, 'dtype') and prop_array.dtype.kind not in ['f', 'i']:
            try:
                prop_array = np.asarray(prop_array, dtype=np.float64)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert property array to numeric format for {prop_name}: {str(e)}") from e
        temp_array, prop_array = ensure_ascending_order(temp_array, prop_array)
        validate_monotonic_energy_density(prop_name, temp_array, prop_array)
        # Apply regression
        from materforge.parsing.config.yaml_keys import SIMPLIFY_KEY, PRE_KEY
        # Temporarily modify config to force PRE_KEY regression
        original_simplify_type = prop_config[REGRESSION_KEY][SIMPLIFY_KEY]
        prop_config[REGRESSION_KEY][SIMPLIFY_KEY] = PRE_KEY
        try:
            piecewise_func = PiecewiseBuilder.build_from_data(temp_array, prop_array, T, prop_config, prop_name)
            logger.debug(f"Successfully created piecewise function for {prop_name} with post-regression")
        finally:  # Restore original config
            prop_config[REGRESSION_KEY][SIMPLIFY_KEY] = original_simplify_type
        setattr(material, prop_name, piecewise_func)
        logger.debug(f"Successfully applied post-regression to property: {prop_name}")
