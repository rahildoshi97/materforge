# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import Dict, List, Set, Tuple, Union, Any
import numpy as np
import sympy as sp
from materforge.core.materials import Material
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.algorithms.interpolation import ensure_ascending_order
from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.algorithms.regression_processor import RegressionProcessor
from materforge.parsing.validation.property_type_detector import PropertyType
from materforge.parsing.config.yaml_keys import REGRESSION_KEY, POST_KEY, SIMPLIFY_KEY, PRE_KEY

logger = logging.getLogger(__name__)

class PropertyPostProcessor:
    """Handles post-processing of properties after initial processing."""

    def post_process_properties(self, material: Material, dependency: Union[float, sp.Symbol],
                                properties: Dict[str, Any],
                                categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                processed_properties: Set[str]) -> None:
        """Applies post-processing regression to all eligible properties.

        Only runs in symbolic mode - skipped entirely when dependency is numeric.

        Args:
            material:                 Material instance to update.
            dependency:               SymPy symbol (symbolic mode) or float (numeric mode).
            properties:               Full property configuration dict.
            categorized_properties:   Properties grouped by PropertyType.
            processed_properties:     Set of already-processed property names.
        """
        logger.debug("Starting post-processing of properties")
        if not isinstance(dependency, sp.Symbol):
            logger.debug("Skipping post-processing for numeric dependency")
            return
        errors = []
        processed_count = len(processed_properties)
        total_count = sum(len(prop_list) for prop_list in categorized_properties.values())
        logger.info("Post-processing: %d/%d properties processed", processed_count, total_count)
        if processed_count < total_count:
            unprocessed = [
                prop_name
                for prop_list in categorized_properties.values()
                for prop_name, _ in prop_list
                if prop_name not in processed_properties
            ]
            logger.warning("Some properties were not processed: %s", unprocessed)
        for prop_name, prop_config in properties.items():
            try:
                if not isinstance(prop_config, dict) or REGRESSION_KEY not in prop_config:
                    continue
                dep_array = DependencyResolver.extract_from_config(prop_config, material)
                has_regression, simplify_type, degree, segments = RegressionProcessor.process_regression_params(
                    prop_config, prop_name, len(dep_array)
                )
                if not has_regression or simplify_type != POST_KEY:
                    continue
                if not hasattr(material, prop_name):
                    logger.warning("Property '%s' not found on material during post-processing", prop_name)
                    continue
                prop_value = getattr(material, prop_name)
                if isinstance(prop_value, sp.Integral):
                    logger.warning("Property '%s' is an unevaluated integral - cannot post-process", prop_name)
                    continue
                if not isinstance(prop_value, sp.Expr):
                    logger.debug("Skipping non-symbolic property: '%s'", prop_name)
                    continue
                logger.debug("Applying post-processing regression to property: '%s'", prop_name)
                self._apply_post_regression(material, prop_name, prop_config, dependency)
                logger.debug("Successfully post-processed property: '%s'", prop_name)
            except Exception as e:
                error_msg = f"Failed to post-process '{prop_name}': {str(e)}"
                logger.error(error_msg, exc_info=True)
                errors.append(error_msg)
        if errors:
            error_summary = "\n".join(errors)
            raise ValueError(f"Post-processing errors occurred:\n{error_summary}")
        logger.debug("Post-processing completed successfully")

    @staticmethod
    def _apply_post_regression(material: Material, prop_name: str,
                               prop_config: Dict, dependency: sp.Symbol) -> None:
        """Evaluates a symbolic property over its dependency range and fits a
        piecewise regression, replacing the original expression on the material.

        Args:
            material:    Material instance to update.
            prop_name:   Name of the property to regress.
            prop_config: Property configuration dict (modified temporarily during call).
            dependency:  SymPy symbol used as the dependency variable.
        """
        logger.debug("Applying post-regression to property: '%s'", prop_name)
        prop_value = getattr(material, prop_name)
        try:
            dep_array = DependencyResolver.extract_from_config(prop_config, material)
        except Exception as e:
            raise ValueError(f"Failed to extract dependency array for '{prop_name}': {str(e)}") from e
        # Ensure dep_array is a float64 numpy array
        if isinstance(dep_array, str):
            raise ValueError(f"Dependency array for '{prop_name}' is a string: '{dep_array}'. Expected numpy array.")
        if not isinstance(dep_array, np.ndarray):
            try:
                dep_array = np.array(dep_array, dtype=np.float64)
            except Exception as e:
                raise ValueError(f"Cannot convert dependency data to numpy array for '{prop_name}': {str(e)}") from e
        if dep_array.dtype.kind not in ['f', 'i']:
            logger.warning("Dependency array for '%s' has non-numeric dtype %s. Converting to float64.",
                prop_name, dep_array.dtype)
            try:
                dep_array = np.asarray(dep_array, dtype=np.float64)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert dependency array to numeric format for '{prop_name}': {str(e)}") from e
        # Evaluate symbolic property over the dependency array
        f_prop = sp.lambdify(dependency, prop_value, 'numpy')
        prop_array = f_prop(dep_array)
        if hasattr(prop_array, 'dtype') and prop_array.dtype.kind not in ['f', 'i']:
            try:
                prop_array = np.asarray(prop_array, dtype=np.float64)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert property array to numeric format for '{prop_name}': {str(e)}") from e
        dep_array, prop_array = ensure_ascending_order(dep_array, prop_array)
        # Temporarily set simplify type to PRE so PiecewiseBuilder runs regression directly
        original_simplify_type = prop_config[REGRESSION_KEY][SIMPLIFY_KEY]
        prop_config[REGRESSION_KEY][SIMPLIFY_KEY] = PRE_KEY
        try:
            piecewise_func = PiecewiseBuilder.build_from_data(dep_array, prop_array, dependency, prop_config, prop_name)
        finally:
            prop_config[REGRESSION_KEY][SIMPLIFY_KEY] = original_simplify_type
        setattr(material, prop_name, piecewise_func)
        logger.debug("Successfully applied post-regression to property: '%s'", prop_name)
