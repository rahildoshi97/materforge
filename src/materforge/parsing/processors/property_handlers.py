# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import Any, Dict, Union
import numpy as np
import sympy as sp
from materforge.core.materials import Material
from materforge.parsing.processors.property_processor_base import PropertyProcessorBase
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.io.data_handler import load_property_data
from materforge.parsing.utils.utilities import create_step_visualization_data
from materforge.algorithms.interpolation import ensure_ascending_order
from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.parsing.config.yaml_keys import (
    DEPENDENCY_KEY, VALUE_KEY, BOUNDS_KEY, FILE_PATH_KEY, EQUATION_KEY
)
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


class BasePropertyHandler(PropertyProcessorBase):
    """Base class for property handlers with common functionality.

    All specialised handlers inherit from this class, which provides shared
    processing logic via PropertyProcessorBase.
    """

    def __init__(self):
        super().__init__()
        logger.debug("BasePropertyHandler initialized")


class ConstantValuePropertyHandler(BasePropertyHandler):
    """Handler for constant (dependency-independent) properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Union[float, str],
                         dependency: sp.Symbol) -> None:
        """Processes a constant numeric property.

        Args:
            material:    Material instance to assign the value to.
            prop_name:   Name of the property.
            prop_config: Raw config value - expected to be a numeric scalar.
            dependency:  SymPy symbol
        """
        try:
            value = float(prop_config)
            setattr(material, prop_name, sp.Float(value))
            logger.debug("Set constant property '%s' = %s", prop_name, value)
            self._visualize_if_enabled(material=material, prop_name=prop_name, dependency=dependency,
                prop_type='CONSTANT_VALUE', x_data=None, y_data=None)

            self.processed_properties.add(prop_name)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Failed to process constant property '{prop_name}'\n -> {str(e)}") from e


class StepFunctionPropertyHandler(BasePropertyHandler):
    """Handler for step function properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any],
                         dependency: sp.Symbol) -> None:
        """Processes a step function property.

        Builds a SymPy Piecewise with two constant segments around a transition
        point resolved from a scalar property reference.

        Args:
            material:    Material instance.
            prop_name:   Name of the property.
            prop_config: Property configuration dict.
            dependency:  SymPy symbol.
        """
        try:
            dep_key = prop_config[DEPENDENCY_KEY]
            val_array = prop_config[VALUE_KEY]
            transition_point = DependencyResolver.resolve_dependency_reference(dep_key, material)
            # Build piecewise using the YAML placeholder, then substitute caller's symbol
            step_function = sp.Piecewise(
                (val_array[0], dependency < transition_point),
                (val_array[1], True)
            )
            # Build visualization data around the transition point
            offset = ProcessingConstants.STEP_FUNCTION_OFFSET
            dep_range = np.array([transition_point - offset, transition_point, transition_point + offset])
            x_data, y_data = create_step_visualization_data(transition_point, val_array, dep_range)
            self.finalize_with_piecewise_function(material=material, prop_name=prop_name, piecewise_func=step_function,
                dependency=dependency, config=prop_config, prop_type='STEP_FUNCTION',
                x_data=x_data, y_data=y_data)
        except Exception as e:
            raise ValueError(f"Failed to process step function property '{prop_name}'\n -> {str(e)}") from e


class FileImportPropertyHandler(BasePropertyHandler):
    """Handler for file-based properties."""

    def process_property(self, material: Material, prop_name: str,
                         file_config: Dict[str, Any],
                         dependency: sp.Symbol) -> None:
        """Processes a property loaded from an external file.

        Args:
            material:    Material instance.
            prop_name:   Name of the property.
            file_config: Property configuration dict containing file_path and column keys.
            dependency:  SymPy symbol.
        """
        try:
            file_path = self.base_dir / file_config[FILE_PATH_KEY]
            file_config[FILE_PATH_KEY] = str(file_path)
            logger.debug("Loading property '%s' from file: %s", prop_name, file_path)
            dep_array, prop_array = load_property_data(file_config)
            logger.debug("Loaded %d data points for property '%s' (range: %s – %s)",
                len(dep_array), prop_name, dep_array.min(), dep_array.max())
            self.finalize_with_data_arrays(material=material, prop_name=prop_name, dep_array=dep_array,
                prop_array=prop_array, dependency=dependency,
                config=file_config, prop_type='FILE_IMPORT')
        except FileNotFoundError as e:
            raise ValueError(f"File not found for property '{prop_name}': {file_path}") from e
        except Exception as e:
            raise ValueError(f"Failed to process file property '{prop_name}'\n -> {str(e)}") from e


class TabularDataPropertyHandler(BasePropertyHandler):
    """Handler for tabular (explicit dependency–value pair) properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any],
                         dependency: sp.Symbol) -> None:
        """Processes a property defined by explicit dependency–value pairs.

        Args:
            material:    Material instance.
            prop_name:   Name of the property.
            prop_config: Property configuration dict.
            dependency:  SymPy symbol.
        """
        try:
            dep_def = prop_config[DEPENDENCY_KEY]
            val_array = prop_config[VALUE_KEY]
            dep_array = DependencyResolver.resolve_dependency_definition(dep_def, len(val_array), material)
            if len(dep_array) != len(val_array):
                raise ValueError(f"Length mismatch in '{prop_name}': "
                    f"dependency array ({len(dep_array)}) and value array ({len(val_array)}) must match")
            dep_array, val_array = ensure_ascending_order(dep_array, val_array)
            self.finalize_with_data_arrays(material=material, prop_name=prop_name, dep_array=dep_array,
                prop_array=val_array, dependency=dependency,
                config=prop_config, prop_type='TABULAR_DATA')
        except Exception as e:
            raise ValueError(f"Failed to process tabular data property '{prop_name}'\n -> {str(e)}") from e


class PiecewiseEquationPropertyHandler(BasePropertyHandler):
    """Handler for piecewise equation properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any],
                         dependency: sp.Symbol) -> None:
        """Processes a property defined by piecewise symbolic equations.

        n breakpoints in dependency define n-1 equations. The YAML placeholder
        symbol is substituted with the caller's dependency symbol.

        Args:
            material:    Material instance.
            prop_name:   Name of the property.
            prop_config: Property configuration dict.
            dependency:  SymPy symbol.
        """
        try:
            eqn_strings = prop_config[EQUATION_KEY]
            dep_def = prop_config[DEPENDENCY_KEY]
            dep_points = DependencyResolver.resolve_dependency_definition(dep_def, len(eqn_strings) + 1)
            # Validate that equations only reference the YAML placeholder
            for eqn in eqn_strings:
                expr = sp.sympify(eqn)
                unexpected = [s for s in expr.free_symbols if s != dependency]
                if unexpected:
                    raise ValueError(f"Unexpected symbol(s) {unexpected} in equation '{eqn}' for property '{prop_name}'. "
                        f"Equations may only reference the dependency symbol '{dependency}'.")
            lower_bound_type, upper_bound_type = prop_config[BOUNDS_KEY]
            dep_points, eqn_strings = ensure_ascending_order(dep_points, eqn_strings)
            # Build piecewise using placeholder, then substitute caller's symbol
            piecewise_func = PiecewiseBuilder.build_from_formulas(dep_points, list(eqn_strings), dependency,
                lower_bound_type, upper_bound_type)

            # Build dense array for visualization using the placeholder expression
            diff = max(np.min(np.diff(np.sort(dep_points))) / 10.0, 1.0)
            dep_dense = np.arange(dep_points[0], dep_points[-1] + diff / 2, diff)
            f_pw = sp.lambdify(dependency, piecewise_func, 'numpy')
            y_dense = f_pw(dep_dense)
            self.finalize_with_piecewise_function(material=material, prop_name=prop_name, piecewise_func=piecewise_func,
                dependency=dependency, config=prop_config, prop_type='PIECEWISE_EQUATION',
                x_data=dep_dense, y_data=y_dense)
        except Exception as e:
            raise ValueError(f"Failed to process piecewise equation property '{prop_name}'\n -> {str(e)}") from e


class ComputedPropertyHandler(BasePropertyHandler):
    """Handler for computed (derived) properties."""

    def __init__(self):
        super().__init__()
        self.computed_property_processor = None

    def set_computed_property_processor(self, properties: Dict[str, Any]):
        """Initialises the ComputedPropertyProcessor with access to all properties."""
        from materforge.parsing.processors.computed_property_processor import ComputedPropertyProcessor
        self.computed_property_processor = ComputedPropertyProcessor(properties, self.processed_properties)
        self.computed_property_processor.set_property_handler(self)

    def process_property(self, material: Material, prop_name: str,
                         config: Dict[str, Any],
                         dependency: sp.Symbol) -> None:
        """Processes a computed property via ComputedPropertyProcessor.

        Args:
            material:    Material instance.
            prop_name:   Name of the property.
            config:      Property configuration dict (passed through to processor).
            dependency:  SymPy symbol.
        """
        if self.computed_property_processor is None:
            raise ValueError("ComputedPropertyProcessor not initialised - call set_computed_property_processor first")
        self.computed_property_processor.process_computed_property(material, prop_name, dependency)

    def finalize_computed_property(self, material: Material, prop_name: str,
                                   dep_array: np.ndarray, prop_array: np.ndarray,
                                   dependency: sp.Symbol,
                                   config: Dict[str, Any]) -> None:
        """Finalises a computed property using the standard data-array path.

        Called by ComputedPropertyProcessor after evaluating the expression
        over the dependency range, ensuring consistent regression and
        visualization handling with all other property types.

        Args:
            material:    Material instance.
            prop_name:   Name of the property.
            dep_array:   Dependency values the expression was evaluated at.
            prop_array:  Corresponding evaluated property values.
            dependency:  SymPy symbol.
            config:      Property configuration dict.
        """
        self.finalize_with_data_arrays(material=material, prop_name=prop_name, dep_array=dep_array,
            prop_array=prop_array, dependency=dependency, config=config, prop_type='COMPUTED_PROPERTY')
