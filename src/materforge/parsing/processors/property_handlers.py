import logging
from typing import Any, Dict, Union
import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.processors.property_processor_base import PropertyProcessorBase
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.io.data_handler import load_property_data
from materforge.parsing.utils.utilities import create_step_visualization_data
from materforge.parsing.validation.property_validator import validate_monotonic_energy_density
from materforge.algorithms.interpolation import ensure_ascending_order
from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.parsing.config.yaml_keys import (
    DEPENDENCY_KEY, VALUE_KEY, BOUNDS_KEY, FILE_PATH_KEY, EQUATION_KEY
)
from materforge.data.constants import PhysicalConstants, ProcessingConstants

logger = logging.getLogger(__name__)


class BasePropertyHandler(PropertyProcessorBase):
    """
    Base class for property handlers with common functionality.

    This class inherits from PropertyProcessorBase to provide shared functionality
    for processing material properties. All specialized handlers inherit from this class.
    """

    def __init__(self):
        super().__init__()
        logger.debug("BasePropertyHandler initialized")


class ConstantValuePropertyHandler(BasePropertyHandler):
    """Handler for constant properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Union[float, str], T: Union[float, sp.Symbol]) -> None:
        """Process constant float property."""
        try:
            value = float(prop_config)
            prop_value = sp.Float(value)
            setattr(material, prop_name, prop_value)
            logger.debug(f"Set constant property {prop_name} = {value}")
            # Only visualize for symbolic temperature
            if isinstance(T, sp.Symbol):
                self._visualize_if_enabled(material=material, prop_name=prop_name, T=T,
                                           prop_type='CONSTANT_VALUE', x_data=None, y_data=None)
            else:
                logger.debug(f"Skipping visualization for constant property '{prop_name}' - numeric temperature")
            self.processed_properties.add(prop_name)
        except (ValueError, TypeError) as e:
            logger.error(f"Failed to process constant property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process constant property \n -> {str(e)}") from e


class StepFunctionPropertyHandler(BasePropertyHandler):
    """Handler for step function properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any], T: Union[float, sp.Symbol]) -> None:
        """Process step function with unified symbol handling."""
        try:
            temp_key = prop_config[DEPENDENCY_KEY]
            val_array = prop_config[VALUE_KEY]
            transition_temp = DependencyResolver.resolve_dependency_reference(temp_key, material)
            T_standard = sp.Symbol('T')
            step_function = sp.Piecewise((val_array[0], T_standard < transition_temp), (val_array[1], True))
            if str(T) != 'T':
                step_function = step_function.subs(T_standard, T)
            # Create visualization data
            offset = ProcessingConstants.STEP_FUNCTION_OFFSET
            val1 = max(transition_temp - offset, PhysicalConstants.ABSOLUTE_ZERO)
            val2 = transition_temp + offset
            step_temp_array = np.array([val1, transition_temp, val2])
            x_data, y_data = create_step_visualization_data(transition_temp, val_array, step_temp_array)
            # Use piecewise finalization
            self.finalize_with_piecewise_function(material=material, prop_name=prop_name, piecewise_func=step_function,
                                                  T=T, config=prop_config, prop_type='STEP_FUNCTION',
                                                  x_data=x_data, y_data=y_data)
        except Exception as e:
            raise ValueError(f"Failed to process step function property '{prop_name}'\n -> {str(e)}") from e


class FileImportPropertyHandler(BasePropertyHandler):
    """Handler for file-based properties."""

    def process_property(self, material: Material, prop_name: str,
                         file_config: Dict[str, Any], T: Union[float, sp.Symbol]) -> None:
        """Process property data from a file configuration."""
        try:
            file_path = self.base_dir / file_config[FILE_PATH_KEY]
            file_config[FILE_PATH_KEY] = str(file_path)
            logger.debug(f"Loading property '{prop_name}' from file: {file_path}")
            temp_array, prop_array = load_property_data(file_config)
            logger.debug(f"Loaded {len(temp_array)} data points for property '{prop_name}' "
                         f"(T range: {np.min(temp_array):.1f}K - {np.max(temp_array):.1f}K)")
            validate_monotonic_energy_density(prop_name, temp_array, prop_array)
            # Use data array finalization
            self.finalize_with_data_arrays(material=material, prop_name=prop_name, temp_array=temp_array,
                                           prop_array=prop_array, T=T, config=file_config, prop_type='FILE_IMPORT')
        except FileNotFoundError as e:
            logger.error(f"File not found for property '{prop_name}': {file_path}", exc_info=True)
            raise ValueError(f"File not found for property '{prop_name}': {file_path}") from e
        except Exception as e:
            logger.error(f"Failed to process file property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process file property {prop_name} \n -> {str(e)}") from e


class TabularDataPropertyHandler(BasePropertyHandler):
    """Handler for key-value properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any], T: Union[float, sp.Symbol]) -> None:
        """Process property defined with key-val pairs."""
        try:
            temp_def = prop_config[DEPENDENCY_KEY]
            val_array = prop_config[VALUE_KEY]
            key_array = DependencyResolver.resolve_dependency_definition(temp_def, len(val_array), material)
            if len(key_array) != len(val_array):
                raise ValueError(f"Length mismatch in {prop_name}: key and val arrays must have same length")
            key_array, val_array = ensure_ascending_order(key_array, val_array)
            validate_monotonic_energy_density(prop_name, key_array, val_array)
            # Use data array finalization
            self.finalize_with_data_arrays(material=material, prop_name=prop_name, temp_array=key_array,
                                           prop_array=val_array, T=T, config=prop_config, prop_type='TABULAR_DATA')
        except Exception as e:
            raise ValueError(f"Failed to process key-val property '{prop_name}' \n -> {str(e)}") from e


class PiecewiseEquationPropertyHandler(BasePropertyHandler):
    """Handler for piecewise equation properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any], T: Union[float, sp.Symbol]) -> None:
        """Process piecewise equation property."""
        try:
            eqn_strings = prop_config[EQUATION_KEY]
            temp_def = prop_config[DEPENDENCY_KEY]
            temp_points = DependencyResolver.resolve_dependency_definition(temp_def, len(eqn_strings) + 1)
            # Validate equations
            for eqn in eqn_strings:
                expr = sp.sympify(eqn)
                for symbol in expr.free_symbols:
                    if str(symbol) != 'T':
                        raise ValueError(
                            f"Unsupported symbol '{symbol}' in equation '{eqn}' for property '{prop_name}'. "
                            f"Only 'T' is allowed.")
            lower_bound_type, upper_bound_type = prop_config[BOUNDS_KEY]
            temp_points, eqn_strings = ensure_ascending_order(temp_points, eqn_strings)
            # Create piecewise function from formulas
            T_standard = sp.Symbol('T')
            piecewise_standard = PiecewiseBuilder.build_from_formulas(temp_points, list(eqn_strings), T_standard,
                                                                      lower_bound_type, upper_bound_type)
            # Substitute the actual temperature symbol if different
            if str(T) != 'T':
                piecewise_func = piecewise_standard.subs(T_standard, T)
            else:
                piecewise_func = piecewise_standard
            # Create dense temperature array for visualization
            diff = max(np.min(np.diff(np.sort(temp_points))) / 10.0, 1.0)
            temp_dense = np.arange(temp_points[0], temp_points[-1] + diff / 2, diff)
            # Evaluate for visualization
            f_pw = sp.lambdify(T_standard, piecewise_standard, 'numpy')
            y_dense = f_pw(temp_dense)
            validate_monotonic_energy_density(prop_name, temp_dense, y_dense)
            # Use piecewise finalization
            self.finalize_with_piecewise_function(material=material, prop_name=prop_name, piecewise_func=piecewise_func,
                                                  T=T, config=prop_config, prop_type='PIECEWISE_EQUATION',
                                                  x_data=temp_dense, y_data=y_dense)
        except Exception as e:
            raise ValueError(f"Failed to process piecewise equation property '{prop_name}' \n -> {str(e)}") from e


class ComputedPropertyHandler(BasePropertyHandler):
    """Handler for computed properties."""

    def __init__(self):
        super().__init__()
        self.computed_property_processor = None

    def set_computed_property_processor(self, properties: Dict[str, Any]):
        """Set the computed property processor with access to all properties."""
        from materforge.parsing.processors.computed_property_processor import ComputedPropertyProcessor
        self.computed_property_processor = ComputedPropertyProcessor(properties, self.processed_properties)
        # Pass reference to this handler for finalization
        self.computed_property_processor.set_property_handler(self)

    def process_property(self, material: Material, prop_name: str,
                         config: Dict[str, Any], T: Union[float, sp.Symbol]) -> None:
        """Process computed properties using dependency processor."""
        if self.computed_property_processor is None:
            raise ValueError("Dependency processor not initialized")
        # Pass the config to the dependency processor if needed, or ignore it
        self.computed_property_processor.process_computed_property(material, prop_name, T)

    def finalize_computed_property(self, material: Material, prop_name: str,
                                   temp_array: np.ndarray, prop_array: np.ndarray,
                                   T: Union[float, sp.Symbol], config: Dict[str, Any]) -> None:
        """
        Public method to finalize computed property processing.

        This method provides a public interface to the protected _finalize_property_processing
        method, allowing the ComputedPropertyProcessor to properly finalize computed properties
        while maintaining consistency with other property handlers.
        """
        # Use data array finalization
        self.finalize_with_data_arrays(material=material, prop_name=prop_name, temp_array=temp_array,
                                       prop_array=prop_array, T=T, config=config, prop_type='COMPUTED_PROPERTY')
