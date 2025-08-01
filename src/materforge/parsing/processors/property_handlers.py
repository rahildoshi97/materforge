import logging
from typing import Any, Dict, Union
import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.processors.property_processor_base import PropertyProcessorBase
from materforge.parsing.processors.temperature_resolver import TemperatureResolver
from materforge.parsing.io.data_handler import load_property_data
from materforge.parsing.utils.utilities import create_step_visualization_data
from materforge.parsing.validation.property_validator import validate_monotonic_energy_density
from materforge.algorithms.interpolation import ensure_ascending_order
from materforge.algorithms.piecewise_builder import PiecewiseBuilder
from materforge.parsing.config.yaml_keys import (
    DEPENDENCIES_KEY, VALUE_KEY, BOUNDS_KEY, FILE_PATH_KEY, EQUATION_KEY,
    RANGES_KEY, COLUMNS_KEY
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
                         prop_config: Union[float, str],
                         dependency_symbols: Dict[str, sp.Symbol]) -> None:
        """Process constant float property."""
        try:
            value = float(prop_config)
            prop_value = sp.Float(value)
            setattr(material, prop_name, prop_value)
            logger.debug(f"Set constant property {prop_name} = {value}")
            # Only visualize for symbolic dependencies
            if self._has_symbolic_dependencies(dependency_symbols):
                primary_symbol = list(dependency_symbols.values())[0]  # Use first symbol for visualization
                self._visualize_if_enabled(material=material, prop_name=prop_name,
                                           primary_symbol=primary_symbol,
                                           prop_type='CONSTANT_VALUE', x_data=None, y_data=None)
            else:
                logger.debug(f"Skipping visualization for constant property '{prop_name}' - numeric dependencies")
            self.processed_properties.add(prop_name)
        except (ValueError, TypeError) as e:
            logger.error(f"Failed to process constant property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process constant property \n -> {str(e)}") from e


class StepFunctionPropertyHandler(BasePropertyHandler):
    """Handler for step function properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any],
                         dependency_symbols: Dict[str, sp.Symbol]) -> None:
        """Process step function with multi-dependency support."""
        try:
            dependencies = prop_config[DEPENDENCIES_KEY]
            primary_dependency = dependencies[0]  # Step functions have single dependency
            # Get ranges and values
            ranges = prop_config[RANGES_KEY]
            temp_key = ranges[primary_dependency]
            val_array = prop_config[VALUE_KEY]
            # Resolve transition temperature
            transition_temp = TemperatureResolver.resolve_temperature_reference(temp_key, material)
            # Get the symbol for this dependency
            primary_symbol = dependency_symbols[primary_dependency]
            # Create step function using standard symbol
            T_standard = sp.Symbol('T')
            step_function = sp.Piecewise(
                (val_array[0], T_standard < transition_temp),
                (val_array[1], True)
            )
            # Replace with actual symbol
            if isinstance(primary_symbol, sp.Symbol):
                step_function = step_function.subs(T_standard, primary_symbol)
            # Create visualization data
            offset = ProcessingConstants.STEP_FUNCTION_OFFSET
            val1 = max(transition_temp - offset, PhysicalConstants.ABSOLUTE_ZERO)
            val2 = transition_temp + offset
            step_temp_array = np.array([val1, transition_temp, val2])
            x_data, y_data = create_step_visualization_data(transition_temp, val_array, step_temp_array)
            # Use piecewise finalization
            self.finalize_with_piecewise_function(
                material=material, prop_name=prop_name, piecewise_func=step_function,
                dependency_symbols=dependency_symbols, config=prop_config,
                prop_type='STEP_FUNCTION', x_data=x_data, y_data=y_data
            )
        except Exception as e:
            raise ValueError(f"Failed to process step function property '{prop_name}'\n -> {str(e)}") from e


class FileImportPropertyHandler(BasePropertyHandler):
    """Handler for file-based properties."""

    def process_property(self, material: Material, prop_name: str,
                         file_config: Dict[str, Any],
                         dependency_symbols: Dict[str, sp.Symbol]) -> None:
        """Process property data from a file configuration."""
        try:
            dependencies = file_config[DEPENDENCIES_KEY]
            primary_dependency = dependencies[0]  # File imports have single dependency
            # Setup file configuration
            file_path = self.base_dir / file_config[FILE_PATH_KEY]
            columns = file_config[COLUMNS_KEY]
            # Create load configuration
            load_config = {
                FILE_PATH_KEY: str(file_path),
                'dependency_column': columns[primary_dependency],
                'property_column': columns['property']
            }
            logger.debug(f"Loading property '{prop_name}' from file: {file_path}")
            temp_array, prop_array = load_property_data(load_config)
            logger.debug(f"Loaded {len(temp_array)} data points for property '{prop_name}' "
                         f"(range: {np.min(temp_array):.1f} - {np.max(temp_array):.1f})")
            validate_monotonic_energy_density(prop_name, temp_array, prop_array)
            # Use data array finalization
            self.finalize_with_data_arrays(
                material=material, prop_name=prop_name, temp_array=temp_array,
                prop_array=prop_array, dependency_symbols=dependency_symbols,
                config=file_config, prop_type='FILE_IMPORT'
            )
        except FileNotFoundError as e:
            logger.error(f"File not found for property '{prop_name}': {file_path}", exc_info=True)
            raise ValueError(f"File not found for property '{prop_name}': {file_path}") from e
        except Exception as e:
            logger.error(f"Failed to process file property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process file property {prop_name} \n -> {str(e)}") from e


class TabularDataPropertyHandler(BasePropertyHandler):
    """Handler for tabular data properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any],
                         dependency_symbols: Dict[str, sp.Symbol]) -> None:
        """Process property defined with tabular data."""
        try:
            dependencies = prop_config[DEPENDENCIES_KEY]
            primary_dependency = dependencies[0]  # Tabular data has single dependency
            # Get ranges and values
            ranges = prop_config[RANGES_KEY]
            temp_def = ranges[primary_dependency]
            val_array = prop_config[VALUE_KEY]
            # Resolve temperature array
            key_array = TemperatureResolver.resolve_temperature_definition(
                temp_def, len(val_array), material
            )
            if len(key_array) != len(val_array):
                raise ValueError(f"Length mismatch in {prop_name}: key and val arrays must have same length")
            key_array, val_array = ensure_ascending_order(key_array, val_array)
            validate_monotonic_energy_density(prop_name, key_array, val_array)
            # Use data array finalization
            self.finalize_with_data_arrays(
                material=material, prop_name=prop_name, temp_array=key_array,
                prop_array=val_array, dependency_symbols=dependency_symbols,
                config=prop_config, prop_type='TABULAR_DATA'
            )
        except Exception as e:
            raise ValueError(f"Failed to process tabular data property '{prop_name}' \n -> {str(e)}") from e


class PiecewiseEquationPropertyHandler(BasePropertyHandler):
    """Handler for piecewise equation properties."""

    def process_property(self, material: Material, prop_name: str,
                         prop_config: Dict[str, Any],
                         dependency_symbols: Dict[str, sp.Symbol]) -> None:
        """Process piecewise equation property."""
        try:
            dependencies = prop_config[DEPENDENCIES_KEY]
            primary_dependency = dependencies[0]  # Piecewise equations have single dependency
            eqn_strings = prop_config[EQUATION_KEY]
            ranges = prop_config[RANGES_KEY]
            temp_def = ranges[primary_dependency]
            # Resolve temperature points
            temp_points = TemperatureResolver.resolve_temperature_definition(temp_def, len(eqn_strings) + 1)
            # Get dependency symbol mapping from config
            from materforge.parsing.config.material_yaml_parser import MaterialYAMLParser
            parser_config = getattr(material, '_parser_config', None)
            if parser_config and 'independent_variables' in parser_config:
                symbol_mapping = parser_config['independent_variables']
            else:
                # Fallback: assume standard mapping
                symbol_mapping = {primary_dependency: 'T'}
            # Validate equations and replace symbols
            validated_equations = []
            for eqn in eqn_strings:
                expr = sp.sympify(eqn)

                # Check that only mapped symbols are used
                for symbol in expr.free_symbols:
                    symbol_str = str(symbol)
                    if symbol_str not in symbol_mapping.values():
                        raise ValueError(
                            f"Unsupported symbol '{symbol}' in equation '{eqn}' for property '{prop_name}'. "
                            f"Only symbols from independent_variables are allowed: {list(symbol_mapping.values())}"
                        )
                validated_equations.append(eqn)

            # Get bounds for primary dependency
            bounds_config = prop_config[BOUNDS_KEY]
            primary_bounds = bounds_config[primary_dependency]
            lower_bound_type, upper_bound_type = primary_bounds
            temp_points, validated_equations = ensure_ascending_order(temp_points, validated_equations)
            # Create piecewise function from formulas using the mapped symbol
            yaml_symbol_name = symbol_mapping[primary_dependency]
            T_yaml = sp.Symbol(yaml_symbol_name)
            piecewise_yaml = PiecewiseBuilder.build_from_formulas(
                temp_points, validated_equations, T_yaml, lower_bound_type, upper_bound_type
            )
            # Replace YAML symbol with actual dependency symbol
            actual_symbol = dependency_symbols[primary_dependency]
            if isinstance(actual_symbol, sp.Symbol):
                piecewise_func = piecewise_yaml.subs(T_yaml, actual_symbol)
            else:
                piecewise_func = piecewise_yaml
            # Create dense temperature array for visualization
            diff = max(np.min(np.diff(np.sort(temp_points))) / 10.0, 1.0)
            temp_dense = np.arange(temp_points[0], temp_points[-1] + diff / 2, diff)
            # Evaluate for visualization
            f_pw = sp.lambdify(T_yaml, piecewise_yaml, 'numpy')
            y_dense = f_pw(temp_dense)
            validate_monotonic_energy_density(prop_name, temp_dense, y_dense)
            # Use piecewise finalization
            self.finalize_with_piecewise_function(
                material=material, prop_name=prop_name, piecewise_func=piecewise_func,
                dependency_symbols=dependency_symbols, config=prop_config,
                prop_type='PIECEWISE_EQUATION', x_data=temp_dense, y_data=y_dense
            )
        except Exception as e:
            raise ValueError(f"Failed to process piecewise equation property '{prop_name}' \n -> {str(e)}") from e


class ComputedPropertyHandler(BasePropertyHandler):
    """Handler for computed properties."""

    def __init__(self):
        super().__init__()
        self.dependency_processor = None

    def set_dependency_processor(self, properties: Dict[str, Any]):
        """Set the dependency processor with access to all properties."""
        from materforge.parsing.processors.dependency_processor import DependencyProcessor
        self.dependency_processor = DependencyProcessor(properties, self.processed_properties)
        # Pass reference to this handler for finalization
        self.dependency_processor.set_property_handler(self)

    def process_property(self, material: Material, prop_name: str,
                         config: Dict[str, Any],
                         dependency_symbols: Dict[str, sp.Symbol]) -> None:
        """Process computed properties using dependency processor."""
        if self.dependency_processor is None:
            raise ValueError("Dependency processor not initialized")
        # Pass the dependency symbols to the dependency processor
        self.dependency_processor.process_computed_property(material, prop_name, dependency_symbols)

    def finalize_computed_property(self, material: Material, prop_name: str,
                                   temp_array: np.ndarray, prop_array: np.ndarray,
                                   dependency_symbols: Dict[str, sp.Symbol],
                                   config: Dict[str, Any]) -> None:
        """
        Public method to finalize computed property processing.

        This method provides a public interface to the property finalization,
        allowing the DependencyProcessor to properly finalize computed properties
        while maintaining consistency with other property handlers.
        """
        # Use data array finalization
        self.finalize_with_data_arrays(
            material=material, prop_name=prop_name, temp_array=temp_array,
            prop_array=prop_array, dependency_symbols=dependency_symbols,
            config=config, prop_type='COMPUTED_PROPERTY'
        )
