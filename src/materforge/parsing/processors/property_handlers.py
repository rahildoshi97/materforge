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
    DEPENDENCIES_KEY, RANGES_KEY, VALUE_KEY, BOUNDS_KEY, FILE_PATH_KEY, EQUATION_KEY,
    DEPENDENCY_KEY, COLUMNS_KEY, PROPERTY_KEY  # Legacy support
)
from materforge.data.constants import PhysicalConstants, ProcessingConstants

logger = logging.getLogger(__name__)


class BasePropertyHandler(PropertyProcessorBase):
    """Base class for property handlers with common functionality."""

    def __init__(self):
        super().__init__()
        logger.debug("BasePropertyHandler initialized")


class ConstantValuePropertyHandler(BasePropertyHandler):
    """Handler for constant properties."""

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Union[float, str]) -> None:
        """Process constant float property with multi-dependency support."""
        try:
            value = float(config)
            prop_value = sp.Float(value)
            setattr(material, prop_name, prop_value)
            logger.debug(f"Set constant property {prop_name} = {value}")

            # Only visualize for symbolic dependencies
            if hasattr(self, 'symbol_mapping') and self.symbol_mapping:
                primary_symbol = list(self.symbol_mapping.values())[0]
                if isinstance(primary_symbol, sp.Symbol):
                    self._visualize_if_enabled(material=material, prop_name=prop_name, T=primary_symbol,
                                               prop_type='CONSTANT_VALUE', x_data=None, y_data=None)

            self.processed_properties.add(prop_name)
        except (ValueError, TypeError) as e:
            logger.error(f"Failed to process constant property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process constant property -> {str(e)}") from e


class StepFunctionPropertyHandler(BasePropertyHandler):
    """Handler for step function properties."""

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Dict[str, Any]) -> None:
        """Process step function with multi-dependency support."""
        try:
            # Extract configuration for new format
            if DEPENDENCIES_KEY in config:
                dependencies = config[DEPENDENCIES_KEY]
                ranges_config = config.get(RANGES_KEY, {})

                if len(dependencies) != 1:
                    raise ValueError("Step functions must have exactly one dependency")

                dep_name = dependencies[0]
                if dep_name not in ranges_config:
                    raise ValueError(f"Missing range for dependency '{dep_name}'")

                # Resolve the dependency range (single transition point)
                resolved_ranges = self.dependency_resolver.resolve_dependency_ranges(config, material)
                temp_array = resolved_ranges[dep_name]

                # For step functions, expect a single transition point
                if len(temp_array) == 1:
                    transition_temp = temp_array[0]
                else:
                    # Take the first value if multiple provided
                    transition_temp = temp_array[0]

                # Get the symbol for this dependency
                dep_symbols = self.dependency_resolver.get_dependency_symbols(dependencies, self.symbol_mapping)
                yaml_symbol = self.dependency_resolver.independent_vars[dep_name]
                actual_symbol = dep_symbols[yaml_symbol]

            else:
                # Legacy format support
                temp_key = config[DEPENDENCY_KEY]
                transition_temp = self.dependency_resolver.resolve_dependency_reference(temp_key, material)
                actual_symbol = list(self.symbol_mapping.values())[0] if self.symbol_mapping else sp.Symbol('T')

            val_array = config[VALUE_KEY]

            # Create step function
            T_standard = sp.Symbol('T')
            step_function = sp.Piecewise((val_array[0], T_standard < transition_temp), (val_array[1], True))

            if str(actual_symbol) != 'T':
                step_function = step_function.subs(T_standard, actual_symbol)

            # Create visualization data
            offset = ProcessingConstants.STEP_FUNCTION_OFFSET
            val1 = max(transition_temp - offset, PhysicalConstants.ABSOLUTE_ZERO)
            val2 = transition_temp + offset
            step_temp_array = np.array([val1, transition_temp, val2])
            x_data, y_data = create_step_visualization_data(transition_temp, val_array, step_temp_array)

            # Set the property directly and visualize
            setattr(material, prop_name, step_function)

            # Generate visualization
            bounds = (np.min(x_data), np.max(x_data)) if x_data is not None else None
            self._visualize_if_enabled(
                material=material, prop_name=prop_name, T=actual_symbol,
                prop_type='STEP_FUNCTION', x_data=x_data, y_data=y_data,
                config=config, bounds=bounds
            )

            self.processed_properties.add(prop_name)

        except Exception as e:
            raise ValueError(f"Failed to process step function property '{prop_name}' -> {str(e)}") from e


class FileImportPropertyHandler(BasePropertyHandler):
    """Handler for file-based properties."""

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Dict[str, Any]) -> None:
        """Process file-based property with multi-dependency support."""
        try:
            file_path = self.base_dir / config[FILE_PATH_KEY]
            config[FILE_PATH_KEY] = str(file_path)
            logger.debug(f"Loading property '{prop_name}' from file: {file_path}")

            # Handle new multi-dependency file format
            if COLUMNS_KEY in config:
                # Multi-dependency file format
                dependencies = config.get(DEPENDENCIES_KEY, ['temperature'])
                columns = config[COLUMNS_KEY]

                # For now, handle single dependency case
                if len(dependencies) == 1:
                    dep_name = dependencies[0]
                    if dep_name not in columns:
                        raise ValueError(f"Column mapping missing for dependency '{dep_name}'")

                    # Create temporary config for load_property_data
                    temp_config = {
                        **config,
                        'dependency_column': columns[dep_name],
                        'property_column': columns.get(PROPERTY_KEY, columns.get('property'))
                    }
                    temp_array, prop_array = load_property_data(temp_config)
                else:
                    raise NotImplementedError("Multi-dependency file import not yet fully implemented")
            else:
                # Legacy format
                temp_array, prop_array = load_property_data(config)

            logger.debug(f"Loaded {len(temp_array)} data points for property '{prop_name}' "
                         f"(range: {np.min(temp_array):.1f} - {np.max(temp_array):.1f})")

            validate_monotonic_energy_density(prop_name, temp_array, prop_array)

            # Get primary symbol
            primary_symbol = list(self.symbol_mapping.values())[0] if self.symbol_mapping else sp.Symbol('T')

            # Use data array finalization
            self.finalize_with_data_arrays(
                material=material, prop_name=prop_name, temp_array=temp_array,
                prop_array=prop_array, T=primary_symbol, config=config, prop_type='FILE_IMPORT'
            )

        except FileNotFoundError as e:
            logger.error(f"File not found for property '{prop_name}': {file_path}", exc_info=True)
            raise ValueError(f"File not found for property '{prop_name}': {file_path}") from e
        except Exception as e:
            logger.error(f"Failed to process file property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process file property {prop_name} -> {str(e)}") from e


class TabularDataPropertyHandler(BasePropertyHandler):
    """Handler for tabular data properties."""

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Dict[str, Any]) -> None:
        """Process tabular data property with multi-dependency support."""
        try:
            # Handle new multi-dependency format
            if DEPENDENCIES_KEY in config:
                dependencies = config[DEPENDENCIES_KEY]
                resolved_ranges = self.dependency_resolver.resolve_dependency_ranges(config, material)

                # For now, handle single dependency case
                if len(dependencies) == 1:
                    dep_name = dependencies[0]
                    key_array = resolved_ranges[dep_name]
                    val_array = config[VALUE_KEY]

                    # Get the symbol for this dependency
                    dep_symbols = self.dependency_resolver.get_dependency_symbols(dependencies, self.symbol_mapping)
                    yaml_symbol = self.dependency_resolver.independent_vars[dep_name]
                    actual_symbol = dep_symbols[yaml_symbol]
                else:
                    raise NotImplementedError("Multi-dependency tabular data not yet fully implemented")
            else:
                # Legacy format
                temp_def = config[DEPENDENCY_KEY]
                val_array = config[VALUE_KEY]
                key_array = self.dependency_resolver.resolve_dependency_definition(temp_def, len(val_array), material)
                actual_symbol = list(self.symbol_mapping.values())[0] if self.symbol_mapping else sp.Symbol('T')

            if len(key_array) != len(val_array):
                raise ValueError(f"Length mismatch in {prop_name}: key and val arrays must have same length")

            key_array, val_array = ensure_ascending_order(key_array, val_array)
            validate_monotonic_energy_density(prop_name, key_array, val_array)

            # Use data array finalization
            self.finalize_with_data_arrays(
                material=material, prop_name=prop_name, temp_array=key_array,
                prop_array=val_array, T=actual_symbol, config=config, prop_type='TABULAR_DATA'
            )

        except Exception as e:
            raise ValueError(f"Failed to process tabular data property '{prop_name}' -> {str(e)}") from e


class PiecewiseEquationPropertyHandler(BasePropertyHandler):
    """Handler for piecewise equation properties."""

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Dict[str, Any]) -> None:
        """Process piecewise equation property with multi-dependency support."""
        try:
            eqn_strings = config[EQUATION_KEY]

            # Handle new multi-dependency format
            if DEPENDENCIES_KEY in config:
                dependencies = config[DEPENDENCIES_KEY]
                resolved_ranges = self.dependency_resolver.resolve_dependency_ranges(config, material)

                # For now, handle single dependency case
                if len(dependencies) == 1:
                    dep_name = dependencies[0]
                    temp_points = resolved_ranges[dep_name]

                    # Get the symbol for this dependency
                    dep_symbols = self.dependency_resolver.get_dependency_symbols(dependencies, self.symbol_mapping)
                    yaml_symbol = self.dependency_resolver.independent_vars[dep_name]
                    actual_symbol = dep_symbols[yaml_symbol]
                else:
                    raise NotImplementedError("Multi-dependency piecewise equations not yet fully implemented")
            else:
                # Legacy format
                temp_def = config[DEPENDENCY_KEY]
                temp_points = self.dependency_resolver.resolve_dependency_definition(temp_def, len(eqn_strings) + 1, material)
                actual_symbol = list(self.symbol_mapping.values())[0] if self.symbol_mapping else sp.Symbol('T')

            # Validate equations
            for eqn in eqn_strings:
                expr = sp.sympify(eqn)
                for symbol in expr.free_symbols:
                    if str(symbol) not in [str(actual_symbol), 'T']:  # Allow both T and actual symbol
                        raise ValueError(
                            f"Unsupported symbol '{symbol}' in equation '{eqn}' for property '{prop_name}'. "
                            f"Only 'T' is allowed.")

            # Extract bounds
            bounds_config = config.get(BOUNDS_KEY, ['constant', 'constant'])
            if isinstance(bounds_config, dict):
                # Multi-dependency bounds format
                dep_name = dependencies[0]
                lower_bound_type, upper_bound_type = bounds_config.get(dep_name, ['constant', 'constant'])
            else:
                # Legacy format
                lower_bound_type, upper_bound_type = bounds_config

            temp_points, eqn_strings = ensure_ascending_order(temp_points, eqn_strings)

            # Create piecewise function from formulas
            T_standard = sp.Symbol('T')
            piecewise_standard = PiecewiseBuilder.build_from_formulas(
                temp_points, list(eqn_strings), T_standard, lower_bound_type, upper_bound_type
            )

            # Substitute the actual symbol if different
            if str(actual_symbol) != 'T':
                piecewise_func = piecewise_standard.subs(T_standard, actual_symbol)
            else:
                piecewise_func = piecewise_standard

            # Create dense temperature array for visualization
            diff = max(np.min(np.diff(np.sort(temp_points))) / 10.0, 1.0)
            temp_dense = np.arange(temp_points[0], temp_points[-1] + diff / 2, diff)

            # Evaluate for visualization
            f_pw = sp.lambdify(T_standard, piecewise_standard, 'numpy')
            y_dense = f_pw(temp_dense)
            validate_monotonic_energy_density(prop_name, temp_dense, y_dense)

            # Set property and visualize
            setattr(material, prop_name, piecewise_func)

            # Generate visualization
            bounds = (np.min(temp_dense), np.max(temp_dense))
            self._visualize_if_enabled(
                material=material, prop_name=prop_name, T=actual_symbol,
                prop_type='PIECEWISE_EQUATION', x_data=temp_dense, y_data=y_dense,
                config=config, bounds=bounds
            )

            self.processed_properties.add(prop_name)

        except Exception as e:
            raise ValueError(f"Failed to process piecewise equation property '{prop_name}' -> {str(e)}") from e


class ComputedPropertyHandler(BasePropertyHandler):
    """Handler for computed properties."""

    def __init__(self):
        super().__init__()
        self.dependency_processor = None

    def set_dependency_processor(self, properties: Dict[str, Any], dependency_resolver: DependencyResolver):
        """Set the dependency processor with access to all properties and dependency resolver."""
        from materforge.parsing.processors.dependency_processor import DependencyProcessor
        self.dependency_processor = DependencyProcessor(properties, self.processed_properties)
        # Pass reference to this handler for finalization
        self.dependency_processor.set_property_handler(self)
        # Store dependency resolver for multi-dependency support
        self.dependency_processor.set_dependency_resolver(dependency_resolver)

    def process_multi_dependency_property(self, material: Material, prop_name: str, config: Dict[str, Any]) -> None:
        """Process computed properties using dependency processor with multi-dependency support."""
        if self.dependency_processor is None:
            raise ValueError("Dependency processor not initialized")

        # Process the computed property
        self.dependency_processor.process_computed_property(material, prop_name, config)

    def finalize_computed_property(self, material: Material, prop_name: str,
                                   temp_array: np.ndarray, prop_array: np.ndarray,
                                   config: Dict[str, Any]) -> None:
        """
        Public method to finalize computed property processing.
        """
        # Get primary symbol
        primary_symbol = list(self.symbol_mapping.values())[0] if self.symbol_mapping else sp.Symbol('T')

        # Use data array finalization
        self.finalize_with_data_arrays(
            material=material,
            prop_name=prop_name,
            temp_array=temp_array,
            prop_array=prop_array,
            T=primary_symbol,
            config=config,
            prop_type='COMPUTED_PROPERTY'
        )
