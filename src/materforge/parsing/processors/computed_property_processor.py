import logging
from typing import Dict, List, Set, Union, Any

import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.core.symbol_registry import SymbolRegistry
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.validation.property_validator import validate_monotonic_energy_density
from materforge.parsing.validation.errors import DependencyError, CircularDependencyError
from materforge.parsing.config.yaml_keys import EQUATION_KEY, DEPENDENCY_KEY
from materforge.parsing.utils.utilities import handle_numeric_temperature

logger = logging.getLogger(__name__)


class ComputedPropertyProcessor:
    """Handles dependency resolution and computed property processing."""

    def __init__(self, properties: Dict[str, Any], processed_properties: Set[str]):
        self.properties = properties
        self.processed_properties = processed_properties
        # Store reference to the property handler for finalization
        self.property_handler = None

    def set_property_handler(self, property_handler):
        """Set reference to the property handler for finalization."""
        self.property_handler = property_handler

    def process_computed_property(self, material: Material, prop_name: str,
                                  T: Union[float, sp.Symbol]) -> None:
        """Process computed properties using predefined models with dependency checking."""
        if prop_name in self.processed_properties:
            logger.debug(f"Property '{prop_name}' already processed, skipping")
            return
        try:
            logger.debug(f"Processing computed property: {prop_name}")
            prop_config = self.properties[prop_name]
            if not isinstance(prop_config, dict) or EQUATION_KEY not in prop_config:
                raise ValueError(f"Invalid COMPUTE property configuration for {prop_name}")
            temp_def = prop_config[DEPENDENCY_KEY]
            temp_array = DependencyResolver.resolve_dependency_definition(temp_def, material=material)
            expression = prop_config[EQUATION_KEY]
            logger.debug(f"Computing property '{prop_name}' with expression: {expression}")
            try:
                material_property = self._parse_and_process_expression(expression, material, T, prop_name)
                logger.debug(f"Successfully parsed expression for property '{prop_name}'")
            except CircularDependencyError:
                raise  # Re-raise without wrapping
            except Exception as e:
                logger.error(f"Failed to parse expression for property '{prop_name}': {e}", exc_info=True)
                raise ValueError(f"Failed to process computed property '{prop_name}' \n -> {str(e)}") from e
            # Handle numeric temperature case using utility function
            if handle_numeric_temperature(self, material, prop_name, material_property, T):
                return
            # Evaluate the expression for symbolic temperature
            T_standard = sp.Symbol('T')
            if isinstance(T, sp.Symbol) and str(T) != 'T':
                standard_expr = material_property.subs(T, T_standard)
                f_pw = sp.lambdify(T_standard, standard_expr, 'numpy')
            else:
                f_pw = sp.lambdify(T_standard, material_property, 'numpy')  # Use T_standard instead of T
            try:
                y_dense = f_pw(temp_array)
                if not np.all(np.isfinite(y_dense)):
                    invalid_count = np.sum(~np.isfinite(y_dense))
                    logger.warning(f"Property '{prop_name}' has {invalid_count} non-finite values. "
                                   f"This may indicate issues with the expression: {expression}")
                validate_monotonic_energy_density(prop_name, temp_array, y_dense)
                if self.property_handler is not None:
                    # Use the property processor's finalization method for consistent handling
                    self.property_handler.finalize_computed_property(
                        material, prop_name, temp_array, y_dense, T, prop_config
                    )
                else:
                    # Fallback: Set property directly (no visualization)
                    setattr(material, prop_name, material_property)
                    self.processed_properties.add(prop_name)
                    logger.warning(f"Property processor not available for '{prop_name}' - skipping visualization")
                logger.debug(f"Successfully computed property '{prop_name}' over {len(temp_array)} temperature points")
            except Exception as e:
                logger.error(f"Error evaluating expression for property '{prop_name}': {e}", exc_info=True)
                raise ValueError(f"Error evaluating expression for '{prop_name}' \n -> {str(e)}") from e
        except (DependencyError, CircularDependencyError):
            raise  # Re-raise without wrapping
        except Exception as e:
            logger.error(f"Failed to process computed property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process computed property '{prop_name}' \n -> {str(e)}") from e

    def _parse_and_process_expression(self, expression: str, material: Material,
                                      T: Union[float, sp.Symbol], prop_name: str) -> sp.Expr:
        """Parse and process a mathematical expression string into a SymPy expression."""
        try:
            logger.debug(f"Parsing expression for '{prop_name}': {expression}")
            T_standard = sp.Symbol('T')
            sympy_expr = sp.sympify(expression, evaluate=False)
            # Extract dependencies
            dependencies = [str(symbol) for symbol in sympy_expr.free_symbols if str(symbol) != 'T']
            if dependencies:
                logger.debug(f"Property '{prop_name}' depends on: {dependencies}")
                # Check for missing dependencies
                missing_deps = []
                for dep in dependencies:
                    if not hasattr(material, dep) and dep not in self.properties:
                        missing_deps.append(dep)
                if missing_deps:
                    available_props = sorted(list(self.properties.keys()))
                    logger.error(f"Missing dependencies for '{prop_name}': {missing_deps}. "
                                 f"Available: {available_props}")
                    raise DependencyError(expression=expression, missing_deps=missing_deps,
                                          available_props=available_props)
                # Check for circular dependencies
                self._validate_circular_dependencies(prop_name, dependencies, set())
                # Process dependencies first
                for dep in dependencies:
                    if not hasattr(material, dep) or getattr(material, dep) is None:
                        if dep in self.properties:
                            logger.debug(f"Processing dependency '{dep}' for property '{prop_name}'")
                            self.process_computed_property(material, dep, T)
                        else:
                            available_props = sorted(list(self.properties.keys()))
                            raise DependencyError(expression=expression, missing_deps=[dep],
                                                  available_props=available_props)
            # Verify all dependencies are now available
            missing_deps = [dep for dep in dependencies if not hasattr(material, dep) or getattr(material, dep) is None]
            if missing_deps:
                raise ValueError(f"Cannot compute expression. Missing dependencies: {missing_deps}")
            # Create substitution dictionary
            substitutions = {}
            for dep in dependencies:
                dep_value = getattr(material, dep, None)
                if dep_value is None:
                    logger.error(f"Dependency '{dep}' was processed but is still not available")
                    raise ValueError(f"Dependency '{dep}' was processed but is still not available on the material")
                dep_symbol = SymbolRegistry.get(dep)
                if dep_symbol is None:
                    raise ValueError(f"Symbol '{dep}' not found in symbol registry")
                substitutions[dep_symbol] = dep_value
            # Handle temperature substitution based on type
            if isinstance(T, sp.Symbol):  # If T is a symbolic variable, substitute the standard 'T' with it
                substitutions[T_standard] = T
            else:  # For numeric T, substitute with the value directly
                substitutions[T_standard] = T
            # Perform substitution and evaluate integrals
            result_expr = sympy_expr.subs(substitutions)
            if isinstance(result_expr, sp.Integral):
                logger.debug(f"Evaluating integral in expression for '{prop_name}'")
                result_expr = result_expr.doit()
            logger.debug(f"Successfully processed expression for '{prop_name}'")
            return result_expr
        except CircularDependencyError:
            raise  # Re-raise without wrapping
        except DependencyError:
            raise  # Re-raise without wrapping
        except Exception as e:
            logger.error(f"Failed to parse expression '{expression}' for property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to parse and process expression: {expression}") from e

    def _validate_circular_dependencies(self, prop_name: str, current_deps: List[str],
                                        visited: Set[str], path: List[str] = None) -> None:
        """
        Check for circular dependencies in property definitions.
        Args:
            prop_name: Current property being checked
            current_deps: Dependencies of the current property
            visited: Set of already visited properties
            path: Current dependency path for error reporting
        Raises:
            CircularDependencyError: If circular dependency is detected
        """
        if path is None:
            path = []
        # Filter out 'T' from dependencies
        current_deps = [dep for dep in current_deps if dep != 'T']
        if prop_name is not None:
            if prop_name in visited:
                cycle_path = path + [prop_name]
                logger.error(f"Circular dependency detected: {' -> '.join(cycle_path)}")
                raise CircularDependencyError(dependency_path=cycle_path)
            visited.add(prop_name)
            path = path + [prop_name]
        for dep in current_deps:
            if dep in self.properties:
                dep_config = self.properties[dep]
                if isinstance(dep_config, dict) and EQUATION_KEY in dep_config:
                    dep_deps = self._extract_equation_dependencies(dep_config[EQUATION_KEY])
                    if dep_deps:
                        self._validate_circular_dependencies(dep, dep_deps, visited.copy(), path)

    @staticmethod
    def _extract_equation_dependencies(equation_data) -> List[str]:
        """Extract dependencies from equation data."""
        symbols = set()
        if isinstance(equation_data, list):
            for eq in equation_data:
                expr = sp.sympify(eq)
                symbols.update(expr.free_symbols)
        else:
            expr = sp.sympify(equation_data)
            symbols.update(expr.free_symbols)
        return [str(symbol) for symbol in symbols if str(symbol) != 'T']
