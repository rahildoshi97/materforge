import logging
from typing import Dict, List, Optional, Set, Union, Any
import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.core.symbol_registry import SymbolRegistry
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.validation.property_validator import validate_monotonic_energy_density
from materforge.parsing.validation.errors import DependencyError, CircularDependencyError
from materforge.parsing.config.yaml_keys import EQUATION_KEY, DEPENDENCIES_KEY, RANGES_KEY
from materforge.parsing.utils.utilities import handle_numeric_temperature

logger = logging.getLogger(__name__)

class DependencyProcessor:
    """Handles dependency resolution and computed property processing with multi-dependency support."""

    def __init__(self, properties: Dict[str, Any], processed_properties: Set[str]):
        self.properties = properties
        self.processed_properties = processed_properties
        self.dependency_resolver: Optional[DependencyResolver] = None
        # Store reference to the property handler for finalization
        self.property_handler = None

    def set_property_handler(self, property_handler):
        """Set reference to the property handler for finalization."""
        self.property_handler = property_handler

    def set_dependency_resolver(self, dependency_resolver: DependencyResolver):
        """Set the dependency resolver for multi-dependency support."""
        self.dependency_resolver = dependency_resolver

    def process_computed_property(self, material: Material, prop_name: str, config: Dict[str, Any]) -> None:
        """Process computed properties using predefined models with multi-dependency support."""
        if prop_name in self.processed_properties:
            logger.debug(f"Property '{prop_name}' already processed, skipping")
            return

        try:
            logger.debug(f"Processing computed property: {prop_name}")

            if not isinstance(config, dict) or EQUATION_KEY not in config:
                raise ValueError(f"Invalid COMPUTE property configuration for {prop_name}")

            # Handle multi-dependency format
            if DEPENDENCIES_KEY in config:
                # New multi-dependency format
                dependencies = config[DEPENDENCIES_KEY]
                ranges_config = config.get(RANGES_KEY, {})

                # Resolve dependency ranges
                resolved_ranges = self.dependency_resolver.resolve_dependency_ranges(config, material)

                # Get dependency symbols
                dep_symbols = self.dependency_resolver.get_dependency_symbols(dependencies, self.property_handler.symbol_mapping)

            else:
                raise ValueError(f"Missing '{DEPENDENCIES_KEY}' in configuration for property '{prop_name}'")

            expression = config[EQUATION_KEY]
            logger.debug(f"Computing property '{prop_name}' with expression: {expression}")

            try:
                material_property = self._parse_and_process_expression(
                    expression, material, dep_symbols, prop_name, dependencies
                )
                logger.debug(f"Successfully parsed expression for property '{prop_name}'")

            except CircularDependencyError:
                raise # Re-raise without wrapping
            except Exception as e:
                logger.error(f"Failed to parse expression for property '{prop_name}': {e}", exc_info=True)
                raise ValueError(f"Failed to process computed property '{prop_name}' \n -> {str(e)}") from e

            # For now, handle single dependency case
            if len(dependencies) == 1:
                primary_dep = dependencies[0]
                primary_array = resolved_ranges[primary_dep]
                primary_symbol = list(dep_symbols.values())[0]

                # Evaluate the expression
                if isinstance(primary_symbol, sp.Symbol):
                    T_standard = sp.Symbol('T')
                    if str(primary_symbol) != 'T':
                        standard_expr = material_property.subs(primary_symbol, T_standard)
                        f_pw = sp.lambdify(T_standard, standard_expr, 'numpy')
                    else:
                        f_pw = sp.lambdify(T_standard, material_property, 'numpy')

                    try:
                        y_dense = f_pw(primary_array)
                        if not np.all(np.isfinite(y_dense)):
                            invalid_count = np.sum(~np.isfinite(y_dense))
                            logger.warning(f"Property '{prop_name}' has {invalid_count} non-finite values. "
                                           f"This may indicate issues with the expression: {expression}")

                        validate_monotonic_energy_density(prop_name, primary_array, y_dense)

                        if self.property_handler is not None:
                            # Use the property processor's finalization method for consistent handling
                            self.property_handler.finalize_computed_property(
                                material, prop_name, primary_array, y_dense, config
                            )
                        else:
                            # Fallback: Set property directly (no visualization)
                            setattr(material, prop_name, material_property)
                            self.processed_properties.add(prop_name)
                            logger.warning(f"Property processor not available for '{prop_name}' - skipping visualization")

                        logger.debug(f"Successfully computed property '{prop_name}' over {len(primary_array)} points")

                    except Exception as e:
                        logger.error(f"Error evaluating expression for property '{prop_name}': {e}", exc_info=True)
                        raise ValueError(f"Error evaluating expression for '{prop_name}' \n -> {str(e)}") from e
                else:
                    # Handle numeric case
                    if self.property_handler and hasattr(self.property_handler, 'handle_numeric_temperature'):
                        if self.property_handler.handle_numeric_temperature(material, prop_name, material_property, primary_symbol):
                            return

            else:
                # Multi-dependency case - for now, set the symbolic expression directly
                setattr(material, prop_name, material_property)
                self.processed_properties.add(prop_name)
                logger.debug(f"Set multi-dependency property '{prop_name}' as symbolic expression")

        except (DependencyError, CircularDependencyError):
            raise # Re-raise without wrapping
        except Exception as e:
            logger.error(f"Failed to process computed property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to process computed property '{prop_name}' \n -> {str(e)}") from e

    def _parse_and_process_expression(self, expression: str, material: Material,
                                      dep_symbols: Dict[str, sp.Symbol], prop_name: str,
                                      dependencies: List[str]) -> sp.Expr:
        """Parse and process a mathematical expression string into a SymPy expression."""
        try:
            logger.debug(f"Parsing expression for '{prop_name}': {expression}")

            sympy_expr = sp.sympify(expression, evaluate=False)

            # Extract dependencies
            expr_symbols = [str(symbol) for symbol in sympy_expr.free_symbols]
            property_deps = [sym for sym in expr_symbols if sym not in dep_symbols.keys() if sym != 'T']

            if property_deps:
                logger.debug(f"Property '{prop_name}' depends on: {property_deps}")

                # Check for missing dependencies
                missing_deps = []
                for dep in property_deps:
                    if not hasattr(material, dep) and dep not in self.properties:
                        missing_deps.append(dep)

                if missing_deps:
                    available_props = sorted(list(self.properties.keys()))
                    logger.error(f"Missing dependencies for '{prop_name}': {missing_deps}. "
                                 f"Available: {available_props}")
                    raise DependencyError(expression=expression, missing_deps=missing_deps,
                                          available_props=available_props)

                # Check for circular dependencies
                self._validate_circular_dependencies(prop_name, property_deps, set())

                # Process dependencies first
                for dep in property_deps:
                    if not hasattr(material, dep) or getattr(material, dep) is None:
                        if dep in self.properties:
                            logger.debug(f"Processing dependency '{dep}' for property '{prop_name}'")
                            self.process_computed_property(material, dep, self.properties[dep])
                        else:
                            available_props = sorted(list(self.properties.keys()))
                            raise DependencyError(expression=expression, missing_deps=[dep],
                                                  available_props=available_props)

                # Verify all dependencies are now available
                missing_deps = [dep for dep in property_deps if not hasattr(material, dep) or getattr(material, dep) is None]
                if missing_deps:
                    raise ValueError(f"Cannot compute expression. Missing dependencies: {missing_deps}")

                # Create substitution dictionary
                substitutions = {}
                for dep in property_deps:
                    dep_value = getattr(material, dep, None)
                    if dep_value is None:
                        logger.error(f"Dependency '{dep}' was processed but is still not available")
                        raise ValueError(f"Dependency '{dep}' was processed but is still not available on the material")

                    dep_symbol = SymbolRegistry.get(dep)
                    if dep_symbol is None:
                        raise ValueError(f"Symbol '{dep}' not found in symbol registry")

                    substitutions[dep_symbol] = dep_value

                # Perform substitution and evaluate integrals
                result_expr = sympy_expr.subs(substitutions)
            else:
                result_expr = sympy_expr

            if isinstance(result_expr, sp.Integral):
                logger.debug(f"Evaluating integral in expression for '{prop_name}'")
                result_expr = result_expr.doit()

            logger.debug(f"Successfully processed expression for '{prop_name}'")
            return result_expr

        except CircularDependencyError:
            raise # Re-raise without wrapping
        except DependencyError:
            raise # Re-raise without wrapping
        except Exception as e:
            logger.error(f"Failed to parse expression '{expression}' for property '{prop_name}': {e}", exc_info=True)
            raise ValueError(f"Failed to parse and process expression: {expression}") from e

    def _validate_circular_dependencies(self, prop_name: str, current_deps: List[str],
                                        visited: Set[str], path: List[str] = None) -> None:
        """Check for circular dependencies in property definitions."""
        if path is None:
            path = []

        # Filter out dependency symbols from dependencies
        current_deps = [dep for dep in current_deps if dep != 'T' and dep not in ['P', 'S']]

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

        return [str(symbol) for symbol in symbols if str(symbol) not in ['T', 'P', 'S']]
