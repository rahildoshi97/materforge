# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from typing import Dict, List, Set, Any
import numpy as np
import sympy as sp
from materforge.core.materials import Material
from materforge.core.symbol_registry import SymbolRegistry
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.parsing.validation.errors import DependencyError, CircularDependencyError
from materforge.parsing.config.yaml_keys import EQUATION_KEY, DEPENDENCY_KEY, YAML_PLACEHOLDER

logger = logging.getLogger(__name__)


class ComputedPropertyProcessor:
    """Handles dependency resolution and computed property processing."""

    def __init__(self, properties: Dict[str, Any], processed_properties: Set[str]):
        self.properties = properties
        self.processed_properties = processed_properties
        self.property_handler = None

    def set_property_handler(self, property_handler):
        """Sets reference to the property handler for finalization."""
        self.property_handler = property_handler

    def process_computed_property(self, material: Material, prop_name: str,
                                  dependency: sp.Symbol) -> None:
        """Processes a computed property by evaluating its equation over the
        resolved dependency range.

        Args:
            material:    Material instance to assign the result to.
            prop_name:   Name of the property to compute.
            dependency:  SymPy symbol (symbolic mode).
        """
        if not isinstance(dependency, sp.Basic):
            raise TypeError(f"dependency must be a SymPy expression, got {type(dependency).__name__}")
        if prop_name in self.processed_properties:
            logger.debug("Property '%s' already processed, skipping", prop_name)
            return
        try:
            logger.debug("Processing computed property: '%s'", prop_name)
            prop_config = self.properties[prop_name]
            if not isinstance(prop_config, dict) or EQUATION_KEY not in prop_config:
                raise ValueError(f"Invalid COMPUTED_PROPERTY configuration for '{prop_name}'")
            dep_def = prop_config[DEPENDENCY_KEY]
            dep_array = DependencyResolver.resolve_dependency_definition(dep_def, material=material)
            expression = prop_config[EQUATION_KEY]
            logger.debug("Computing property '%s' with expression: %s", prop_name, expression)
            try:
                material_property = self._parse_and_process_expression(expression, material, dependency, prop_name)
            except (DependencyError, CircularDependencyError):
                raise
            except Exception as e:
                raise ValueError(f"Failed to process computed property '{prop_name}' \n -> {str(e)}") from e
            # Symbolic dependency: lambdify over the dependency array
            # material_property is already expressed in terms of `dependency` (the caller's symbol).
            f_pw = sp.lambdify(dependency, material_property, 'numpy')
            try:
                y_dense = f_pw(dep_array)
                if not np.all(np.isfinite(y_dense)):
                    invalid_count = np.sum(~np.isfinite(y_dense))
                    logger.warning("Property '%s' has %d non-finite values. "
                                   "Check expression: %s", prop_name, invalid_count, expression)
                if self.property_handler is not None:
                    self.property_handler.finalize_computed_property(
                        material, prop_name, dep_array, y_dense, dependency, prop_config
                    )
                else:
                    setattr(material, prop_name, material_property)
                    self.processed_properties.add(prop_name)
                    logger.warning("Property handler not available for '%s' - skipping visualization", prop_name)
                logger.debug("Successfully computed property '%s' over %d dependency points", prop_name, len(dep_array))
            except (DependencyError, CircularDependencyError):
                raise
            except Exception as e:
                raise ValueError(f"Error evaluating expression for '{prop_name}' \n -> {str(e)}") from e
        except (DependencyError, CircularDependencyError):
            raise
        except Exception as e:
            raise ValueError(f"Failed to process computed property '{prop_name}' \n -> {str(e)}") from e

    def _parse_and_process_expression(self, expression: str, material: Material,
                                      dependency: sp.Symbol,
                                      prop_name: str) -> sp.Expr:
        """Parses a mathematical expression string into a SymPy expression,
        substitutes all property dependencies and the placeholder symbol,
        and returns the final expression in terms of `dependency`.

        The YAML placeholder symbol `T` is substituted with `dependency`
        so the returned expression is expressed in the caller's symbol.
        """
        try:
            logger.debug("Parsing expression for '%s': %s", prop_name, expression)
            sympy_expr = sp.sympify(expression, evaluate=False)
            # Identify non-placeholder free symbols - these are property dependencies
            dependencies = [str(s) for s in sympy_expr.free_symbols if s != YAML_PLACEHOLDER]
            if dependencies:
                logger.debug("Property '%s' depends on: %s", prop_name, dependencies)
                missing_deps = [dep for dep in dependencies
                                if not hasattr(material, dep) and dep not in self.properties]
                if missing_deps:
                    available_props = sorted(self.properties.keys())
                    raise DependencyError(expression=expression, missing_deps=missing_deps,
                                          available_props=available_props)
                self._validate_circular_dependencies(prop_name, dependencies, set())
                # Ensure all dependencies are processed before this property
                for dep in dependencies:
                    if not hasattr(material, dep) or getattr(material, dep) is None:
                        if dep in self.properties:
                            logger.debug("Processing dependency '%s' for '%s'", dep, prop_name)
                            self.process_computed_property(material, dep, dependency)
                        else:
                            raise DependencyError(expression=expression, missing_deps=[dep],
                                                  available_props=sorted(self.properties.keys()))
            # Verify all dependencies are now available
            still_missing = [dep for dep in dependencies
                             if not hasattr(material, dep) or getattr(material, dep) is None]
            if still_missing:
                raise ValueError(f"Cannot compute expression. Missing dependencies: {still_missing}")
            # Build substitution dict: property symbols -> their material values
            substitutions = {}
            for dep in dependencies:
                dep_value = getattr(material, dep, None)
                if dep_value is None:
                    raise ValueError(f"Dependency '{dep}' was processed but is still not available on the material")
                dep_symbol = SymbolRegistry.get(dep)
                if dep_symbol is None:
                    raise ValueError(f"Symbol '{dep}' not found in symbol registry")
                substitutions[dep_symbol] = dep_value
            # Substitute the YAML placeholder with the caller's dependency symbol (or value)
            substitutions[YAML_PLACEHOLDER] = dependency
            result_expr = sympy_expr.subs(substitutions)
            # Evaluate any integrals (e.g. Integral(heat_capacity, T))
            if isinstance(result_expr, sp.Integral):
                logger.debug("Evaluating integral in expression for '%s'", prop_name)
                result_expr = result_expr.doit()
            logger.debug("Successfully processed expression for '%s'", prop_name)
            return result_expr
        except (CircularDependencyError, DependencyError):
            raise
        except Exception as e:
            raise ValueError(f"Failed to parse and process expression: {expression}") from e

    def _validate_circular_dependencies(self, prop_name: str, current_deps: List[str],
                                        visited: Set[str], path: List[str] = None) -> None:
        """Checks for circular dependencies in property definitions.

        Args:
            prop_name:    Current property being checked.
            current_deps: Dependencies of the current property.
            visited:      Set of already visited properties.
            path:         Current dependency path for error reporting.
        Raises:
            CircularDependencyError: If a circular dependency is detected.
        """
        if path is None:
            path = []
        if prop_name is not None:
            if prop_name in visited:
                cycle_path = path + [prop_name]
                # Origin point - log here before raising
                logger.error("Circular dependency detected: %s", ' -> '.join(cycle_path))
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
        """Extracts non-placeholder free symbols from equation data."""
        symbols = set()
        if isinstance(equation_data, list):
            for eq in equation_data:
                symbols.update(sp.sympify(eq).free_symbols)
        else:
            symbols.update(sp.sympify(equation_data).free_symbols)
        return [str(s) for s in symbols if s != YAML_PLACEHOLDER]
