import logging
from typing import List, Tuple, Union, Dict
import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


# --- Core Utility Functions ---

def handle_numeric_dependency(processor_instance, material: Material,
                              prop_name: str, piecewise_expr: sp.Expr,
                              dependency_symbols: Dict[str, sp.Symbol]) -> bool:
    """
    Handle numeric dependency evaluation with consistent error handling.

    Returns:
        bool: True if numeric evaluation was performed, False if symbolic
    """
    # Check if all dependencies are numeric
    all_numeric = all(not isinstance(symbol, sp.Symbol) for symbol in dependency_symbols.values())
    if not all_numeric:
        return False
    try:
        # For single dependency case, use the standard approach
        if len(dependency_symbols) == 1:
            dependency_value = list(dependency_symbols.values())[0]
            T_standard = sp.Symbol('T')
            value = float(piecewise_expr.subs(T_standard, dependency_value).evalf())
            setattr(material, prop_name, sp.Float(value))
            processor_instance.processed_properties.add(prop_name)
            logger.debug(f"Numeric evaluation completed for '{prop_name}': {value}")
            return True
        # For multi-dependency case, substitute all symbols
        else:
            # Get symbol mapping from material config for YAML symbols
            if hasattr(material, '_parser_config') and 'independent_variables' in material._parser_config:
                symbol_mapping = material._parser_config['independent_variables']
            else:
                # Fallback: assume temperature mapping
                symbol_mapping = {'temperature': 'T'}
            # Create substitutions for YAML symbols
            substitutions = {}
            for dep_name, dep_value in dependency_symbols.items():
                if dep_name in symbol_mapping:
                    yaml_symbol_name = symbol_mapping[dep_name]
                    yaml_symbol = sp.Symbol(yaml_symbol_name)
                    substitutions[yaml_symbol] = dep_value
            # Apply substitutions and evaluate
            value = float(piecewise_expr.subs(substitutions).evalf())
            setattr(material, prop_name, sp.Float(value))
            processor_instance.processed_properties.add(prop_name)
            logger.debug(f"Numeric evaluation completed for '{prop_name}': {value}")
            return True
    except Exception as e:
        dep_values_str = {k: float(v) if not isinstance(v, sp.Symbol) else str(v)
                          for k, v in dependency_symbols.items()}
        raise ValueError(f"Failed to evaluate {prop_name} at dependency values {dep_values_str}: {str(e)}")


def create_step_visualization_data(transition_temp: float, val_array: List[float],
                                   temp_range: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Create proper step function visualization data."""
    margin = (np.max(temp_range) - np.min(temp_range)) * ProcessingConstants.TEMPERATURE_PADDING_FACTOR
    epsilon = ProcessingConstants.TEMPERATURE_EPSILON
    x_data = np.array([
        np.min(temp_range) - margin,
        transition_temp - epsilon,
        transition_temp,
        transition_temp + epsilon,
        np.max(temp_range) + margin
    ])
    y_data = np.array([
        val_array[0],
        val_array[0],
        val_array[0],
        val_array[1],
        val_array[1]
    ])
    return x_data, y_data


def ensure_sympy_compatible(value):
    """Ensure value is compatible with SymPy operations."""
    if hasattr(value, 'item'):  # NumPy scalar
        return float(value.item())
    elif isinstance(value, (np.float64, np.int64, np.float32, np.int32, np.number)):
        return float(value)
    elif isinstance(value, (list, np.ndarray)):
        return [float(x) for x in value]
    else:
        return float(value)


def extract_dependency_from_config(config: Dict, dependency_name: str):
    """Extract dependency configuration for a specific dependency name."""
    dependencies = config.get('dependencies', [])
    if dependency_name not in dependencies:
        raise ValueError(f"Dependency '{dependency_name}' not found in configuration")
    ranges = config.get('ranges', {})
    if dependency_name not in ranges:
        raise ValueError(f"Range for dependency '{dependency_name}' not found in configuration")
    return ranges[dependency_name]


def validate_dependency_symbols(dependency_symbols: Dict[str, sp.Symbol],
                                required_dependencies: List[str]) -> None:
    """Validate that provided dependency symbols match required dependencies."""
    provided_deps = set(dependency_symbols.keys())
    required_deps = set(required_dependencies)
    missing_deps = required_deps - provided_deps
    if missing_deps:
        raise ValueError(f"Missing dependency symbols: {sorted(list(missing_deps))}. "
                         f"Required: {sorted(list(required_deps))}")
    extra_deps = provided_deps - required_deps
    if extra_deps:
        raise ValueError(f"Extra dependency symbols: {sorted(list(extra_deps))}. "
                         f"Expected: {sorted(list(required_deps))}")
    # Validate that all symbols are SymPy symbols or numeric values
    for dep_name, symbol in dependency_symbols.items():
        if not isinstance(symbol, (sp.Symbol, int, float)):
            raise ValueError(f"Invalid symbol type for dependency '{dep_name}': {type(symbol)}. "
                             "Must be sympy.Symbol, int, or float.")


def get_primary_dependency(config: Dict) -> str:
    """Get the primary (first) dependency from configuration."""
    dependencies = config.get('dependencies', [])
    if not dependencies:
        return 'temperature'  # Default fallback
    return dependencies[0]


def create_symbol_mapping(independent_variables: Dict[str, str],
                          dependency_symbols: Dict[str, sp.Symbol]) -> Dict[sp.Symbol, sp.Symbol]:
    """Create mapping from YAML symbols to actual dependency symbols."""
    mapping = {}
    for dep_name, yaml_symbol_name in independent_variables.items():
        if dep_name in dependency_symbols:
            yaml_symbol = sp.Symbol(yaml_symbol_name)
            actual_symbol = dependency_symbols[dep_name]
            mapping[yaml_symbol] = actual_symbol
    return mapping


def substitute_dependency_symbols(expression: sp.Expr,
                                  dependency_symbols: Dict[str, sp.Symbol],
                                  symbol_mapping: Dict[str, str]) -> sp.Expr:
    """Replace YAML symbols in expression with actual dependency symbols."""
    substitutions = {}
    for dep_name, actual_symbol in dependency_symbols.items():
        if dep_name in symbol_mapping:
            yaml_symbol_name = symbol_mapping[dep_name]
            yaml_symbol = sp.Symbol(yaml_symbol_name)
            substitutions[yaml_symbol] = actual_symbol
    return expression.subs(substitutions)

def validate_dependency_symbols1(dependency_symbols: Dict[str, sp.Symbol],
                                config_dependencies: Dict[str, str]) -> None:
    """
    Validate that provided dependency symbols match the configuration.

    Args:
        dependency_symbols: Dictionary mapping dependency names to SymPy symbols
        config_dependencies: Dictionary from YAML independent_variables section

    Raises:
        ValueError: If dependencies don't match
    """
    config_deps = set(config_dependencies.keys())
    provided_deps = set(dependency_symbols.keys())

    missing_deps = config_deps - provided_deps
    if missing_deps:
        logger.error("Missing dependency symbols: %s", missing_deps)
        raise ValueError(f"Missing dependency symbols: {sorted(list(missing_deps))}. "
                         f"Required: {sorted(list(config_deps))}")

    extra_deps = provided_deps - config_deps
    if extra_deps:
        logger.error("Extra dependency symbols: %s", extra_deps)
        raise ValueError(f"Extra dependency symbols: {sorted(list(extra_deps))}. "
                         f"Expected: {sorted(list(config_deps))}")

    logger.debug("Dependency symbols validated successfully")


def create_dependency_mapping(independent_variables: Dict[str, str],
                              dependency_symbols: Dict[str, sp.Symbol]) -> Dict[str, sp.Symbol]:
    """
    Create mapping from YAML symbols to actual symbols.

    Args:
        independent_variables: Mapping from dependency names to YAML symbol names
        dependency_symbols: Mapping from dependency names to actual symbols

    Returns:
        Dictionary mapping YAML symbol names to actual symbols
    """
    mapping = {}
    for dep_name, yaml_symbol in independent_variables.items():
        if dep_name in dependency_symbols:
            actual_symbol = dependency_symbols[dep_name]
            mapping[yaml_symbol] = actual_symbol

    return mapping


def substitute_dependency_symbols1(expression: sp.Expr,
                                  symbol_mapping: Dict[str, sp.Symbol]) -> sp.Expr:
    """
    Replace YAML placeholder symbols with actual dependency symbols.

    Args:
        expression: SymPy expression with YAML symbols
        symbol_mapping: Mapping from YAML symbol names to actual symbols

    Returns:
        Expression with substituted symbols
    """
    substitutions = {}
    for yaml_symbol_name, actual_symbol in symbol_mapping.items():
        yaml_symbol = sp.Symbol(yaml_symbol_name)
        substitutions[yaml_symbol] = actual_symbol

    return expression.subs(substitutions)


def extract_primary_dependency(config: Dict, default: str = 'temperature') -> str:
    """
    Extract the primary dependency from property configuration.

    Args:
        config: Property configuration dictionary
        default: Default dependency name if none found

    Returns:
        Primary dependency name
    """
    dependencies = config.get('dependencies', [])
    if dependencies:
        return dependencies[0]
    return default


def format_dependency_error(prop_name: str, missing_deps: List[str],
                            available_deps: List[str]) -> str:
    """
    Format a user-friendly dependency error message.

    Args:
        prop_name: Name of the property with missing dependencies
        missing_deps: List of missing dependency names
        available_deps: List of available dependency names

    Returns:
        Formatted error message
    """
    message = f"Property '{prop_name}' has missing dependencies: {', '.join(missing_deps)}"
    if available_deps:
        message += f"\nAvailable dependencies: {', '.join(available_deps)}"
    return message
