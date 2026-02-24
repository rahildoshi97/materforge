# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import pandas as pd
import sympy as sp

from materforge.algorithms.piecewise_inverter import PiecewiseInverter
from materforge.parsing.api import (
    create_material,
    evaluate_material_properties,
    get_material_info,
    get_material_property_names,
    validate_yaml_file,
)


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s %(levelname)s %(name)s -> %(message)s",
    )
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('fontTools').setLevel(logging.WARNING)


def _scalar_properties(mat) -> dict:
    """Return all constant (non-symbolic) properties on a material."""
    return {
        n: getattr(mat, n)
        for n in mat.property_names()
        if not (hasattr(getattr(mat, n), 'free_symbols')
                and getattr(mat, n).free_symbols)
        and isinstance(getattr(mat, n), (int, float, sp.Float, sp.Integer))
    }

def demonstrate_material_properties() -> None:
    """Comprehensive demonstration of material property evaluation."""
    setup_logging()

    T1 = sp.Symbol('T_Al')
    T2 = sp.Symbol('T_SS')

    current_file = Path(__file__)
    yaml_path_Al = (current_file.parent.parent / "src" / "materforge" / "data"
                    / "materials" / "pure_metals" / "Al" / "Al.yaml")
    yaml_path_SS304L = (current_file.parent.parent / "src" / "materforge" / "data"
                        / "materials" / "alloys" / "1.4301" / "1.4301.yaml")

    materials = []
    material_symbols = {}

    print(f"\n{'=' * 80}")
    print("MaterForge MATERIAL PROPERTY DEMONSTRATION")
    print(f"{'=' * 80}")

    # ===================================================================
    # 1. YAML FILE VALIDATION
    # ===================================================================
    print(f"\n{'1. YAML FILE VALIDATION':<50}")
    print(f"{'-' * 50}")

    for path, name in [(yaml_path_Al, "Aluminum"), (yaml_path_SS304L, "Steel 1.4301")]:
        if path.exists():
            try:
                is_valid = validate_yaml_file(path)
                print(f"{name} YAML validation: {'PASSED' if is_valid else 'FAILED'}")
            except Exception as e:
                print(f"{name} YAML validation: FAILED - {e}")
        else:
            print(f"{name} YAML file not found: {path}")

    # ===================================================================
    # 2. MATERIAL INFO EXTRACTION
    # ===================================================================
    print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
    print(f"{'-' * 50}")

    for path, name in [(yaml_path_Al, "Aluminum"), (yaml_path_SS304L, "SS304L")]:
        if not path.exists():
            print(f"{name} YAML file not found: {path}")
            continue
        try:
            info = get_material_info(path)
            print(f"\n{name} Information:")
            print(f"  Name:             {info['name']}")
            print(f"  Total Properties: {info['total_properties']}")
            print(f"  Properties:       {info['properties']}")
            if 'property_types' in info:
                print("  Property Types:")
                for pt, count in info['property_types'].items():
                    print(f"    {pt:<30}: {count}")
        except Exception as e:
            print(f"Failed to get {name} info: {e}")

    # ===================================================================
    # 3. MATERIAL CREATION
    # ===================================================================
    print(f"\n{'3. MATERIAL CREATION':<50}")
    print(f"{'-' * 50}")

    for path, symbol, label in [
        (yaml_path_Al,    T1, "Aluminum"),
        (yaml_path_SS304L, T2, "Steel 1.4301"),
    ]:
        if not path.exists():
            print(f"{label}: YAML not found")
            continue
        try:
            mat = create_material(yaml_path=path, dependency=symbol, enable_plotting=True)
            materials.append(mat)
            material_symbols[mat.name] = symbol
            print(f"Successfully created: {mat.name} (symbol: {symbol})")
        except Exception as e:
            print(f"Failed to create {label}: {e}")

    # ===================================================================
    # 4. PROPERTY TYPES OVERVIEW
    # ===================================================================
    print(f"\n{'4. PROPERTY TYPES OVERVIEW':<50}")
    print(f"{'-' * 50}")
    print("MaterForge supports 6 property definition types (any property name is valid):")
    for i, (ptype, desc) in enumerate([
        ("CONSTANT_VALUE",     "Single numeric value, e.g. density: 7000.0"),
        ("STEP_FUNCTION",      "Two values split at a scalar-property reference"),
        ("TABULAR_DATA",       "Paired dependency / value lists"),
        ("FILE_IMPORT",        "Load data from .csv / .xlsx / .txt"),
        ("PIECEWISE_EQUATION", "Symbolic equations over temperature breakpoints"),
        ("COMPUTED_PROPERTY",  "Derived from other properties via expression"),
    ], 1):
        print(f"  {i}. {ptype:<30} - {desc}")

    # ===================================================================
    # 5. MATERIAL PROPERTY ANALYSIS
    # ===================================================================
    print(f"\n{'5. MATERIAL PROPERTY ANALYSIS':<50}")
    print(f"{'-' * 50}")

    for mat in materials:
        T_mat = material_symbols[mat.name]

        print(f"\n{'=' * 80}")
        print(f"MATERIAL: {mat.name}  (temperature symbol: {T_mat})")
        print(f"{'=' * 80}")
        print(f"Name: {mat.name}")

        # Show scalar constants (temperature references, fixed densities, etc.)
        scalars = _scalar_properties(mat)
        if scalars:
            print("Scalar constants:")
            for prop_name, val in sorted(scalars.items()):
                print(f"  {prop_name:<35}: {val}")

        # -------------------------------------------------------------------
        # 5.1 Available Properties
        # -------------------------------------------------------------------
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")

        available_props = get_material_property_names(mat)
        print(f"Material has {len(available_props)} processed properties:")

        for prop_name in sorted(available_props):
            prop_value = getattr(mat, prop_name)
            print(f"  {prop_name:<30}: {prop_value}")
            print(f"  {'':30}  Type: {type(prop_value).__name__}")
            allowed_types = (sp.Piecewise, sp.Float, sp.Expr, float, int, type(None))
            marker = ("Valid type" if isinstance(prop_value, allowed_types) else "! WARNING: Unexpected type")
            print(f"  {'':30}  {marker}")

        # -------------------------------------------------------------------
        # 5.2 Manual Property Evaluation
        # -------------------------------------------------------------------
        test_temp = 500.15
        print(f"\n{'PROPERTY VALUES AT ' + str(test_temp) + 'K:':<50}")
        print(f"{'-' * 50}")

        print("Method 1: Manual Property Evaluation")
        for prop_name in sorted(available_props):
            try:
                prop_value = getattr(mat, prop_name)
                if isinstance(prop_value, sp.Expr) and prop_value.free_symbols:
                    numerical_value = prop_value.subs(T_mat, test_temp).evalf()
                    print(f"  {prop_name:<30}: {numerical_value} (symbolic)")
                else:
                    print(f"  {prop_name:<30}: {prop_value} (constant)")
            except Exception as e:
                print(f"  {prop_name:<30}: Error - {e}")

        # -------------------------------------------------------------------
        # 5.3 API Evaluation Methods
        # -------------------------------------------------------------------
        print(f"\n{'API METHODS:':<50}")
        print(f"{'-' * 50}")

        print("Method 2: material.evaluate_properties_at_temperature()")
        try:
            all_values = mat.evaluate_properties_at_temperature(test_temp)
            print(f"All properties at {test_temp} K:")
            for prop, value in sorted(all_values.items()):
                print(f"  {prop:<30}: {value:.6e}")
        except Exception as e:
            print(f"Error in Method 2: {e}")

        print("\nMethod 3: evaluate_material_properties()")
        try:
            all_values_func = evaluate_material_properties(mat, test_temp)
            print(f"  Results match Method 2: {all_values == all_values_func}")
        except Exception as e:
            print(f"Error in Method 3: {e}")

        print("\nMethod 4: Specific Properties Only")
        try:
            specific_props = sorted(available_props)[:3]
            if specific_props:
                specific_values = mat.evaluate_properties_at_temperature(
                    test_temp, properties=specific_props)
                print(f"  Requested: {specific_props}")
                for prop, value in sorted(specific_values.items()):
                    print(f"  {prop:<30}: {value:.6e}")
            else:
                print("  No properties available")
        except Exception as e:
            print(f"Error in Method 4: {e}")

        print("\nMethod 5: Temperature-Dependent Properties Only")
        try:
            temp_dependent = mat.evaluate_properties_at_temperature(
                test_temp, include_constants=False)
            if temp_dependent:
                for prop, value in sorted(temp_dependent.items()):
                    print(f"  {prop:<30}: {value:.6e}")
            else:
                print("  No temperature-dependent properties found")
        except Exception as e:
            print(f"Error in Method 5: {e}")

        # -------------------------------------------------------------------
        # 5.4 Batch Temperature Evaluation
        # -------------------------------------------------------------------
        print(f"\n{'BATCH TEMPERATURE EVALUATION:':<50}")
        print(f"{'-' * 50}")

        temperatures = [300, 400, 500, 600, 700]
        batch_results = {}
        for temp in temperatures:
            try:
                batch_results[temp] = mat.evaluate_properties_at_temperature(temp)
            except Exception as e:
                print(f"  Error at {temp}K: {e}")
                batch_results[temp] = {}

        if batch_results:
            try:
                df = pd.DataFrame(batch_results).T
                display_cols = [c for c in sorted(available_props) if c in df.columns][:3]
                if display_cols:
                    header = f"{'Temperature':<12}" + "".join(
                        f"{c:<18}" for c in display_cols)
                    print(header)
                    print("-" * (12 + 18 * len(display_cols)))
                    for temp, row in df.iterrows():
                        line = f"{temp:<12.0f}"
                        for col in display_cols:
                            val = row.get(col)
                            line += f"{val:.4e}  " if val is not None else f"{'N/A':<18}"
                        print(line)
                else:
                    print("  No numeric properties available for batch display")
            except Exception as e:
                print(f"  DataFrame creation failed: {e}")

    # ===================================================================
    # 6. INVERSE FUNCTION TESTING
    # ===================================================================
    print(f"\n{'6. INVERSE FUNCTION TESTING':<50}")
    print(f"{'-' * 50}")
    test_inverse_functions(materials, material_symbols)


def test_inverse_functions(materials: list, material_symbols: dict) -> None:
    """Tests inverse function creation and accuracy for materials with energy_density."""
    print(f"\n{'=' * 80}")
    print("INVERSE FUNCTION TESTING")
    print(f"{'=' * 80}")

    for mat in materials:
        T_mat = material_symbols[mat.name]
        print(f"\n--- Testing Inverse Functions for {mat.name} (symbol: {T_mat}) ---")

        if 'energy_density' not in mat.property_names():
            print(f"  {mat.name} has no energy_density property - skipping")
            continue

        try:
            print("  Method 2: Direct approach...")
            try:
                energy_density = getattr(mat, 'energy_density')
                energy_symbols = energy_density.free_symbols
                if len(energy_symbols) == 1:
                    temp_symbol = next(iter(energy_symbols))
                    E_symbol = sp.Symbol('E')
                    inverse_func2 = PiecewiseInverter.create_inverse(
                        energy_density, temp_symbol, E_symbol)
                    print("  Method 2 succeeded!")
                    test_round_trip_accuracy(
                        mat, inverse_func2, temp_symbol, E_symbol, method="Method 2")
                else:
                    print(f"  Unexpected symbols in energy_density: {energy_symbols}")
            except Exception as e:
                print(f"  Method 2 failed: {e}")
        except Exception as e:
            print(f"  Inverse testing failed for {mat.name}: {e}")


def test_round_trip_accuracy(material, inverse_func, temp_symbol,
                              energy_symbol, method: str = "Unknown") -> None:
    """Tests round-trip accuracy: T -> E -> T."""
    print(f"  Round-trip accuracy test ({method}):")

    test_temperatures = [300, 350, 400, 450, 500,
                         600, 700, 800, 900, 1000,
                         1200, 1500, 1800, 2000, 2500, 3000]
    max_error = 0.0
    successful = 0
    energy_density = getattr(material, 'energy_density')

    for temp in test_temperatures:
        try:
            energy_val = float(energy_density.subs(temp_symbol, temp))
            recovered = float(inverse_func.subs(energy_symbol, energy_val))
            error = abs(temp - recovered)
            max_error = max(max_error, error)
            successful += 1
            print(f"    T={temp:4.0f}K -> E={energy_val:.2e} -> "
                  f"T={recovered:6.1f}K  (error: {error:.2e})")
        except Exception as e:
            print(f"    T={temp:4.0f}K -> Error: {e}")

    print(f"  Completed: {successful}/{len(test_temperatures)}")
    print(f"  Max error: {max_error:.2e} K")
    if max_error < 1e-8:
        print("  Excellent accuracy")
    elif max_error < 1e-4:
        print("  Good accuracy")
    else:
        print("  ! Consider reviewing inverse function accuracy")


def demonstrate_advanced_usage() -> None:
    """Demonstrates advanced usage patterns."""
    print(f"\n{'=' * 80}")
    print("ADVANCED USAGE PATTERNS")
    print(f"{'=' * 80}")

    print("\n7. TEMPERATURE SWEEP ANALYSIS")
    print("-" * 50)

    current_file = Path(__file__)
    yaml_path_Al = (current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml")
    yaml_path_SS304L = (current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml")

    T1 = sp.Symbol('T_Al')
    T2 = sp.Symbol('T_SS')
    materials_advanced = []

    for path, symbol, label in [
        (yaml_path_Al,    T1, "Aluminum"),
        (yaml_path_SS304L, T2, "Steel 1.4301"),
    ]:
        if not path.exists():
            print(f"{label}: YAML not found")
            continue
        try:
            mat = create_material(path, dependency=symbol, enable_plotting=False)
            materials_advanced.append((mat, symbol, label))
            print(f"Created {mat.name} with symbol {symbol}")
        except Exception as e:
            print(f"Failed to create {label}: {e}")

    temp_range = range(300, 3001, 100)

    for mat, T_symbol, mat_name in materials_advanced:
        print(f"\n--- Temperature Sweep for {mat_name} ---")

        # Pick up to 3 temperature-dependent (symbolic) properties
        sweep_props = [
            p for p in sorted(mat.property_names())
            if (hasattr(getattr(mat, p, None), 'free_symbols')
                and getattr(mat, p).free_symbols)
        ][:3]

        if not sweep_props:
            print(f"  No temperature-dependent properties found for {mat_name}")
            continue

        sweep_data = []
        for temp in temp_range:
            try:
                props = mat.evaluate_properties_at_temperature(
                    temp, properties=sweep_props, include_constants=False)
                props['temperature'] = temp
                sweep_data.append(props)
            except Exception as e:
                print(f"  Error at {temp}K: {e}")

        if sweep_data:
            df = pd.DataFrame(sweep_data).set_index('temperature')
            print(f"Sweep for {mat_name} (properties: {sweep_props}):")
            print(df.head(5).to_string())
            for col in sweep_props:
                if col in df.columns and len(df[col].dropna()) >= 2:
                    first, last = df[col].iloc[0], df[col].iloc[-1]
                    if first != 0:
                        change = ((last - first) / first) * 100
                        print(f"  {col} change "
                              f"{list(temp_range)[0]}K -> {list(temp_range)[-1]}K: "
                              f"{change:.2f}%")


if __name__ == "__main__":
    demonstrate_material_properties()
    demonstrate_advanced_usage()
