# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import sympy as sp

from materforge.parsing.api import (
    create_material,
    evaluate_material_properties,
    get_material_property_names,
    validate_yaml_file,
    get_material_info,
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
    """Return all constant (non-symbolic) properties on the material."""
    return {
        n: getattr(mat, n)
        for n in mat.property_names()
        if not (hasattr(getattr(mat, n), 'free_symbols')
                and getattr(mat, n).free_symbols)
        and isinstance(getattr(mat, n), (int, float, sp.Float, sp.Integer))
    }

def demonstrate_material_properties() -> None:
    setup_logging()

    T = sp.Symbol('T')

    current_file = Path(__file__)
    yaml_path = current_file.parent / "myAlloy.yaml"
    # yaml_path = (current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml")
    # yaml_path = (current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml")

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

    if yaml_path.exists():
        try:
            is_valid = validate_yaml_file(yaml_path)
            print(f"myAlloy YAML validation: {'PASSED' if is_valid else 'FAILED'}")
        except Exception as e:
            raise ValueError(f"myAlloy YAML validation: FAILED - {e}")
    else:
        raise ValueError(f"myAlloy YAML file not found: {yaml_path}")

    # ===================================================================
    # 2. MATERIAL INFO EXTRACTION
    # ===================================================================
    print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
    print(f"{'-' * 50}")

    if yaml_path.exists():
        try:
            info = get_material_info(yaml_path)
            print("\nmyAlloy Information:")
            print(f"  Name:             {info['name']}")
            print(f"  Total Properties: {info['total_properties']}")
            print(f"  Properties:       {info['properties']}")
            if 'property_types' in info:
                print("  Property Types:")
                for ptype, count in info['property_types'].items():
                    print(f"    {ptype:<30}: {count}")
        except Exception as e:
            raise ValueError(f"Failed to get myAlloy info: {e}")

    # ===================================================================
    # 3. MATERIAL CREATION
    # ===================================================================
    print(f"\n{'3. MATERIAL CREATION':<50}")
    print(f"{'-' * 50}")

    if yaml_path.exists():
        try:
            myAlloy = create_material(yaml_path=yaml_path, dependency=T, enable_plotting=True)
            materials.append(myAlloy)
            material_symbols[myAlloy.name] = T
            print(f"Successfully created: {myAlloy.name} (symbol: {T})")
            print(f"Representation: {repr(myAlloy)}")
            print(f"density1: {myAlloy.density1}  (type: {type(myAlloy.density1).__name__})")
        except Exception as e:
            raise ValueError(f"Failed to create myAlloy: {e}")

    # ===================================================================
    # 4. PROPERTY TYPES OVERVIEW
    # ===================================================================
    print(f"\n{'4. PROPERTY TYPES OVERVIEW':<50}")
    print(f"{'-' * 50}")
    print("MaterForge supports 6 property definition types "
          "(any property name is valid):")
    for i, (ptype, desc) in enumerate([
        ("CONSTANT_VALUE",     "Single numeric value, e.g. density: 7000."),
        ("STEP_FUNCTION",      "Two values split at a scalar-property reference"),
        ("TABULAR_DATA",       "Paired dependency / value lists"),
        ("FILE_IMPORT",        "Load from .csv / .xlsx / .txt"),
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
        # Scalar constants (temperature references, fixed values) - dynamic lookup
        scalars = _scalar_properties(mat)
        if scalars:
            print("Scalar constants:")
            for prop_name, val in sorted(scalars.items()):
                print(f"  {prop_name:<35}: {val}")
        available_props = get_material_property_names(mat)

        # -------------------------------------------------------------------
        # 5.1 Available Properties (names only)
        # -------------------------------------------------------------------
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")
        print(f"{mat.name} has {len(available_props)} processed properties:")
        for prop in sorted(available_props):
            print(f"  - {prop}")

        # -------------------------------------------------------------------
        # 5.2 Full Symbolic Expressions
        # Shows the complete SymPy expression for every property - constants
        # print their value, symbolic properties print the full Piecewise/Expr.
        # -------------------------------------------------------------------
        print(f"\n{'PROPERTY EXPRESSIONS:':<50}")
        print(f"{'-' * 50}")
        for prop_name in sorted(available_props):
            prop_value = getattr(mat, prop_name)
            if isinstance(prop_value, sp.Expr) and prop_value.free_symbols:
                print(f"  {prop_name:<30}: {prop_value}")
            else:
                print(f"  {prop_name:<30}: {prop_value}  (constant)")

        # -------------------------------------------------------------------
        # 5.3 Numerical Evaluation at a Single Temperature
        # -------------------------------------------------------------------
        test_temp = 500.15
        print(f"\n{'PROPERTY VALUES AT ' + str(test_temp) + 'K:':<50}")
        print(f"{'-' * 50}")

        print("Method 1: Manual substitution")
        for prop_name in sorted(available_props):
            try:
                prop_value = getattr(mat, prop_name)
                if isinstance(prop_value, sp.Expr) and prop_value.free_symbols:
                    numerical = prop_value.subs(T_mat, test_temp).evalf()
                    print(f"  {prop_name:<30}: {numerical}  (symbolic)")
                else:
                    print(f"  {prop_name:<30}: {prop_value}  (constant)")
            except Exception as e:
                raise ValueError(f"  {prop_name:<30}: Error - {e}")

        # -------------------------------------------------------------------
        # 5.4 API Evaluation Methods
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
            raise ValueError(f"Error in Method 2: {e}")

        print("\nMethod 3: evaluate_material_properties()")
        try:
            all_values_func = evaluate_material_properties(mat, test_temp)
            print(f"  Results match Method 2: {all_values == all_values_func}")
        except Exception as e:
            raise ValueError(f"Error in Method 3: {e}")

        print("\nMethod 4: Specific Properties Only")
        try:
            specific_props = sorted(available_props)[:2]
            if specific_props:
                specific_values = mat.evaluate_properties_at_temperature(
                    test_temp, properties=specific_props)
                print(f"  Requested: {specific_props}")
                for prop, value in sorted(specific_values.items()):
                    print(f"  {prop:<30}: {value:.6e}")
            else:
                print("  No properties available")
        except Exception as e:
            raise ValueError(f"Error in Method 4: {e}")

        print("\nMethod 5: Temperature-Dependent Properties Only")
        try:
            temp_dependent = mat.evaluate_properties_at_temperature(
                test_temp, include_constants=False)
            if temp_dependent:
                print(f"  Temperature-dependent properties at {test_temp} K:")
                for prop, value in sorted(temp_dependent.items()):
                    print(f"  {prop:<30}: {value:.6e}")
            else:
                print("  No temperature-dependent properties found")
        except Exception as e:
            raise ValueError(f"Error in Method 5: {e}")


if __name__ == "__main__":
    demonstrate_material_properties()
