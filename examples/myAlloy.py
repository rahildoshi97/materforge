# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import sympy as sp

from materforge import (
    create_material,
    evaluate_material_properties,
    get_material_property_names,
    validate_yaml_file,
    get_material_info,
)

MATERIAL_CONFIGS = [
    #("src/materforge/data/materials/Al.yaml", "T_Al"),
    #("src/materforge/data/materials/1.4301.yaml", "T_SS"),
    ("examples/myAlloy.yaml", "u_C"),
    #("src/materforge/data/materials/Al2O3.yaml", "T_Al2O3"),
]

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
        if not (hasattr(getattr(mat, n), 'free_symbols') and getattr(mat, n).free_symbols)
        and isinstance(getattr(mat, n), (int, float, sp.Float, sp.Integer))
    }

def demonstrate_material_properties() -> None:
    repo_root = Path(__file__).parent.parent

    targets = []
    for rel_path, sym_name in MATERIAL_CONFIGS:
        yaml_path = repo_root / rel_path
        if not yaml_path.exists():
            raise FileNotFoundError(f"YAML file not found: {yaml_path} (configured in {Path(__file__)})")
        targets.append((yaml_path, sp.Symbol(sym_name)))

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

    for yaml_path, _ in targets:
        is_valid = validate_yaml_file(yaml_path)
        print(f"{yaml_path.stem} YAML validation: {'PASSED' if is_valid else 'FAILED'}")

    # ===================================================================
    # 2. MATERIAL INFO EXTRACTION
    # ===================================================================
    print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
    print(f"{'-' * 50}")

    for yaml_path, _ in targets:
        info = get_material_info(yaml_path)
        print(f"\n{info['name']} Information:")
        print(f"  Name:             {info['name']}")
        print(f"  Total Properties: {info['total_properties']}")
        print(f"  Properties:       {info['properties']}")
        if 'property_types' in info:
            print("  Property Types:")
            for ptype, count in info['property_types'].items():
                print(f"    {ptype:<30}: {count}")

    # ===================================================================
    # 3. MATERIAL CREATION
    # ===================================================================
    print(f"\n{'3. MATERIAL CREATION':<50}")
    print(f"{'-' * 50}")

    for yaml_path, T in targets:
        mat = create_material(yaml_path=yaml_path, dependency=T, enable_plotting=True)
        materials.append(mat)
        material_symbols[mat.name] = T
        print(f"Successfully created: {mat.name} (symbol: {T})")
        print(f"Representation: {repr(mat)}")
        scalars = _scalar_properties(mat)
        if scalars:
            first_name, first_val = next(iter(sorted(scalars.items())))
            print(f"{first_name}: {first_val}  (type: {type(first_val).__name__})")

    # ===================================================================
    # 4. PROPERTY TYPES OVERVIEW
    # ===================================================================
    print(f"\n{'4. PROPERTY TYPES OVERVIEW':<50}")
    print(f"{'-' * 50}")
    print("MaterForge supports 6 property definition types (any property name is valid):")
    for i, (ptype, desc) in enumerate([
        ("CONSTANT_VALUE",     "Single numeric value, e.g. density: 7000."),
        ("STEP_FUNCTION",      "Two values split at a scalar-property reference"),
        ("TABULAR_DATA",       "Paired dependency / value lists"),
        ("FILE_IMPORT",        "Load from .csv / .xlsx / .txt"),
        ("PIECEWISE_EQUATION", "Symbolic equations over dependency breakpoints"),
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
        print(f"MATERIAL: {mat.name}  (dependency symbol: {T_mat})")
        print(f"{'=' * 80}")

        scalars = _scalar_properties(mat)
        if scalars:
            print("Scalar constants:")
            for prop_name, val in sorted(scalars.items()):
                print(f"  {prop_name:<35}: {val}")

        available_props = get_material_property_names(mat)

        # -------------------------------------------------------------------
        # 5.1 Available Properties
        # -------------------------------------------------------------------
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")
        print(f"{mat.name} has {len(available_props)} processed properties:")
        for prop in sorted(available_props):
            print(f"  - {prop}")

        # -------------------------------------------------------------------
        # 5.2 Full Symbolic Expressions
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
        print(f"\n{'PROPERTY VALUES AT ' + str(test_temp) + ' K:':<50}")
        print(f"{'-' * 50}")

        print("Method 1: Manual substitution")
        for prop_name in sorted(available_props):
            prop_value = getattr(mat, prop_name)
            if isinstance(prop_value, sp.Expr) and prop_value.free_symbols:
                numerical = prop_value.subs(T_mat, test_temp).evalf()
                print(f"  {prop_name:<30}: {numerical}  (symbolic)")
            else:
                print(f"  {prop_name:<30}: {prop_value}  (constant)")

        # -------------------------------------------------------------------
        # 5.4 API Evaluation Methods
        # -------------------------------------------------------------------
        print(f"\n{'API METHODS:':<50}")
        print(f"{'-' * 50}")

        print("Method 2: material.evaluate()")
        evaluated_mat = mat.evaluate(T_mat, test_temp)
        print(f"All properties at {test_temp} K:")
        for prop in sorted(evaluated_mat.property_names()):
            value = getattr(evaluated_mat, prop)
            print(f"  {prop:<30}: {type(value).__name__:<15}  {float(value):.6e}")

        print("\nMethod 3: evaluate_material_properties()")
        all_values = evaluate_material_properties(mat, T_mat, test_temp)
        print(f"  Results type: {type(all_values).__name__}")
        for prop in sorted(all_values.property_names()):
            value = getattr(all_values, prop)
            print(f"  {prop:<30}: {float(value):.6e}")


if __name__ == "__main__":
    setup_logging()
    demonstrate_material_properties()
