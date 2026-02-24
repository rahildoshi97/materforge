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


def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.WARNING,  # DEBUG/INFO/WARNING/ERROR/CRITICAL
        format="%(asctime)s %(levelname)s %(name)s -> %(message)s",
    )
    # Silence noisy libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('fontTools').setLevel(logging.WARNING)


def demonstrate_material_properties():
    """Comprehensive demonstration of material property evaluation."""
    setup_logging()

    T = sp.Symbol('T')

    current_file = Path(__file__)
    yaml_path = current_file.parent / "myAlloy.yaml"

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

    for path, name in [(yaml_path, "myAlloy")]:
        if path.exists():
            try:
                is_valid = validate_yaml_file(path)
                print(f"✓ {name} YAML validation: {'PASSED' if is_valid else 'FAILED'}")
            except Exception as e:
                raise ValueError(f"✗ {name} YAML validation: FAILED - {e}")
        else:
            raise ValueError(f"✗ {name} YAML file not found: {path}")

    # ===================================================================
    # 2. MATERIAL INFO EXTRACTION
    # ===================================================================
    print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
    print(f"{'-' * 50}")

    for path, name in [(yaml_path, "myAlloy")]:
        if path.exists():
            try:
                info = get_material_info(path)
                print(f"\n{name} Information:")
                print(f"  Name:             {info['name']}")
                print(f"  Type:             {info['material_type']}")
                print(f"  Total Properties: {info['total_properties']}")
                print(f"  Composition:      {info['composition']}")

                if info['material_type'] == 'pure_metal':
                    print(f"  Melting Temperature: {info.get('melting_temperature', 'N/A')} K")
                    print(f"  Boiling Temperature: {info.get('boiling_temperature', 'N/A')} K")
                else:
                    print(f"  Solidus Temperature:          {info.get('solidus_temperature', 'N/A')} K")
                    print(f"  Liquidus Temperature:         {info.get('liquidus_temperature', 'N/A')} K")
                    print(f"  Initial Boiling Temperature:  {info.get('initial_boiling_temperature', 'N/A')} K")
                    print(f"  Final Boiling Temperature:    {info.get('final_boiling_temperature', 'N/A')} K")

                if 'property_types' in info:
                    print(f"  Property Types:")
                    for ptype, count in info['property_types'].items():
                        print(f"    {ptype:<30}: {count}")
            except Exception as e:
                raise ValueError(f"✗ Failed to get {name} info: {e}")

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
            print(f"✓ Successfully created: {myAlloy.name} (symbol: {T})")
        except Exception as e:
            raise ValueError(f"Failed to create myAlloy: {e}")

    # ===================================================================
    # 4. PROPERTY TYPES OVERVIEW
    # ===================================================================
    # get_supported_properties() has been removed - property names are no longer
    # restricted to a fixed list. The YAML parser accepts any property name.
    # Use get_material_info() to inspect what a specific YAML defines,
    # or material.property_names() on a created instance.
    print(f"\n{'4. PROPERTY TYPES OVERVIEW':<50}")
    print(f"{'-' * 50}")
    print("MaterForge supports 6 property types (any name is valid):")
    property_types = [
        ("CONSTANT_VALUE",     "Single numeric value: `density: 7000.`"),
        ("STEP_FUNCTION",      "Two values split at a temperature reference"),
        ("TABULAR_DATA",       "Key-value pairs (dependency list + value list)"),
        ("FILE_IMPORT",        "Load from .csv / .xlsx / .txt file"),
        ("PIECEWISE_EQUATION", "Symbolic equations over temperature breakpoints"),
        ("COMPUTED_PROPERTY",  "Derived from other properties via expression"),
    ]
    for i, (ptype, description) in enumerate(property_types):
        print(f"  {i:2d}. {ptype:<30} — {description}")

    # ===================================================================
    # 5. MATERIAL PROPERTY ANALYSIS
    # ===================================================================
    print(f"\n{'5. MATERIAL PROPERTY ANALYSIS':<50}")
    print(f"{'-' * 50}")

    for mat in materials:
        T_mat = material_symbols[mat.name]

        print(f"\n{'=' * 80}")
        print(f"MATERIAL: {mat.name} (Temperature symbol: {T_mat})")
        print(f"{'=' * 80}")

        print(f"Name:     {mat.name}")
        print(f"Type:     {mat.material_type}")
        print(f"Elements: {[elem.name for elem in mat.elements]}")
        for elem, frac in zip(mat.elements, mat.composition):
            print(f"  {elem.name}: {frac}")

        # Temperature fields — guard with is not None, not hasattr
        # (all temperature fields always exist on the dataclass, even when None)
        if mat.solidus_temperature is not None:
            print(f"Solidus Temperature:         {mat.solidus_temperature} K")
        if mat.liquidus_temperature is not None:
            print(f"Liquidus Temperature:        {mat.liquidus_temperature} K")
        if mat.melting_temperature is not None:
            print(f"Melting Temperature:         {mat.melting_temperature} K")
        if mat.boiling_temperature is not None:
            print(f"Boiling Temperature:         {mat.boiling_temperature} K")
        if mat.initial_boiling_temperature is not None:
            print(f"Initial Boiling Temperature: {mat.initial_boiling_temperature} K")
        if mat.final_boiling_temperature is not None:
            print(f"Final Boiling Temperature:   {mat.final_boiling_temperature} K")

        # ===================================================================
        # 5.1 Available Properties
        # ===================================================================
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")

        available_props = get_material_property_names(mat)
        print(f"{mat.name} has {len(available_props)} processed properties:")
        for prop in sorted(available_props):
            print(f"  - {prop}")

        # ===================================================================
        # 5.2 Manual Property Evaluation
        # ===================================================================
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
                raise ValueError(f"  {prop_name:<30}: Error - {str(e)}")

        # ===================================================================
        # 5.3 API Evaluation Methods
        # ===================================================================
        print(f"\n{'NEW API METHODS:':<50}")
        print(f"{'-' * 50}")

        # Method 2: OO approach
        print("Method 2: material.evaluate_properties_at_temperature()")
        try:
            all_values = mat.evaluate_properties_at_temperature(test_temp)
            print(f"All properties at {test_temp} K:")
            for prop, value in sorted(all_values.items()):
                print(f"  {prop:<30}: {value:.6e}")
        except Exception as e:
            raise ValueError(f"Error in Method 2: {e}")

        # Method 3: Functional approach
        print("\nMethod 3: evaluate_material_properties()")
        try:
            all_values_func = evaluate_material_properties(mat, test_temp)
            print(f"  Results match Method 2: {all_values == all_values_func}")
        except Exception as e:
            raise ValueError(f"Error in Method 3: {e}")

        # Method 4: Specific properties — uses actual names from this material
        print("\nMethod 4: Specific Properties Only")
        try:
            # Request the first 2 available properties as a concrete example
            specific_props = sorted(available_props)[:2] if len(available_props) >= 2 else available_props
            if specific_props:
                specific_values = mat.evaluate_properties_at_temperature(
                    test_temp, properties=specific_props)
                print(f"  Requested: {specific_props}")
                for prop, value in sorted(specific_values.items()):
                    print(f"  {prop:<30}: {value:.6e}")
            else:
                print("  No properties available to demonstrate")
        except Exception as e:
            raise ValueError(f"Error in Method 4: {e}")

        # Method 5: Temperature-dependent only (exclude constants)
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
