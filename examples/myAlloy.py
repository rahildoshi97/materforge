# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import sympy as sp

from materforge.parsing.api import (
    create_material,
    get_supported_properties,
    evaluate_material_properties,
    get_material_property_names,
    validate_yaml_file,
    get_material_info
)


def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.WARNING,  # DEBUG/INFO/WARNING/ERROR/CRITICAL
        format="%(asctime)s %(levelname)s %(name)s -> %(message)s"
    )
    # Silence noisy libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('fontTools').setLevel(logging.WARNING)


def demonstrate_material_properties():
    """Comprehensive demonstration of material property evaluation."""
    setup_logging()

    # Create symbolic temperature variables for different materials
    T = sp.Symbol('T')  # Temperature symbol for myAlloy

    # Set up file paths
    current_file = Path(__file__)
    yaml_path = current_file.parent / "myAlloy.yaml"

    materials = []
    material_symbols = {}  # Dictionary to track which symbol goes with which material

    print(f"\n{'=' * 80}")
    print("MaterForge MATERIAL PROPERTY DEMONSTRATION")
    print(f"{'=' * 80}")

    # ===================================================================
    # 1. YAML FILE VALIDATION
    # ===================================================================
    print(f"\n{'1. YAML FILE VALIDATION':<50}")
    print(f"{'-' * 50}")

    for yaml_path, name in [(yaml_path, "myAlloy")]:
        if yaml_path.exists():
            try:
                is_valid = validate_yaml_file(yaml_path)
                print(f"✓ {name} YAML validation: {'PASSED' if is_valid else 'FAILED'}")
            except Exception as e:
                print(f"✗ {name} YAML validation: FAILED - {e}")
        else:
            print(f"✗ {name} YAML file not found: {yaml_path}")

    # ===================================================================
    # 2. MATERIAL INFO EXTRACTION
    # ===================================================================
    print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
    print(f"{'-' * 50}")

    for yaml_path, name in [(yaml_path, "myAlloy")]:
        if yaml_path.exists():
            try:
                info = get_material_info(yaml_path)
                print(f"\n{name} Information:")
                print(f"  Name: {info['name']}")
                print(f"  Type: {info['material_type']}")
                print(f"  Total Properties: {info['total_properties']}")
                print(f"  Composition: {info['composition']}")

                # Temperature properties
                if info['material_type'] == 'pure_metal':
                    print(f"  Melting Temperature: {info.get('melting_temperature', 'N/A')} K")
                    print(f"  Boiling Temperature: {info.get('boiling_temperature', 'N/A')} K")
                else:
                    print(f"  Solidus Temperature: {info.get('solidus_temperature', 'N/A')} K")
                    print(f"  Liquidus Temperature: {info.get('liquidus_temperature', 'N/A')} K")

            except Exception as e:
                print(f"✗ Failed to get {name} info: {e}")

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
            print(f"Successfully created: {myAlloy.name} (using symbol {T})")
        except Exception as e:
            print(f"Failed to create myAlloy: {e}")

    # ===================================================================
    # 4. SUPPORTED PROPERTIES OVERVIEW
    # ===================================================================
    print(f"\n{'4. SUPPORTED PROPERTIES OVERVIEW':<50}")
    print(f"{'-' * 50}")

    supported_props = get_supported_properties()
    print(f"MaterForge supports {len(supported_props)} property types:")
    for i, prop in enumerate(sorted(supported_props), 0):
        print(f"  {i:2d}. {prop}")

    # ===================================================================
    # 5. MATERIAL PROPERTY ANALYSIS
    # ===================================================================
    print(f"\n{'5. MATERIAL PROPERTY ANALYSIS':<50}")
    print(f"{'-' * 50}")

    for mat in materials:
        # Get the appropriate temperature symbol for this material
        T_mat = material_symbols[mat.name]

        print(f"\n{'=' * 80}")
        print(f"MATERIAL: {mat.name} (Temperature symbol: {T_mat})")
        print(f"{'=' * 80}")

        # Basic material info
        print(f"Name: {mat.name}")
        print(f"Type: {mat.material_type}")
        print(f"Elements: {[elem.name for elem in mat.elements]}")
        print(f"Composition: {mat.composition}")
        for i in range(len(mat.composition)):
            print(f"  {mat.elements[i].name}: {mat.composition[i]}")

        # Temperature properties
        if hasattr(mat, 'solidus_temperature') and hasattr(mat, 'liquidus_temperature'):
            print(f"Solidus Temperature: {mat.solidus_temperature}")
            print(f"Liquidus Temperature: {mat.liquidus_temperature}")
        if hasattr(mat, 'melting_temperature') and hasattr(mat, 'boiling_temperature'):
            print(f"Melting Temperature: {mat.melting_temperature}")
            print(f"Boiling Temperature: {mat.boiling_temperature}")
        if hasattr(mat, 'initial_boiling_temperature') and hasattr(mat, 'final_boiling_temperature'):
            print(f"Initial Boiling Temperature: {mat.initial_boiling_temperature} K")
            print(f"Final Boiling Temperature: {mat.final_boiling_temperature} K")

        # ===================================================================
        # 5.1 Available Properties on Material
        # ===================================================================
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")

        available_props = get_material_property_names(mat)
        print(f"{myAlloy.name} has {len(available_props)} properties:")
        print(available_props)

        # ===================================================================
        # 5.2 Property Evaluation at Specific Temperature
        # ===================================================================
        test_temp = 500.15  # Kelvin
        print(f"\n{'PROPERTY VALUES AT ' + str(test_temp) + 'K:':<50}")
        print(f"{'-' * 50}")

        # Method 1: Manual evaluation (old way)
        print("Method 1: Manual Property Evaluation")
        for prop_name in sorted(available_props):
            try:
                prop_value = getattr(mat, prop_name)
                if isinstance(prop_value, sp.Expr):
                    numerical_value = prop_value.subs(T_mat, test_temp).evalf() # type: ignore
                    print(f"{prop_name:<30}: {numerical_value} (symbolic)")
                else:
                    print(f"{prop_name:<30}: {prop_value} (constant)")
            except Exception as e:
                print(f"{prop_name:<30}: Error - {str(e)}")

        # ===================================================================
        # 5.3 Using New Property Evaluation APIs
        # ===================================================================
        print(f"\n{'NEW API METHODS:':<50}")
        print(f"{'-' * 50}")

        # Method 2: Object-oriented approach
        print("Method 2: Material.evaluate_properties_at_temperature()")
        try:
            all_values = mat.evaluate_properties_at_temperature(test_temp)
            print(f"All properties at {test_temp}K:")
            for prop, value in sorted(all_values.items()):
                print(f"  {prop:<28}: {value:.6e}")
        except Exception as e:
            print(f"Error: {e}")

        # Method 3: Functional approach
        print("\nMethod 3: evaluate_material_properties()")
        try:
            all_values_func = evaluate_material_properties(mat, test_temp)
            print(f"Functional API results match: {all_values == all_values_func}") # type: ignore
        except Exception as e:
            print(f"Error: {e}")

        # Method 4: Specific properties only
        print("\nMethod 4: Specific Properties Only")
        try:
            specific_props = ['density', 'heat_capacity', 'energy_density']
            available_specific = [p for p in specific_props if p in available_props]

            if available_specific:
                specific_values = mat.evaluate_properties_at_temperature(
                    test_temp,
                    properties=available_specific
                )
                print(f"Specific properties at {test_temp}K:")
                for prop, value in sorted(specific_values.items()):
                    print(f"  {prop:<28}: {value:.6e}")
            else:
                print("None of the requested specific properties available")
        except Exception as e:
            print(f"Error: {e}")

        # Method 5: Temperature-dependent properties only (exclude constants)
        print("\nMethod 5: Temperature-Dependent Properties Only")
        try:
            temp_dependent = mat.evaluate_properties_at_temperature(
                test_temp,
                include_constants=False
            )
            print(f"Temperature-dependent properties at {test_temp}K:")
            for prop, value in sorted(temp_dependent.items()):
                print(f"  {prop:<28}: {value:.6e}")
        except Exception as e:
            print(f"Error: {e}")


if __name__ == "__main__":
    demonstrate_material_properties()
