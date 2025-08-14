"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import sympy as sp
import pandas as pd

from materforge.parsing.api import (
    create_material,
    get_supported_properties,
    evaluate_material_properties,
    get_material_property_names,
    validate_yaml_file,
    get_material_info
)
from materforge.algorithms.piecewise_inverter import PiecewiseInverter


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
    T1 = sp.Symbol('T_Al')  # Temperature symbol for Aluminum
    T2 = sp.Symbol('T_SS')  # Temperature symbol for Steel 1.4301

    # Set up file paths
    current_file = Path(__file__)
    yaml_path_Al = current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"
    yaml_path_SS304L = current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

    materials = []
    material_symbols = {}  # Dictionary to track which symbol goes with which material

    print(f"\n{'=' * 80}")
    print("PYMATLIB MATERIAL PROPERTY DEMONSTRATION")
    print(f"{'=' * 80}")

    # ===================================================================
    # 1. YAML FILE VALIDATION
    # ===================================================================
    print(f"\n{'1. YAML FILE VALIDATION':<50}")
    print(f"{'-' * 50}")

    for yaml_path, name in [(yaml_path_Al, "Aluminum"), (yaml_path_SS304L, "Steel 1.4301")]:
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

    for yaml_path, name in [(yaml_path_Al, "Aluminum"), (yaml_path_SS304L, "SS304L")]:
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

    if yaml_path_Al.exists():
        try:
            mat_Al = create_material(yaml_path=yaml_path_Al, dependency=T1, enable_plotting=True)
            materials.append(mat_Al)
            material_symbols[mat_Al.name] = T1
            print(f"✓ Successfully created: {mat_Al.name} (using symbol {T1})")
        except Exception as e:
            print(f"✗ Failed to create Aluminum: {e}")

    if yaml_path_SS304L.exists():
        try:
            mat_SS304L = create_material(yaml_path=yaml_path_SS304L, dependency=T2, enable_plotting=True)
            materials.append(mat_SS304L)
            material_symbols[mat_SS304L.name] = T2
            print(f"✓ Successfully created: {mat_SS304L.name} (using symbol {T2})")
        except Exception as e:
            print(f"✗ Failed to create SS304L: {e}")

    # ===================================================================
    # 4. SUPPORTED PROPERTIES OVERVIEW
    # ===================================================================
    print(f"\n{'4. SUPPORTED PROPERTIES OVERVIEW':<50}")
    print(f"{'-' * 50}")

    supported_props = get_supported_properties()
    print(f"PyMatLib supports {len(supported_props)} property types:")
    for i, prop in enumerate(sorted(supported_props), 1):
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
        if hasattr(mat, 'solidus_temperature'):
            print(f"Solidus Temperature: {mat.solidus_temperature}")
        if hasattr(mat, 'liquidus_temperature'):
            print(f"Liquidus Temperature: {mat.liquidus_temperature}")
        if hasattr(mat, 'melting_temperature'):
            print(f"Melting Temperature: {mat.melting_temperature}")
        if hasattr(mat, 'boiling_temperature'):
            print(f"Boiling Temperature: {mat.boiling_temperature}")

        # ===================================================================
        # 5.1 Available Properties on Material
        # ===================================================================
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")

        available_props = get_material_property_names(mat)
        print(f"Material has {len(available_props)} properties:")

        for prop_name in sorted(available_props):
            prop_value = getattr(mat, prop_name)
            print(f"{prop_name:<30}: {prop_value}")
            print(f"{' ' * 24} Type : {type(prop_value)}")

            # Validate property types
            allowed_types = (sp.Piecewise, sp.Float, type(None))
            if not isinstance(prop_value, allowed_types):
                print(f"{' ' * 30}  !  WARNING: Unexpected type")
            else:
                print(f"{' ' * 30}  ✓ Valid type")

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
        print(f"\nMethod 3: evaluate_material_properties()")
        try:
            all_values_func = evaluate_material_properties(mat, test_temp)
            print(f"Functional API results match: {all_values == all_values_func}") # type: ignore
        except Exception as e:
            print(f"Error: {e}")

        # Method 4: Specific properties only
        print(f"\nMethod 4: Specific Properties Only")
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
        print(f"\nMethod 5: Temperature-Dependent Properties Only")
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

        # ===================================================================
        # 5.4 Batch Temperature Evaluation
        # ===================================================================
        print(f"\n{'BATCH TEMPERATURE EVALUATION:':<50}")
        print(f"{'-' * 50}")

        temperatures = [300, 400, 500, 600, 700]
        batch_results = {}

        for temp in temperatures:
            try:
                batch_results[temp] = mat.evaluate_properties_at_temperature(temp)
            except Exception as e:
                print(f"Error at {temp}K: {e}")
                batch_results[temp] = {}

        # Create DataFrame for analysis
        if batch_results:
            try:
                df = pd.DataFrame(batch_results).T  # Transpose so temperatures are rows
                print("Temperature variation analysis:")
                print(f"{'Temperature':<12} {'density':<12} {'heat_capacity':<15} {'energy_density':<15}")
                print("-" * 55)

                for temp, row in df.iterrows():
                    density_val = row.get('density', 'N/A')
                    heat_cap_val = row.get('heat_capacity', 'N/A')
                    energy_val = row.get('energy_density', 'N/A')

                    density_str = f"{density_val:.2e}" if density_val != 'N/A' else 'N/A'
                    heat_str = f"{heat_cap_val:.2e}" if heat_cap_val != 'N/A' else 'N/A'
                    energy_str = f"{energy_val:.2e}" if energy_val != 'N/A' else 'N/A'

                    print(f"{temp:<12.0f} {density_str:<12} {heat_str:<15} {energy_str:<15}")

            except Exception as e:
                print(f"DataFrame creation failed: {e}")

    # ===================================================================
    # 6. INVERSE FUNCTION TESTING
    # ===================================================================
    print(f"\n{'6. INVERSE FUNCTION TESTING':<50}")
    print(f"{'-' * 50}")

    test_inverse_functions(materials, material_symbols)


def test_inverse_functions(materials, material_symbols):
    """Test inverse function creation and accuracy for materials with energy density."""
    print(f"\n{'=' * 80}")
    print("INVERSE FUNCTION TESTING")
    print(f"{'=' * 80}")

    for mat in materials:
        # Get the appropriate temperature symbol for this material
        T_mat = material_symbols[mat.name]

        print(f"\n--- Testing Inverse Functions for {mat.name} (symbol: {T_mat}) ---")

        if not hasattr(mat, 'energy_density') or mat.energy_density is None:
            print(f"!  {mat.name} has no energy_density property - skipping inverse test")
            continue

        try:
            # Method 1: Try convenience function (may fail for non-piecewise)
            print("Method 1: Convenience function...")
            try:
                E = sp.Symbol('E')
                inverse_func1 = PiecewiseInverter.create_energy_density_inverse(mat, 'E')
                print("✓ Method 1 succeeded!")
                test_round_trip_accuracy(mat, inverse_func1, T_mat, E, method="Method 1")
            except Exception as e:
                print(f"✗ Method 1 failed: {e}")

            # Method 2: Direct approach (more likely to work)
            print("Method 2: Direct approach...")
            try:
                energy_symbols = mat.energy_density.free_symbols
                if len(energy_symbols) == 1:
                    temp_symbol = list(energy_symbols)[0]
                    E_symbol = sp.Symbol('E')
                    inverse_func2 = PiecewiseInverter.create_inverse(
                        mat.energy_density, temp_symbol, E_symbol
                    )
                    print(f"✓ Method 2 succeeded!")
                    test_round_trip_accuracy(mat, inverse_func2, temp_symbol, E_symbol, method="Method 2")
                else:
                    print(f"✗ Unexpected symbols in energy density: {energy_symbols}")
            except Exception as e:
                print(f"✗ Method 2 failed: {e}")

        except Exception as e:
            print(f"✗ Inverse testing failed for {mat.name}: {e}")


def test_round_trip_accuracy(material, inverse_func, temp_symbol, energy_symbol, method="Unknown"):
    """Test round-trip accuracy: T -> E -> T."""
    print(f"  Round-trip accuracy test ({method}):")

    test_temperatures = [
        # Low temperature region
        300, 350, 400, 450, 500,
        # Medium temperature region
        600, 700, 800, 900, 1000,
        # High temperature region
        1200, 1500, 1800, 2000, 2500,
        # Very high temperature
        3000
    ]

    max_error = 0.0
    successful_tests = 0

    for temp in test_temperatures:
        try:
            # Forward: T -> E
            energy_val = float(material.energy_density.subs(temp_symbol, temp))

            # Backward: E -> T
            recovered_temp = float(inverse_func.subs(energy_symbol, energy_val))

            # Calculate error
            error = abs(temp - recovered_temp)
            max_error = max(max_error, error)
            successful_tests += 1

            print(f"    T={temp:4.0f}K -> E={energy_val:.2e} -> T={recovered_temp:6.1f}K (Error: {error:.2e})")

        except Exception as e:
            print(f"    T={temp:4.0f}K -> Error: {e}")

    print(f"  Tests completed: {successful_tests}/{len(test_temperatures)}")
    print(f"  Maximum error: {max_error:.2e} K")

    if max_error < 1e-8:
        print("  ✓ Excellent accuracy!")
    elif max_error < 1e-4:
        print("  ✓ Good accuracy")
    else:
        print("  !  Consider reviewing inverse function accuracy")


def demonstrate_advanced_usage():
    """Demonstrate advanced usage patterns."""
    print(f"\n{'=' * 80}")
    print("ADVANCED USAGE PATTERNS")
    print(f"{'=' * 80}")

    # Temperature sweep analysis
    print("\n7. TEMPERATURE SWEEP ANALYSIS")
    print("-" * 50)

    # Create materials with different symbols
    current_file = Path(__file__)
    yaml_path_Al = current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"
    yaml_path_SS304L = current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

    T1 = sp.Symbol('T_Al')  # Temperature symbol for Aluminum
    T2 = sp.Symbol('T_SS')  # Temperature symbol for Steel

    materials_advanced = []

    if yaml_path_Al.exists():
        try:
            mat_Al = create_material(yaml_path_Al, dependency=T1, enable_plotting=True)
            materials_advanced.append((mat_Al, T1, "Aluminum"))
            print(f"✓ Created {mat_Al.name} with symbol {T1}")
        except Exception as e:
            print(f"✗ Failed to create Aluminum: {e}")

    if yaml_path_SS304L.exists():
        try:
            mat_SS = create_material(yaml_path_SS304L, dependency=T2, enable_plotting=True)
            materials_advanced.append((mat_SS, T2, "Steel 1.4301"))
            print(f"✓ Created {mat_SS.name} with symbol {T2}")
        except Exception as e:
            print(f"✗ Failed to create Steel: {e}")

    # Temperature sweep for each material
    temp_range = range(300, 3001, 100)  # 300K to 3000K in 100K steps

    for mat, T_symbol, mat_name in materials_advanced:
        print(f"\n--- Temperature Sweep for {mat_name} ---")

        sweep_data = []
        for temp in temp_range:
            try:
                props = mat.evaluate_properties_at_temperature(
                    temp,
                    properties=['density', 'heat_capacity', 'thermal_diffusivity'],
                    include_constants=False
                )
                props['temperature'] = temp
                sweep_data.append(props)
            except Exception as e:
                print(f"Error at {temp}K: {e}")

        if sweep_data:
            df = pd.DataFrame(sweep_data)
            print(f"Temperature sweep results for {mat_name}:")
            print(df.head(5))

            # Simple analysis
            if 'density' in df.columns:
                density_change = ((df['density'].iloc[-1] - df['density'].iloc[0]) / df['density'].iloc[0]) * 100
                print(f"Density change from {temp_range[0]}K to {temp_range[-1]}K: {density_change:.2f}%")


if __name__ == "__main__":
    demonstrate_material_properties()
    demonstrate_advanced_usage()
