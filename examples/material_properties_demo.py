"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import sympy as sp

from pymatlib.parsing.api import create_material, get_material_info, get_supported_properties
from pymatlib.algorithms.piecewise_inverter import PiecewiseInverter


def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,  # DEBUG/INFO/WARNING/ERROR/CRITICAL
        format="%(asctime)s %(levelname)s %(name)s -> %(message)s"
    )
    # Silence noisy libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('fontTools').setLevel(logging.WARNING)


def demonstrate_material_properties():
    """Demonstrate material property evaluation."""
    setup_logging()
    T = sp.Symbol('T')
    T1 = sp.Symbol('u_C')
    # T = 300.15
    current_file = Path(__file__)
    yaml_path_Al = current_file.parent.parent / "src" / "pymatlib" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"
    yaml_path_SS304L = current_file.parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"
    materials = []
    print(f"\n{'=' * 80}")
    if yaml_path_Al.exists():
        mat_Al = create_material(yaml_path=yaml_path_Al, T=T, enable_plotting=True)
        materials.append(mat_Al)
    else:
        raise FileNotFoundError(f"Aluminum YAML file not found: {yaml_path_Al}")
    if yaml_path_SS304L.exists():
        mat_SS304L = create_material(yaml_path=yaml_path_SS304L, T=T1, enable_plotting=True)
        materials.append(mat_SS304L)
    else:
        raise FileNotFoundError(f"1.4301 YAML file not found: {yaml_path_SS304L}")
    print(f"{'=' * 80}")
    print(f"get_material_info_Al:\n{get_material_info(yaml_path=yaml_path_Al)}")
    print(f"{'=' * 80}")
    print(f"get_material_info_SS304L:\n{get_material_info(yaml_path=yaml_path_SS304L)}")
    print(f"{'=' * 80}")
    for mat in materials:
        print(f"\n{'=' * 80}")
        print(f"MATERIAL: {mat.name}")
        print(f"{'=' * 80}")
        print(f"Name: {mat.name}")
        print(f"Type: {mat.material_type}")
        print(f"Elements: {[elem.name for elem in mat.elements]}")
        print(f"Composition: {mat.composition}")
        for i in range(len(mat.composition)):
            print(f"  {mat.elements[i].name}: {mat.composition[i]}")
        if hasattr(mat, 'solidus_temperature'):
            print(f"Solidus Temperature: {mat.solidus_temperature}")
        if hasattr(mat, 'liquidus_temperature'):
            print(f"Liquidus Temperature: {mat.liquidus_temperature}")
        if hasattr(mat, 'melting_temperature'):
            print(f"Melting Temperature: {mat.melting_temperature}")
        if hasattr(mat, 'boiling_temperature'):
            print(f"Boiling Temperature: {mat.boiling_temperature}")
        print(f"{'=' * 80}")
        # Test computed properties at specific temperature
        test_temp = 300.15  # Kelvin
        valid_properties = get_supported_properties()
        print(f"\n{'=' * 80}")
        print(f"PROPERTIES FOR '{mat.name}'")
        print(f"\n{'=' * 80}")
        for valid_prop_name in sorted(valid_properties):
            try:
                if hasattr(mat, valid_prop_name):
                    valid_prop_value = getattr(mat, valid_prop_name)
                    print(f"{valid_prop_name:<30}: {valid_prop_value}")
                    print(f"{' ' * 24} Type : {type(valid_prop_value)}")
                    allowed_types = (sp.Piecewise, sp.Float, type(None))
                    if not isinstance(valid_prop_value, allowed_types):
                        raise TypeError(
                            f"{' ' * 24}  !  WARNING: Unexpected type for property '{valid_prop_name}': {type(valid_prop_value)}."
                            f"Expected: sp.Piecewise, sp.Float, or None")
                    else:
                        print(f"{' ' * 30}  ✓ Valid type")
                else:
                    print(f"{valid_prop_name:<30}: Not defined in YAML")
                    print(f"{' ' * 30} Type: Not defined")
            except Exception as e:
                print(f"{valid_prop_name:<30}: Error - {str(e)}")
        print(f"\n{'=' * 80}")
        print(f"PROPERTY VALUES FOR '{mat.name}' AT {test_temp}K")
        print(f"{'=' * 80}")
        for valid_prop_name in sorted(valid_properties):
            try:
                if hasattr(mat, valid_prop_name):
                    valid_prop_value = getattr(mat, valid_prop_name)
                    if isinstance(valid_prop_value, sp.Expr):
                        numerical_value = valid_prop_value.subs(T, test_temp).evalf()
                        print(f"{valid_prop_name:<30}: {numerical_value} (symbolic)")
                    else:
                        print(f"{valid_prop_name:<30}: {valid_prop_value} (constant)")
                else:
                    print(f"{valid_prop_name:<30}: Not defined in YAML")
            except Exception as e:
                print(f"{valid_prop_name:<30}: Error - {str(e)}")
    test_inverse_functions(materials, T)


def test_inverse_functions(materials, T):
    """Test inverse function creation and accuracy for materials with energy density."""
    print(f"\n{'=' * 80}")
    print("INVERSE FUNCTION TESTING")
    print(f"{'=' * 80}")
    for mat in materials:
        print(f"\n--- Testing Inverse Functions for {mat.name} ---")
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
                test_round_trip_accuracy(mat, inverse_func1, T, E, method="Method 1")
            except ValueError as e:
                print(f"x Method 1 failed: {e}")
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
                    print(f"✓ Method 2 succeeded!)")
                    test_round_trip_accuracy(mat, inverse_func2, temp_symbol, E_symbol, method="Method 2")
                else:
                    print(f"x Unexpected symbols in energy density: {energy_symbols}")
            except Exception as e:
                print(f"x Method 2 failed: {e}")
        except Exception as e:
            print(f"x Inverse testing failed for {mat.name}: {e}")


def test_round_trip_accuracy(material, inverse_func, temp_symbol, energy_symbol, method="Unknown"):
    """Test round-trip accuracy: T -> E -> T."""
    print(f"  Round-trip accuracy test ({method}):")
    test_temperatures = [
        # Low temperature region
        0, 50, 150, 250, 299,
        # Boundary transitions
        300, 301, 350, 500,
        # Heating phase
        800, 1000, 1200, 1500, 1602.147, 1605, 1606.5, 1667, 1667.5,
        # Phase transition region
        1668, 1668.5, 1669, 1670, 1685, 1700, 1725, 1734, 1734.5, 1735,
        # High temperature liquid phase
        1800, 2000, 2500, 2999,
        # Ultra-high temperature
        3000, 3500, 4000
    ]
    max_error = 0.0
    for temp in test_temperatures:
        try:
            # Forward: T -> E
            energy_val = float(material.energy_density.subs(temp_symbol, temp))
            # Backward: E -> T
            recovered_temp = float(inverse_func.subs(energy_symbol, energy_val))
            # Calculate error
            error = abs(temp - recovered_temp)
            max_error = max(max_error, error)
            print(f"    T={temp:4.0f}K -> E={energy_val:.2e} -> T={recovered_temp:6.1f}K (Error: {error:.2e})")
        except Exception as e:
            print(f"    T={temp:4.0f}K -> Error: {e}")
    print(f"  Maximum error: {max_error:.2e} K")
    if max_error < 1e-8:
        print("  ✓ Excellent accuracy!")
    elif max_error < 1e-4:
        print("  ✓ Good accuracy")
    else:
        print("  !  Consider reviewing inverse function accuracy")


if __name__ == "__main__":
    demonstrate_material_properties()
