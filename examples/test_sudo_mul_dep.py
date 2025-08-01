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

    for yaml_path, name in [(yaml_path_Al, "Aluminum")]:
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

    for yaml_path, name in [(yaml_path_Al, "Aluminum")]:
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
            mat_Al = create_material(yaml_path=yaml_path_Al, dependency_symbols={'temperature': sp.Symbol('u_C'), 'pressure': sp.Symbol('p')}, enable_plotting=True)
            materials.append(mat_Al)
            material_symbols[mat_Al.name] = T1
            print(f"✓ Successfully created: {mat_Al.name} (using symbol {T1})")
        except Exception as e:
            print(f"✗ Failed to create Aluminum: {e}")

    """if yaml_path_SS304L.exists():
        try:
            mat_SS304L = create_material(yaml_path=yaml_path_SS304L, dependency_symbols=T2, enable_plotting=True)
            materials.append(mat_SS304L)
            material_symbols[mat_SS304L.name] = T2
            print(f"✓ Successfully created: {mat_SS304L.name} (using symbol {T2})")
        except Exception as e:
            print(f"✗ Failed to create SS304L: {e}")"""

if __name__ == "__main__":
    demonstrate_material_properties()