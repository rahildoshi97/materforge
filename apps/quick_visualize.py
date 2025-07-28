"""Enhanced visualization script with better material property handling."""
import sys
from pathlib import Path
from pycallgraph2 import PyCallGraph, Config
from pycallgraph2.output import GraphvizOutput
from pycallgraph2.globbing_filter import GlobbingFilter

# Add pymatlib to path
sys.path.append(str(Path(__file__).parent.parent))


def quick_visualize(target_function, output_name, include_patterns):
    """Quick visualization for any target function."""
    image_folder = Path(__file__).parent / 'callgraph_images'
    image_folder.mkdir(exist_ok=True)

    config = Config(max_depth=10)
    config.trace_filter = GlobbingFilter(
        include=include_patterns,
        exclude=['pycallgraph2.*', 'logging.*', 'matplotlib.*']
    )

    output_file = image_folder / f'{output_name}_callgraph.svg'
    graphviz = GraphvizOutput(
        output_file=str(output_file),
        font_size=8,
        output_type='svg',
        dpi=300,
    )

    with PyCallGraph(config=config, output=graphviz):
        try:
            result = target_function()
            print(f"Function executed successfully: {type(result)}")
        except Exception as e:
            print(f"Function failed with error: {e}")

    print(f"Visualization saved to: {output_file}")


def analyze_material_properties(mat):
    """Analyze and report material properties."""
    print(f"\n--- Material Analysis: {mat.name} ---")
    print(f"Type: {mat.material_type}")
    print(f"Elements: {[elem.name for elem in mat.elements]}")

    # Check for energy density
    if hasattr(mat, 'energy_density') and mat.energy_density is not None:
        print(f"✓ Energy density: {type(mat.energy_density)}")
        print(f"  Symbols: {mat.energy_density.free_symbols}")
        return True
    else:
        print(f"✗ Energy density: None (not defined in YAML)")
        return False


def test_material_with_energy_density():
    """Test material creation focusing on materials with energy density."""
    import sympy as sp
    from pymatlib.parsing.api import create_material

    T = sp.Symbol('T')

    # Test multiple materials to find ones with energy density
    yaml_paths = [
        Path(__file__).parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml",
        Path(__file__).parent.parent / "src" / "pymatlib" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml",
    ]

    materials_with_energy = []

    for yaml_path in yaml_paths:
        if yaml_path.exists():
            try:
                mat = create_material(yaml_path=yaml_path, T=T, enable_plotting=False)
                has_energy = analyze_material_properties(mat)
                if has_energy:
                    materials_with_energy.append(mat)
            except Exception as e:
                print(f"Failed to create material from {yaml_path}: {e}")

    return materials_with_energy


def test_inverse_with_working_material():
    """Test inverse function with materials that have energy density."""
    import sympy as sp
    from pymatlib.algorithms.piecewise_inverter import PiecewiseInverter

    materials = test_material_with_energy_density()

    for mat in materials:
        if hasattr(mat, 'energy_density') and mat.energy_density is not None:
            try:
                print(f"\n--- Testing inverse for {mat.name} ---")

                # Method 2: Direct approach
                energy_symbols = mat.energy_density.free_symbols
                if len(energy_symbols) == 1:
                    temp_symbol = list(energy_symbols)[0]
                    E_symbol = sp.Symbol('E')

                    inverse_func = PiecewiseInverter.create_inverse(mat.energy_density, temp_symbol, E_symbol)
                    print(f"✓ Inverse function created successfully")

                    # Test evaluation
                    test_temp = 300.0
                    energy_val = float(mat.energy_density.subs(temp_symbol, test_temp))
                    recovered_temp = float(inverse_func.subs(E_symbol, energy_val))
                    error = abs(test_temp - recovered_temp)
                    print(
                        f"✓ Round-trip test: T={test_temp} -> E={energy_val:.2e} -> T={recovered_temp:.1f}, Error={error:.2e}")

                    return inverse_func
                else:
                    print(f"✗ Unexpected symbols in energy density: {energy_symbols}")
            except Exception as e:
                print(f"✗ Inverse creation failed: {e}")

    return None


def test_heat_equation_workflow():
    """Test the complete heat equation workflow."""
    import sympy as sp
    import pystencils as ps
    from pymatlib.parsing.api import create_material
    from pymatlib.algorithms.piecewise_inverter import PiecewiseInverter

    print("\n--- Heat Equation Workflow Test ---")

    # Create fields
    data_type = "float64"
    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')

    # Test with 1.4301 (known to have energy density)
    yaml_path = Path(
        __file__).parent.parent / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

    if yaml_path.exists():
        try:
            mat = create_material(yaml_path=yaml_path, T=u.center(), enable_plotting=False)
            print(f"✓ Created material: {mat.name}")

            if hasattr(mat, 'energy_density') and mat.energy_density is not None:
                # Create inverse function
                energy_symbols = mat.energy_density.free_symbols
                if len(energy_symbols) == 1:
                    temp_symbol = list(energy_symbols)[0]
                    E_symbol = sp.Symbol('E')

                    inverse_func = PiecewiseInverter.create_inverse(mat.energy_density, temp_symbol, E_symbol)
                    print(f"✓ Inverse function created for heat equation")

                    # Test thermal diffusivity
                    if hasattr(mat, 'thermal_diffusivity'):
                        print(f"✓ Thermal diffusivity available: {type(mat.thermal_diffusivity)}")

                    return mat
            else:
                print(f"✗ Material {mat.name} has no energy density")
        except Exception as e:
            print(f"✗ Heat equation workflow failed: {e}")

    return None


def create_comprehensive_demo():
    """Create a comprehensive demonstration of working functionality."""
    from examples.material_properties_demo import demonstrate_material_properties

    print("\n--- Running Comprehensive Demo ---")
    try:
        demonstrate_material_properties()
        return True
    except Exception as e:
        print(f"Demo failed: {e}")
        return False


if __name__ == "__main__":
    print("Creating enhanced material visualization call graphs...")

    # Test material creation with detailed analysis
    quick_visualize(
        test_material_with_energy_density,
        'enhanced_material_creation',
        ['pymatlib.parsing.*', 'pymatlib.core.materials.*', 'create_material']
    )

    # Test inverse function with working materials
    quick_visualize(
        test_inverse_with_working_material,
        'working_inverse_function',
        ['pymatlib.algorithms.piecewise_inverter.*', 'PiecewiseInverter.*', 'create_inverse']
    )

    # Test heat equation workflow
    quick_visualize(
        test_heat_equation_workflow,
        'heat_equation_workflow',
        ['pymatlib.*', 'pystencils.*', '!pymatlib.*.logging.*']
    )

    # Test comprehensive demo
    quick_visualize(
        create_comprehensive_demo,
        'comprehensive_demo',
        ['pymatlib.*', '!pymatlib.*.logging.*', '!pymatlib.*.matplotlib.*']
    )

    print("All enhanced visualizations completed!")
