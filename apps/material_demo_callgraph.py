"""Visualize material properties demo functionality with robust error handling."""
import os
import sys
from pathlib import Path
from pycallgraph2 import PyCallGraph, Config
from pycallgraph2.output import GraphvizOutput
from pycallgraph2.globbing_filter import GlobbingFilter

# Add materforge to path
sys.path.append(str(Path(__file__).parent.parent))


def setup_visualization_environment():
    """Setup directories and paths for visualization."""
    apps_directory = Path(__file__).parent
    image_folder = apps_directory / 'callgraph_images'
    image_folder.mkdir(exist_ok=True)
    return image_folder


def create_robust_inverse_visualization():
    """Create visualization that handles both successful and failed inverse creation."""
    image_folder = setup_visualization_environment()

    config = Config(max_depth=12, min_depth=1)
    config.trace_filter = GlobbingFilter(
        include=[
            'materforge.algorithms.piecewise_inverter.*',
            'materforge.parsing.api.*',
            'materforge.core.materials.*',
            'PiecewiseInverter.*',
            'create_material',
            'create_inverse',
            'create_energy_density_inverse',
        ],
        exclude=[
            'pycallgraph2.*',
            'logging.*',
            'matplotlib.*',
        ]
    )

    output_file = image_folder / 'robust_inverse_function_callgraph.svg'
    graphviz = GraphvizOutput(
        output_file=str(output_file),
        font_name='Verdana',
        font_size=8,
        output_type='svg',
        dpi=300,
    )

    with PyCallGraph(config=config, output=graphviz):
        import sympy as sp
        from materforge.parsing.api import create_material
        from materforge.algorithms.piecewise_inverter import PiecewiseInverter

        T = sp.Symbol('T')
        current_file = Path(__file__)
        yaml_paths = [
            current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml",
            current_file.parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"
        ]

        for yaml_path in yaml_paths:
            if yaml_path.exists():
                try:
                    mat = create_material(yaml_path=yaml_path, dependency=T, enable_plotting=True)
                    print(f"Created material: {mat.name}")

                    if hasattr(mat, 'energy_density'):
                        # Method 1: Try convenience function (may fail)
                        try:
                            E = sp.Symbol('E')
                            inverse_func1 = PiecewiseInverter.create_energy_density_inverse(mat, 'E')
                            print(f"Method 1 succeeded for {mat.name}")
                        except ValueError as e:
                            print(f"Method 1 failed for {mat.name}: {e}")

                        # Method 2: Try direct approach (more likely to work)
                        try:
                            energy_symbols = mat.energy_density.free_symbols # type: ignore
                            if len(energy_symbols) == 1:
                                temp_symbol = list(energy_symbols)[0]
                                E_symbol = sp.Symbol('E')
                                inverse_func2 = PiecewiseInverter.create_inverse(mat.energy_density, temp_symbol, # type: ignore
                                                                                 E_symbol)
                                print(f"Method 2 succeeded for {mat.name}")

                                # Test a few evaluations
                                test_temps = [300, 500, 1000]
                                for temp in test_temps:
                                    try:
                                        energy_val = float(mat.energy_density.subs(temp_symbol, temp)) # type: ignore
                                        recovered_temp = float(inverse_func2.subs(E_symbol, energy_val)) # type: ignore
                                        print(f"T={temp} -> E={energy_val:.2e} -> T={recovered_temp:.1f}")
                                    except:
                                        pass
                        except Exception as e:
                            print(f"Method 2 failed for {mat.name}: {e}")

                except Exception as e:
                    print(f"Failed to create material from {yaml_path}: {e}")

    print(f"Robust inverse visualization saved to: {output_file}")


def create_heat_equation_workflow_visualization():
    """Create visualization for the heat equation workflow from your existing code."""
    image_folder = setup_visualization_environment()

    config = Config(max_depth=10)
    config.trace_filter = GlobbingFilter(
        include=[
            'materforge.*',
            'pystencils.*',
            'create_material',
            'PiecewiseInverter.*',
        ],
        exclude=[
            'pycallgraph2.*',
            'logging.*',
            'matplotlib.*',
            'walberla.*',  # Exclude walberla internals
        ]
    )

    output_file = image_folder / 'heat_equation_workflow_callgraph.svg'
    graphviz = GraphvizOutput(
        output_file=str(output_file),
        font_name='Verdana',
        font_size=8,
        output_type='svg',
        dpi=300,
    )

    with PyCallGraph(config=config, output=graphviz):
        # Simplified version of your heat equation code
        import sympy as sp
        import pystencils as ps
        from materforge.parsing.api import create_material
        from materforge.algorithms.piecewise_inverter import PiecewiseInverter

        # Create fields
        data_type = "float64"
        u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')

        # Load material
        yaml_path = Path(__file__).parent / '1.4301_HeatEquationKernelWithMaterial.yaml'
        if not yaml_path.exists():
            yaml_path = Path(
                __file__).parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

        if yaml_path.exists():
            try:
                mat = create_material(yaml_path=yaml_path, dependency=u.center(), enable_plotting=True) # type: ignore
                print(f"Created material for heat equation: {mat.name}")

                # Test inverse function creation (Method 2)
                if hasattr(mat, 'energy_density'):
                    try:
                        energy_symbols = mat.energy_density.free_symbols # type: ignore
                        if len(energy_symbols) == 1:
                            temp_symbol = list(energy_symbols)[0]
                            E_symbol = sp.Symbol('E')
                            inverse_func = PiecewiseInverter.create_inverse(mat.energy_density, temp_symbol, E_symbol) # type: ignore
                            print("Inverse function created for heat equation")
                    except Exception as e:
                        print(f"Inverse creation failed: {e}")
            except Exception as e:
                print(f"Material creation failed: {e}")

    print(f"Heat equation workflow visualization saved to: {output_file}")


if __name__ == "__main__":
    print("Creating robust material and inverse function visualizations...")

    # Create visualizations that handle both success and failure cases
    create_robust_inverse_visualization()
    create_heat_equation_workflow_visualization()

    print("All visualizations completed!")
