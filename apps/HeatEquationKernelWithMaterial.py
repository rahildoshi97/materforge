# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This application depends on waLBerla and pystencils (GPLv3), requiring GPL licensing.

import logging
import sympy as sp
import pystencils as ps
from pathlib import Path
from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep

from materforge import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

logging.basicConfig(
    level=logging.WARNING,  # DEBUG/INFO/WARNING/ERROR/CRITICAL
    format="%(asctime)s %(levelname)s %(name)s -> %(message)s"
)

# Silence matplotlib and other noisy libraries
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)

with SourceFileGenerator() as sfg:
    data_type = "float64"

    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    thermal_diffusivity_symbol = sp.Symbol("thermal_diffusivity")
    thermal_diffusivity_field = ps.fields(f"thermal_diffusivity_field: {data_type}[2D]", layout='fzyx')
    dx, dt = sp.Symbol("dx"), sp.Symbol("dt")

    heat_pde = ps.fd.transient(u) - thermal_diffusivity_symbol * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1)) # type: ignore

    discretize = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify() # type: ignore

    yaml_path = Path(__file__).parent / '1.4301_HeatEquationKernelWithMaterial.yaml'
    #yaml_path_Al = Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials" / "Al.yaml"
    #yaml_path_SS304L = Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials" / "1.4301.yaml"

    S = sp.Symbol('S')
    mat = create_material(yaml_path=yaml_path, dependency=S, enable_plotting=True)
    #mat_Al = create_material(yaml_path=yaml_path_Al, dependency=u.center(), enable_plotting=True) # type: ignore
    #mat_SS304L = create_material(yaml_path=yaml_path_SS304L, dependency=u.center(), enable_plotting=True) # type: ignore

    print(f"Energy density function: {mat.energy_density}")
    print(f"Type: {type(mat.energy_density)}")
    print("=" * 80)

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

    print("Create inverse energy density function Using PiecewiseInverter")
    print("-" * 60)

    if hasattr(mat, 'energy_density'):
        try:
            energy_symbols = mat.energy_density.free_symbols # type: ignore
            if len(energy_symbols) != 1:
                raise ValueError(f"Energy density function must have exactly one symbol, found: {energy_symbols}")

            temp_symbol = list(energy_symbols)[0]
            E_symbol = sp.Symbol('E')

            # Create inverter with custom tolerance
            inverse_func = PiecewiseInverter.create_inverse(mat.energy_density, temp_symbol, E_symbol) # type: ignore

            print("Inverse created successfully!")
            print(f"Temperature symbol used: {temp_symbol}")
            print(f"Inverse function: {inverse_func}")

            print("\nRound-trip accuracy test:")
            errors = []
            passed = 0
            failed = 0

            for temp in test_temperatures:
                try:
                    # Forward: T -> E
                    energy_val = float(mat.energy_density.subs(temp_symbol, temp).evalf()) # type: ignore
                    # Backward: E -> T
                    recovered_temp = float(inverse_func.subs(E_symbol, energy_val)) # type: ignore
                    error = abs(temp - recovered_temp)
                    errors.append(error)

                    status = "✅" if error < 1e-6 else "⚠️" if error < 1e-3 else "❌"
                    if error < 1e-3:
                        passed += 1
                    else:
                        failed += 1

                    print(f"{status} T={temp:6.1f}K -> E={energy_val:12.2e} -> T={recovered_temp:6.1f}K, Error={error:.2e}")

                except Exception as e:
                    failed += 1
                    print(f"Error at T={temp}K: {e}")

            max_error = max(errors) if errors else float('inf')
            print("\nSummary:")
            print(f"  Passed: {passed}/{len(test_temperatures)}")
            print(f"  Failed: {failed}/{len(test_temperatures)}")
            print(f"  Maximum error: {max_error:.2e}")

        except ValueError as e:
            if "degree" in str(e).lower():
                print(f"Expected error: {e}")
                print("  This is expected if the material has non-linear energy density.")
                print("  The simplified inverter only supports linear piecewise functions.")
                passed = 0
                failed = len(test_temperatures)
            else:
                print(f"Failed: {e}")
                passed = 0
                failed = len(test_temperatures)
        except Exception as e:
            print(f"Failed with unexpected error: {e}")
            passed = 0
            failed = len(test_temperatures)
    else:
        print("Material does not have energy_density property")
        passed = 0
        failed = len(test_temperatures)

    print("\n" + "=" * 80)

    subexp = [
        ps.Assignment(thermal_diffusivity_symbol, mat.thermal_diffusivity.subs(S, u.center())),
    ]

    ac = ps.AssignmentCollection(
        subexpressions=subexp,
        main_assignments=[
            ps.Assignment(u_tmp.center(), heat_pde_discretized), # type: ignore
            ps.Assignment(thermal_diffusivity_field.center(), thermal_diffusivity_symbol) # type: ignore
        ])

    print(f"ac\n{ac}, type = {type(ac)}")

    sweep = Sweep("HeatEquationKernelWithMaterial", ac)
    sfg.generate(sweep)
