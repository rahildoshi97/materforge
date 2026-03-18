# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This application depends on waLBerla and pystencils (GPLv3), requiring GPL licensing.

"""Visualization script for materforge call graphs using pycallgraph2."""
import sys
from pathlib import Path

import sympy as sp
import pystencils as ps
from pycallgraph2 import PyCallGraph, Config
from pycallgraph2.output import GraphvizOutput
from pycallgraph2.globbing_filter import GlobbingFilter

from materforge import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

sys.path.append(str(Path(__file__).parent.parent))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_APPS_DIR   = Path(__file__).parent
_SRC_DIR    = _APPS_DIR.parent / "src" / "materforge" / "data" / "materials"
_IMAGE_DIR  = _APPS_DIR / "pycallgraph2_images"

YAML_AL    = _SRC_DIR / "Al.yaml"
YAML_SS    = _SRC_DIR / "1.4301.yaml"
YAML_SS_HE = _APPS_DIR / "1.4301_HeatEquationKernelWithMaterial.yaml"

# Any plain sp.Symbol is valid here.
# Do NOT pass u.center() or any pystencils Field.Access directly -
# substitute it after creation at the pystencils boundary.
S = sp.Symbol('S')

# ---------------------------------------------------------------------------
# Infrastructure
# ---------------------------------------------------------------------------

def _image_dir() -> Path:
    _IMAGE_DIR.mkdir(exist_ok=True)
    return _IMAGE_DIR


def visualize(target_fn, output_name: str, include_patterns: list,
              max_depth: int = 10) -> None:
    """Run `target_fn` inside a PyCallGraph context and save the SVG."""
    config = Config(max_depth=max_depth)
    config.trace_filter = GlobbingFilter(
        include=include_patterns,
        exclude=['pycallgraph2.*', 'logging.*', 'matplotlib.*', 'PIL.*', 'fontTools.*'],
    )
    output_file = _image_dir() / f"{output_name}_callgraph.svg"
    graphviz = GraphvizOutput(
        output_file=str(output_file),
        font_name='Verdana',
        font_size=8,
        output_type='svg',
        dpi=300,
    )
    with PyCallGraph(config=config, output=graphviz):
        try:
            target_fn()
        except Exception as e:
            print(f"  [!] {target_fn.__name__} raised: {e}")
    print(f"  Saved: {output_file}")

# ---------------------------------------------------------------------------
# Material helpers
# ---------------------------------------------------------------------------

def _create(yaml_path: Path, enable_plotting: bool = False):
    if not yaml_path.exists():
        print(f"  [!] YAML not found: {yaml_path}")
        return None
    try:
        mat = create_material(yaml_path=yaml_path, dependency=S, enable_plotting=enable_plotting)
        print(f"  Created: {mat.name}  ({len(mat.property_names())} properties)")
        return mat
    except Exception as e:
        print(f"  [!] Failed to create material from {yaml_path.name}: {e}")
        return None

def _analyze(mat) -> bool:
    """Print a summary of the material. Returns True if energy_density is present."""
    print(f"\n  Material : {mat.name}")
    print(f"  Properties: {sorted(mat.property_names())}")
    has_ed = hasattr(mat, 'energy_density') and mat.energy_density is not None
    print(f"  energy_density: {'present' if has_ed else 'absent'}")
    if has_ed:
        print(f"    free_symbols: {mat.energy_density.free_symbols}")
    return has_ed

def _test_inverse(mat) -> sp.Expr | None:
    """Create and round-trip-test the inverse of energy_density. Returns inverse or None."""
    if not (hasattr(mat, 'energy_density') and mat.energy_density is not None):
        print(f"  [!] {mat.name} has no energy_density - skipping inverse")
        return None
    try:
        E = sp.Symbol('E')
        inv = PiecewiseInverter.create_inverse(mat.energy_density, S, E)
        print(f"  Inverse created for {mat.name}")
        for temp in [300.0, 500.0, 1000.0, 1500.0]:
            try:
                e_val    = float(mat.energy_density.subs(S, temp).evalf())
                t_rec    = float(inv.subs(E, e_val))
                error    = abs(temp - t_rec)
                status   = "✓" if error < 1e-6 else "!" if error < 1e-3 else "✗"
                print(f"    {status} S={temp:7.1f}K -> E={e_val:.3e} -> S={t_rec:7.1f}K  err={error:.2e}")
            except Exception as e:
                print(f"    [!] S={temp}K: {e}")
        return inv
    except Exception as e:
        print(f"  [!] Inverse creation failed for {mat.name}: {e}")
        return None

# ---------------------------------------------------------------------------
# Workflow functions (called inside visualize())
# ---------------------------------------------------------------------------

def workflow_material_creation() -> None:
    """Create all standard materials and analyse their properties."""
    print("\n=== Material Creation ===")
    for path in (YAML_AL, YAML_SS):
        mat = _create(path, enable_plotting=True)
        if mat:
            _analyze(mat)

def workflow_inverse_functions() -> None:
    """Create materials and test inverse energy density functions."""
    print("\n=== Inverse Function Testing ===")
    for path in (YAML_AL, YAML_SS):
        mat = _create(path)
        if mat:
            _test_inverse(mat)

def workflow_heat_equation() -> None:
    """Full heat equation kernel workflow: material -> pystencils assignment collection.

    materforge produces expressions in S (plain SymPy symbol).
    The single .subs(S, u.center()) call is the explicit coupling point
    between materforge and pystencils - it must NOT move into create_material().
    """
    print("\n=== Heat Equation Workflow ===")

    yaml_path = YAML_SS_HE if YAML_SS_HE.exists() else YAML_SS
    mat = _create(yaml_path, enable_plotting=True)
    if mat is None:
        return

    _analyze(mat)
    _test_inverse(mat)

    # pystencils setup
    data_type = "float64"
    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    thermal_diffusivity_sym   = sp.Symbol("thermal_diffusivity")
    thermal_diffusivity_field = ps.fields(f"thermal_diffusivity_field: {data_type}[2D]", layout='fzyx')
    dx, dt = sp.Symbol("dx"), sp.Symbol("dt")

    heat_pde         = ps.fd.transient(u) - thermal_diffusivity_sym * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1))
    discretize       = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
    heat_pde_disc    = discretize(heat_pde)
    heat_pde_disc    = heat_pde_disc.args[1] + heat_pde_disc.args[0].simplify()

    def to_ps(expr: sp.Expr) -> sp.Expr:
        """Substitute materforge placeholder S -> pystencils field accessor."""
        return expr.subs(S, u.center())

    ac = ps.AssignmentCollection(
        subexpressions=[
            ps.Assignment(thermal_diffusivity_sym, to_ps(mat.thermal_diffusivity)),
        ],
        main_assignments=[
            ps.Assignment(u_tmp.center(), heat_pde_disc),
            ps.Assignment(thermal_diffusivity_field.center(), thermal_diffusivity_sym),
        ],
    )
    print(f"\n  AssignmentCollection:\n{ac}")

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("materforge pycallgraph2 visualizations\n" + "=" * 60)

    visualize(
        workflow_material_creation,
        output_name="material_creation",
        include_patterns=[
            'materforge.parsing.*',
            'materforge.core.materials.*',
            'create_material',
        ],
        max_depth=12,
    )

    visualize(
        workflow_inverse_functions,
        output_name="inverse_functions",
        include_patterns=[
            'materforge.algorithms.piecewise_inverter.*',
            'PiecewiseInverter.*',
            'create_inverse',
        ],
    )

    visualize(
        workflow_heat_equation,
        output_name="heat_equation_workflow",
        include_patterns=[
            'materforge.*',
            'pystencils.*',
        ],
        max_depth=10,
    )

    print("\nAll visualizations completed.")
