import logging
# import sympy as sp
import pystencils as ps
from pathlib import Path
from pystencilssfg import SourceFileGenerator
from pystencils import SymbolCreator
from walberla.codegen import Sweep

from materforge.parsing.api import create_material

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

    f_u, f_u_tmp, f_k, f_rho, f_cp, f_alpha = ps.fields(f"f_u, f_u_tmp, f_k, f_rho, f_cp, f_alpha: {data_type}[2D]", layout='fzyx')
    s = SymbolCreator()
    s_dx, s_dt, s_k, s_rho, s_cp, s_alpha = s.s_dx, s.s_dt, s.s_k, s.s_rho, s.s_cp, s.s_alpha

    alpha = s_k / (s_rho * s_cp)
    heat_pde = ps.fd.transient(f_u) - alpha * (ps.fd.diff(f_u, 0, 0) + ps.fd.diff(f_u, 1, 1))

    discretize = ps.fd.Discretization2ndOrder(dx=s_dx, dt=s_dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify() # type: ignore

    yaml_path = Path(__file__).parent / '1.4301_HeatEquationKernelWithMaterial.yaml'
    yaml_path_Al = Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"
    yaml_path_SS304L = Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"

    mat = create_material(yaml_path=yaml_path, dependency=f_u.center(), enable_plotting=True)
    mat_Al = create_material(yaml_path=yaml_path_Al, dependency=f_u.center(), enable_plotting=True)
    mat_SS304L = create_material(yaml_path=yaml_path_SS304L, dependency=f_u.center(), enable_plotting=True)

    subexp = [
        ps.Assignment(s_k, mat.heat_conductivity),
        ps.Assignment(s_rho, mat.density),
        ps.Assignment(s_cp, mat.heat_capacity),
        ps.Assignment(s_alpha, mat.thermal_diffusivity),
    ]

    ac = ps.AssignmentCollection(
        subexpressions=subexp,
        main_assignments=[
            ps.Assignment(f_u_tmp.center(), heat_pde_discretized),
            # ps.Assignment(f_k.center(), s_k),
            # ps.Assignment(f_rho.center(), s_rho),
            # ps.Assignment(f_cp.center(), s_cp),
            ps.Assignment(f_alpha.center(), s_alpha)
        ])

    print(f"ac\n{ac}")

    sweep = Sweep("HeatEquationKernelWithMaterial", ac)
    sfg.generate(sweep)
