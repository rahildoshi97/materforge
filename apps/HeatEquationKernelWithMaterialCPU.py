import logging
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

material_library = 0 # Set to 1 to use material library, 0 to use hardcoded values

with SourceFileGenerator() as sfg:
    data_type = "float64"

    f_u, f_u_tmp, f_k, f_rho, f_cp, f_alpha = ps.fields(f"f_u, f_u_tmp, f_k, f_rho, f_cp, f_alpha: {data_type}[3D]", layout='fzyx')
    s = SymbolCreator()
    s_dx, s_dt, s_k, s_rho, s_cp, s_alpha = s.s_dx, s.s_dt, s.s_k, s.s_rho, s.s_cp, s.s_alpha

    alpha = s_k / (s_rho * s_cp)

    heat_pde = ps.fd.transient(f_u) - alpha * (ps.fd.diff(f_u, 0, 0) + ps.fd.diff(f_u, 1, 1))

    discretize = ps.fd.Discretization2ndOrder(dx=s_dx, dt=s_dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify()

    yaml_path = Path(__file__).parent / '1.4301_HeatEquationKernelWithMaterialBM.yaml'
    mat = create_material(yaml_path=yaml_path, dependency=f_u.center(), enable_plotting=True)

    if material_library:
        subexp = [
            ps.Assignment(s_k, mat.heat_conductivity),
            ps.Assignment(s_rho, mat.density),
            ps.Assignment(s_cp, mat.heat_capacity),
        ]
    else:
        subexp = [
            ps.Assignment(s_k, 35.3375739419170),
            ps.Assignment(s_rho, 6824.46293263393),
            ps.Assignment(s_cp, 835.533555786031),
        ]

    ac = ps.AssignmentCollection(
        subexpressions=subexp,
        main_assignments=[
            ps.Assignment(f_u_tmp.center(), heat_pde_discretized),
            ps.Assignment(f_alpha.center(), alpha)
        ])

    print(f"ac:\n{ac}")

    cpu_config = ps.CreateKernelConfig(target=ps.Target.CPU)

    sweep = Sweep("HeatEquationKernelWithMaterialCPU", ac, config=cpu_config)
    sfg.generate(sweep)
