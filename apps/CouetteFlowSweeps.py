#!/usr/bin/env python3

"""
CouetteFlowKernel.py - 3D Thermal Couette Flow
Code generation script for Couette flow simulation with temperature-dependent viscosity
"""

import logging
from pathlib import Path

from argparse import ArgumentParser
from dataclasses import replace

from lbmpy import (
    LBStencil,
    Stencil,
    LBMConfig,
    LBMOptimisation,
    Method,
    relaxation_rate_from_lattice_viscosity,
    create_lb_method,
    create_lb_update_rule
)
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

import sympy as sp
import pystencils as ps

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep, get_build_config
from sweepgen.boundaries import GenericHBB
from sweepgen.symbolic import cell, domain
from sweepgen.prefabs import LbmBulk
from sweepgen.build_config import DEBUG

DEBUG.use_cpu_default()
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(levelname)s %(name)s -> %(message)s")

print(f"Starting code generation at {Path(__file__).resolve()}")

use_materforge = False

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    sfg.namespace("CouetteFlow::gen")
    parser = ArgumentParser()
    parser.add_argument("-t", "--target", choices=["cpu", "gpu"], default="cpu")
    args = parser.parse_args(sfg.context.argv)

    build_cfg = get_build_config(sfg)
    build_cfg.target = ps.Target[args.target.upper()]
    
    data_type = "float64"
    stencil = LBStencil(Stencil.D3Q19)
    D, Q = stencil.D, stencil.Q
    # print(f"LBM Configuration: Stencil={stencil}, Dimensions={D}, Q={Q}, Data type={data_type}")

    # Symbolic parameters
    nu, u_max = sp.symbols("nu, u_max")
    T_top, T_bottom = sp.symbols("T_top, T_bottom")  # Temperature bounds

    f_pdfs, f_pdfs_tmp, f_density, f_temperature, f_viscosity, f_velocity = ps.fields(
        f"f_pdfs({Q}), f_pdfs_tmp({Q}), f_density, f_temperature, f_viscosity, f_velocity({D}): {data_type}[{D}D]",
        layout='fzyx'
    )
    
    const_nu = 1./6.  # Constant viscosity if MaterForge is not used

    # Viscosity handling
    if use_materforge:
        print("Code generation: Using temperature-dependent viscosity from MaterForge")
        try:
            from materforge.parsing.api import create_material
            from pathlib import Path

            yaml_path = Path(__file__).parent / 'CouetteFlowMaterial_Updated.yaml'
            
            if yaml_path.exists():
                mat = create_material(yaml_path=yaml_path, dependency=f_temperature.center(), enable_plotting=False)
                nu_expr = mat.dynamic_viscosity if hasattr(mat, 'dynamic_viscosity') else const_nu
                print(f"MaterForge viscosity expression loaded: {nu_expr}")
            else:
                print(f"Warning: Material file {yaml_path} not found, using constant viscosity")
                nu_expr = const_nu
        except Exception as e:
            print(f"MaterForge error: {e}, using constant viscosity")
            nu_expr = const_nu
    else:
        print("Code generation: Using constant viscosity")
        nu_expr = const_nu

    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.TRT,
        compressible=True,
        relaxation_rate=relaxation_rate_from_lattice_viscosity(nu_expr),
    )
    print(f"relaxation rate: {lbm_config.relaxation_rate}")
    
    # Generate LBM bulk sweeps
    with sfg.namespace("Couette"):
        # lbm_bulk = LbmBulk(sfg, "LBM", lbm_config)
        # sfg.generate(lbm_bulk)
        # rho, u = lbm_bulk.rho, lbm_bulk.u

        lbm_opt = LBMOptimisation(symbolic_field=f_pdfs, 
                                  symbolic_temporary_field=f_pdfs_tmp, 
                                  field_layout='fzyx',)
        
        output_fields = dict({"density": f_density, "velocity": f_velocity})
        collision_rule = create_lb_update_rule(lbm_config=lbm_config, 
                                               lbm_optimisation=lbm_opt, 
                                               output=output_fields)
        
        combined_assignments = ps.AssignmentCollection(main_assignments=collision_rule.main_assignments + [ps.Assignment(f_viscosity.center(), nu)], 
                                                       subexpressions=collision_rule.subexpressions + [ps.Assignment(nu, nu_expr)])

        stream_collide = Sweep("StreamCollide", combined_assignments)
        stream_collide.swap_fields(f_pdfs, f_pdfs_tmp)
        sfg.generate(stream_collide)

        # lb_method = collision_rule.method
        lb_method = create_lb_method(lbm_config)
        init_rule = macroscopic_values_setter(
            lb_method=lb_method,
            density=f_density.center,
            velocity=f_velocity.center_vector,
            pdfs=f_pdfs,
            set_pre_collision_pdfs=True,
        )
        sfg.generate(Sweep("InitPdfs", init_rule))

        # Analytical solution for Couette flow: u(z) = u_max * z / H (linear velocity profile)
        couette_analytical = [
            ps.Assignment(f_density(), 1),
            ps.Assignment(f_velocity(0), u_max * cell.z() / domain.z_max()),
            ps.Assignment(f_velocity(1), 0),
            ps.Assignment(f_velocity(2), 0),
        ]
        sfg.generate(Sweep("SetAnalyticalSolution", couette_analytical))

        # Temperature initialization - linear gradient from bottom to top
        # T(z) = T_bottom + (T_top - T_bottom) * z / H
        temperature_init = [
            ps.Assignment(f_temperature.center(), T_bottom + (T_top - T_bottom) * cell.z() / domain.z_max())
        ]
        sfg.generate(Sweep("InitializeTemperature", temperature_init))

        # Error computation: L-infinity norm
        ux = sp.Symbol("ux")
        error_ux = ps.TypedSymbol("error_ux", ps.DynamicType.NUMERIC_TYPE)

        error_calc = [
            ps.Assignment(ux, u_max * cell.z() / domain.z_max()),
            ps.MaxReductionAssignment(error_ux, sp.Abs(f_velocity(0) - ux))
        ]
        sfg.generate(Sweep("VelocityErrorLmax", error_calc))

    # Boundary conditions
    noSlip = GenericHBB(NoSlip(name="NoSlip"), lb_method, f_pdfs)
    sfg.generate(noSlip)

    wall_velocity = (u_max, 0, 0)
    ubb = GenericHBB(UBB(wall_velocity, name="UBB"), lb_method, f_pdfs)
    sfg.generate(ubb)
    
print("Code generation complete!")
