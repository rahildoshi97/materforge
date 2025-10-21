#!/usr/bin/env python3

"""
CouetteFlowSweepsGPU.py - 3D Thermal Couette Flow with integrated analytical solution
"""

import logging
from pathlib import Path
from argparse import ArgumentParser

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
import numpy as np
from scipy.integrate import cumulative_trapezoid

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep, get_build_config
from sweepgen.boundaries import GenericBoundary
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

    # Symbolic parameters
    nu = sp.Symbol("nu")
    u_max = sp.Symbol("u_max")
    T_top, T_bottom = sp.symbols("T_top, T_bottom")

    # Define fields
    f_pdfs, f_pdfs_tmp, f_density, f_temperature, f_viscosity, f_velocity = ps.fields(
        f"f_pdfs({Q}), f_pdfs_tmp({Q}), f_density, f_temperature, f_viscosity, f_velocity({D}): {data_type}[{D}D]",
        layout='fzyx'
    )
    
    const_nu = 1./6.  # Constant viscosity if MaterForge is not used

    analytical_velocity_expr = None
    
    if use_materforge:
        print("Code generation: Using temperature-dependent viscosity from MaterForge")
        try:
            from materforge.parsing.api import create_material
            
            yaml_path = Path(__file__).parent / 'CouetteFlowMaterial_Updated.yaml'
            
            if yaml_path.exists():
                # Create MaterForge material with temperature field dependency
                mat = create_material(yaml_path=yaml_path, dependency=f_temperature.center(), enable_plotting=False)
                
                if hasattr(mat, 'dynamic_viscosity'):
                    nu_expr = mat.dynamic_viscosity
                    print(f"MaterForge viscosity expression loaded: {nu_expr}")
                    # ========================================================
                    # COMPUTE ANALYTICAL SOLUTION
                    # ========================================================
                    print("Computing analytical solution for error calculation...")
                    
                    # Extract the symbolic temperature variable from nu_expr
                    temp_symbols = [s for s in nu_expr.free_symbols if 'temperature' in str(s)]
                    if temp_symbols:
                        temp_var = temp_symbols[0]
                        print(f"Found temperature symbol: {temp_var}")
                        
                        # Create a simple evaluation symbol
                        T_eval = sp.Symbol('T_eval')
                        nu_expr_eval = nu_expr.subs(temp_var, T_eval)
                        
                        # Lambdify for numerical evaluation
                        nu_func = sp.lambdify(T_eval, nu_expr_eval, modules='numpy')
                        
                        # High-resolution grid for analytical solution
                        z_analytical = np.linspace(0, 1, 1000)
                        
                        # Temperature profile (300K to 600K)
                        T_analytical = 300.0 + 300.0 * z_analytical
                        
                        # Evaluate viscosity
                        nu_analytical = nu_func(T_analytical)
                        
                        # Numerical integration
                        integrand = 1.0 / nu_analytical
                        integral_z = cumulative_trapezoid(integrand, z_analytical, initial=0)
                        integral_total = integral_z[-1]
                        u_normalized = integral_z / integral_total
                        
                        # Polynomial fit (degree 5)
                        from numpy.polynomial import polynomial as P
                        coeffs = P.polyfit(z_analytical, u_normalized, 5)
                        
                        # Create symbolic polynomial
                        z_norm = cell.z() / domain.z_max()
                        analytical_velocity_expr = u_max * sum(
                            float(coeffs[i]) * z_norm**i for i in range(len(coeffs))
                        ) # type: ignore
                        
                        fit_error = np.max(np.abs(u_normalized - P.polyval(z_analytical, coeffs)))
                        print(f"Analytical solution: polynomial fit max error = {fit_error:.6e}")
                        
                    else:
                        print("Warning: Could not extract temperature symbol, using linear solution")
                        nu_expr = const_nu
                        analytical_velocity_expr = u_max * cell.z() / domain.z_max()
                else:
                    print("Warning: No dynamic_viscosity attribute, using constant")
                    nu_expr = const_nu
                    analytical_velocity_expr = u_max * cell.z() / domain.z_max()
            else:
                print(f"Warning: Material file not found")
                nu_expr = const_nu
                analytical_velocity_expr = u_max * cell.z() / domain.z_max()
                
        except Exception as e:
            print(f"MaterForge error: {e}")
            import traceback
            traceback.print_exc()
            nu_expr = const_nu
            analytical_velocity_expr = u_max * cell.z() / domain.z_max()
    else:
        print("Code generation: Using constant viscosity")
        nu_expr = const_nu
        analytical_velocity_expr = u_max * cell.z() / domain.z_max()

    # LBM configuration
    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.TRT,
        compressible=True,
        relaxation_rate=relaxation_rate_from_lattice_viscosity(nu_expr),
    )
    print(f"relaxation rate: {lbm_config.relaxation_rate}")
    
    # Generate LBM sweeps
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
        
        combined_assignments = ps.AssignmentCollection(main_assignments=collision_rule.main_assignments + [ps.Assignment(f_viscosity.center(), nu_expr)], 
                                                       subexpressions=collision_rule.subexpressions)

        stream_collide = Sweep("StreamCollide", combined_assignments)
        stream_collide.swap_fields(f_pdfs, f_pdfs_tmp)
        sfg.generate(stream_collide)

        # Initialize PDFs
        lb_method = create_lb_method(lbm_config)
        init_rule = macroscopic_values_setter(
            lb_method=lb_method,
            density=f_density.center,
            velocity=f_velocity.center_vector,
            pdfs=f_pdfs,
            set_pre_collision_pdfs=True,
        )
        sfg.generate(Sweep("InitPdfs", init_rule))

        # Set analytical solution
        couette_analytical = [
            ps.Assignment(f_density(), 1),
            ps.Assignment(f_velocity(0), analytical_velocity_expr),
            ps.Assignment(f_velocity(1), 0),
            ps.Assignment(f_velocity(2), 0),
        ]
        sfg.generate(Sweep("SetAnalyticalSolution", couette_analytical))

        # Temperature initialization
        temperature_init = [
            ps.Assignment(
                f_temperature.center(), 
                T_bottom + (T_top - T_bottom) * cell.z() / domain.z_max()
            )
        ]
        sfg.generate(Sweep("InitializeTemperature", temperature_init))

        # Error computation with analytical solution
        # ux_analytical = sp.Symbol("ux_analytical")
        # error_ux = ps.TypedSymbol("error_ux", ps.DynamicType.NUMERIC_TYPE)

        """error_calc = [
            ps.Assignment(ux_analytical, analytical_velocity_expr),
            ps.MaxReductionAssignment(error_ux, sp.Abs(f_velocity(0) - ux_analytical))
        ]"""
        # sfg.generate(Sweep("VelocityErrorLmax", error_calc))

    # Boundary conditions
    noSlip = GenericBoundary(NoSlip(name="NoSlip"), lb_method, f_pdfs)
    sfg.generate(noSlip)

    wall_velocity = (u_max, 0, 0)
    ubb = GenericBoundary(UBB(wall_velocity, name="UBB"), lb_method, f_pdfs)
    sfg.generate(ubb)
    
print("Code generation complete!")
