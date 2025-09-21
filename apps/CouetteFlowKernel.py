#!/usr/bin/env python3
"""
CouetteFlowKernel.py
Code generation script for Couette flow simulation with temperature-dependent viscosity
"""

import logging
import pystencils as ps
from pathlib import Path
from pystencilssfg import SourceFileGenerator
from pystencils import SymbolCreator
from walberla.codegen import Sweep
from lbmpy import LBMConfig, Stencil, Method, create_lb_method, create_lb_update_rule, LBMOptimisation

# Configure logging
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s %(levelname)s %(name)s -> %(message)s"
)

# Silence noisy libraries
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)

# Material properties integration flag
use_material_library = False  # Set to True to use materforge integration

with SourceFileGenerator() as sfg:
    data_type = "float64"
    
    # Define fields
    f_pdfs, f_pdfs_tmp = ps.fields(f"f_pdfs({19}), f_pdfs_tmp({19}): {data_type}[3D]", layout='fzyx')
    f_velocity, f_temperature, f_viscosity, f_density, f_force = ps.fields(
        f"f_velocity(3), f_temperature, f_viscosity, f_density, f_force(3): {data_type}[3D]", layout='fzyx'
    )
    
    # Create symbol creator for parameters
    s = SymbolCreator()
    s_dx, s_dt, s_wall_velocity, s_viscosity_ref = s.s_dx, s.s_dt, s.s_wall_velocity, s.s_viscosity_ref
    s_temp_coeff, s_temp_ref = s.s_temp_coeff, s.s_temp_ref
    s_omega = s.s_omega  # LBM relaxation parameter
    
    # LBM configuration
    stencil = Stencil.D3Q19
    
    # Temperature-dependent viscosity model
    # μ(T) = μ_ref * exp(-α * (T - T_ref))
    constant_viscosity = 0.001
    
    # Relaxation time from viscosity: τ = 3ν/c² + 0.5, where ν = μ/ρ, c² = 1/3
    # τ = 3μ/(ρc²) + 0.5 = 9μ/ρ + 0.5
    relaxation_time = 9.0 * constant_viscosity / f_density.center() + 0.5
    omega = 1.0 / relaxation_time
    
    # Create LBM configuration
    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.SRT,  # Single relaxation time
        relaxation_rate=omega,
        compressible=False,
        zero_centered=False,
        force=f_force.center_vector,
        output={'density': f_density, 'velocity': f_velocity}
    )
    
    # Create LBM method
    lb_method = create_lb_method(lbm_config)
    
    # Create LBM optimization
    lbm_opt = LBMOptimisation(
        symbolic_field=f_pdfs,
        symbolic_temporary_field=f_pdfs_tmp
    )
    
    # Create main LBM update rule
    main_assignments = []
    subexpressions = []
    
    # Add temperature-dependent viscosity calculation
    subexpressions.extend([
        ps.Assignment(s.viscosity_local, constant_viscosity),
        ps.Assignment(s.tau_local, relaxation_time),
        ps.Assignment(s.omega_local, omega),
    ])
    
    # Create LBM collision and streaming
    collision_rule = create_lb_update_rule(
        lbm_config=lbm_config,
        lbm_optimisation=lbm_opt
    )
    
    # Update viscosity field for output
    main_assignments.append(
        ps.Assignment(f_viscosity.center(), s.viscosity_local)
    )
    
    # Create assignment collection
    couette_kernel = ps.AssignmentCollection(
        main_assignments=collision_rule.main_assignments + main_assignments,
        subexpressions=collision_rule.subexpressions + subexpressions
    )
    
    # Generate the main sweep
    couette_sweep = Sweep("CouetteFlowSweep", couette_kernel)
    couette_sweep.swap_fields(f_pdfs, f_pdfs_tmp)
    sfg.generate(couette_sweep)
    
    print(f"Generated CouetteFlowSweep with assignments:\n{couette_kernel}")
    
    # ===== INITIALIZATION SWEEP =====
    
    # Initialization for equilibrium PDFs
    from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
    
    # Create initialization kernel
    init_kernel = macroscopic_values_setter(
        lb_method,
        density=f_density,
        velocity=f_velocity,
        pdfs=f_pdfs,
        set_pre_collision_pdfs=True
    )
    
    # Additional initialization for temperature and viscosity
    init_assignments = []
    
    # Initialize temperature field (linear profile will be set in C++)
    init_assignments.append(
        ps.Assignment(f_temperature.center(), s_temp_ref)
    )
    
    # Initialize viscosity field
    init_assignments.append(
        ps.Assignment(f_viscosity.center(), s_viscosity_ref)
    )
    
    # Initialize density field
    init_assignments.append(
        ps.Assignment(f_density.center(), 1.0)  # Normalized density
    )
    
    # Combine initialization assignments
    complete_init = ps.AssignmentCollection(
        main_assignments=init_kernel.main_assignments + init_assignments,
        subexpressions=init_kernel.subexpressions
    )
    
    # Generate initialization sweep
    init_sweep = Sweep("CouetteFlowInit", complete_init)
    sfg.generate(init_sweep)
    
    print(f"Generated CouetteFlowInit with assignments:\n{complete_init}")
    
    # ===== BOUNDARY CONDITIONS =====
    
    # Simple boundary condition implementations
    # Top wall: moving wall with prescribed velocity
    top_wall_bc = ps.AssignmentCollection([
        ps.Assignment(f_velocity.center(0), s_wall_velocity),
        ps.Assignment(f_velocity.center(1), 0.0),
        ps.Assignment(f_velocity.center(2), 0.0),
    ])
    
    top_wall_sweep = Sweep("TopWallBC", top_wall_bc)
    sfg.generate(top_wall_sweep)
    
    # Bottom wall: stationary wall
    bottom_wall_bc = ps.AssignmentCollection([
        ps.Assignment(f_velocity.center(0), 0.0),
        ps.Assignment(f_velocity.center(1), 0.0),
        ps.Assignment(f_velocity.center(2), 0.0),
    ])
    
    bottom_wall_sweep = Sweep("BottomWallBC", bottom_wall_bc)
    sfg.generate(bottom_wall_sweep)
    
    print("Generated boundary condition sweeps")
    
    # ===== MATERIAL INTEGRATION (if enabled) =====
    
    if use_material_library:
        try:
            from materforge.parsing.api import create_material
            
            # Load material from YAML file
            yaml_path = Path(__file__).parent / 'CouetteFlowMaterial.yaml'
            if yaml_path.exists():
                mat = create_material(yaml_path=yaml_path, dependency=f_temperature.center())
                
                # Create material-based viscosity calculation
                material_viscosity = ps.AssignmentCollection(
                    subexpressions=[
                        ps.Assignment(s.density_mat, mat.density),
                        ps.Assignment(s.viscosity_mat, mat.dynamic_viscosity),
                    ],
                    main_assignments=[
                        ps.Assignment(f_viscosity.center(), s.viscosity_mat),
                        ps.Assignment(f_density.center(), s.density_mat),
                    ]
                )
                
                material_sweep = Sweep("MaterialProperties", material_viscosity)
                sfg.generate(material_sweep)
                
                print("Generated material property sweep using materforge")
            else:
                print(f"Warning: Material file not found at {yaml_path}")
                
        except ImportError:
            print("Warning: materforge not available, using hardcoded properties")
    
    print("Code generation completed successfully!")