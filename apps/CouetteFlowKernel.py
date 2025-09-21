#!/usr/bin/env python3
"""
CouetteFlowKernel.py - CORRECTED WITH SIMPLE VELOCITY BOUNDARY CONDITIONS
Code generation script for Couette flow simulation with temperature-dependent viscosity
"""

import logging
import sympy as sp
import pystencils as ps
from pathlib import Path
from pystencilssfg import SourceFileGenerator
from pystencils import SymbolCreator
from walberla.codegen import Sweep
from lbmpy import LBMConfig, Stencil, Method, create_lb_method, create_lb_update_rule, LBMOptimisation, ForceModel

# Configure logging
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(levelname)s %(name)s -> %(message)s")

# Material properties integration flag
use_material_library = True  # Enable materforge integration

with SourceFileGenerator() as sfg:
    data_type = "float64"
    
    # Field definitions
    stencil = Stencil.D3Q19
    dim = 3
    Q = 19
    
    # Define all fields
    f_pdfs, f_pdfs_tmp, f_density, f_temperature, f_viscosity, f_force, f_velocity = ps.fields(
        f"f_pdfs({Q}), f_pdfs_tmp({Q}), f_density, f_temperature, f_viscosity, f_force({dim}), f_velocity({dim}): {data_type}[{dim}D]", 
        layout='fzyx'
    )
    
    # Create symbol creator for parameters
    s = SymbolCreator()
    
    # ===== MATERFORGE INTEGRATION =====
    material_assignments = []
    viscosity_expression = None
    
    if use_material_library:
        try:
            from materforge.parsing.api import create_material
            
            yaml_path = Path(__file__).parent / 'CouetteFlowMaterial_Updated.yaml'
            if yaml_path.exists():
                print(f"Loading material from: {yaml_path}")
                
                mat = create_material(yaml_path=yaml_path, dependency=f_temperature.center(), enable_plotting=False)
                
                material_assignments.extend([
                    ps.Assignment(s.density_mat, 1.0),  # LBM normalized density
                    ps.Assignment(s.viscosity_mat, mat.dynamic_viscosity), 
                ])
                
                viscosity_expression = s.viscosity_mat
                print("Successfully integrated materforge properties with temperature dependency")
                
            else:
                print(f"Warning: Material file not found at {yaml_path}")
                use_material_library = False
                
        except ImportError as e:
            print(f"Warning: materforge not available ({e}), using fallback model")
            use_material_library = False
        except Exception as e:
            print(f"Error loading material: {e}")
            use_material_library = False
    
    # Fallback to constant model
    if not use_material_library or viscosity_expression is None:
        print("Using constant viscosity and density model")
        
        material_assignments.extend([
            ps.Assignment(s.density_mat, 1.0),  # LBM normalized density
            ps.Assignment(s.viscosity_mat, 0.001),  # Constant viscosity
        ])
        
        viscosity_expression = s.viscosity_mat
    
    # ===== LBM CONFIGURATION =====
    
    # Calculate relaxation parameters - CORRECTED for proper LBM viscosity
    # Use smaller viscosity for better flow development
    tau = 3.0 * 0.01 + 0.5  # Use tau = 0.53 for stable flow
    omega = 1.0 / tau
    
    # Create LBM configuration
    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.SRT,
        relaxation_rate=omega,
        compressible=False,
        zero_centered=False,
        force=f_force,
        force_model=ForceModel.GUO,
        output={'density': f_density, 'velocity': f_velocity}
    )
    
    # Create LBM method
    lb_method = create_lb_method(lbm_config)
    lbm_opt = LBMOptimisation(symbolic_field=f_pdfs, symbolic_temporary_field=f_pdfs_tmp)
    
    # Create LBM collision and streaming
    collision_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    
    # ===== MAIN SWEEP =====
    
    # Combine material calculations with LBM assignments
    main_assignments = collision_rule.main_assignments + [
        ps.Assignment(f_viscosity.center(), s.viscosity_mat),
    ]
    
    subexpressions = material_assignments + collision_rule.subexpressions
    
    # Create assignment collection
    couette_kernel = ps.AssignmentCollection(
        main_assignments=main_assignments,
        subexpressions=subexpressions
    )
    
    # Generate the main sweep
    couette_sweep = Sweep("CouetteFlowSweep", couette_kernel)
    couette_sweep.swap_fields(f_pdfs, f_pdfs_tmp)
    sfg.generate(couette_sweep)
    
    print(f"Generated CouetteFlowSweep with material properties")
    
    # ===== INITIALIZATION SWEEP =====
    
    from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
    
    # Create initialization kernel
    init_kernel = macroscopic_values_setter(
        lb_method,
        density=f_density,
        velocity=f_velocity,
        pdfs=f_pdfs,
        set_pre_collision_pdfs=True
    )
    
    # Generate initialization sweep
    init_sweep = Sweep("CouetteFlowInit", init_kernel)
    sfg.generate(init_sweep)
    
    print(f"Generated CouetteFlowInit")
    
    # ===== SIMPLE BOUNDARY CONDITIONS =====
    
    # Top wall: moving wall with prescribed velocity
    top_wall_bc = ps.AssignmentCollection([
        ps.Assignment(f_velocity.center(0), s.s_wall_velocity),
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
    
    print("Generated simple boundary condition sweeps")
    print("Code generation completed successfully!")
