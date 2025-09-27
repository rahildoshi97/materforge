#!/usr/bin/env python3
"""
CouetteFlowKernel.py - 3D Thermal Couette Flow
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
use_material_library = False  # Enable materforge integration

print("=== COUETTE FLOW CODE GENERATION DEBUG LOG ===")
print(f"Starting code generation at {Path(__file__).resolve()}")
print(f"Material library integration: {'ENABLED' if use_material_library else 'DISABLED'}")

with SourceFileGenerator() as sfg:
    data_type = "float64"
    
    # Field definitions
    stencil = Stencil.D3Q19
    dim = 3
    Q = 19
    
    print(f"LBM Configuration: Stencil={stencil}, Dimensions={dim}, Q={Q}, Data type={data_type}")
    
    # Define all fields
    f_pdfs, f_pdfs_tmp, f_density, f_temperature, f_viscosity, f_force, f_velocity = ps.fields(
        f"f_pdfs({Q}), f_pdfs_tmp({Q}), f_density, f_temperature, f_viscosity, f_force({dim}), f_velocity({dim}): {data_type}[{dim}D]", 
        layout='fzyx'
    )
    
    print("DEBUG: Field definitions created successfully")
    print(f"  - PDFs: {f_pdfs.name} ({Q} components)")
    print(f"  - Velocity: {f_velocity.name} ({dim} components)")
    print(f"  - Force: {f_force.name} ({dim} components)")
    print(f"  - Scalars: density, temperature, viscosity (1 component each)")
    
    # Create symbol creator for parameters
    s = SymbolCreator()
    
    # ===== MATERFORGE INTEGRATION =====
    material_assignments = []
    viscosity_expression = None
    
    print("\n=== MATERFORGE INTEGRATION DEBUG ===")
    
    if use_material_library:
        try:
            from materforge.parsing.api import create_material
            
            yaml_path = Path(__file__).parent / 'CouetteFlowMaterial_Updated.yaml'
            print(f"DEBUG: Looking for YAML at: {yaml_path.resolve()}")
            
            if yaml_path.exists():
                print(f"Loading material from: {yaml_path}")
                
                mat = create_material(yaml_path=yaml_path, dependency=f_temperature.center(), enable_plotting=False)
                print(f"DEBUG: Material has dynamic_viscosity: {hasattr(mat, 'dynamic_viscosity')}")
                
                if hasattr(mat, 'dynamic_viscosity'):
                    print(f"DEBUG: Dynamic viscosity expression: {mat.dynamic_viscosity}")
                
                material_assignments.extend([
                    ps.Assignment(s.density_mat, 1.0),  # Material density
                    ps.Assignment(s.viscosity_mat, mat.dynamic_viscosity), 
                ])
                
                viscosity_expression = s.viscosity_mat
                print("DEBUG: Successfully integrated materforge properties with temperature dependency")
                print(f"DEBUG: Viscosity expression symbol: {viscosity_expression}")
                
            else:
                print(f"ERROR: Material file not found at {yaml_path}")
                print("DEBUG: Falling back to constant viscosity model")
                use_material_library = False
                
        except ImportError as e:
            print(f"WARNING: materforge not available ({e}), using fallback model")
            use_material_library = False
        except Exception as e:
            print(f"ERROR: Error loading material: {e}")
            print(f"ERROR: Exception type: {type(e).__name__}")
            import traceback
            traceback.print_exc()
            use_material_library = False
    
    # Fallback to simple constant model
    if not use_material_library or viscosity_expression is None:
        print("DEBUG: Using constant viscosity and density model")
        print("DEBUG: This means materforge integration failed or is disabled")
        
        material_assignments.extend([
            ps.Assignment(s.density_mat, 1.0),  # LBM normalized density
            ps.Assignment(s.viscosity_mat, 0.1667),
        ])
        
        viscosity_expression = s.viscosity_mat
        print(f"DEBUG: Constant viscosity value: 0.1667")
        print(f"DEBUG: Fallback viscosity expression: {viscosity_expression}")
    
    print(f"DEBUG: Final material assignments count: {len(material_assignments)}")
    for i, assignment in enumerate(material_assignments):
        print(f"  Assignment {i}: {assignment}")
    
    # ===== LBM CONFIGURATION DEBUG =====
    print("\n=== LBM CONFIGURATION DEBUG ===")
    
    # Calculate relaxation parameters - temperature dependent
    relaxation_time = 3.0 * viscosity_expression + 0.5
    omega = 1.0 / relaxation_time
    
    print(f"DEBUG: Relaxation time expression: {relaxation_time}")
    print(f"DEBUG: Omega expression: {omega}")
    print(f"DEBUG: For viscosity=0.1667: tau={3.0*0.1667+0.5}={3.0*0.1667+0.5:.3f}, omega={1.0/(3.0*0.1667+0.5):.3f}")
    
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
    
    print("DEBUG: LBM Configuration created successfully")
    print(f"  - Method: {lbm_config.method}")
    print(f"  - Stencil: {lbm_config.stencil}")
    print(f"  - Force model: {lbm_config.force_model}")
    print(f"  - Compressible: {lbm_config.compressible}")
    print(f"  - Zero centered: {lbm_config.zero_centered}")
    print(f"  - Output fields: {lbm_config.output}")
    
    # Create LBM method
    lb_method = create_lb_method(lbm_config)
    print(f"DEBUG: LBM method created with type: {type(lb_method)}")
    
    lbm_opt = LBMOptimisation(symbolic_field=f_pdfs, symbolic_temporary_field=f_pdfs_tmp)
    print("DEBUG: LBM optimization settings created")
    
    # Create LBM collision and streaming
    collision_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    print("DEBUG: LBM collision rule created successfully")
    print(f"DEBUG: Main assignments count: {len(collision_rule.main_assignments)}")
    print(f"DEBUG: Subexpressions count: {len(collision_rule.subexpressions)}")
    
    # ===== MAIN SWEEP GENERATION DEBUG =====
    print("\n=== MAIN SWEEP GENERATION DEBUG ===")
    
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
    
    print("DEBUG: CouetteFlowSweep generated successfully")
    print(f"DEBUG: Sweep class name: CouetteFlowSweep")
    print(f"DEBUG: PDF field swapping enabled: {f_pdfs.name} <-> {f_pdfs_tmp.name}")
    
    # ===== INITIALIZATION SWEEP DEBUG =====
    print("\n=== INITIALIZATION SWEEP DEBUG ===")
    
    from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
    
    # Create initialization kernel
    init_kernel = macroscopic_values_setter(
        lb_method,
        density=f_density,
        velocity=f_velocity,
        pdfs=f_pdfs,
        set_pre_collision_pdfs=True
    )
    
    print("DEBUG: Initialization kernel created successfully")
    print(f"DEBUG: Pre-collision PDFs: True")
    
    # Generate initialization sweep
    init_sweep = Sweep("CouetteFlowInit", init_kernel)
    sfg.generate(init_sweep)
    
    print("DEBUG: CouetteFlowInit generated successfully")
    
    # ===== BOUNDARY CONDITIONS DEBUG =====
    print("\n=== BOUNDARY CONDITIONS DEBUG ===")
    
    # Top wall: moving wall with prescribed velocity
    top_wall_bc = ps.AssignmentCollection([
        ps.Assignment(f_velocity.center(0), s.s_wall_velocity),
        ps.Assignment(f_velocity.center(1), 0.0),
        ps.Assignment(f_velocity.center(2), 0.0),
    ])
    
    top_wall_sweep = Sweep("TopWallBC", top_wall_bc)
    sfg.generate(top_wall_sweep)
    print("DEBUG: TopWallBC generated - moving wall boundary condition")
    
    # Bottom wall: stationary wall
    bottom_wall_bc = ps.AssignmentCollection([
        ps.Assignment(f_velocity.center(0), 0.0),
        ps.Assignment(f_velocity.center(1), 0.0),
        ps.Assignment(f_velocity.center(2), 0.0),
    ])
    
    bottom_wall_sweep = Sweep("BottomWallBC", bottom_wall_bc)
    sfg.generate(bottom_wall_sweep)
    print("DEBUG: BottomWallBC generated - stationary wall boundary condition")
    
    # ===== TEMPERATURE FIELD UPDATE DEBUG =====
    print("\n=== TEMPERATURE FIELD UPDATE DEBUG ===")
    
    # Keep temperature profile linear and constant throughout simulation
    temp_update_kernel = ps.AssignmentCollection([
        ps.Assignment(f_temperature.center(), 
                     s.T_bottom + (s.T_top - s.T_bottom) * s.y_normalized)
    ])
    
    temp_update_sweep = Sweep("TemperatureUpdate", temp_update_kernel)
    sfg.generate(temp_update_sweep)
    
    print("Generated TemperatureUpdate sweep for linear profile")
    
    print("Code generation completed successfully!")
