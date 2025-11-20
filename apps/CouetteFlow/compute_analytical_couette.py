#!/usr/bin/env python3
"""
compute_analytical_couette.py - Compute analytical solution for variable viscosity Couette flow
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import sympy as sp
from materforge.parsing.api import create_material
from pathlib import Path

# Parameters
T_bottom = 300.0  # K
T_top = 600.0     # K
u_max = 0.025
H = 32  # Number of cells in z-direction

# Temperature profile (linear): Create temperature array
z = np.linspace(0, H, 1000)  # 1000 points from bottom to top
T = T_bottom + (T_top - T_bottom) * z / H  # Linear temperature

# Load viscosity model from MaterForge
print("Loading viscosity model from MaterForge...")
yaml_path = Path(__file__).parent / 'CouetteFlowMaterial_Updated.yaml'

if yaml_path.exists():
    # Create symbolic variable for temperature
    T_sym = sp.Symbol('T')
    
    # Load material properties
    mat = create_material(yaml_path=yaml_path, dependency=T_sym, enable_plotting=False)
    
    if hasattr(mat, 'dynamic_viscosity'):
        nu_expr = mat.dynamic_viscosity
        print(f"Symbolic viscosity expression:\n{nu_expr}\n")
        
        # Convert symbolic expression to numerical function using lambdify
        # This allows evaluation on NumPy arrays
        nu_func = sp.lambdify(T_sym, nu_expr, modules='numpy')
    else:
        print("Warning: No dynamic_viscosity attribute found, using constant viscosity")
        nu_func = lambda T: 1.0/6.0
else:
    print(f"Warning: Material file {yaml_path} not found, using constant viscosity")
    nu_func = lambda T: 1.0/6.0

# Evaluate viscosity at all temperature points
nu_z = nu_func(T)

print(f"Viscosity evaluation:")
print(f"  At T={T_bottom}K: nu = {nu_func(T_bottom):.8f}")
print(f"  At T={T_top}K: nu = {nu_func(T_top):.8f}")
print(f"  Viscosity variation: {(nu_func(T_bottom)/nu_func(T_top) - 1)*100:.2f}%\n")

# Compute velocity profile: u(z) = U_wall * int_0^z (1/nu(z')) dz' / int_0^H (1/nu(z')) dz'
integrand = 1.0 / nu_z  # Compute integrand (1/viscosity)
# Numerical integration (cumulative trapezoid)
integral_z = cumulative_trapezoid(integrand, z, initial=0)
integral_H = integral_z[-1]  # Total integral from 0 to H

# Velocity profile
u_analytical = u_max * integral_z / integral_H

# Linear profile (constant viscosity for comparison)
u_linear = u_max * z / H

# Compute deviation
max_deviation = np.max(np.abs(u_analytical - u_linear))
relative_deviation = max_deviation / u_max * 100

print(f"Analytical solution statistics:")
print(f"  Maximum deviation from linear: {max_deviation:.6e}")
print(f"  Relative deviation: {relative_deviation:.2f}%")

# Plot comparison
plt.figure(figsize=(12, 5))

# Subplot 1: Velocity profiles
plt.subplot(1, 2, 1)
plt.plot(u_analytical, z, 'b-', linewidth=2, label='Variable viscosity (analytical)')
plt.plot(u_linear, z, 'r--', linewidth=2, label='Constant viscosity (linear)')
plt.xlabel('Velocity u(z)', fontsize=12)
plt.ylabel('Height z [cells]', fontsize=12)
plt.title('Couette Flow: Velocity Profiles', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)

# Subplot 2: Viscosity profile
plt.subplot(1, 2, 2)
plt.plot(nu_z, z, 'g-', linewidth=2)
plt.xlabel('Viscosity ν(z)', fontsize=12)
plt.ylabel('Height z [cells]', fontsize=12)
plt.title('Temperature-Dependent Viscosity', fontsize=14)
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('couette_analytical_comparison.png', dpi=300, bbox_inches='tight')
print(f"\nPlot saved to 'couette_analytical_comparison.png'")
plt.show()

# Plot deviation
plt.figure(figsize=(8, 6))
deviation = u_analytical - u_linear
plt.plot(deviation * 1000, z, 'k-', linewidth=2)  # Convert to mm/s
plt.xlabel('Velocity Deviation (×10⁻³)', fontsize=12)
plt.ylabel('Height z [cells]', fontsize=12)
plt.title('Deviation from Linear Profile', fontsize=14)
plt.grid(True, alpha=0.3)
plt.axvline(x=0, color='r', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('couette_deviation.png', dpi=300, bbox_inches='tight')
print(f"Deviation plot saved to 'couette_deviation.png'")

# Save analytical solution for comparison with simulation
output_data = np.column_stack([z, u_analytical, u_linear, nu_z, T])
np.savetxt('couette_analytical.dat', 
           output_data,
           header='z[cells]  u_variable[LBM]  u_linear[LBM]  nu[LBM]  T[K]',
           fmt='%.10e')
print(f"Analytical solution saved to 'couette_analytical.dat'\n")

# Summary
print("="*60)
print("SUMMARY")
print("="*60)
print(f"Domain height: {H} cells")
print(f"Temperature range: {T_bottom}K to {T_top}K")
print(f"Wall velocity: {u_max}")
print(f"Viscosity range: {nu_func(T_bottom):.6f} to {nu_func(T_top):.6f}")
print(f"Maximum velocity deviation: {max_deviation:.6e} ({relative_deviation:.2f}%)")
print(f"\nConclusion: {'Small' if relative_deviation < 5 else 'Significant'} effect of variable viscosity")
print("="*60)
