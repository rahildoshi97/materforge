# Converting Between Energy Density and Temperature

This guide explains how to perform bilateral conversions between energy density and temperature in MaterForge using the current architecture.

## Why Energy-Temperature Conversion Matters

In many material simulations, particularly those involving heat transfer and phase changes, you need to convert between:

- **Temperature**: The conventional measure of thermal state
- **Energy Density**: The amount of thermal energy per unit volume

This conversion is essential because:
- Energy is conserved in physical processes
- Phase transitions involve latent heat where temperature remains constant
- Many numerical methods work better with energy as the primary variable

## Current Implementation Status

**Important Note**: The energy-temperature conversion functionality is currently under development in MaterForge v0.3.0. The following features are planned but not yet implemented:

- InterpolationArrayContainer class
- Automatic C++ code generation for interpolation
- Binary search and double lookup methods

## Basic Temperature-Energy Conversion (Current)

### 1. Create Material with Energy Density
```python
import sympy as sp
from materforge.parsing.api import create_material

# Create symbolic temperature variable
T = sp.Symbol('T')

# Load material with energy density property
material = create_material("1.4301.yaml", T)

# Check if energy density is available
if hasattr(material, 'energy_density'):
    print("Energy density property is available")
else:
    print("Energy density property not defined")
```

### 2. Temperature to Energy Conversion
```python
# Temperature to energy conversion using symbolic evaluation

temperature = 1500.0 # Kelvin
if hasattr(material, 'energy_density'):
    energy = material.energy_density.evalf(T, temperature)
    print(f"Energy density at {temperature}K: {energy} J/m³")
```

### 3. Energy to Temperature Conversion

The inverse conversion (energy to temperature) will be available through:
```python
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

# Create inverse function (when implemented)
if hasattr(material, 'energy_density'):
    E = sp.Symbol('E')
    inverse_func = PiecewiseInverter.create_energy_density_inverse(material, 'E')
```
# Use inverse function
```python
energy_value = 1.5e9  # J/m³
temperature = float(inverse_func.subs(E, energy_value))
```
