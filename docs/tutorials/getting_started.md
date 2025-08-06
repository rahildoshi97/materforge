# Getting Started with materforge

This tutorial will guide you through the basics of using materforge to model material properties for simulations.

## Prerequisites

Before starting, ensure you have:

- Python 3.10 or newer
- Basic knowledge of material properties
- Basic familiarity with Python

## Installation

Install materforge using pip:

```bash
pip install "git+https://i10git.cs.fau.de/rahil.doshi/materforge.git"
```


For development, clone the repository and install in development mode:

```bash
git clone https://i10git.cs.fau.de/rahil.doshi/materforge.git
cd materforge
pip install -e .
```


## Basic Concepts

MaterForge organizes material data around these key concepts:

1. **Materials**: Materials with variable-dependent properties
2. **Property Types**: Six different ways to define material properties
3. **Symbolic Processing**: SymPy-based symbolic mathematics for property relationships
4. **Dependencies**: Automatic handling of variable-dependent evaluations

## Your First Material

Let's create a simple steel alloy definition:

1. Create a file named `simple_steel.yaml`:
```yaml
name: SimpleSteel

composition:
    Fe: 0.98
    C: 0.02

solidus_temperature: 1450
liquidus_temperature: 1520
initial_boiling_temperature: 3000
final_boiling_temperature: 3100

properties:
    density:
        dependency: [300, 800, 1300, 1800]
        value: [7850, 7800, 7750, 7700]
        bounds: [constant, constant]

    heat_conductivity:
        dependency: [300, 800, 1300, 1800]
        value: [18.5, 25, 32, 36.5]
        bounds: [constant, constant]

    heat_capacity:
        dependency: [300, 800, 1300, 1800]
        value: [450, 500, 550, 600]
        bounds: [constant, constant]

    thermal_diffusivity:
        dependency: (300, 3000, 5.0)
        equation: heat_conductivity / (density * heat_capacity)
        bounds: [extrapolate, extrapolate]
```

2. Load the material in Python:
```python
import sympy as sp
from materforge.parsing.api import create_material

# Create symbolic temperature variable
T = sp.Symbol('T')

# Load the material definition
material = create_material("simple_steel.yaml", T)

# Print basic information
print(f"Material: {material.name}")
print(f"Composition: {material.composition}")
print(f"Temperature range: {material.solidus_temperature}K - {material.liquidus_temperature}K")

# Evaluate properties at specific temperatures
temperature = 500 # Kelvin
density = material.density.evalf(T, temperature)
conductivity = material.heat_conductivity.evalf(T, temperature)
heat_capacity = material.heat_capacity.evalf(T, temperature)
diffusivity = material.thermal_diffusivity.evalf(T, temperature)

print(f"At {temperature}K:")
print(f" Density: {density:.2f} kg/m³")
print(f" Thermal Conductivity: {conductivity:.2f} W/(m·K)")
print(f" Heat Capacity: {heat_capacity:.2f} J/(kg·K)")
print(f" Thermal Diffusivity: {diffusivity:.2e} m²/s")
```


## Next Steps

Now that you've created your first material, you can:

- Learn how to [create your first simulation](first_simulation.md)
- Explore [defining custom material properties](../how-to/define_materials.md)
- Understand the [YAML schema](../reference/yaml_schema.md)
- Read about [energy-temperature conversion](../how-to/energy_temperature_conversion.md)
