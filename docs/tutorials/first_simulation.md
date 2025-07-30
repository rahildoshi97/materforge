# Creating Your First Material Simulation

This tutorial will guide you through creating a basic heat equation simulation using [MaterForge](https://i10git.cs.fau.de/rahil.doshi/materforge) and [pystencils](https://pycodegen.pages.i10git.cs.fau.de/pystencils/).
It builds upon the existing [waLBerla tutorial for code generation](https://walberla.net/doxygen/tutorial_codegen01.html), adding material property handling with materforge.

## Prerequisites

Before starting, ensure you have:
- Completed the [Getting Started](getting_started.md) tutorial
- Installed pystencils and [pystencilssfg](https://pycodegen.pages.i10git.cs.fau.de/pystencils-sfg/)
- Basic understanding of heat transfer and the heat equation
- Created the `simple_steel.yaml` file from the [Getting Started](getting_started.md) tutorial

## Overview

We'll create a 2D heat equation simulation with temperature-dependent material properties using MaterForge's architecture.

## Step 1: Set Up the Simulation Framework
```python
import sympy as sp
import pystencils as ps
from pystencilssfg import SourceFileGenerator
from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

with SourceFileGenerator() as sfg:
    data_type = "float64"

    # Define temperature fields and output field for thermal diffusivity
    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    thermal_diffusivity_field = ps.fields(f"thermal_diffusivity_field: {data_type}[2D]", layout='fzyx')
```

## Step 2: Create Material with Temperature Dependencies
```python
# Create symbolic temperature variable
T = sp.Symbol('T')

# Load material with temperature dependency
material = create_material("simple_steel.yaml", u.center())

# Access material properties
thermal_diffusivity = material.thermal_diffusivity
```

## Step 3: Set Up the Heat Equation
```python
# Create symbolic variables for the equation
dx = sp.Symbol("dx")
dt = sp.Symbol("dt")
thermal_diffusivity_sym = sp.Symbol("thermal_diffusivity")

# Define the heat equation using finite differences
heat_pde = ps.fd.transient(u) - thermal_diffusivity_sym * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1))

# Discretize the PDE
discretize = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
heat_pde_discretized = discretize(heat_pde)
heat_pde_discretized = heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify()
```

## Step 4: Generate the Simulation Code
```python
from sfg_walberla import Sweep

# Create assignment collection with material property assignments
ac = ps.AssignmentCollection(
    subexpressions=subexp,
    main_assignments=[
        ps.Assignment(u_tmp.center(), heat_pde_discretized),
        ps.Assignment(thermal_diffusivity_field.center(), thermal_diffusivity_symbol)
    ])

# Generate the sweep
sweep = Sweep("HeatEquationKernelWithMaterial", ac)
sfg.generate(sweep)
```

## Step 5: Energy-Temperature Conversion
For simulations requiring energy-temperature conversion:
```python
# Create inverse function for temperature from energy density

if hasattr(material, 'energy_density'):
    E = sp.Symbol('E')
    inverse_energy_density = PiecewiseInverter.create_energy_density_inverse(material, 'E')
```

## Complete Example

Here's the complete example based on your current codebase:
```python
import sympy as sp
import pystencils as ps
from pystencilssfg import SourceFileGenerator
from walberla.codegen import Sweep

from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

with SourceFileGenerator() as sfg:
    data_type = "float64"

    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    thermal_diffusivity_symbol = sp.Symbol("thermal_diffusivity")
    thermal_diffusivity_field = ps.fields(f"thermal_diffusivity_field: {data_type}[2D]", layout='fzyx')
    dx, dt = sp.Symbol("dx"), sp.Symbol("dt")

    heat_pde = ps.fd.transient(u) - thermal_diffusivity_symbol * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1))

    discretize = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify()

    # Create material with temperature-dependent properties
    material = create_material("simple_steel.yaml", u.center())

    if hasattr(material, 'energy_density'):
        E = sp.Symbol('E')
        inverse_energy_density = PiecewiseInverter.create_energy_density_inverse(material, 'E')

    subexp = [
        ps.Assignment(thermal_diffusivity_symbol, material.thermal_diffusivity),
    ]

    ac = ps.AssignmentCollection(
        subexpressions=subexp,
        main_assignments=[
            ps.Assignment(u_tmp.center(), heat_pde_discretized),
            ps.Assignment(thermal_diffusivity_field.center(), thermal_diffusivity_symbol)
        ])

    sweep = Sweep("HeatEquationKernelWithMaterial", ac)
    sfg.generate(sweep)
```

## Next Steps

Now that you've created your first simulation with temperature-dependent material properties, you can:
- Learn about [energy-temperature conversion](../how-to/energy_temperature_conversion.md)
- Explore more complex [material properties](../explanation/material_properties.md)
- Understand the [API reference](../reference/api/material.md) for advanced usage
- Check the original [waLBerla tutorial](https://walberla.net/doxygen/tutorial_codegen01.html) for details on running the simulation
