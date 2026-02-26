# Creating Your First Material Simulation

This tutorial guides you through creating a basic heat equation simulation using
[MaterForge](https://i10git.cs.fau.de/rahil.doshi/materforge) and
[pystencils](https://pycodegen.pages.i10git.cs.fau.de/pystencils/).
It builds upon the existing [waLBerla tutorial for code generation](https://walberla.net/doxygen/tutorial_codegen01.html),
adding dependency-driven material property handling with MaterForge.

## Prerequisites

Before starting, ensure you have:
- Completed the [Getting Started](getting_started.md) tutorial
- Installed pystencils and [pystencilssfg](https://pycodegen.pages.i10git.cs.fau.de/pystencils-sfg/)
- Basic understanding of heat transfer and the heat equation
- Created the `simple_steel.yaml` file from the [Getting Started](getting_started.md) tutorial

---

## Overview

We will create a 2D heat equation simulation where material properties are driven by the
field variable `u` (representing temperature). MaterForge expresses all dependency-driven
properties as SymPy expressions - these plug directly into pystencils assignment collections.

---

## Step 1: Set Up the Simulation Framework

```python
import sympy as sp
import pystencils as ps
from pystencilssfg import SourceFileGenerator
from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

with SourceFileGenerator() as sfg:
    data_type = "float64"

    # Define the field variable and output field
    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    thermal_diffusivity_field = ps.fields(
        f"thermal_diffusivity_field: {data_type}[2D]", layout='fzyx'
    )
```

---

## Step 2: Create the Material

Pass `u.center()` as the dependency so MaterForge substitutes the field variable
directly into all symbolic property expressions.

```python
    # u.center() is a pystencils field access - MaterForge treats it as the
    # dependency variable, replacing the T placeholder in all YAML equations
    material = create_material("simple_steel.yaml", dependency=u.center())

    # Access the symbolic expression for thermal diffusivity
    # This is a SymPy Piecewise expression in u.center()
    thermal_diffusivity = material.thermal_diffusivity
```

---

## Step 3: Set Up the Heat Equation

```python
    dx = sp.Symbol("dx")
    dt = sp.Symbol("dt")
    thermal_diffusivity_sym = sp.Symbol("thermal_diffusivity")

    heat_pde = (
        ps.fd.transient(u)
        - thermal_diffusivity_sym * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1))
    )

    discretize = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = (
        heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify()
    )
```

---

## Step 4: Build the Assignment Collection and Generate Code

```python
    from sfg_walberla import Sweep

    subexp = [
        ps.Assignment(thermal_diffusivity_sym, material.thermal_diffusivity),
    ]

    ac = ps.AssignmentCollection(
        subexpressions=subexp,
        main_assignments=[
            ps.Assignment(u_tmp.center(), heat_pde_discretized),
            ps.Assignment(thermal_diffusivity_field.center(), thermal_diffusivity_sym),
        ],
    )

    sweep = Sweep("HeatEquationKernelWithMaterial", ac)
    sfg.generate(sweep)
```

---

## Step 5: Property Inversion (Optional)

If your simulation requires recovering the dependency value from a property value
(e.g. recovering temperature from energy density), use `PiecewiseInverter`.
This works for any piecewise-linear property - not just energy density.

```python
import sympy as sp
from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

T = sp.Symbol('T')
material = create_material("simple_steel.yaml", dependency=T)

if hasattr(material, 'energy_density'):
    E = sp.Symbol('E')
    # Creates inverse: given energy density E, returns T
    inverse = PiecewiseInverter.create_inverse(material.energy_density, 'T', 'E')

    # Evaluate
    energy_value = 1.5e9   # J/m^3
    recovered_T = float(inverse.subs(E, energy_value))
    print(f"Recovered temperature: {recovered_T:.2f}")
```

Note: `PiecewiseInverter` requires the property to use `degree: 1` regression.
Higher-degree segments cannot be inverted.

---

## Complete Example

```python
import sympy as sp
import pystencils as ps
from pystencilssfg import SourceFileGenerator
from sfg_walberla import Sweep

from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

with SourceFileGenerator() as sfg:
    data_type = "float64"

    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    thermal_diffusivity_field = ps.fields(
        f"thermal_diffusivity_field: {data_type}[2D]", layout='fzyx'
    )
    thermal_diffusivity_sym = sp.Symbol("thermal_diffusivity")
    dx, dt = sp.Symbol("dx"), sp.Symbol("dt")

    # Discretise heat equation
    heat_pde = (
        ps.fd.transient(u)
        - thermal_diffusivity_sym * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1))
    )
    discretize = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = (
        heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify()
    )

    # Load material - u.center() is the dependency variable
    material = create_material("simple_steel.yaml", dependency=u.center())

    # Optional: build inverse for energy-to-temperature recovery
    if hasattr(material, 'energy_density'):
        T = sp.Symbol('T')
        E = sp.Symbol('E')
        # Rebuild with symbolic T for inversion (inversion needs a sp.Symbol)
        mat_sym = create_material("simple_steel.yaml", dependency=T)
        inverse = PiecewiseInverter.create_inverse(mat_sym.energy_density, 'T', 'E')

    subexp = [
        ps.Assignment(thermal_diffusivity_sym, material.thermal_diffusivity),
    ]

    ac = ps.AssignmentCollection(
        subexpressions=subexp,
        main_assignments=[
            ps.Assignment(u_tmp.center(), heat_pde_discretized),
            ps.Assignment(thermal_diffusivity_field.center(), thermal_diffusivity_sym),
        ],
    )

    sweep = Sweep("HeatEquationKernelWithMaterial", ac)
    sfg.generate(sweep)
```

---

## Next Steps

- Learn about [property inversion](../how-to/property_inversion.md)
- Explore more complex [material properties](../explanation/material_properties.md)
- Read the [API reference](../reference/api/material.md) for advanced usage
- See the original [waLBerla tutorial](https://walberla.net/doxygen/tutorial_codegen01.html)
  for details on compiling and running the generated code
