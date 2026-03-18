# Getting Started with MaterForge

This tutorial guides you through the basics of using MaterForge to define and evaluate
material properties for simulations.

## Prerequisites

- Python 3.10 or newer
- Basic familiarity with Python and YAML

---

## Installation

Install MaterForge from PyPI:

```bash
pip install materforge
```

For development, clone the repository and install in editable mode:

```bash
git clone https://i10git.cs.fau.de/rahil.doshi/materforge.git
cd materforge
pip install -e .   # changes to source files take effect immediately
```

---

## Core Concepts

MaterForge is built around a small set of ideas:

1. **Schema-agnostic materials**: Any material kind, any property name - the YAML file
   drives everything. No required fields beyond `name` and `properties`.
2. **Six property types**: Constants, step functions, file imports, tabular data,
   piecewise equations, and computed properties.
3. **Dependency-driven expressions**: Properties are symbolic SymPy expressions driven
   by any independent variable you choose - temperature, pressure, composition, etc.
4. **Automatic dependency resolution**: Computed properties are always processed after
   their dependencies, regardless of order in the YAML file.

---

## Your First Material

### 1. Create `simple_steel.yaml`

```yaml
name: SimpleSteel

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
        bounds: [linear, linear]
```

### 2. Load and Evaluate in Python

```python
import sympy as sp
from materforge import create_material

# Any SymPy symbol works as the dependency variable
T = sp.Symbol('T')

# Load the material - all properties become SymPy expressions in T
mat = create_material("simple_steel.yaml", dependency=T)

print(mat)        # Material: SimpleSteel (4 properties)
print(repr(mat))  # Material(name='SimpleSteel', properties=[...])

# Access symbolic expressions directly
print(mat.density)             # SymPy Piecewise in T
print(mat.thermal_diffusivity) # SymPy Piecewise in T

# Evaluate all properties at a specific value - returns a new Material
evaluated = mat.evaluate(T, 500.0)
print("At 500 K:")
for prop in sorted(evaluated.property_names()):
    print(f"  {prop:<30}: {float(getattr(evaluated, prop)):.6e}")
```

Expected output:

```
At 500 K:
  density                       : 7.825000e+03
  heat_capacity                 : 4.625000e+02
  heat_conductivity             : 2.018750e+01
  thermal_diffusivity           : 5.590e-06
```

### 3. Using a Different Dependency Symbol

The symbol name in the YAML (`T` in equations) is a placeholder only - MaterForge
replaces it with whatever symbol you pass at runtime:

```python
P = sp.Symbol('P')   # pressure-driven material
mat2 = create_material("simple_steel.yaml", dependency=P)
print(mat2.density)  # Piecewise expression in P, not T
```

---

## Validating a YAML File

Check a file for structural correctness without creating the material:

```python
from materforge import validate_yaml_file

is_valid = validate_yaml_file("simple_steel.yaml")
print(f"Valid: {is_valid}")
```

---

## Inspecting a Material Without Creating It

```python
from materforge import get_material_info

info = get_material_info("simple_steel.yaml")
print(info['name'])              # SimpleSteel
print(info['total_properties'])  # 4
print(info['property_types'])    # {'TABULAR_DATA': 3, 'COMPUTED_PROPERTY': 1}
```

---

## Next Steps

- [Create your first simulation](first_simulation.md) with pystencils
- [Define custom material properties](../how-to/define_materials.md) - all six property types
- [YAML schema reference](../reference/yaml_schema.md)
- [Property inversion](../how-to/property_inversion.md) - recovering the dependency value from a property value
