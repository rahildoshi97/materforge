# Defining Custom Material Properties

This guide explains how to define custom material properties in MaterForge using
different property types.

## Overview

MaterForge uses a schema-agnostic YAML format - any material kind (metals, alloys,
polymers, ceramics, composites, or hypothetical materials) and any property name are
valid. Properties can be constants or expressions driven by any independent variable
(temperature, pressure, composition, strain, or any SymPy symbol). SI units are
recommended throughout.

## YAML Configuration Options

### 1. Constant Value

For properties that do not vary with the dependency variable:

```yaml
properties:
    thermal_expansion_coefficient: 16.3e-6
    density: 7000.0
```

### 2. Step Functions

For properties that change abruptly at a transition point:

```yaml
properties:
    latent_heat_of_fusion:
        dependency: solidus_temp + 5   # reference to another scalar property
        value: [0.0, 171401.0]         # [before_transition, after_transition]
        bounds: [constant, constant]
```

`dependency` here is a single scalar reference or arithmetic expression that defines
the transition point.

### 3. Importing from External Files

For properties defined in spreadsheets or data files:

```yaml
properties:
    # Excel file import
    density:
        file_path: ./material_data.xlsx
        dependency_column: T (K)
        property_column: Density (kg/m^3)
        bounds: [constant, constant]

    # CSV file import
    heat_capacity:
        file_path: ./heat_capacity_data.csv
        dependency_column: Temperature
        property_column: Cp
        bounds: [constant, constant]

    # Text file import (space/tab separated, headerless - use column index)
    thermal_conductivity:
        file_path: ./conductivity_data.txt
        dependency_column: 0
        property_column: 1
        bounds: [constant, constant]
```

Supported formats: `.txt` (space/tab separated), `.csv`, `.xlsx`.

### 4. Tabular Data

For properties defined by paired dependency-value lists:

```yaml
properties:
    heat_conductivity:
        dependency: [500, 1000, 1600, 1700, 1750, 2000, 2500]
        value: [19.25, 25.47, 32.94, 33.52, 31.53, 35.33, 42.95]
        bounds: [constant, constant]

    # Using scalar property references in the dependency list
    latent_heat_of_fusion:
        dependency: [solidus_temp - 1, liquidus_temp + 1]
        value: [0, 171401.0]
        bounds: [constant, constant]

    # Using tuple notation
    density:
        dependency: (1735.0, -5)   # start=1735, decrement by 5 per value
        value: [7037.47, 7060.15, 7088.80, 7110.46, 7127.68]
        bounds: [constant, constant]
```

### 5. Piecewise Equations

For properties with different expressions over different ranges:

```yaml
properties:
    viscosity:
        dependency: [300, 1660, 1736, 3000]        # n breakpoints -> n-1 equations
        equation: [7877.39-0.37*T, 11816.63-2.74*T, 8596.40-0.88*T]
        bounds: [constant, constant]
```

The symbol in the equation (`T` here) is just a placeholder.
It is replaced by the symbol passed to `create_material(..., dependency=X)` at runtime.
The symbol used in equations (`T` in the example above) is a placeholder only.

So this is valid:

```python
P = sp.Symbol('P')
mat = create_material('myAlloy.yaml', dependency=P)
print(mat.viscosity)  # Piecewise in P, not T
```

The YAML file never needs to change regardless of what symbol the caller uses.

### 6. Computed Properties

For properties derived from other already-defined properties:

```yaml
properties:
    thermal_diffusivity:
        dependency: (3000, 300, -5.0)
        equation: heat_conductivity / (density * heat_capacity)
        bounds: [linear, linear]
        regression:
            simplify: post
            degree: 2
            segments: 3
```

MaterForge resolves property dependencies automatically - `thermal_diffusivity` will
always be processed after `heat_conductivity`, `density`, and `heat_capacity`
regardless of their order in the YAML file.

---

## Dependency Definition Formats

The `dependency` key accepts five formats:

```yaml
# 1. Explicit list of values
dependency: [300, 500, 800, 1000]

# 2. Tuple (start, increment) - length inferred from value list
dependency: (300, 50)            # 300, 350, 400, ...

# 3. Tuple (start, stop, step) - step is float
dependency: (300, 1000, 10.0)    # 300, 310, 320, ..., 1000

# 4. Tuple (start, stop, points) - points is integer
dependency: (300, 1000, 71)      # 71 evenly spaced values

# 5. Decreasing range
dependency: (1000, 300, -5.0)    # 1000, 995, 990, ..., 300
```

Scalar property references and arithmetic expressions are also valid inside lists:

```yaml
dependency: [solidus_temp - 1, liquidus_temp + 1]
dependency: solidus_temp + 5     # single reference for step functions
```

---

## Bounds Options

Controls behaviour outside the defined range:

```yaml
bounds: [constant, linear]
```

| Option      | Behaviour                                        |
|-------------|--------------------------------------------------|
| constant    | Clamp to the boundary value outside the range    |
| linear | Linear extrapolation beyond the range            |

---

## Regression Configuration

Reduces expression complexity and smooths noisy data:

```yaml
regression:
    simplify: pre     # 'pre': applied to raw data before symbolic processing
                      # 'post': applied after symbolic expressions are evaluated
    degree: 1         # polynomial degree per segment (1=linear, 2=quadratic, ...)
    segments: 3       # number of piecewise segments (<=8 recommended)
```

Note: High segment counts (>6) risk overfitting. MaterForge will warn if this
threshold is exceeded.

---

## Complete Example

```yaml
name: "myAlloy"

properties:
    density: 7000.0

    solidus_temp: 1605.
    liquidus_temp: 1735.

    latent_heat_of_fusion:
        dependency: solidus_temp + 5
        value: [0.0, 171401.0]
        bounds: [constant, constant]

    heat_conductivity:
        dependency: [500, 1000, 1600, 1700, 1750, 2000, 2500]
        value: [19.25, 25.47, 32.94, 33.52, 31.53, 35.33, 42.95]
        bounds: [linear, linear]

    heat_capacity:
        file_path: ./myAlloy.csv
        dependency_column: T (K)
        property_column: Specific heat (J/(Kg K))
        bounds: [constant, constant]
        regression:
            simplify: pre
            degree: 3
            segments: 6

    viscosity:
        dependency: [300, 1660, 1736, 3000]
        equation: [7877.39-0.37*T, 11816.63-2.74*T, 8596.40-0.88*T]
        bounds: [constant, constant]

    thermal_diffusivity:
        dependency: (3000, 300, -5.0)
        equation: heat_conductivity / (density * heat_capacity)
        bounds: [linear, linear]
        regression:
            simplify: post
            degree: 3
            segments: 7
```

Load it in Python:

```python
import sympy as sp
from materforge import create_material

T = sp.Symbol('T')   # any symbol works: sp.Symbol('P'), sp.Symbol('x'), etc.
mat = create_material('myAlloy.yaml', dependency=T, enable_plotting=True)

# Access symbolic expressions directly
print(mat.heat_conductivity)     # SymPy Piecewise expression in T
print(mat.density)               # 7000.0 (constant float)

# Evaluate all properties at a specific value
evaluated = mat.evaluate(T, 500.0)
print(evaluated.heat_conductivity)           # sp.Float, prints as number e.g. 25.47
print(float(evaluated.heat_conductivity))    # Python float if downstream code requires i
```

---

## Best Practices

- Use SI units throughout and document units in comments
- Cover the full range of your dependency variable in tabular and piecewise data
- Validate property values against experimental data where possible
- All numeric values must use '.' as the decimal separator, not ','
- Interpolation between data points is automatic for tabular and file-import properties
- Keep regression segments <= 8; prefer simplify: post for computed properties
