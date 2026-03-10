# YAML Schema for Material Definition

This document defines the schema for material definition YAML files in MaterForge.

---

## Schema Overview

A valid material definition requires exactly two top-level fields:

```yaml
name: "myMaterial"   # string identifier - required
properties:          # map of property definitions - required
  ...
```

No other top-level fields are required. Any property name is valid inside `properties`.
The schema is fully driven by what you put there - no hardcoded material types,
compositions, or required property names.

---

## Property Definition Types

### 1. Constant Value

A single numeric value. Not dependency-driven.

```yaml
properties:
    density: 7000.0
    thermal_expansion_coefficient: 16.3e-6
```

Note: Use float notation (`7000.0` not `7000`). Scientific notation is supported (`1.71401e5`).

---

### 2. Step Function

A property that changes discontinuously at a single transition point.
`dependency` must be a scalar property reference or arithmetic expression.

```yaml
properties:
    solidus_temp: 1605.

    latent_heat_of_fusion:
        dependency: solidus_temp + 5
        value: [0.0, 171401.0]   # [before_transition, after_transition]
        bounds: [constant, constant]
```

---

### 3. File Import

Data loaded from an external file. Relative paths are resolved from the YAML file location.

```yaml
properties:
    # Excel
    heat_capacity:
        file_path: ./material_data.xlsx
        dependency_column: T (K)
        property_column: Specific heat (J/(Kg K))
        bounds: [constant, constant]

    # CSV
    density:
        file_path: ./density.csv
        dependency_column: Temperature
        property_column: Density
        bounds: [constant, constant]

    # Text file - headerless, use column index
    thermal_conductivity:
        file_path: ./conductivity.txt
        dependency_column: 0
        property_column: 1
        bounds: [constant, constant]
```

Supported formats: `.xlsx`, `.csv`, `.txt` (space/tab separated).

---

### 4. Tabular Data

Explicit paired dependency-value lists. Both lists must have the same length.

```yaml
properties:
    heat_conductivity:
        dependency: [500, 1000, 1600, 1700, 1750, 2000, 2500]
        value: [19.25, 25.47, 32.94, 33.52, 31.53, 35.33, 42.95]
        bounds: [constant, constant]

    # Scalar property references in dependency
    latent_heat_of_fusion:
        dependency: [solidus_temp - 1, liquidus_temp + 1]
        value: [0.0, 171401.0]
        bounds: [constant, constant]

    # Tuple (start, increment) - length inferred from value list
    heat_capacity:
        dependency: (273.15, 100.0)
        value: [897, 921, 950, 980, 1010, 1040, 1070, 1084]
        bounds: [constant, constant]

    # Tuple with negative increment
    density:
        dependency: (1735.0, -5)
        value: [7037.47, 7060.15, 7088.80, 7110.46, 7127.68]
        bounds: [constant, constant]
```

---

### 5. Piecewise Equation

`n` breakpoints in `dependency` define `n-1` equations.

```yaml
properties:
    viscosity:
        dependency: [300, 1660, 1736, 3000]
        equation: [7877.39-0.37*T, 11816.63-2.74*T, 8596.40-0.88*T]
        bounds: [constant, constant]
```

`T` in equations is a **placeholder only**. MaterForge substitutes it with the symbol
passed to `create_material(..., dependency=symbol)` at runtime. The final symbolic
expressions use your symbol, not `T`.

---

### 6. Computed Property

Derived from other properties using a symbolic expression.
MaterForge resolves processing order automatically via topological sort -
no manual ordering needed.

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

    energy_density:
        dependency: (300, 3000, 5.0)
        equation: density * specific_enthalpy
        bounds: [linear, linear]
```

---

## Dependency Definition Formats

```yaml
# 1. Explicit list
dependency: [300, 500, 800, 1000]

# 2. Tuple (start, increment) - length inferred from value list
dependency: (300, 50)

# 3. Tuple (start, stop, step) - step is float
dependency: (300, 1000, 10.0)     # 300, 310, ..., 1000

# 4. Tuple (start, stop, points) - points is integer
dependency: (300, 1000, 71)       # 71 evenly spaced values

# 5. Decreasing range
dependency: (1000, 300, -5.0)     # 1000, 995, ..., 300
```

Scalar references and arithmetic are valid inside lists:

```yaml
dependency: [solidus_temp - 1, liquidus_temp + 1]
dependency: solidus_temp + 5      # single reference for step functions
```

---

## Bounds Options

Required for all non-constant property types.

```yaml
bounds: [lower_bound, upper_bound]
```

| Option      | Behaviour                                      |
|-------------|------------------------------------------------|
| constant    | Clamp to the boundary value outside the range  |
| linear | Linear extrapolation beyond the range          |

---

## Regression Configuration

Optional. Reduces expression complexity and smooths noisy input data.

```yaml
regression:
    simplify: pre    # 'pre': on raw data before symbolic processing
                     # 'post': after symbolic expressions are evaluated
    degree: 1        # polynomial degree per segment (1=linear, 2=quadratic, ...)
    segments: 3      # number of piecewise segments (<=8 recommended)
```

Note: segments > 6 risk overfitting - MaterForge logs a warning.
Use `degree: 1` for any property you intend to invert with `PiecewiseInverter`.

---

## Visualization

Plots are generated automatically when `dependency` is a SymPy symbol:

```python
import sympy as sp
from materforge.parsing.api import create_material

T = sp.Symbol('T')   # any symbol - not limited to temperature
mat = create_material('myAlloy.yaml', dependency=T, enable_plotting=True)
# Plots saved automatically
```

Pass a float instead to skip plotting and evaluate immediately:

```python
mat = create_material('myAlloy.yaml', dependency=500.0, enable_plotting=False)
```

---

## Property Inversion

For any piecewise-linear property, create its inverse with `PiecewiseInverter`:

```python
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

T = sp.Symbol('T')
E = sp.Symbol('E')

mat = create_material('myAlloy.yaml', dependency=T)
inverse = PiecewiseInverter.create_inverse(mat.energy_density, 'T', 'E')

# Given a property value, recover the dependency value
recovered = float(inverse.subs(E, 1.5e9))
```

Limitation: linear piecewise only (`degree: 1`). Higher-degree segments are not supported.

---

## Validation Rules

1. `name` and `properties` are the only required top-level fields
2. A property cannot be defined more than once
3. All dependencies of a computed property must be defined in the same `properties` block
4. Dependency arrays must be monotonic for interpolation
5. For tabular data, `dependency` and `value` lists must have equal length
6. File paths must exist and be readable; columns must match exactly
7. Tuple increment must be non-zero
8. `bounds` is required for all non-constant property types
9. Regression `segments` must be a positive integer

---

## Important Notes

- All numeric values must use `.` as decimal separator, not `,`
- Interpolation between data points is automatic for tabular and file-import properties
- Property processing order is resolved automatically - position in the file does not matter
- The `T` placeholder in equations carries no special meaning; it is substituted at runtime
