# Inverting Material Properties

This guide explains how to evaluate dependency-driven properties and create their
inverse functions in MaterForge.

## Why Property Inversion Matters

Material properties in MaterForge are expressed as dependency-driven symbolic functions.
Sometimes you need to go in the reverse direction - given a property value, find the
corresponding dependency value. Common examples:

- Given energy density, find temperature
- Given density, find the pressure at which it occurs
- Given a measured property value, recover the driving variable

MaterForge supports this through `PiecewiseInverter`, which creates an inverse symbolic
expression from any piecewise property.

---

## Forward Evaluation

Given a material with a dependency-driven property, evaluate it at a specific value:

```python
import sympy as sp
from materforge.parsing.api import create_material

T = sp.Symbol('T')   # any symbol works
mat = create_material('myAlloy.yaml', dependency=T)

# Evaluate all properties at T=1500
evaluated = mat.evaluate(T, 1500.0)
print(evaluated.energy_density)        # sp.Float
print(float(evaluated.energy_density)) # Python float if needed

# Or access the symbolic expression and substitute manually
expr = mat.energy_density              # SymPy Piecewise in T
value = float(expr.subs(T, 1500.0).evalf())
```

---

## Inverse: From Property Value Back to Dependency

If your property is piecewise linear, `PiecewiseInverter` creates the inverse function
as a new SymPy Piecewise expression:

```python
import sympy as sp
from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

T = sp.Symbol('T')
E = sp.Symbol('E')
mat = create_material('myAlloy.yaml', dependency=T)

# Create inverse: energy_density(T) -> T(E)
inverse = PiecewiseInverter.create_inverse(mat.energy_density, T, E)

# Evaluate: given energy density, recover temperature
energy_value = 1.5e9   # J/m^3
recovered_T = float(inverse.subs(E, energy_value))
print(f"Temperature at E={energy_value:.2e}: {recovered_T:.2f}")
```

`input_symbol` and `output_symbol` are the actual SymPy symbols used in the
original expression and its inverse respectively
---

## Round-Trip Validation

Always validate the inverse by checking the round-trip error:

```python
test_values = [500.0, 800.0, 1200.0, 1500.0, 1800.0]

print(f"{'T_original':>12}  {'E':>14}  {'T_recovered':>12}  {'error':>10}")
print("-" * 55)
for t_orig in test_values:
    e = float(mat.energy_density.subs(T, t_orig).evalf())
    t_rec = float(inverse.subs(E, e))
    err = abs(t_rec - t_orig)
    print(f"{t_orig:>12.2f}  {e:>14.4e}  {t_rec:>12.4f}  {err:>10.2e}")
```

---

## Using a Different Dependency Symbol

The workflow is identical regardless of what symbol drives the property:

```python
import sympy as sp
from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

P   = sp.Symbol('P')    # pressure-driven material
rho = sp.Symbol('rho')
mat = create_material('myMaterial.yaml', dependency=P)

inverse = PiecewiseInverter.create_inverse(mat.density, P, rho)

# Given a density value, recover pressure
target_density = 7200.0
pressure = float(inverse.subs(rho, target_density))
```

---

## Known Limitation

`PiecewiseInverter` currently supports **linear piecewise functions only**. Properties
with higher-degree polynomial segments (e.g. `degree: 2` or higher in regression
configuration) cannot be inverted automatically. Use `degree: 1` in the regression
block for any property you intend to invert:

```yaml
energy_density:
    dependency: (300, 3000, 5.0)
    equation: density * specific_enthalpy
    bounds: [linear, linear]
    regression:
        simplify: pre
        degree: 1      # required for inversion
        segments: 6
```

---

## Best Practices

- Validate round-trip accuracy before using an inverse in production
- Use `degree: 1` regression for properties you intend to invert
- The input property must be monotonic over its domain for the inverse to be
  well-defined - verify this by inspecting the forward plot
- Prefer `bounds: [linear, linear]` for properties used in inversion
  so out-of-range queries degrade gracefully rather than clamping silently
