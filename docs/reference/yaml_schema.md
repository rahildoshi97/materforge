# YAML Schema for Material Definition
This document defines the schema for material definition YAML files in materforge v0.5.5.

## Schema Overview

### A valid `pure metal` definition must include:
- `name`: String identifier for the material
- `composition`: Map of element symbols to their mass fractions
- `melting_temperature`: Numeric value in Kelvin
- `boiling_temperature`: Numeric value in Kelvin
- `properties`: Map of property names to their definitions
### A valid `alloy` definition must include:
- `name`: String identifier for the material
- `composition`: Map of element symbols to their mass fractions
- `solidus_temperature`: Numeric value in Kelvin
- `liquidus_temperature`: Numeric value in Kelvin
- `initial_boiling_temperature`: Numeric value in Kelvin
- `final_boiling_temperature`: Numeric value in Kelvin
- `properties`: Map of property names to their definitions

## 🔧 Property Definition Types
Properties can be defined in six different ways:
### 1. Constant Value
Simple numeric values, for properties that don't vary with any dependency:
```yaml
thermal_expansion_coefficient: 16.3e-6 # Single numeric value
```
### 2. Step Function
Phase transitions with before/after values:
```yaml
latent_heat_of_fusion:
  dependency: melting_temperature - 10
  value: [0.0, 10790.0]
  bounds: [constant, constant]
```
### 3. File Import
Import data from Excel, CSV, or text files:
```yaml
heat_capacity:
  file_path: ./data.xlsx
  dependency_column: T (K)
  property_column: Cp (J/kg·K)
  bounds: [constant, constant]
```
### 4. Tabular Data
Explicit dependency-property pairs:
```yaml
thermal_expansion_coefficient:
  dependency: [300, 400, 500, 600]
  value: [1.2e-5, 1.4e-5, 1.6e-5, 1.8e-5]
  bounds: [constant, constant]
```
### 5. Piecewise Equation
Multiple equations for different dependency ranges:
```yaml
heat_conductivity:
  dependency: [500, 1700, 3000]
  equation: ["0.012*T + 13", "0.015*T + 5"]
  bounds: [constant, constant]
```
### 6. Computed Property
Properties calculated from other properties:
```yaml
thermal_diffusivity:
  dependency: (300, 3000, 5.0)
  equation: heat_conductivity / (density * heat_capacity)
  bounds: [extrapolate, extrapolate]
```

## 📊 Dependency Definition Formats
### Explicit Lists
```yaml
dependency: [300, 400, 500, 600]  # Explicit values
```
### Tuple Formats
```yaml
# (start, increment) - requires matching value list length
dependency: (300, 50)  # 300, 350, 400, ... (length from values)

# (start, stop, step) - step size
dependency: (300, 1000, 10.0)  # 300, 310, 320, ..., 1000

# (start, stop, points) - number of points
dependency: (300, 1000, 71)  # 71 evenly spaced points

# Decreasing dependency
dependency: (1000, 300, -5.0)  # 1000, 995, 990, ..., 300
```
### Temperature References
```yaml
# Direct references
temperature: melting_temperature
temperature: solidus_temperature

# Arithmetic expressions
temperature: melting_temperature + 50
temperature: liquidus_temperature - 10
```

## 🎯 Advanced Features
### Regression Configuration
Control data simplification and memory usage:
```yaml
regression:
  simplify: pre    # Apply before symbolic processing
  degree: 1        # Linear regression
  segments: 3      # Number of piecewise segments
```
- `simplify: pre`: Apply regression to raw data before processing
- `simplify: post`: Apply regression after symbolic expressions are evaluated
- `degree`: Polynomial degree (1=linear, 2=quadratic, etc.)
- `segments`: Number of segments for piecewise functions
### Boundary Behavior
Control extrapolation outside data range:
```yaml
bounds: [constant, extrapolate]
```
- `constant`: Use boundary values as constants outside range
- `extrapolate`: Linear extrapolation outside range
### Dependency Resolution
Properties are automatically processed in correct order:
```yaml
# These will be processed in dependency order automatically
specific_enthalpy:
  equation: Integral(heat_capacity, T)

energy_density:
  equation: density * specific_enthalpy  # Depends on specific_enthalpy

thermal_diffusivity:
  equation: heat_conductivity / (density * heat_capacity)  # Multiple dependencies
```

## 📈 Visualization
Automatic plot generation when using symbolic temperature:

```python
import sympy as sp
from materforge.parsing.api import create_material

T = sp.Symbol('T')
material = create_material('1.4301.yaml', T, enable_plotting=True)
# Plots automatically saved to 'materforge_plots/' directory
```

## 🧪 Energy-Temperature Inversion
For applications requiring temperature from energy density:

```python
import sympy as sp
from materforge.parsing.api import create_material
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

# Create inverse function T = f_inv(E)
T = sp.Symbol('T')
E = sp.Symbol('E')
material = create_material('1.4301.yaml', T)

# Create inverse (only for linear piecewise functions)
inverse_func = PiecewiseInverter.create_energy_density_inverse(material, 'E')

# Use inverse function
energy_value = 1.5e9  # J/m³
temperature = float(inverse_func.subs(E, energy_value))
```

## 🔍 Supported Properties

```python
from materforge.parsing.api import get_supported_properties

print(get_supported_properties())
```

## Validation Rules

1. All required top-level fields must be present
2. Material type must be "pure_metal" or "alloy"
3. Composition fractions must sum to approximately 1.0
4. Pure metals must have exactly one element with composition 1.0
5. Alloys must have at least 2 elements with non-zero composition
6. Liquidus temperature must be greater than or equal to solidus temperature
7. Properties cannot be defined in multiple ways or multiple times
8. Required dependencies for computed properties must be present
9. Dependency arrays must be monotonic
10. Energy density arrays must be monotonic with respect to the dependency
11. File paths must be valid and files must exist
12. For tabular data, dependency and value arrays must have the same length
13. When using tuple notation for dependency arrays, the increment must be non-zero

## Important Notes

1. All numerical values must use period (.) as decimal separator, not comma
2. Interpolation between data points is performed automatically for file-import and tabular data properties
3. Properties will be computed in the correct order regardless of their position in the file
4. To retrieve temperature from energy_density, use the interpolation methods from the generated InterpolationArrayContainer class
5. The new architecture provides enhanced error messages and validation compared to previous versions
