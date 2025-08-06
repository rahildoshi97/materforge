# Defining Custom Material Properties

This guide explains how to define custom material properties in MaterForge using different methods.

## YAML Configuration Options

MaterForge supports several ways to define material properties in YAML files, following SI units (m, s, kg, A, V, K, etc.).

### 1. Constant Value

For properties that don't vary with temperature:

```yaml
properties:
    thermal_expansion_coefficient: 16.3e-6
```

### 2. Step Functions

For properties that change abruptly at phase transitions:
```YAML
properties:
    latent_heat_of_fusion:
        dependency: solidus_temperature
        value: [0.0, 171401.0]
        bounds: [constant, constant]
```

### 3. Importing from External Files

For properties defined in spreadsheets:

```yaml
properties:
    # Excel file import
    density:
        file_path: ./1.4301_data.xlsx
        dependency_column: T (K)
        property_column: Density (kg/(m)^3)
        bounds: [constant, constant]
    
    # CSV file import
    heat_capacity:
        file_path: ./heat_capacity_data.csv
        dependency_column: Temperature
        property_column: Cp
        bounds: [constant, constant]
    
    # Text file import (space/tab separated)
    thermal_conductivity:
        file_path: ./conductivity_data.txt
        dependency_column: 0  # Column index for headerless files
        property_column: 1
        bounds: [constant, constant]
```
Supported file imports include .txt (space/tab separated), .csv, and .xlsx.

### 4. Tabular Data for Interpolation

For properties that vary with temperature:

```yaml
properties:
    # Using explicit temperature list
    heat_conductivity:
        dependency: [1200, 1800, 2200, 2400]  # Temperatures in Kelvin
        value: [25, 30, 33, 35]  # Property values
        bounds: [constant, constant]
    
    # Using references to defined temperatures
    latent_heat_of_fusion:
        dependency: [solidus_temperature, liquidus_temperature]
        value: [171401]
        bounds: [constant, constant]
    
    # Using tuple for temperature generation
    heat_capacity:
        dependency: (1000, 200)  # Start at 1000K and increment by 200K for each value
        value:
        bounds: [constant, constant]
    
    # Using tuple with negative increment
    density:
        dependency: (1735.00, -5)  # Start at 1735.00K and decrement by 5K for each value
        value: [7037.470, 7060.150, 7088.800, 7110.460, 7127.680, 7141.620, 7156.800, 7172.590, 7184.010, 7192.780]
        bounds: [constant, constant]
```

### 5. Piecewise Equations

For properties with different equations in different temperature ranges:
```yaml
properties:
    heat_conductivity:
        dependency: [1700][3000]
        equation: ["0.012T + 13", "0.015T + 5"]
        bounds: [constant, constant]
```

### 6. Computed Properties

For properties that can be derived from others:

```yaml
properties:
    # Thermal diffusivity computation
    thermal_diffusivity:
        dependency: (300, 3000, 5.0)
        equation: heat_conductivity / (density * heat_capacity)
        bounds: [extrapolate, extrapolate]
    
    # Energy density with different models
    energy_density:
        dependency: (300, 3000, 541)  # 541 evenly spaced points
        equation: density * specific_enthalpy
        bounds: [extrapolate, extrapolate]
```

## Dependency Definition Formats

### Explicit Lists
```yaml
dependency: # Explicit values
```

### Tuple Formats
```yaml
# (start, increment) - requires matching value list length
dependency: (300, 50) # 300, 350, 400, ... (length from values)

# (start, stop, step) - step size
dependency: (300, 1000, 10.0) # From 300K to 3000K in steps of 10K

# (start, stop, points) - number of points
dependency: (300, 1000, 71) # 71 evenly spaced points

# Decreasing temperature
dependency: (1000, 300, -5.0) # 1000, 995, 990, ..., 300
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

## Advanced Features

### Regression Configuration
Control data simplification and memory usage:
```yaml
properties:
    heat_conductivity:
        dependency: [1200, 1800, 2200, 2400]  # Temperature values in Kelvin
        value: [25, 30, 33, 35]  # Property values
        bounds: [constant, constant]
        regression:
            simplify: pre # Apply before symbolic processing
            degree: 1 # Linear regression
            segments: 3 # Number of piecewise segments
```

### Boundary Behavior
Control extrapolation outside data range:
```yaml
bounds: [constant, extrapolate]
```
- `constant`: Use boundary values as constants outside range
- `extrapolate`: Linear extrapolation outside range

## Creating a Complete Material Definition

Here's a complete example for Aluminum (Pure metal):
```yaml
# ====================================================================================================
# MATERFORGE MATERIAL CONFIGURATION FILE - PURE METAL
# ====================================================================================================
# This file defines material properties for pure Aluminum using the materforge format.
# MaterForge supports 6 property types: CONSTANT_VALUE, STEP_FUNCTION, FILE_IMPORT, TABULAR_DATA, PIECEWISE_EQUATION, and COMPUTED_PROPERTY
# Pure metals require 'melting_temperature' and 'boiling_temperature' instead of solidus/liquidus
# ====================================================================================================

name: Aluminum
material_type: pure_metal  # Must be 'pure_metal' for single-element materials

# Composition must sum to 1.0 (for pure metals, single element = 1.0)
composition:
  Al: 1.0  # Aluminum purity (99.95%)

# Required temperature properties for pure metals
melting_temperature: 933.47   # Temperature where solid becomes liquid (K)
boiling_temperature: 2743.0   # Temperature where liquid becomes gas (K)

properties:
  
  # ====================================================================================================
  # STEP_FUNCTION EXAMPLES
  # Step functions use a single dependency transition point with before/after values
  # ====================================================================================================
  
  latent_heat_of_fusion:
    dependency: melting_temperature - 1    # Transition 1K before melting point
    value: [0.0, 10790.0]                  # [before_transition, after_transition] in J/kg
    bounds: [constant, constant]           # Keep constant values outside range
  
  latent_heat_of_vaporization:
    dependency: boiling_temperature        # Transition at boiling point
    value: [0.0, 294000.0]                 # [before_boiling, after_boiling] in J/kg
    bounds: [constant, constant]
  
  # ====================================================================================================
  # TABULAR_DATA EXAMPLES
  # Tabular data pairs with explicit dependency-property relationships
  # ====================================================================================================
  
  heat_capacity:
    dependency: (273.15, 100.0)            # Tuple: (start=273.15K, increment=100K)
    # Generates: [273.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
    value: [897, 921, 950, 980, 1010, 1040, 1070, 1084]  # J/kg·K values
    bounds: [constant, constant]
    regression:
      simplify: post                       # Apply linear regression after processing
      degree: 1                            # Linear fit
      segments: 1                          # Single segment (no breakpoints)

  thermal_expansion_coefficient:
    dependency: [373.15, 473.15, 573.15, 673.15, 773.15, 873.15]  # Explicit temperature list
    value: [24.56e-06, 26.54e-06, 28.51e-06, 30.49e-06, 32.47e-06, 34.45e-06]  # 1/K values
    bounds: [constant, constant]
    regression:
      simplify: post
      degree: 1
      segments: 1

  # ====================================================================================================
  # PIECEWISE_EQUATION EXAMPLE
  # Multiple equations for different dependency ranges
  # ====================================================================================================
  
  heat_conductivity:
    temperature: [500, 1700, 3000]         # Three breakpoints define two ranges
    equation: [0.0124137215440647*T + 13.0532171803243, 0.0124137215440647*T + 13.0532171803243]
    # Equation 1: 500-1700K, Equation 2: 1700-3000K
    bounds: [constant, constant]
    regression:
      simplify: post
      degree: 1
      segments: 1
  
  # ====================================================================================================
  # COMPUTE EXAMPLES
  # Properties calculated from other properties using symbolic expressions
  # ====================================================================================================
  
  density:
    dependency: (300, 3000, 541)           # (start, stop, points) - 541 evenly spaced points
    equation: 2700 / (1 + 3 * thermal_expansion_coefficient * (T - 293.15))  # Thermal expansion model
    bounds: [constant, constant]
    regression:
      simplify: post
      degree: 1
      segments: 1
  
  thermal_diffusivity:
    dependency: (3000, 300, -5.)           # (start, stop, negative_step) - decreasing by 5K
    equation: heat_conductivity / (density * heat_capacity)  # Standard thermal diffusivity formula
    bounds: [extrapolate, extrapolate]     # Allow extrapolation outside range
    regression:
      simplify: post
      degree: 2
      segments: 1

  energy_density:
    dependency: (300, 3000, 5.0)           # (start, stop, step) - 5K increments
    equation: density * specific_enthalpy  # Property dependencies automatically resolved
    bounds: [extrapolate, extrapolate]
    regression:
      simplify: pre
      degree: 1
      segments: 6

  specific_enthalpy:
    dependency: (300, 3000, 541)           # (start, stop, num_points) - 541 evenly spaced points
    equation: Integral(heat_capacity, T)   # Symbolic integration
    bounds: [ constant, constant ]
    regression:
      simplify: post
      degree: 1
      segments: 2

# ====================================================================================================
# PURE METAL SPECIFIC NOTES:
# 
# 1. REQUIRED FIELDS:
#    - name, material_type: pure_metal, composition, melting_temperature, boiling_temperature
# 
# 2. TEMPERATURE REFERENCES:
#    - melting_temperature, boiling_temperature (instead of solidus/liquidus for alloys)
#    - Arithmetic expressions: melting_temperature + 50, boiling_temperature - 100
# 
# 3. PROPERTY DEPENDENCIES:
#    - Properties are processed in dependency order automatically
#    - thermal_diffusivity depends on heat_conductivity, density, heat_capacity
#    - density depends on thermal_expansion_coefficient
# 
# 4. DEPENDENCY TUPLE FORMATS:
#    - (start, increment): Requires n_values from 'value' list length
#    - (start, stop, step): Step size (float) or number of points (integer)
#    - Negative steps create decreasing temperature arrays
# 
# 5. REGRESSION BENEFITS:
#    - Reduces memory usage and computational complexity
#    - Smooths noisy data while preserving essential behavior
#    - 'pre': Applied to raw data before symbolic processing
#    - 'post': Applied after symbolic expressions are evaluated
# ====================================================================================================
```

Here's a complete example for 1.4301 (Alloy):

```yaml
# ====================================================================================================
# MATERFORGE MATERIAL CONFIGURATION FILE
# ====================================================================================================
# This file defines material properties for 1.4301 (Steel) using the new materforge format.
# MaterForge supports 6 property types: CONSTANT_VALUE, STEP_FUNCTION, FILE_IMPORT, TABULAR_DATA, PIECEWISE_EQUATION, and COMPUTED_PROPERTY
#
# IMPORTANT: All property configurations must include 'bounds' parameter (except CONSTANT properties)
# ====================================================================================================

name: "1.4301"
material_type: alloy  # Must be 'alloy' or 'pure_metal'

# Composition fractions must sum to 1.0
composition:
  Fe: 0.675
  Cr: 0.170
  Ni: 0.120
  Mo: 0.025
  Mn: 0.01

# Required temperature properties for alloys
solidus_temperature: 1605.          # Temperature where melting begins (K)
liquidus_temperature: 1735.         # Temperature where material is completely melted (K)
initial_boiling_temperature: 3090.  # Temperature where boiling begins (K)
final_boiling_temperature: 3200.    # Temperature where material is completely vaporized (K)

properties:

  # ====================================================================================================
  # PROPERTY TYPE 1: CONSTANT
  # Format: property_name: numeric_value
  # Note: Must be float, not integer (use 1.0 instead of 1)
  # ====================================================================================================

  latent_heat_of_vaporization: 1.71401E5  # J/kg - Scientific notation supported

  # ====================================================================================================
  # PROPERTY TYPE 2: STEP_FUNCTION
  # Format: Uses 'dependency' (single reference) and 'value' (list of 2 values)
  # Supports temperature references and arithmetic expressions
  # ====================================================================================================

  # Example with temperature reference arithmetic:
  # latent_heat_of_fusion:
  #   dependency: solidus_temperature - 1  # Arithmetic expressions supported
  #   value: [0.0, 171401.0]               # [before_transition, after_transition]
  #   bounds: [constant, constant]         # Required for all non-constant properties

  # ====================================================================================================
  # PROPERTY TYPE 3: TABULAR_DATA
  # Format: Uses 'dependency' (list/tuple) and 'value' (list) with same length
  # Supports temperature references, explicit lists, and tuple notation
  # ====================================================================================================

  latent_heat_of_fusion:
    dependency: [solidus_temperature - 1, liquidus_temperature + 1]  # Temperature references with arithmetic
    value: [0, 171401.]                    # Corresponding property values
    bounds: [constant, constant]           # Boundary behavior: 'constant' or 'extrapolate'

  thermal_expansion_coefficient:
    dependency: [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1500, 2000, 2500, 3000]
    value: [1.2632e-5, 1.468e-5, 1.524e-5, 1.581e-5, 1.639e-5, 1.699e-5, 1.759e-5, 1.821e-5, 1.885e-5, 2.1e-5, 2.3e-5, 2.5e-5, 2.7e-5]
    bounds: [constant, constant]
    regression:                            # Optional regression configuration
      simplify: pre                        # 'pre' (before processing) or 'post' (after processing)
      degree: 1                            # Polynomial degree for regression
      segments: 2                          # Number of piecewise segments

  # ====================================================================================================
  # PROPERTY TYPE 4: FILE_IMPORT
  # Format: Uses 'file_path', 'dependency_column', 'property_column'
  # Supports .txt, .csv, .xlsx files with automatic missing value handling
  # ====================================================================================================

  heat_capacity:
    file_path: ./1.4301.xlsx               # Relative path from YAML file location
    dependency_column: T (K)               # Column name for temperature data
    property_column: Specific heat (J/(Kg K)) # Column name for property data
    bounds: [constant, constant]           # Required boundary behavior
    regression:                            # Optional regression for data simplification
      simplify: pre                        # Apply regression before processing
      degree: 1                            # Linear regression
      segments: 4                          # Divide into 4 piecewise segments

  density:
    file_path: ./1.4301.xlsx
    dependency_column: T (K)
    property_column: Density (kg/(m)^3)
    bounds: [constant, constant]
    regression:
      simplify: post                       # Apply regression after processing
      degree: 2
      segments: 3

  # ====================================================================================================
  # PROPERTY TYPE 5: PIECEWISE_EQUATION
  # Format: Uses 'dependency' (breakpoints) and 'equation' (list of expressions)
  # Number of equations = number of breakpoints - 1
  # ====================================================================================================

  heat_conductivity:
    temperature: [500, 1700, 3000]         # Temperature breakpoints (K)
    equation: [0.0124137215440647*T + 13.0532171803243, 0.0124137215440647*T + 13.0532171803243]
    # Two equations for three breakpoints: [500-1700K] and [1700-3000K]
    bounds: [constant, constant]           # Boundary behavior outside range
    regression:
      simplify: post                       # Apply regression after symbolic processing
      degree: 1
      segments: 2

  # ====================================================================================================
  # PROPERTY TYPE 6: COMPUTE
  # Format: Uses 'equation' (single expression) and 'dependency' (range definition)
  # Supports dependency resolution and automatic ordering
  # ====================================================================================================

  specific_enthalpy:
    dependency: (300, 3000, 541)           # (start, stop, num_points) - 541 evenly spaced points
    equation: Integral(heat_capacity, T)   # Symbolic integration
    bounds: [constant, constant]
    regression:
      simplify: post
      degree: 1
      segments: 2

  energy_density:
    dependency: (300, 3000, 5.0)           # (start, stop, step) - 5K increments
    equation: density * specific_enthalpy  # Property dependencies automatically resolved
    bounds: [extrapolate, extrapolate]
    regression:
      simplify: pre
      degree: 1
      segments: 6

  thermal_diffusivity:
    dependency: (3000, 300, -5.0)          # (start, stop, negative_step) - decreasing temperature
    equation: heat_conductivity /(density * heat_capacity)  # Multi-property dependency
    bounds: [constant, constant]
    regression:
      simplify: post
      degree: 1
      segments: 3

# ====================================================================================================
# DEPENDENCY/TEMPERATURE DEFINITION FORMATS:
#
# 1. Single value: 1000.0
# 2. List: [300, 400, 500, 600]
# 3. Temperature references: solidus_temperature, liquidus_temperature, etc.
# 4. Arithmetic expressions: melting_temperature + 50, solidus_temperature - 10
# 5. Tuple formats:
#    - (start, increment): (300, 50) with n_values parameter
#    - (start, stop, step): (300, 1000, 10.0) for 10K increments
#    - (start, stop, points): (300, 1000, 71) for 71 evenly spaced points
#
# BOUNDS OPTIONS:
# - 'constant': Use boundary values as constants outside range
# - 'extrapolate': Linear extrapolation outside range
#
# REGRESSION OPTIONS:
# - simplify: 'pre' (before symbolic processing) or 'post' (after symbolic processing)
# - degree: Polynomial degree (1=linear, 2=quadratic, etc.)
# - segments: Number of piecewise segments for regression
# ====================================================================================================
```

## Best Practices

- Use consistent units throughout your definitions (SI units recommended)
- Document the expected units for each property
- For temperature-dependent properties, cover the full range of temperatures you expect in your simulation
- Validate your property data against experimental values when possible
- All numerical values must use period (.) as decimal separator, not comma
- Interpolation between data points is performed automatically for file-import and tabular data properties

## Validation and Error Handling

MaterForge's new architecture includes comprehensive validation:

- **Composition validation**: Ensures fractions sum to 1.0
- **Temperature validation**: Checks for monotonicity and physical validity
- **Property dependency validation**: Ensures computed properties have required dependencies
- **File validation**: Validates file existence and format
- **Type validation**: Ensures proper data types throughout

## Important Notes

- Properties cannot be defined in multiple ways or multiple times
- Required dependencies for computed properties must be present
- Properties will be computed in the correct order regardless of their position in the file
- Temperature arrays must be monotonic
- Energy density arrays must be monotonic with respect to temperature
