# Material Properties in MaterForge

This document explains the conceptual framework behind temperature-dependent material properties in MaterForge, how they are represented internally, and the mathematical models used for computed properties.

## Conceptual Framework

Material properties in MaterForge are designed around these key principles:

1. **Temperature Dependence**: Most material properties vary with temperature, especially during phase transitions
2. **Symbolic Representation**: Properties are represented as symbolic expressions for mathematical manipulation
3. **Flexible Definition**: Properties can be defined through various methods (constants, data points, files, or computation)
4. **Physical Consistency**: Computed properties follow established physical relationships

## Internal Representation

### Material Class

At the core of MaterForge's property system is the `Material` class, which contains:

- Basic material information (name, type, composition)
- Temperature properties (melting/boiling points, solidus/liquidus temperatures)
- Optional material properties as SymPy expressions
- Validation and calculation methods

```python
@dataclass
class Material:
    name: str
    material_type: str # 'alloy' or 'pure_metal'
    elements: List[ChemicalElement]
    composition: Union[np.ndarray, List, Tuple]
    # Temperature properties vary by material type
    # Optional properties as SymPy expressions
    density: sp.Expr = None
    heat_capacity: sp.Expr = None
    # ... other properties
```

### Property Processing Pipeline

Properties are processed through a sophisticated pipeline:

1. **Type Detection**: `PropertyConfigAnalyzer` automatically determines property definition type
2. **Validation**: Ensures configuration is valid for the detected type
3. **Processing**: `PropertyManager` converts configuration to SymPy expressions
4. **Dependency Resolution**: Handles property interdependencies automatically

## Property Definition Methods

materforge supports multiple ways to define material properties:

1. Constant Values

Properties that don't vary with temperature are defined as simple numeric values:

```yaml
thermal_expansion_coefficient: 16.3e-6
```

Internally converted to `sp.Float(16.3e-6)`.

### 2. Step Functions

Properties with discontinuous changes at phase transitions:
```python
latent_heat_of_fusion:
    temperature: solidus_temperature
    value: [0.0, 171401.0]
    bounds: [constant, constant]
```

Represented as `sp.Piecewise` expressions with temperature-dependent conditions.

### 3. File Import Properties

Properties loaded from external data files:
```python
density:
    file_path: ./material_data.xlsx
    temperature_column: T (K)
    property_column: Density (kg/(m)^3)
    bounds: [constant, constant]
```

Data is loaded via `read_data_from_file` and converted to piecewise interpolation functions.

### 4. Tabular Data

Explicit temperature-property relationships:
```yaml
heat_conductivity:
    temperature: [1200, 1800, 2200, 2400]  # Temperatures in Kelvin
    value: [25, 30, 33, 35]          # Property values
    bounds: [constant, constant]
```

Converted to piecewise linear interpolation functions through `PiecewiseBuilder`.

### 5. Piecewise Equations

Multiple equations for different temperature ranges:
```python
heat_conductivity:
    temperature: [1700][3000]
    equation: ["0.012T + 13", "0.015T + 5"]
    bounds: [constant, constant]
```

Each equation is parsed as a SymPy expression and combined into a piecewise function.

### 6. Computed Properties

Properties calculated from other properties:
```python
thermal_diffusivity:
    temperature: (300, 3000, 5.0)
    equation: heat_conductivity / (density * heat_capacity)
    bounds: [extrapolate, extrapolate]
```

Symbolic expressions that reference other material properties with automatic dependency resolution.

## Temperature Processing

MaterForge provides sophisticated dependency definition processing through `DependencyResolver`:

### Temperature Definition Formats

1. **Explicit Lists**: `[300, 400, 500, 600]`
2. **Tuple Formats**:
    - `(300, 50)` - start and increment
    - `(300, 1000, 10.0)` - start, stop, step
    - `(300, 1000, 71)` - start, stop, points
3. **Temperature References**: `solidus_temperature`, `melting_temperature + 50`

### Temperature Resolution

The `DependencyResolver` class handles:
- Reference resolution to material properties
- Arithmetic expression evaluation
- Validation of temperature ranges
- Conversion to numpy arrays

## Interpolation and Evaluation

### Piecewise Functions

Properties are represented as piecewise functions that:
- Perform linear interpolation between data points
- Handle boundary conditions (constant or extrapolation)
- Support symbolic evaluation with SymPy
- Can be evaluated at specific temperatures using `.evalf()`

### Boundary Handling

Two boundary types are supported:
- **Constant**: Use boundary values outside the defined range
- **Extrapolate**: Linear extrapolation beyond the data range

## Dependency Management

### Dependency Detection

The system automatically:
- Extracts symbols from mathematical expressions using SymPy
- Identifies required properties for computed properties
- Validates that all dependencies are available

### Circular Dependency Prevention

Sophisticated checking prevents circular dependencies:
- Tracks dependency chains during processing in `PropertyManager`
- Detects cycles before they cause infinite loops
- Provides clear error messages through `CircularDependencyError`

### Processing Order

Properties are processed in dependency order:
- Independent properties first
- Dependent properties after their dependencies
- Automatic topological sorting of the dependency graph

## Validation and Quality Assurance

### Data Validation

Comprehensive validation includes:
- Temperature monotonicity checking through `is_monotonic`
- Energy density monotonicity validation via `validate_energy_density_monotonicity`
- Physical reasonableness checks
- Data quality assessment in file processing

### Error Handling

Clear, actionable error messages for:
- Invalid property configurations through `PropertyConfigAnalyzer`
- Missing dependencies via `DependencyError`
- Data quality issues in file processing
- Physical inconsistencies in material properties

## Integration with Simulations

### SymPy Integration

Properties as SymPy expressions enable:
- Symbolic differentiation and integration
- Algebraic manipulation and simplification
- Direct evaluation at specific temperatures

### Simulation Framework Integration

Properties can be used in:
- pystencils-based simulations through symbolic expressions
- Custom simulation frameworks via `.evalf()` method
- Scientific computing workflows through NumPy integration

## Best Practices

### Property Definition

1. **Use Consistent Units**: All properties should use SI units
2. **Cover Full Temperature Range**: Define properties across the entire simulation range
3. **Add Points Around Transitions**: Use more data points near phase transitions
4. **Validate Against Experiments**: Compare with experimental data when possible

### Performance Optimization

1. **Use Appropriate Property Types**: Choose the most efficient definition method
2. **Consider Regression**: Use regression for large datasets via `RegressionManager`
3. **Optimize Temperature Arrays**: Balance accuracy and performance

### Maintainability

1. **Document Property Sources**: Keep track of data origins
2. **Use Descriptive Names**: Clear property and file naming
3. **Validate Configurations**: Use `validate_yaml_file` function
4. **Version Control**: Track changes to material definitions

This framework provides a robust, flexible foundation for modeling complex material behavior in scientific simulations while maintaining ease of use and extensibility.

