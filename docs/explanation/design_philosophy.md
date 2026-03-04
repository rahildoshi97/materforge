# Design Philosophy of MaterForge

This document explains the core design principles, architectural decisions,
and the rationale behind MaterForge's structure and implementation.

## Core Principles

MaterForge is built upon several core principles:

- **Modularity**: Clearly separated components for ease of maintenance, testing, and extensibility
- **Flexibility**: Any material kind, any property name, any dependency variable - the YAML
  schema drives everything with no hardcoded constraints
- **Performance**: Leverage symbolic computation for high-performance simulations
- **Transparency and Reproducibility**: Clearly document material property definitions and
  computations to ensure reproducibility

## Layered Architecture

MaterForge follows a layered architecture to separate concerns clearly:

### 1. User Interface Layer (YAML Configuration)

- Provides a simple, human-readable format for defining any material and its properties
- Allows users to specify properties using multiple intuitive methods (constants, step functions,
  interpolation points, file-based data, piecewise equations, computed properties)
- Schema-agnostic: no required fields beyond `name` and `properties`
- Ensures clarity, readability, and ease of use

### 2. Parsing Layer (Python)

- **Configuration Processing**: Handles YAML parsing and validation through `MaterialYAMLParser`
- **Property Type Detection**: Automatically determines property definition types using
  `PropertyTypeDetector`
- **Dependency Resolution**: Processes various dependency definition formats via `DependencyResolver`
- **Data Handling**: Manages file I/O for external data sources through `load_property_data`

### 3. Symbolic Representation Layer (SymPy)

- Uses symbolic mathematics (via SymPy) internally to represent material properties
- The dependency symbol is supplied at runtime via `create_material(..., dependency=symbol)` -
  any SymPy symbol is valid (temperature, pressure, composition fraction, strain, etc.)
- Equations in YAML files use `T` as a placeholder; MaterForge substitutes it with the
  caller-supplied symbol automatically
- Enables symbolic manipulation, simplification, and validation of property definitions
- Facilitates automatic computation of derived properties through symbolic expressions

### 4. Algorithms Layer (Python)

- **Interpolation**: Robust interpolation methods for evaluating dependency-driven properties
- **Regression**: Data simplification and piecewise function generation via `RegressionProcessor`
- **Piecewise Functions**: Creating piecewise expressions through `PiecewiseBuilder`
- **Inversion**: Creating inverse functions for specialized applications

### 5. Visualization Layer (Python)

- Automatic plot generation for material properties through `PropertyVisualizer`
- Property visualization and validation
- Integration with matplotlib for scientific plotting

### 6. Core Layer (Python)

- **Material**: Schema-agnostic material container with fully dynamic property tracking
- **Interfaces**: Abstract base classes for extensibility

## Modular Architecture

MaterForge's architecture is organized into distinct modules:

### Core Module (`materforge.core`)
- **Material**: Schema-agnostic material container - properties are assigned dynamically
  via `setattr` and tracked automatically; no material type or composition required
- **Interfaces**: Abstract base classes for extensibility (`PropertyProcessor`,
  `DependencyResolver`, etc.)

### Parsing Module (`materforge.parsing`)
- **API**: Main entry points (`create_material`, `validate_yaml_file`, `get_material_info`)
- **Configuration**: YAML parsing and validation through `MaterialYAMLParser`
- **Processors**: Property and dependency processing (`PropertyProcessor`, `DependencyResolver`)
- **I/O**: File handling for external data (`load_property_data`)
- **Validation**: Type detection and error handling (`PropertyTypeDetector`)

### Algorithms Module (`materforge.algorithms`)
- **Interpolation**: Dependency-driven property evaluation
- **Regression**: Data simplification and fitting through `RegressionProcessor`
- **Piecewise**: Piecewise function construction via `PiecewiseBuilder`
- **Inversion**: Inverse function creation for specialized applications

### Visualization Module (`materforge.visualization`)
- **Property Plots**: Automatic visualization generation through `PropertyVisualizer`
- **Scientific Plotting**: Integration with matplotlib for high-quality plots

### Data Module (`materforge.data`)
- **Materials**: Example material configurations (Steel 1.4301, Aluminum, etc.)

## Why YAML?

YAML was chosen as the primary configuration format because:

- **Human-readable**: Easy to edit manually and understand
- **Structured**: Naturally supports nested structures required by complex material definitions
- **Ecosystem Integration**: Seamless integration with Python via ruamel.yaml
- **Reference Support**: Allows referencing previously defined scalar properties within the file
  (e.g., `solidus_temp`, `liquidus_temp`)
- **Version Control Friendly**: Text-based format works well with git and other VCS

## Integration with pystencils

MaterForge integrates with [pystencils](https://pycodegen.pages.i10git.cs.fau.de/pystencils/)
through the following workflow:

1. **Symbolic Definition**: Material properties are defined symbolically in MaterForge using
   YAML configurations
2. **Property Processing**: The parsing system converts YAML definitions into SymPy expressions
3. **Symbolic Evaluation**: Material properties can be evaluated at specific dependency values
   or kept as symbolic expressions for direct use in simulations
4. **Simulation Integration**: Symbolic expressions can be used directly in pystencils-based
   simulations

This integration allows MaterForge to leverage:
- Symbolic mathematics from SymPy for property relationships
- Dependency-driven material properties in numerical simulations
- Flexible property definitions that adapt to simulation needs

## Property Type System

MaterForge uses a property type detection system with six distinct types:

### Six Property Types

1. **CONSTANT_VALUE**: Simple numeric values for dependency-independent properties
2. **STEP_FUNCTION**: Discontinuous changes at a scalar transition point
3. **FILE_IMPORT**: Data loaded from external files (Excel, CSV, text)
4. **TABULAR_DATA**: Explicit dependency-property value pairs
5. **PIECEWISE_EQUATION**: Multiple equations for different dependency variable ranges
6. **COMPUTED_PROPERTY**: Properties calculated from other already-defined properties

### Automatic Type Detection

The `PropertyTypeDetector` automatically detects property types based on configuration structure:
- Rule-based detection with priority ordering
- Comprehensive validation for each type
- Clear error messages for invalid configurations

## Dependency Resolution

MaterForge automatically handles property dependencies:

### Dependency Analysis
- Extracts dependencies from symbolic expressions using SymPy
- Validates that all required properties are available
- Detects and prevents circular dependencies through graph analysis

### Processing Order
- Automatically determines correct processing order using topological sorting
- Processes dependencies before dependent properties
- Handles complex dependency chains transparently

## Schema-Agnostic Material Model

The `Material` class imposes no structural constraints:

- Any property name is valid - `density`, `viscosity`, `band_gap`, `yield_strength`, etc.
- All properties are assigned dynamically and tracked automatically
- No material type, composition, or element list required
- The YAML schema drives everything; the Python class is a generic container

## Error Handling Philosophy

MaterForge emphasizes clear, actionable error messages:

### Validation at Every Layer
- YAML syntax and structure validation
- Property configuration validation through `PropertyTypeDetector`
- Data quality validation in file processing
- Dependency validation with circular dependency detection

### Contextual Error Messages
- Specific error locations and suggestions for fixes
- Clear explanation of what went wrong and how to fix it
- Symbol mismatch detection: if the wrong dependency symbol is passed to `evaluate()`,
  the error reports exactly which symbol the expression requires

## Performance Optimization

### Symbolic Computation
- Efficient SymPy expression handling
- Optimized property evaluation via `Material.evaluate(symbol, value)`
- Minimal symbolic overhead for numeric evaluations

### Memory Efficiency
- Efficient data structures for large datasets
- Optional regression for memory reduction and expression simplification
- Streaming processing for large files

## Testing and Validation

### Unit Testing
- Individual component testing for all modules
- Property type validation testing
- Algorithm correctness verification

### Integration Testing
- End-to-end workflow testing
- File format compatibility testing
- Property dependency resolution testing

## Scientific Computing Integration

### SymPy Integration
- Properties as SymPy expressions enable symbolic manipulation
- Any SymPy symbol can serve as the dependency variable
- Mathematical expression validation and simplification

### NumPy Integration
- Efficient array processing for dependency variable and property data
- Vectorized operations for performance
- Seamless conversion between symbolic and numeric representations

### Matplotlib Integration
- Automatic plot generation for property visualization
- Publication-quality scientific plots
- Customizable visualization options

This design philosophy ensures MaterForge is both powerful and user-friendly, suitable for
research applications while maintaining the flexibility needed for diverse materials science
applications - from classical metals and alloys to polymers, ceramics, and beyond.
