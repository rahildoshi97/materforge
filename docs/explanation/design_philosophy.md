# Design Philosophy of MaterForge

This document explains the core design principles, architectural decisions, and the rationale behind MaterForge's structure and implementation.

## Core Principles

MaterForge is built upon several core principles:

- **Modularity**: Clearly separated components for ease of maintenance, testing, and extensibility
- **Flexibility**: Allow users to define material properties in various intuitive ways
- **Performance**: Leverage symbolic computation for high-performance simulations
- **Transparency and Reproducibility**: Clearly document material property definitions and computations to ensure reproducibility

## Layered Architecture

MaterForge follows a layered architecture to separate concerns clearly:

### 1. User Interface Layer (YAML Configuration)

- Provides a simple, human-readable format for defining materials and their properties
- Allows users to specify properties using multiple intuitive methods (constants, interpolation points, file-based data, computed properties)
- Ensures clarity, readability, and ease of use

### 2. Parsing Layer (Python)

- **Configuration Processing**: Handles YAML parsing and validation through `MaterialConfigParser`
- **Property Type Detection**: Automatically determines property definition types using `PropertyConfigAnalyzer`
- **Dependency Resolution**: Processes various dependency definition formats via `DependencyResolver`
- **Data Handling**: Manages file I/O for external data sources through `read_data_from_file`

### 3. Symbolic Representation Layer (SymPy)

- Uses symbolic mathematics (via SymPy) internally to represent material properties
- Enables symbolic manipulation, simplification, and validation of property definitions
- Facilitates automatic computation of derived properties through symbolic expressions

### 4. Algorithms Layer (Python)

- **Interpolation**: Robust interpolation methods for evaluating temperature-dependent properties
- **Regression**: Data simplification and piecewise function generation via `RegressionManager`
- **Piecewise Functions**: Creating piecewise expressions through `PiecewiseBuilder`
- **Inversion**: Creating inverse functions for specialized applications

### 5. Visualization Layer (Python)

- Automatic plot generation for material properties through `PropertyVisualizer`
- Property visualization and validation
- Integration with matplotlib for scientific plotting

### 6. Core Layer (Python)

- **Material**: Fundamental material representation with validation
- **ChemicalElement**: Element data and properties
- **Interfaces**: Abstract base classes for extensibility
- **Symbol Registry**: SymPy symbol management for consistency

## Modular Architecture

MaterForge's architecture is organized into distinct modules:

### Core Module (`materforge.core`)
- **Material**: Fundamental material representation with composition validation
- **ChemicalElement**: Element data and properties from periodic table
- **Interfaces**: Abstract base classes for extensibility (`PropertyProcessor`, `DependencyResolver`, etc.)
- **Symbol Registry**: SymPy symbol management to ensure uniqueness

### Parsing Module (`materforge.parsing`)
- **API**: Main entry points (`create_material`, `validate_yaml_file`, `get_supported_properties`)
- **Configuration**: YAML parsing and validation through `MaterialConfigParser`
- **Processors**: Property and dependency processing (`PropertyManager`, `DependencyResolver`)
- **I/O**: File handling for external data (`read_data_from_file`)
- **Validation**: Type detection and error handling (`PropertyConfigAnalyzer`)

### Algorithms Module (`materforge.algorithms`)
- **Interpolation**: Temperature-dependent property evaluation
- **Regression**: Data simplification and fitting through `RegressionManager`
- **Piecewise**: Piecewise function construction via `PiecewiseBuilder`
- **Inversion**: Inverse function creation for specialized applications

### Visualization Module (`materforge.visualization`)
- **Property Plots**: Automatic visualization generation through `PropertyVisualizer`
- **Scientific Plotting**: Integration with matplotlib for publication-quality plots

### Data Module (`materforge.data`)
- **Elements**: Chemical element database with periodic table data
- **Constants**: Physical and processing constants
- **Materials**: Pre-defined material configurations (Steel 1.4301, aluminum, etc.)

## Why YAML?

YAML was chosen as the primary configuration format because:

- **Human-readable**: Easy to edit manually and understand
- **Structured**: Naturally supports nested structures required by complex material definitions
- **Ecosystem Integration**: Seamless integration with Python via ruamel.yaml
- **Reference Support**: Allows referencing previously defined variables within the file (e.g., `solidus_temperature`, `liquidus_temperature`)
- **Version Control Friendly**: Text-based format works well with git and other VCS

## Integration with pystencils

MaterForge integrates with [pystencils](https://pycodegen.pages.i10git.cs.fau.de/pystencils/) through the following workflow:

1. **Symbolic Definition**: Material properties are defined symbolically in MaterForge using YAML configurations
2. **Property Processing**: The parsing system converts YAML definitions into SymPy expressions
3. **Symbolic Evaluation**: Material properties can be evaluated at specific temperatures or kept as symbolic expressions
4. **Simulation Integration**: Symbolic expressions can be used directly in pystencils-based simulations

This integration allows MaterForge to leverage:
- Symbolic mathematics from SymPy for property relationships
- Temperature-dependent material properties in numerical simulations
- Flexible property definitions that adapt to simulation needs

## Property Type System

MaterForge uses a sophisticated property type detection system with six distinct types:

### Six Property Types

1. **CONSTANT_VALUE**: Simple numeric values for temperature-independent properties
2. **STEP_FUNCTION**: Discontinuous changes at phase transitions
3. **FILE_IMPORT**: Data loaded from external files (Excel, CSV, text)
4. **TABULAR_DATA**: Explicit temperature-property pairs
5. **PIECEWISE_EQUATION**: Multiple equations for different temperature ranges
6. **COMPUTED_PROPERTY**: Properties calculated from other properties

### Automatic Type Detection

The `PropertyConfigAnalyzer` automatically detects property types based on configuration structure:
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

## Extensibility Framework

MaterForge is designed for extensibility through abstract interfaces:

### Abstract Interfaces
- `PropertyProcessor`: For custom property processing logic
- `DependencyResolver`: For custom dependency handling
- `DataHandler`: For custom file formats
- `Visualizer`: For custom visualization approaches

### Plugin Architecture
- Users can extend functionality without modifying core code
- New property types can be added through the type detection system
- Custom algorithms can be integrated through the algorithms module

## Error Handling Philosophy

MaterForge emphasizes clear, actionable error messages:

### Validation at Every Layer
- YAML syntax and structure validation
- Property configuration validation through `PropertyConfigAnalyzer`
- Data quality validation in file processing
- Dependency validation with circular dependency detection

### Contextual Error Messages
- Specific error locations and suggestions for fixes
- Available options and corrections
- Clear explanation of what went wrong and how to fix it

## Performance Optimization

Performance-critical operations are optimized:

### Symbolic Computation
- Efficient SymPy expression handling
- Optimized property evaluation methods
- Minimal symbolic overhead for numeric evaluations

### Memory Efficiency
- Efficient data structures for large datasets
- Optional regression for memory reduction
- Streaming processing for large files

### Caching and Optimization
- Symbol registry prevents duplicate symbol creation
- Efficient temperature array processing
- Optimized interpolation algorithms

## Testing and Validation

MaterForge includes comprehensive testing:

### Unit Testing
- Individual component testing for all modules
- Property type validation testing
- Algorithm correctness verification

### Integration Testing
- End-to-end workflow testing
- File format compatibility testing
- Property dependency resolution testing

### Validation Testing
- Physical property validation against known values
- Monotonicity checking for energy density
- Composition validation for materials

## Scientific Computing Integration

MaterForge is designed to integrate seamlessly with the scientific Python ecosystem:

### SymPy Integration
- Properties as SymPy expressions enable symbolic manipulation
- Automatic differentiation and integration capabilities
- Mathematical expression validation and simplification

### NumPy Integration
- Efficient array processing for temperature and property data
- Vectorized operations for performance
- Seamless conversion between symbolic and numeric representations

### Matplotlib Integration
- Automatic plot generation for property visualization
- Publication-quality scientific plots
- Customizable visualization options

This design philosophy ensures MaterForge is both powerful and user-friendly, suitable for research applications while maintaining the flexibility needed for diverse materials science applications.
