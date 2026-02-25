# Material API Reference

## Core Classes

### Material

A dataclass representing a material with constant or dependent properties.

```python
from materforge.core.materials import Material
```

#### Properties

**Basic Properties:**
- `name`: String identifier for the material
- `properties`: Dict of user-defined material properties

#### Example Usage
```python
import sympy as sp
from materforge.parsing.api import create_material

# Create symbolic temperature
T = sp.Symbol('T')

# Load material from YAML
material = create_material('1.4301.yaml', T)

# Access basic properties
print(f"Material: {material.name}")

# Access temperature-dependent properties
if hasattr(material, 'density'):
    density_at_500K = material.evaluate_properties_at_temperature(500.0)
    print(f"Density at 500K: {density_at_500K} kg/m³")
```

## Main API Functions

### create_material

Create material instance from YAML configuration file.
```python
from materforge.parsing.api import create_material

def create_material(yaml_path: Union[str, Path],
                    T: Union[float, sp.Symbol],
                    enable_plotting: bool = True) -> Material:
```

**Parameters:**
- `yaml_path`: Path to the YAML configuration file
- `T`: Temperature value or symbol for property evaluation
    - Use a float value for a specific temperature
    - Use a symbolic variable (e.g., `sp.Symbol('T')`) for symbolic expressions
- `enable_plotting`: Whether to generate visualization plots (default: True)

**Returns:**
- `Material`: The material instance with all properties initialized

**Example:**
```python
import sympy as sp
from materforge.parsing.api import create_material

# Create material with symbolic temperature
T = sp.Symbol('T')
material = create_material('1.4301.yaml', T)

# Create material with custom temperature symbol
u_C = sp.Symbol('u_C')
material = create_material('copper.yaml', u_C)
```

### validate_yaml_file

Validate a YAML file without creating the material.
```python
from materforge.parsing.api import validate_yaml_file

def validate_yaml_file(yaml_path: Union[str, Path]) -> bool:
```

**Parameters:**
- `yaml_path`: Path to the YAML configuration file to validate

**Returns:**
- `bool`: True if the file is valid

**Raises:**
- `FileNotFoundError`: If the file doesn't exist
- `ValueError`: If the YAML content is invalid

**Example:**
```python
try:
  is_valid = validate_yaml_file('material.yaml')
  print(f"File is valid: {is_valid}")
  except ValueError as e:
  print(f"Validation error: {e}")
```

## Symbol Registry

### SymbolRegistry

Registry for SymPy symbols to ensure uniqueness across the application.
```python
from materforge.core.symbol_registry import SymbolRegistry

# Get or create a symbol
T = SymbolRegistry.get('T')

# Get all registered symbols
all_symbols = SymbolRegistry.get_all()

# Clear all symbols (useful for testing)
SymbolRegistry.clear()
```

## Error Classes

### Material Errors

```python
from materforge.core.materials import MaterialCompositionError, MaterialTemperatureError
```
These are raised automatically during material validation

### Property Errors

```python
from materforge.parsing.validation.errors import (
  PropertyError,
  DependencyError,
  CircularDependencyError
)
```
These are raised during property processing

## Type Definitions

### PropertyType Enum

```python
from materforge.parsing.validation.property_type_detector import PropertyType

# Available property types:
PropertyType.CONSTANT_VALUE
PropertyType.STEP_FUNCTION
PropertyType.FILE_IMPORT
PropertyType.TABULAR_DATA
PropertyType.PIECEWISE_EQUATION
PropertyType.COMPUTED_PROPERTY
```

This API provides a comprehensive interface for working with materials in MaterForge,
from basic material creation to advanced property manipulation and validation.
