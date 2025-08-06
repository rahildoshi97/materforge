# Material API Reference

## Core Classes

### Material

A dataclass representing a material (pure metal or alloy) with temperature-dependent properties.

```python
from materforge.core.materials import Material
```

#### Properties

**Basic Properties:**
- `name`: String identifier for the material
- `material_type`: Either 'pure_metal' or 'alloy'
- `elements`: List of ChemicalElement objects
- `composition`: List/array of element fractions (must sum to 1.0)

**Temperature Properties (Pure Metals):**
- `melting_temperature`: Melting point in Kelvin
- `boiling_temperature`: Boiling point in Kelvin

**Temperature Properties (Alloys):**
- `solidus_temperature`: Solidus temperature in Kelvin
- `liquidus_temperature`: Liquidus temperature in Kelvin
- `initial_boiling_temperature`: Initial boiling temperature in Kelvin
- `final_boiling_temperature`: Final boiling temperature in Kelvin

**Computed Properties:**
- `atomic_mass`: Composition-weighted atomic mass
- `atomic_number`: Composition-weighted atomic number

**Material Properties (Optional):**
- `density`: Density as function of temperature
- `dynamic_viscosity`: Dynamic viscosity as function of temperature
- `energy_density`: Energy density as function of temperature
- `heat_capacity`: Specific heat capacity as function of temperature
- `heat_conductivity`: Thermal conductivity as function of temperature
- `kinematic_viscosity`: Kinematic viscosity as function of temperature
- `latent_heat_of_fusion`: Latent heat of fusion
- `latent_heat_of_vaporization`: Latent heat of vaporization
- `specific_enthalpy`: Specific enthalpy as function of temperature
- `surface_tension`: Surface tension as function of temperature
- `thermal_diffusivity`: Thermal diffusivity as function of temperature
- `thermal_expansion_coefficient`: Thermal expansion coefficient

#### Methods

**solidification_interval()**
```python
def solidification_interval(self) -> Tuple[sp.Float, sp.Float]:
```
Returns the solidification interval (solidus, liquidus) for alloys.

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
print(f"Type: {material.material_type}")
print(f"Composition: {dict(zip([e.name for e in material.elements], material.composition))}")

# Access temperature-dependent properties
if hasattr(material, 'density'):
    density_at_500K = material.evaluate_properties_at_temperature(500.0)
    print(f"Density at 500K: {density_at_500K} kg/m³")

# For alloys, get solidification interval
if material.material_type == 'alloy':
    solidus, liquidus = material.solidification_interval()
    print(f"Solidification range: {solidus}K - {liquidus}K")
```

### ChemicalElement

A dataclass representing a chemical element with its properties.
````python
from materforge.core.elements import ChemicalElement
````

#### Properties

- `name`: Element name (e.g., "Iron")
- `atomic_number`: Atomic number
- `atomic_mass`: Atomic mass in u
- `melting_temperature`: Melting temperature in K
- `boiling_temperature`: Boiling temperature in K
- `latent_heat_of_fusion`: Latent heat of fusion in J/kg
- `latent_heat_of_vaporization`: Latent heat of vaporization in J/kg

#### Example Usage
```python
from materforge.data.elements.element_data import element_map

# Access element data
iron = element_map['Fe']
print(f"Iron melting point: {iron.melting_temperature}K")
print(f"Iron atomic mass: {iron.atomic_mass}u")
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

### get_supported_properties

Returns a list of all supported material properties.
```python
from materforge.parsing.api import get_supported_properties

def get_supported_properties() -> list:
```
**Returns:**
- `list`: List of strings representing valid property names

**Example:**
```python
properties = get_supported_properties()
print("Supported properties:", properties)
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

## Utility Functions

### Element Interpolation Functions
```python
from materforge.core.elements import (
interpolate_atomic_number,
interpolate_atomic_mass,
interpolate_melting_temperature,
interpolate_boiling_temperature
)
```

These functions interpolate element properties based on composition:

```python
# Example usage

elements = [element_map['Fe'], element_map['C']]
composition = [0.98, 0.02]

avg_atomic_mass = interpolate_atomic_mass(elements, composition)
avg_melting_temp = interpolate_melting_temperature(elements, composition)
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
