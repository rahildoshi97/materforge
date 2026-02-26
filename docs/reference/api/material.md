# Material API Reference

## Core Classes

### Material

A dataclass representing a material with fully dynamic property tracking.
All properties - constants and dependency-driven expressions - are assigned
dynamically and tracked automatically. No material type, composition, or
fixed schema required.

```python
from materforge.core.materials import Material
```

#### Attributes

- `name` (`str`): Human-readable material identifier
- `_dynamic_properties` (`set`): Automatically tracked set of all assigned property names

#### Methods

##### `property_names() -> set`

Returns all dynamically assigned property names.

```python
props = material.property_names()
# {'density', 'heat_capacity', 'solidus_temp', ...}
```

##### `evaluate(symbol, value) -> Dict[str, float]`

Evaluates all properties by substituting `symbol = value`.

```python
def evaluate(self, symbol: sp.Symbol, value: Union[float, int]) -> Dict[str, float]:
```

**Parameters:**
- `symbol`: SymPy symbol to substitute (must match the symbol used in `create_material`)
- `value`: Numeric value to substitute - must be a number

**Returns:**
- `Dict[str, float]`: All property names mapped to evaluated float values.
  Properties that fail evaluation are excluded and logged as errors.

**Raises:**
- `ValueError`: If `symbol` is not `sp.Symbol`, or `value` is non-numeric

**Example:**
```python
import sympy as sp
from materforge.parsing.api import create_material

T = sp.Symbol('T')
mat = create_material('myAlloy.yaml', dependency=T)

# Evaluate all properties at a specific value
results = mat.evaluate(T, 500.0)
print(results['heat_capacity'])    # float

# Access symbolic expressions directly via dot notation
print(mat.density)                 # 7000.0 (constant)
print(mat.heat_conductivity)       # SymPy Piecewise expression in T
```

##### `evaluate` - removed

This method was renamed to `evaluate()`. Update any existing callsites:

```python
# Before
mat.evaluate(T, 500.0)

# After
mat.evaluate(T, 500.0)
```

#### Repr

```python
str(mat)   # Material: myAlloy (10 properties)
repr(mat)  # Material(name='myAlloy', properties=['density', 'heat_capacity', ...])
```

---

## Main API Functions

### `create_material`

Creates a `Material` instance from a YAML configuration file.

```python
from materforge.parsing.api import create_material

def create_material(
    yaml_path: Union[str, Path],
    dependency: Union[float, sp.Symbol],
    enable_plotting: bool = True,
) -> Material:
```

**Parameters:**
- `yaml_path`: Path to the YAML configuration file
- `dependency`: SymPy symbol used as the independent variable in property expressions.
  Pass a `float` to evaluate all properties immediately at that value instead.
- `enable_plotting`: Whether to generate and save property plots (default: `True`)

**Returns:**
- `Material`: Fully initialised material with all properties assigned

**Raises:**
- `ValueError`: If configuration is invalid or property processing fails

**Example:**
```python
import sympy as sp
from materforge.parsing.api import create_material

# Symbolic - properties stored as SymPy expressions
T = sp.Symbol('T')
mat = create_material('myAlloy.yaml', dependency=T)

# Any symbol works
P = sp.Symbol('P')
mat2 = create_material('pressureMaterial.yaml', dependency=P)

# Numeric - properties evaluated immediately at that value
mat3 = create_material('myAlloy.yaml', dependency=500.0)
```

---

### `validate_yaml_file`

Validates a YAML file structure without creating the material.

```python
from materforge.parsing.api import validate_yaml_file

def validate_yaml_file(yaml_path: Union[str, Path]) -> bool:
```

**Parameters:**
- `yaml_path`: Path to the YAML configuration file

**Returns:**
- `bool`: `True` if the file is structurally valid

**Raises:**
- `FileNotFoundError`: If the file does not exist
- `ValueError`: If the YAML content is invalid

**Example:**
```python
from materforge.parsing.api import validate_yaml_file

try:
    is_valid = validate_yaml_file('myAlloy.yaml')
    print(f"Valid: {is_valid}")
except ValueError as e:
    print(f"Validation error: {e}")
```

---

### `get_material_info`

Returns metadata about a YAML file without fully processing the material.

```python
from materforge.parsing.api import get_material_info

info = get_material_info('myAlloy.yaml')
print(info['name'])              # 'myAlloy'
print(info['total_properties'])  # 10
print(info['properties'])        # ['density', 'heat_capacity', ...]
print(info['property_types'])    # {'CONSTANT_VALUE': 5, 'TABULAR_DATA': 1, ...}
```

---

### `get_material_property_names`

Returns all property names on an already-created material.

```python
from materforge.parsing.api import get_material_property_names

names = get_material_property_names(mat)
# ['density', 'heat_capacity', 'heat_conductivity', ...]
```

---

### `evaluate_material_properties`

Functional wrapper around `Material.evaluate()`.

```python
from materforge.parsing.api import evaluate_material_properties

results = evaluate_material_properties(mat, T, 500.0)
```

Equivalent to `mat.evaluate(T, 500.0)`. Prefer calling `mat.evaluate()` directly.

---

## Property Type Enum

```python
from materforge.parsing.validation.property_type_detector import PropertyType

PropertyType.CONSTANT_VALUE      # single numeric value
PropertyType.STEP_FUNCTION       # discontinuous transition at a scalar reference
PropertyType.FILE_IMPORT         # data loaded from .csv / .xlsx / .txt
PropertyType.TABULAR_DATA        # explicit dependency-value pairs
PropertyType.PIECEWISE_EQUATION  # symbolic equations over dependency ranges
PropertyType.COMPUTED_PROPERTY   # derived from other properties
```

---

## Error Handling

MaterForge raises standard Python exceptions with descriptive messages:

| Situation | Exception |
|---|---|
| YAML file not found | `FileNotFoundError` |
| Invalid YAML syntax | `ruamel.yaml.scanner.ScannerError` |
| Duplicate key in YAML | `ruamel.yaml.constructor.DuplicateKeyError` |
| Missing `name` or `properties` block | `ValueError` |
| Wrong symbol passed to `evaluate()` | `ValueError` (lists required symbols) |
| Invalid property configuration | `ValueError` |
| Circular property dependency | `ValueError` |
| File import column not found | `ValueError` |

All errors include the property name and config path where the failure occurred.
