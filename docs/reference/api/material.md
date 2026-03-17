# Material API Reference

## Core Classes

### Material

A dataclass representing a material with fully dynamic property tracking.
All properties - constants and dependency-driven expressions - are assigned
dynamically and tracked automatically via a `properties` dict. No material
type, composition, or fixed schema required.

```python
from materforge.core.materials import Material
```

#### Methods

##### `property_names() -> set`

Returns all dynamically assigned property names.

```python
props = material.property_names()
# {'density', 'heat_capacity', 'solidus_temperature', ...}
```

##### `evaluate(symbol, value) -> Material`

Evaluates all properties by substituting `symbol = value`. Returns a **new
`Material`** instance with all expressions reduced to numeric SymPy values.
Properties that still have free symbols after substitution, or that fail
evaluation, are excluded and logged as errors.

```python
def evaluate(self, symbol: sp.Symbol, value: Union[float, int]) -> Material:
```

**Parameters:**
- `symbol`: SymPy symbol to substitute - must match the symbol passed to
  `create_material()`. Must be `sp.Symbol`.
- `value`: Numeric value to substitute - must be convertible to `float`

**Returns:**
- A new `Material` named `"{name}@{symbol}={value}"` with all properties
  evaluated to numeric SymPy scalars (`sp.Float`)

**Raises:**
- `ValueError`: If `symbol` is not `sp.Symbol`
- `ValueError`: If `value` is `None` or not convertible to `float`

**Example:**
```python
import sympy as sp
from materforge import create_material

T = sp.Symbol('T')
mat = create_material('myAlloy.yaml', dependency=T)

# Evaluate all properties at T = 500 K - returns a new Material
mat_at_500 = mat.evaluate(T, 500.0)
print(mat_at_500.heat_capacity)           # sp.Float numeric value
print(float(mat_at_500.heat_capacity))    # Python float if needed

# Access symbolic expressions directly via dot notation on the original
print(mat.density)            # sp.Float(7000.0) for a constant
print(mat.heat_conductivity)  # SymPy Piecewise expression in T
```

> **pystencils users:** Do not pass `u.center()` to `create_material()`.
> Always create the material with a plain `sp.Symbol`, then substitute
> to the field accessor at the pystencils boundary:
> ```python
> T = sp.Symbol('T')
> mat = create_material('myAlloy.yaml', dependency=T)
> expr = mat.thermal_diffusivity.subs(T, u.center())   # one explicit coupling point
> ```

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
from materforge import create_material

def create_material(
    yaml_path: Union[str, Path],
    dependency: sp.Symbol,
    enable_plotting: bool = True,
) -> Material:
```

**Parameters:**
- `yaml_path`: Path to the YAML configuration file
- `dependency`: SymPy symbol used as the independent variable in all property
  expressions. Must be a plain `sp.Symbol`. pystencils field accessors such
  as `u.center()` are not valid here - apply `.subs()` after creation instead.
- `enable_plotting`: Whether to generate and save property plots (default: `True`)

**Returns:**
- `Material`: Fully initialised material with all properties assigned

**Raises:**
- `FileNotFoundError`: If the YAML file does not exist
- `TypeError`: If `dependency` is not a `sp.Symbol` instance
- `MaterialConfigError`: If the YAML content is invalid or property processing fails
- `PropertyConfigError`: If a specific property block is structurally invalid

**Example:**
```python
import sympy as sp
from materforge import create_material

T = sp.Symbol('T')
mat = create_material('myAlloy.yaml', dependency=T)
mat_no_plots = create_material('myAlloy.yaml', dependency=T, enable_plotting=False)
```

---

### `validate_yaml_file`

Validates a YAML file structure without creating the material.

```python
from materforge import validate_yaml_file

def validate_yaml_file(yaml_path: Union[str, Path]) -> bool:
```

**Parameters:**
- `yaml_path`: Path to the YAML configuration file

**Returns:**
- `bool`: `True` if the file is structurally valid

**Raises:**
- `FileNotFoundError`: If the file does not exist
- `MaterialConfigError`: If the top-level YAML structure is invalid
- `PropertyConfigError`: If a specific property block is structurally invalid

**Example:**
```python
from materforge import validate_yaml_file
from materforge.parsing.validation.errors import MaterialConfigError, PropertyConfigError

try:
    is_valid = validate_yaml_file('myAlloy.yaml')
    print(f"Valid: {is_valid}")
except PropertyConfigError as e:
    print(f"Property config error: {e}")
except MaterialConfigError as e:
    print(f"Config error: {e}")
```

---

### `get_material_info`

Returns metadata about a YAML file without fully processing the material.

```python
from materforge import get_material_info

info = get_material_info('myAlloy.yaml')
print(info['name'])              # 'myAlloy'
print(info['total_properties'])  # 10
print(info['properties'])        # ['density', 'heat_capacity', ...]
print(info['property_types'])    # {'CONSTANT_VALUE': 5, 'TABULAR_DATA': 1, ...}
```

---

### `get_material_property_names`

Returns all property names on an already-created material. Equivalent to
`mat.property_names()` - prefer calling `mat.property_names()` directly.

```python
from materforge import get_material_property_names

names = get_material_property_names(mat)   # returns a set
```

**Raises:**
- `TypeError`: If `material` is not a `Material` instance

---

### `evaluate_material_properties`

Functional wrapper around `Material.evaluate()`. Equivalent to
`mat.evaluate(symbol, value)` - prefer calling `mat.evaluate()` directly.

```python
from materforge import evaluate_material_properties

mat_at_500 = evaluate_material_properties(mat, T, 500.0)
```

**Raises:**
- `TypeError`: If `material` is not a `Material` instance

---

## YAML Placeholder Contract

All YAML equation strings **must** use `T` as the dependency variable.
`T` is the fixed placeholder (`YAML_PLACEHOLDER = sp.Symbol('T')`),
hardcoded in materforge. At runtime, `T` is substituted with whatever
`sp.Symbol` the caller passes to `create_material()`.

```yaml
# Correct - always use T in YAML equations
specific_enthalpy:
  dependency: (1773, 293, 541)
  equation: Integral(heat_capacity, T)
  bounds: [constant, constant]

thermal_diffusivity:
  dependency: (1773, 293, 100)
  equation: thermal_conductivity / (density * specific_heat)
  bounds: [constant, constant]
```

```yaml
# Wrong - any symbol other than T in the equation will cause an error
specific_enthalpy:
  equation: Integral(heat_capacity, u_C)   # crashes unconditionally
```

When YAML uses any symbol other than `T`, the exact exception depends on
property type:

- **`COMPUTED_PROPERTY`**: raises `DependencyError` - `sp.sympify()` parses
  the non-`T` symbol as a plain `sp.Symbol`. The dependency filter
  (`if s != YAML_PLACEHOLDER`) does not exclude it, so it lands in the
  dependency list. Since `u_C` (or any other name) is not a material property,
  materforge raises:
  ```
  DependencyError: Missing dependencies ['u_C'] in expression
      'Integral(heat_capacity, u_C)': u_C
  Available properties: density, heat_capacity, ...
  Please check for typos or add the missing properties to your configuration.
  ```

- **`PIECEWISE_EQUATION`**: raises `ValueError` - equations in this type are
  validated upfront and may only reference `T`. Any other symbol is rejected
  immediately:
  ```
  ValueError: Unexpected symbol(s) [u_C] in equation '7877.39-0.37*u_C'
      for property 'viscosity'. PIECEWISE_EQUATION equations must use 'T'
      as the only variable.
  ```
  To reference other material properties in an expression, use
  `COMPUTED_PROPERTY` instead.

Both errors occur regardless of what symbol the caller passes to
`create_material()`. The YAML and the Python caller are fully decoupled -
the YAML always uses `T`, and the caller can use any `sp.Symbol`.

---

## Property Type Enum

```python
from materforge.parsing.validation.property_type_detector import PropertyType

PropertyType.CONSTANT_VALUE      # single numeric value
PropertyType.STEP_FUNCTION       # discontinuous transition at a scalar reference
PropertyType.FILE_IMPORT         # data loaded from .csv / .xlsx / .txt
PropertyType.TABULAR_DATA        # explicit dependency-value pairs
PropertyType.PIECEWISE_EQUATION  # symbolic equations over dependency ranges
PropertyType.COMPUTED_PROPERTY   # derived from other properties via equation
```

Note: `PIECEWISE_EQUATION` equations may only reference `T` (the YAML
placeholder). To combine multiple material properties in an expression, use
`COMPUTED_PROPERTY` instead.

---

## Error Handling

MaterForge raises standard Python exceptions with descriptive messages that
include the property name and config path where the failure occurred.

| Situation | Exception |
|---|---|
| YAML file not found | `FileNotFoundError` |
| Invalid YAML syntax | `ruamel.yaml.scanner.ScannerError` |
| Duplicate key in YAML | `ruamel.yaml.constructor.DuplicateKeyError` |
| Unknown top-level YAML key | `MaterialConfigError` |
| Missing `name` or `properties` block | `MaterialConfigError` |
| `properties` block is empty | `MaterialConfigError` |
| `dependency` not `sp.Symbol` in `create_material()` | `TypeError` |
| `material` wrong type in API helpers | `TypeError` |
| `symbol` not `sp.Symbol` in `evaluate()` | `ValueError` |
| `value` non-numeric in `evaluate()` | `ValueError` |
| Invalid property configuration | `PropertyConfigError` |
| Non-`T` symbol in `PIECEWISE_EQUATION` | `ValueError` |
| Non-`T` symbol in `COMPUTED_PROPERTY` | `DependencyError` |
| Missing property dependency | `DependencyError` |
| Circular property dependency | `CircularDependencyError` |
| File import column not found | `ValueError` |

`MaterialConfigError`, `PropertyConfigError`, `DependencyError`, and
`CircularDependencyError` are importable from
`materforge.parsing.validation.errors`.
