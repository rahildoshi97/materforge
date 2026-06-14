# Load Bundled Example Materials

MaterForge is built around YAML configs that you write yourself - that is the whole
point of the library, and the [material definition guide](define_materials.md) walks
through it. To make trying the API easier, the package also ships a few reference
materials you can load by name without locating the installed file.

These bundled materials are examples, not a curated database. Use them for demos,
quick experiments, and tests; for real work, define your own configs.

---

## When to Use Which

There are two ways to load a material, and they cover different needs:

- **Your own config, anywhere on disk** - pass a path to `create_material`. This
  is the normal workflow; keep the YAML wherever suits your project.

  ```python
  from pathlib import Path
  from materforge import create_material

  yaml_path = Path(__file__).parent / "my_material.yaml"
  mat = create_material(yaml_path, dependency=T)
  ```

- **One of the examples shipped with the package** - pass a name to
  `load_material` (or `get_material_path`). These helpers only see the files
  bundled inside the installed package; they will not find a YAML sitting
  elsewhere on disk, so `load_material("my_material")` raises `ValueError`.

The catalog helpers sit on top of `create_material` as a convenience - they don't
replace it. `validate_yaml_file` and `get_material_info` likewise take a path and
work with any file you point them at.

---

## Listing What Is Available

```python
from materforge import list_materials

print(list_materials())   # ['1.4301', 'Al', 'Al2O3']
```

The names are the YAML file stems: `Al` (aluminium), `Al2O3` (alumina), and
`1.4301` (the EN grade for AISI 304 stainless steel).

---

## Loading a Material by Name

`load_material` is a thin wrapper around `create_material` for the shipped files.
Matching is case-insensitive.

```python
import sympy as sp
from materforge import load_material

T = sp.Symbol('T')
steel = load_material('1.4301', T)

print(steel)              # Material: 1.4301 (13 properties)
print(steel.heat_capacity)  # SymPy expression in T

# Evaluate at a temperature, as with any material
at_800 = steel.evaluate(T, 800.0)
print(float(at_800.heat_capacity))
```

Plotting is off by default here (unlike `create_material`), so loading an example
never writes plot files next to the installed package. Pass `enable_plotting=True`
if you want the plots.

---

## Getting the File Path

When you would rather read or copy the YAML - for instance, to use a bundled file
as the starting point for your own - ask for its path:

```python
from materforge import get_material_path

path = get_material_path('Al')
print(path)               # .../materforge/data/materials/Al.yaml
print(path.read_text())   # the raw YAML
```

An unknown name raises `ValueError` and lists what is available:

```python
get_material_path('inconel')
# ValueError: Unknown material 'inconel'. Available materials: 1.4301, Al, Al2O3
```

---

## Next Steps

- [Define your own material properties](define_materials.md) - all six property types
- [YAML schema reference](../reference/yaml_schema.md)
