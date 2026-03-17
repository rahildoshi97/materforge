# MaterForge - Materials Formulation Engine with Python

A high-performance Python library for materials modeling and simulation. MaterForge enables
efficient definition of any material - metals, alloys, polymers, ceramics, composites, or
hypothetical materials - through YAML configuration files, providing symbolic and numerical
property evaluation as a function of any SymPy symbol.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Latest Release](https://i10git.cs.fau.de/rahil.doshi/materforge/-/badges/release.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/releases)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE)
[![Documentation Status](https://readthedocs.org/projects/materforge/badge/?version=latest)](https://materforge.readthedocs.io/)
[![Pipeline Status](https://i10git.cs.fau.de/rahil.doshi/materforge/badges/master/pipeline.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/pipelines)
[![Code Coverage](https://i10git.cs.fau.de/rahil.doshi/materforge/badges/master/coverage.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/commits/master)

**Documentation:** [https://materforge.readthedocs.io](https://materforge.readthedocs.io)

## Table of Contents
- [Key Features](#-key-features)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [YAML Configuration Format](#-yaml-configuration-format)
- [Documentation](#-documentation)
- [Contributing](#-contributing)
- [Known Limitations](#-known-limitations)
- [License](#-license)
- [Citation](#-citation)
- [Support](#-support)
- [Acknowledgments](#-acknowledgments)

---

## 🚀 Key Features

- **Schema-Agnostic Material Definition**: Any material kind, any property name - the YAML
  schema drives everything with no hardcoded material types or required fields
- **Dependency-Driven Properties**: Properties expressed as symbolic functions of any
  independent variable - temperature, pressure, composition, strain, or any SymPy symbol
- **Modular Architecture**: Clean separation with algorithms, parsing, and visualization modules
- **Symbolic Mathematics**: Built on SymPy for precise mathematical expressions
- **Piecewise Functions**: Advanced piecewise function support with regression capabilities
- **Property Inversion**: Create inverse functions for any piecewise-linear property
- **Visualization**: Automatic plotting of material properties with customizable options
- **Multiple Property Types**: Constants, step functions, file-based data, tabular data,
  piecewise equations, and computed properties
- **Regression Analysis**: Built-in piecewise linear fitting with configurable parameters

---

## 📦 Installation

### Install from PyPI

```bash
pip install materforge
```

All required dependencies are installed automatically.

### For contributors

Clone the repository and install in editable mode:

```bash
git clone https://i10git.cs.fau.de/rahil.doshi/materforge.git
cd materforge
pip install -e .
```

**Note:** The `-e` flag (`--editable`) installs the package by symlinking directly to the
source directory. Changes to source files take effect immediately without reinstalling -
intended for development only. Use `pip install materforge` for regular use.

---

## 🏃 Quick Start

### Basic Material Creation

```python
import sympy as sp
from materforge import create_material

# Any SymPy symbol works as the dependency variable
T = sp.Symbol('T')

# Load a material from YAML (see examples/myAlloy.yaml)
mat = create_material('examples/myAlloy.yaml', dependency=T, enable_plotting=True)

# Access symbolic property expressions directly
print(mat.heat_capacity)   # SymPy Piecewise expression in T
print(mat.density)         # 7000.0 (constant float)

# Evaluate all properties at a specific value
results = mat.evaluate(T, 500.0)
print(results['heat_capacity'])   # float
```

### Property Inversion

```python
from materforge.algorithms.piecewise_inverter import PiecewiseInverter
import sympy as sp
from materforge import create_material

T = sp.Symbol('T')
mat = create_material('examples/myAlloy.yaml', dependency=T)

if hasattr(mat, 'energy_density'):
    E = sp.Symbol('E')
    inverse = PiecewiseInverter.create_inverse(mat.energy_density, 'T', 'E')

    # Validate round-trip accuracy
    test_val = 500.0
    e = float(mat.energy_density.subs(T, test_val))
    recovered = float(inverse.subs(E, e))
    print(f"Round-trip: T={test_val} -> E={e:.2e} -> T={recovered:.2f}")
```

---

## 📋 YAML Configuration Format

### Supported Property Types
- **CONSTANT_VALUE**: Single numeric value - not dependency-driven
- **STEP_FUNCTION**: Discontinuous transition at a scalar reference point
- **FILE_IMPORT**: Data loaded from CSV, Excel, or text files
- **TABULAR_DATA**: Explicit dependency-value pairs
- **PIECEWISE_EQUATION**: Symbolic equations over dependency variable ranges
- **COMPUTED_PROPERTY**: Properties derived from other defined properties

The symbol used in YAML equations (`T`) is a **placeholder only** - MaterForge substitutes
it with whatever symbol you pass to `create_material(..., dependency=symbol)` at runtime.

See the [YAML schema documentation](https://materforge.readthedocs.io/en/latest/reference/yaml_schema.html)
for full configuration options.

**Example YAML files:**
- [Example alloy](https://github.com/rahildoshi97/materforge/blob/main/examples/myAlloy.yaml)
- [Aluminum](https://github.com/rahildoshi97/materforge/blob/main/src/materforge/data/materials/pure_metals/Al/Al.yaml)
- [Steel 1.4301](https://github.com/rahildoshi97/materforge/blob/main/src/materforge/data/materials/alloys/1.4301/1.4301.yaml)

---

## 📚 Documentation

Full documentation is available at **https://materforge.readthedocs.io**

The documentation follows the [Diátaxis](https://diataxis.fr/) framework:

| Type | Content |
|------|---------|
| **Tutorials** | [Getting Started](https://materforge.readthedocs.io/en/latest/tutorials/getting_started.html) · [First Simulation](https://materforge.readthedocs.io/en/latest/tutorials/first_simulation.html) |
| **How-to Guides** | [Defining Material Properties](https://materforge.readthedocs.io/en/latest/how-to/define_materials.html) · [Property Inversion](https://materforge.readthedocs.io/en/latest/how-to/property_inversion.html) |
| **Reference** | [API Reference](https://materforge.readthedocs.io/en/latest/reference/api.html) · [YAML Schema](https://materforge.readthedocs.io/en/latest/reference/yaml_schema.html) |
| **Explanation** | [Design Philosophy](https://materforge.readthedocs.io/en/latest/explanation/design_philosophy.html) · [Material Properties](https://materforge.readthedocs.io/en/latest/explanation/material_properties.html) |

---

## 🤝 Contributing

Contributions are welcome! Please see our
[Contributing Guide](https://github.com/rahildoshi97/materforge/blob/main/CONTRIBUTING.md)
for details on how to get started.

---

## 🐛 Known Limitations

- **Piecewise Inverter**: Currently supports linear piecewise functions only (`degree: 1`)
- **File Formats**: Limited to CSV, Excel, and text files
- **Memory Usage**: Large datasets may require optimization for very high-resolution data
- **Regression**: Maximum 8 segments recommended for stability

---

## 📄 License

### Core Library (BSD-3-Clause)

The MaterForge library itself (`src/materforge/`, `examples/`, `tests/`, `docs/`) is licensed
under the **BSD 3-Clause License**. See the
[LICENSE](https://github.com/rahildoshi97/materforge/blob/main/LICENSE) file for full details.

### Application Examples (GPL-3.0-or-later)

The `apps/` directory contains demonstration applications that integrate MaterForge with
[waLBerla](https://i10git.cs.fau.de/walberla/walberla) and
[pystencils](https://pypi.org/project/pystencils/). Because these dependencies are
GPLv3-licensed, the apps directory is licensed under **GPL-3.0-or-later**. See
[apps/LICENSE](https://github.com/rahildoshi97/materforge/blob/main/apps/LICENSE) for
full details.

### PyPI Distribution

`pip install materforge` includes **only the BSD-3-Clause licensed core library**. The
GPL-licensed apps are excluded from the PyPI distribution.

| Component | Location | License | In PyPI |
|-----------|----------|---------|---------|
| Core library | `src/materforge/` | BSD-3-Clause | ✅ |
| Example scripts | `examples/` | BSD-3-Clause | ❌ |
| Tests | `tests/` | BSD-3-Clause | ❌ |
| Documentation | `docs/` | BSD-3-Clause | ❌ |
| Apps | `apps/` | GPL-3.0-or-later | ❌ |

---

## 📖 Citation

If you use MaterForge in your research, please cite it using the information in our
[CITATION.cff](https://github.com/rahildoshi97/materforge/blob/main/CITATION.cff) file.

---

## 📞 Support

- **Author**: Rahil Doshi
- **Email**: [rahil.doshi@fau.de](mailto:rahil.doshi@fau.de)
- **Documentation**: [materforge.readthedocs.io](https://materforge.readthedocs.io)
- **Bug Tracker**: [GitHub Issues](https://github.com/rahildoshi97/materforge/issues)
- **GitLab**: [i10git.cs.fau.de](https://i10git.cs.fau.de/rahil.doshi/materforge)

---

## 🙏 Acknowledgments

- Built with [SymPy](https://www.sympy.org/) for symbolic mathematics
- Data handling powered by [pandas](https://pandas.pydata.org/)
- Uses [pwlf](https://github.com/cjekel/piecewise_linear_fit_py) for piecewise linear fitting
- Visualization powered by [Matplotlib](https://matplotlib.org/)
- YAML parsing with [ruamel.yaml](https://yaml.dev/doc/ruamel.yaml/)

---

*MaterForge - Empowering materials simulation with Python*
