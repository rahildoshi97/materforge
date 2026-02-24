# MaterForge - Materials Formulation Engine with Python

A high-performance Python library for material simulation and analysis. MaterForge enables efficient modeling of pure metals and alloys through YAML configuration files, providing symbolic and numerical property evaluation for various material properties.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Latest Release](https://i10git.cs.fau.de/rahil.doshi/materforge/-/badges/release.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/releases)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://github.com/rahildoshi97/materforge/blob/main/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/materforge/badge/?version=latest)](https://materforge.readthedocs.io/)
[![Pipeline Status](https://i10git.cs.fau.de/rahil.doshi/materforge/badges/master/pipeline.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/pipelines)
[![Code Coverage](https://i10git.cs.fau.de/rahil.doshi/materforge/badges/master/coverage.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/commits/master)

**Documentation:** https://materforge.readthedocs.io

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

## 🚀 Key Features
- **Modular Architecture**: Clean separation with algorithms, parsing, and visualization modules
- **Flexible Material Definition**: Support for both pure metals and alloys
- **YAML-Driven Configuration**: Define materials using intuitive YAML files
- **Temperature-Dependent Properties**: Support for complex temperature-dependent material properties
- **Symbolic Mathematics**: Built on SymPy for precise mathematical expressions
- **Piecewise Functions**: Advanced piecewise function support with regression capabilities
- **Property Inversion**: Create inverse functions for energy density and other properties
- **Visualization**: Automatic plotting of material properties with customizable options
- **Multiple Property Types**: Constants, step functions, file-based data, key-value pairs,
  piecewise equations, and computed properties
- **Regression Analysis**: Built-in piecewise linear fitting with configurable parameters

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

**Note:** The `-e` flag (`--editable`) installs the package by symlinking directly to the source directory.
Changes to source files take effect immediately without reinstalling - intended for development only. 
Use `pip install materforge` for regular use.

## 🏃 Quick Start

### Basic Material Creation

```python
import sympy as sp
from materforge.parsing.api import create_material

# Create a material with symbolic temperature
T = sp.Symbol('T')

# Load an example material (see examples/myAlloy.yaml)
material = create_material('examples/myAlloy.yaml', T, enable_plotting=True)

# Access symbolic property expressions
print(f"Heat capacity: {material.heat_capacity}")
print(f"Density:       {material.density}")
```

### Working with Piecewise Inverse Functions

```python
from materforge.algorithms.piecewise_inverter import PiecewiseInverter
import sympy as sp

T = sp.Symbol('T')
material = create_material('examples/myAlloy.yaml', T)

if hasattr(material, 'energy_density'):
    E = sp.Symbol('E')
    inverse_func = PiecewiseInverter.create_energy_density_inverse(material, 'E')

    # Test round-trip accuracy
    test_temp = 500.0
    energy_val = float(material.energy_density.subs(T, test_temp))
    recovered_temp = float(inverse_func.subs(E, energy_val))
    print(f"Round-trip: T={test_temp} -> E={energy_val:.2e} -> T={recovered_temp:.2f}")
```

## 📋 YAML Configuration Format

### Supported Property Types
- **CONSTANT_VALUE**: Simple numeric values
- **FILE_IMPORT**: Data loaded from CSV/Excel/text files
- **TABULAR_DATA**: Temperature and corresponding property value pairs
- **STEP_FUNCTION**: Discontinuous transitions
- **PIECEWISE_EQUATION**: Symbolic equations over temperature ranges
- **COMPUTED_PROPERTY**: Properties calculated from other properties

See the [YAML schema documentation](https://materforge.readthedocs.io/en/latest/reference/yaml_schema.html)
for detailed configuration options.

YAML configuration examples:
- [Test alloy (standalone example)](https://github.com/rahildoshi97/materforge/blob/main/examples/myAlloy.yaml)
- [Pure metal - Aluminum](https://github.com/rahildoshi97/materforge/blob/main/src/materforge/data/materials/pure_metals/Al/Al.yaml)
- [Alloy - Steel 1.4301 (full data)](https://github.com/rahildoshi97/materforge/blob/main/src/materforge/data/materials/alloys/1.4301/1.4301.yaml)

## 📚 Documentation

Full documentation is available at **https://materforge.readthedocs.io**

The documentation follows the [Diátaxis](https://diataxis.fr/) framework:

| Type | Content |
|------|---------|
| **Tutorials** | [Getting Started](https://materforge.readthedocs.io/en/latest/tutorials/getting_started.html) · [First Simulation](https://materforge.readthedocs.io/en/latest/tutorials/first_simulation.html) |
| **How-to Guides** | [Defining Material Properties](https://materforge.readthedocs.io/en/latest/how-to/define_materials.html) · [Energy-Temperature Conversion](https://materforge.readthedocs.io/en/latest/how-to/energy_temperature_conversion.html) |
| **Reference** | [API Reference](https://materforge.readthedocs.io/en/latest/reference/api.html) · [YAML Schema](https://materforge.readthedocs.io/en/latest/reference/yaml_schema.html) |
| **Explanation** | [Design Philosophy](https://materforge.readthedocs.io/en/latest/explanation/design_philosophy.html) · [Material Properties](https://materforge.readthedocs.io/en/latest/explanation/material_properties.html) |

## 🤝 Contributing

Contributions are welcome! Please see our
[Contributing Guide](https://github.com/rahildoshi97/materforge/blob/main/CONTRIBUTING.md)
for details on how to get started.

## 🐛 Known Limitations
- **Piecewise Inverter**: Currently supports only linear piecewise functions
- **File Formats**: Limited to CSV, Excel, and text files
- **Memory Usage**: Large datasets may require optimization for very high-resolution data
- **Regression**: Maximum 8 segments recommended for stability

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

## 📖 Citation

If you use MaterForge in your research, please cite it using the information in our
[CITATION.cff](https://github.com/rahildoshi97/materforge/blob/main/CITATION.cff) file.

## 📞 Support
- **Author**: Rahil Doshi
- **Email**: [rahil.doshi@fau.de](mailto:rahil.doshi@fau.de)
- **Documentation**: [materforge.readthedocs.io](https://materforge.readthedocs.io)
- **Bug Tracker**: [GitHub Issues](https://github.com/rahildoshi97/materforge/issues)
- **GitLab**: [i10git.cs.fau.de](https://i10git.cs.fau.de/rahil.doshi/materforge)

## 🙏 Acknowledgments
- Built with [SymPy](https://www.sympy.org/) for symbolic mathematics
- Data handling powered by [pandas](https://pandas.pydata.org/)
- Uses [pwlf](https://github.com/cjekel/piecewise_linear_fit_py) for piecewise linear fitting
- Visualization powered by [Matplotlib](https://matplotlib.org/)
- YAML parsing with [ruamel.yaml](https://yaml.dev/doc/ruamel.yaml/)

---

*MaterForge - Empowering material simulation with Python*
