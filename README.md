# MaterForge - Materials Formulation Engine with Python

A high-performance Python library for material simulation and analysis. MaterForge enables efficient modeling of pure metals and alloys through YAML configuration files, providing symbolic and numerical property evaluation for various material properties.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Latest Release](https://i10git.cs.fau.de/rahil.doshi/materforge/-/badges/release.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/releases)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE)
[![Documentation Status](https://readthedocs.org/projects/materforge/badge/?version=latest)](https://materforge.readthedocs.io/en/latest/?badge=latest)
[![Pipeline Status](https://i10git.cs.fau.de/rahil.doshi/materforge/badges/master/pipeline.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/pipelines)
[![Code Coverage](https://i10git.cs.fau.de/rahil.doshi/materforge/badges/master/coverage.svg)](https://i10git.cs.fau.de/rahil.doshi/materforge/-/commits/master)

## Table of Contents
- [Key Features](#-key-features)
- [Installation](#-installation)
- [Quick Sart](#-quick-start)
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
- **Material Types**: Support for both pure metals and alloys with appropriate phase transition temperatures
- **Multiple Property Types**: Support for constants, step functions, file-based data, key-value pairs, and computed properties
- **Regression Analysis**: Built-in piecewise linear fitting with configurable parameters

## 📦 Installation
### Prerequisites
- Python 3.10 or higher
- Required dependencies: `numpy`, `sympy`, `matplotlib`, `pandas`, `ruamel.yaml`

### Install using pip
```
pip install materforge
```
### Development Installation
```bash
git clone https://i10git.cs.fau.de/rahil.doshi/materforge.git
cd materforge
pip install -e .[dev]
```

## 🏃 Quick Start
### Basic Material Creation

```python
import sympy as sp
from materforge.parsing.api import create_material

# Create a material with symbolic temperature
T = sp.Symbol('T')
material_T = create_material('path/to/material.yaml', T)

# Access properties
print(f"Heat capacity: {material_T.heat_capacity}")

# Evaluate at specific temperature
temp_value = 1500.0  # Kelvin
density_at_temp = material_T.evaluate_properties_at_temperature(temp_value)
print(f"Density at {temp_value}K: {density_at_temp:.2f} kg/m³")

# For symbolic expressions with automatic plotting
material_with_plot = create_material('steel.yaml', T, enable_plotting=True)
```
### Working with Piecewise Inverse Functions

```python
from materforge.algorithms.piecewise_inverter import PiecewiseInverter

# Create inverse energy density function: T = f_inv(E)
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

See [the YAML schema documentation](docs/reference/yaml_schema.md) for detailed configuration options.
YAML configuration examples can be found here:
- [Pure Metals](src/materforge/data/materials/pure_metals/Al/Al.yaml)
- [Alloys](src/materforge/data/materials/alloys/1.4301/1.4301.yaml)

## 📚 Documentation
Our documentation follows the _Diátaxis_ framework with four distinct types:
### Tutorials - Learning-oriented guides
- [Getting Started with materforge](docs/tutorials/getting_started.md)
- [Creating Your First Material Simulation](docs/tutorials/first_simulation.md)
### How-to Guides - Problem-oriented instructions
- [Defining Custom Material Properties](docs/how-to/define_materials.md)
- [Converting Between Energy Density and Temperature](docs/how-to/energy_temperature_conversion.md)
### Reference - Information-oriented documentation
- [API Reference](docs/reference/api)
- [YAML Configuration Schema](docs/reference/yaml_schema.md)
### Explanation - Understanding-oriented discussions
- [Material Properties Concepts](docs/explanation/material_properties.md)
- [Design Philosophy](docs/explanation/design_philosophy.md)

## 🤝 Contributing
Contributions are welcome! Please see our [Contributing Guide](CONTRIBUTING.md) for details on how to get started.

## 🐛 Known Limitations
- **Piecewise Inverter**: Currently supports only linear piecewise functions
- **File Formats**: Limited to CSV, Excel, and text files
- **Memory Usage**: Large datasets may require optimization for very high-resolution data
- **Regression**: Maximum 8 segments recommended for stability

## 📄 License
This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

## 📖 Citation
If you use MaterForge in your research, please cite it using the information in our [CITATION.cff](CITATION.cff) file.

## 📞 Support
- **Author**: Rahil Doshi
- **Email**: rahil.doshi@fau.de
- **Project Homepage**: [materforge](https://i10git.cs.fau.de/rahil.doshi/materforge)
- **Bug Tracker**: [Issues](https://i10git.cs.fau.de/rahil.doshi/materforge/-/issues)

## 🙏 Acknowledgments
- Built with [SymPy](https://www.sympy.org/) for symbolic mathematics
- Data handling powered by [pandas](https://pandas.pydata.org/)
- Uses [pwlf](https://github.com/cjekel/piecewise_linear_fit_py) for piecewise linear fitting
- Visualization powered by [Matplotlib](https://matplotlib.org/)
- YAML parsing with [ruamel.yaml](https://yaml.dev/doc/ruamel.yaml/)

#### MaterForge - Empowering material simulation with Python 🚀

---
