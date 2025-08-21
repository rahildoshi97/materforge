---
title: 'MaterForge: Materials Formulation Engine with Python'
tags:
  - Python
  - materials science
  - symbolic computation
  - scientific computing
  - YAML configuration
  - high-performance computing
authors:
  - name: Rahil Miten Doshi
    orcid: 0009-0008-3570-9841
    affiliation: 1, 2
  - name: Harald Koestler
    orcid: 0000-0002-6992-2690
    affiliation: 1, 3
  - name: Matthias Markl
    orcid: 0000-0002-3000-7081
    affiliation: 2
affiliations: 
  - name: Chair for System Simulation, Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany
    index: 1
  - name: Chair of Materials Science and Engineering for Metals, Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany
    index: 2
  - name: Erlangen National High Performance Computing Center (NHR@FAU), Erlangen, Germany
    index: 3
date: 20 August 2025
bibliography: paper.bib
---

# Summary

MaterForge is an extensible, open-source Python library that streamlines the definition and use of
material properties in numerical simulations.
The library supports complex material behaviors, from simple constants to experimental data
in user-friendly YAML configurations.
These are internally converted into symbolic mathematical expressions for scientific computing frameworks.
MaterForge supports various material types,
provides flexible property definitions,
and automatically resolves dependency order for derived properties while detecting cycles.
It is designed for high-performance computing (HPC) applications
and serves as a bridge between experimental data and numerical simulation.

# Statement of Need

Accurate numerical simulation requires material properties such as thermal conductivity, density, and viscosity
that depend on variables like temperature, pressure, or strain rate [@lewis1996finite].
This challenge is compounded by the wide variation in data availability,
from well-characterized models for established materials to sparse experimental points for novel materials.
Property definitions consequently range from simple constants to complex tabular datasets or sophisticated equations,
creating significant integration hurdles for researchers.

To manage this complexity, researchers often resort to manual interpolation, custom scripting, or proprietary software,
which compromises reproducibility and standardization [@ashby2013materials].
While valuable resources like the NIST WebBook [@linstrom2001nist] and CoolProp [@coolprop] provide valuable raw data,
they lack integrated processing to unify these varied formats.
CALPHAD databases [@calphad] are powerful but often require proprietary software
and do not easily integrate with general-purpose simulation codes.

This leads to ad hoc solutions, hindering workflow efficiency and FAIR data adoption [@wilkinson2016fair].
MaterForge bridges this gap by providing a unified framework that leverages
symbolic mathematics, automatic regression, and dependency resolution
to standardizes and simplify the integration of realistic material behavior into scientific simulations.

# Key Functionality

- **Flexible Input Methods**: The library supports various property definition methods such as
  constant values, step functions, file-based data (.xlsx, .csv, .txt), tabular data, piecewise equations, and computed properties 
  (\autoref{fig:input_methods_new}).
  This versatility allows users to leverage data from diverse sources.

![MaterForge's property definition methods with corresponding YAML examples and automatically generated validation plots.\label{fig:input_methods_new}](figures/input_methods_new.jpg)

- **Extensible Material Support**: The framework supports any material type through its extensible architecture.
  Currently implemented for pure metals and alloys,
  its modular design allows straightforward extension to materials such as
  ceramics, polymers, or composites.

- **Automatic Dependency Resolution**: For dependent properties 
  (e.g., density calculated from thermal expansion coefficient),
  MaterForge automatically determines the correct processing order, resolves mathematical dependencies,
  and detects circular references.

- **Regression and Data Reduction**: The library performs piecewise regression for large datasets, 
  simplifying complex property curves into efficient mathematical representations with configurable polynomial degrees and segments,
  reducing computational overhead while maintaining accuracy.

- **Intelligent Simplification Timing**: MaterForge provides sophisticated control over when data simplification occurs
  via the `simplify` parameter.
  `simplify: pre` optimizes performance by simplifying properties before being used in dependent calculations,
  while `simplify: post` defers simplification until all dependent properties have been computed, maximizing numerical accuracy.

- **Configurable Boundary Behavior**: Users can define how properties behave outside their specified ranges,
  choosing between `constant`-value or `extrapolation` to best match the physical behavior of the material.
  The boundary behavior options work seamlessly with the regression capabilities to provide comprehensive data processing control 
  (\autoref{fig:regression_options_with_boundary_behavior_new}).

```yaml
    bounds: [constant, extrapolate]
    regression:
      simplify: pre
      degree: 2
      segments: 3
```
![MaterForge's data processing capabilities: regression and data reduction showing raw data (points) fitted with different polynomial degrees and segment configurations, and configurable boundary behavior options demonstrating constant versus extrapolate settings for the same density property, illustrating how MaterForge can reduce complexity while maintaining physical accuracy and providing flexible boundary control.\label{fig:regression_options_with_boundary_behavior_new}](figures/regression_options_with_boundary_behavior_new.png)

- **Inverse Property Computation**: The library can generate inverse piecewise-linear functions,
  enabling the determination of independent variables from known property values.
  This capability is essential for energy-based numerical methods [@voller1987fixed],
  where temperature is computed via the inverse function of the enthalpy.

- **Built-in Validation Framework**: A comprehensive validation framework checks YAML configurations for correctness,
  including composition sums, required fields, and valid property names, 
  preventing common configuration errors [@roache1998verification].

- **Integrated Visualization**: An integrated visualization tool
  automatically generates plots to verify property definitions,
  with the option to disable visualization for production workflows.

# Usage

Materials are defined in YAML files and loaded via `create_material`, which returns a fully configured material object.

## YAML Configuration Example: Alloy (`steel.yaml`)
```yaml
name: Steel 1.4301
material_type: alloy
composition: {Fe: 0.675, Cr: 0.170, Ni: 0.120, Mo: 0.025, Mn: 0.010}
solidus_temperature: 1605.0
liquidus_temperature: 1735.0
initial_boiling_temperature: 3090.0
final_boiling_temperature: 3200.0
properties:
  density:
    file_path: ./1.4301.xlsx
    dependency_column: T (K)
    property_column: rho (kg/m^3)
    bounds: [constant, extrapolate]
    regression:
      simplify: pre
      degree: 1
      segments: 3
```

## Python Integration
```python
    import sympy as sp
    from materforge.parsing.api import create_material

    T = sp.Symbol('T')
    steel = create_material('steel.yaml', T, enable_plotting=True)
    steel_density = steel.density
    density_500k = steel.evaluate_properties_at_temperature(500.0, ['density'])
```

# Comparison with Existing Tools

| **Feature**              | **MaterForge** | **CoolProp** | **NIST WebBook** | **CALPHAD Tools** |
|:-------------------------|:---------------|:-------------|:-----------------|:------------------|
| Symbolic Integration     | Yes            | No           | No               | Limited           |
| Dependency Resolution    | Automatic      | No           | No               | No                |
| Input Methods            | 6 types        | 1            | 1                | 1                 |
| Solid Materials          | Yes            | Limited      | Yes              | Yes               |
| Custom Properties        | Any            | No           | No               | Limited           |
| Variable Support         | Any            | T, P only    | Static           | T, P, Comp.       |
| Open Source              | Yes            | Yes          | No               | Mixed             |
| Python Integration       | Native         | Yes          | API only         | Limited           |

**Key Advantage**: MaterForge's native symbolic mathematics via SymPy [@sympy],
automatic dependency resolution, and multiple input methods provide flexibility and integration
not found in existing tools, enabling more reproducible and sophisticated scientific simulations.

# Research Applications

MaterForge is applicable to alloy design [@callister2018materials],
finite element analysis [@hughes2012finite], multiscale modeling [@tadmor2011modeling],
computational fluid dynamics, and heat transfer.
Its architecture promotes reproducible science and is well-suited for HPC environments,
with demonstrated integrations into frameworks like pystencils [@pystencils] and waLBerla [@walberla].

# Availability

MaterForge is distributed under the [BSD-3-Clause License](https://github.com/rahildoshi97/materforge/blob/master/LICENSE). 
The source code is hosted on [GitHub](https://github.com/rahildoshi97/materforge), with
[full documentation](https://materforge.readthedocs.io/en/stable/) 
and [YAML examples](https://materforge.readthedocs.io/en/stable/how-to/define_materials.html).
The package can be installed via [PyPI](https://pypi.org/project/materforge/) using `pip install materforge`.

# Acknowledgements

This work was funded by the European High Performance Computing Joint Undertaking (Grant No. 101093457)
and the Deutsche Forschungsgemeinschaft within Research Unit FOR-5134 (Grant No. 434946896).
We thank Carola Forster for providing the material data for Steel 1.4301 using JMatPro.

# References
