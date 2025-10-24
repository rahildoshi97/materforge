# Changelog

All notable changes to MaterForge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.1] - 2025-10-24

### Added
- SPDX license headers to all Python source files in `src/materforge/`
- SPDX (GPL-3.0-or-later) headers for application files in `apps/`
- `apps/LICENSE` and `apps/README.md` describing separate GPL licensing
- Tag regex support in `pyproject.toml` to handle `v`-prefixed tags
- Expanded license section and summary table in `README.md`

### Changed
- Split BSD‑3‑Clause core license and GPLv3 examples licensing
- Updated installation guide with submodule instructions
- Excluded GPL apps from PyPI distribution

### Fixed
- `setuptools_scm` version detection for `v0.x`-style tags now returns the correct release version (0.6.1)

## [0.6.0] - 2025-08-20

### Added
- Release automation improvements and documentation badge updates
- Improved type annotations for dependency parameters (`Union[float, sp.Symbol]`)
- Excluded `materforge_plots` from package distributions to reduce build size

### Changed
- Updated documentation for stable version release
- Enhanced installation and usage documentation in `README.md`
- Updated `CITATION.cff` metadata for improved scholarly citation support
- Adjusted `pyproject.toml` for refined dependency management
- Improved docstrings and code comments for maintainability

### Fixed
- Resolved redundant file inclusion issues in packaging
- Minor fixes across documentation files and CI pipeline

## [0.5.6] - 2025-08-06

### Changed
- Updated README.md with improved installation instructions and usage examples
- Improved material definition guide in `docs/how-to/define_materials.md`
- Updated API documentation in `docs/reference/api/material.md`
- Refined YAML schema documentation in `docs/reference/yaml_schema.md`
- Enhanced getting started tutorial with clearer examples

## [0.5.5] - 2025-07-30

### Added
- Complete rebranding from PyMatLib to MaterForge
- Materials Analysis & Transformation Engine functionality
- Enhanced symbolic computation with SymPy integration
- Comprehensive material property modeling
- Piecewise function building and regression capabilities
- Property visualization and plotting tools
- YAML-based material configuration system
- Support for pure metals and alloys
- Temperature-dependent property evaluation

### Changed
- Package name: `pymatlib` → `materforge`
- Import statements: `from pymatlib` → `from materforge`
- Repository structure and naming conventions
- Enhanced error handling and validation
- Improved logging and debugging capabilities

### Fixed
- NumPy float64 compatibility issues with SymPy
- Property interpolation edge cases
- Temperature validation and bounds checking
- Dependency resolution for computed properties

### Breaking Changes
- Package name changed from `pymatlib` to `materforge`
- All import statements must be updated
- Installation command: `pip install materforge`

### Migration Guide
```


# Old way

from pymatlib.parsing.api import create_material
material = create_material('steel.yaml', T)

# New way

from materforge.parsing.api import create_material
material = create_material('steel.yaml', T)

```
