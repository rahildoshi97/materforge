# Changelog

All notable changes to MaterForge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.3] - 2025-11-12

### Added
- Material-specific integration tests for piecewise inverter with real material data
- Support for four boundary configuration types in material properties:
  - `[constant, constant]`: domain [T_lower, T_upper]
  - `[constant, extrapolate]`: domain [T_lower, inf)
  - `[extrapolate, constant]`: domain (-inf, T_upper]
  - `[extrapolate, extrapolate]`: domain (-inf, inf)
- Domain validation guidance for piecewise function inversion
- Improved test documentation with domain specifications for each test case

### Changed
- Reorganized test assertions to use `sp.simplify()` for robust symbolic equivalence checking
- Updated test cases to respect material property domain boundaries
- Enhanced test coverage for constant piece boundary handling
- Improved error messages for out-of-domain temperature evaluations

### Fixed
- Test domain validation: removed out-of-domain test temperatures that were incorrectly failing
  - Constant pieces now correctly understood to map all values in their domain to a single output
  - Inverse function behavior validated to correctly return boundary temperatures for constant pieces
  - Tests now only evaluate within valid material property domains
- Test assertions for final constant pieces: now correctly verify boundary temperature returns

### Testing
- 5 new material-specific test methods covering all boundary configurations
- Comprehensive `test_material_specific_enthalpy_all_configurations()` meta-test
- All tests passing with domain-respecting temperature ranges
- Real Aluminum specific_enthalpy piecewise functions used for validation
- Round-trip accuracy <1e-9 K across all valid temperature ranges
- Total: 328 unit tests passing (100% pass rate)

### Documentation
- Added detailed docstrings explaining domain constraints for each boundary configuration
- Clarified mathematical behavior of constant pieces in piecewise functions
- Added domain notation (e.g., `domain: [300K, 4500K]`) to test documentation

## [0.6.2] - 2025-11-10

### Added
- Comprehensive piecewise function inversion algorithm for symbolic material properties
- Support for monotonic piecewise linear functions with automatic slope detection
- Boundary condition validation and energy-to-temperature mapping
- `PiecewiseInverter.create_energy_density_inverse()` static method for material energy density inverse functions
- Robust handling of constant and linear piecewise pieces
- Monotonicity validation to prevent inverting non-monotonic functions
- Comprehensive test suite for piecewise inversion (13+ unit tests)

### Changed
- Enhanced `piecewise_inverter.py` with two-pass boundary collection and inversion algorithm
- Improved constant piece handling: returns SymPy expressions instead of raw floats for consistency
- Updated boundary temperature tracking for proper inverse function domain mapping
- Refined logging with detailed slope and monotonicity diagnostics
- Enhanced error messages for debugging piecewise function issues

### Fixed
- **Critical bug in piecewise inversion**: Constant pieces in inverse functions now correctly map to boundary temperatures instead of constant energy values
  - Previously: `T=300K → E=1.20e+05 → T=119597.6K (Error: 1.19e+05)`
  - Now: `T=300K → E=1.20e+05 → T=300.0K (Error: ~0)`
- Fixed boundary condition direction for decreasing (negative slope) piecewise pieces
- Corrected inverse condition logic: now uses `E > boundary_energy` for decreasing functions and `E < boundary_energy` for increasing functions
- Missing `boundary_temp` parameter in loop iteration (was causing 4500K offset errors)
- SymPy expression type consistency: all inverse pieces now return SymPy objects via `sp.sympify()`
- Test assertions now use `sp.simplify()` for symbolic equivalence comparison
- Round-trip accuracy improved from ~1.19e+05 K error to <1e-10 K precision

### Performance
- Added monotonicity validation to fail fast on invalid piecewise functions
- Efficient two-pass algorithm: first pass collects boundaries, second pass builds inverse
- Minimal computational overhead for piecewise inversion

### Testing
- All 306 unit tests passing (100% pass rate)
- 13 dedicated piecewise inversion tests covering:
  - Linear piecewise functions with multiple segments
  - Constant piece handling
  - Boundary condition extraction
  - Decreasing and increasing slopes
  - Edge cases and non-linear function rejection
  - Material energy density inverse creation
  - Round-trip accuracy validation
- Integration tests confirm near-zero round-trip error (<1e-13) for all temperature ranges

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
