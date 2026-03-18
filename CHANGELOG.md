# Changelog

All notable changes to MaterForge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.6] - 2026-03-16

### Changed
- Allowed int for constant material property type in YAML config 

## [0.6.5] - 2026-03-09

### Added
- `MaterialConfigError` and `PropertyConfigError` typed exception classes in `errors.py`;
  both subclass `ValueError` so existing callers are unaffected
- `UnicodeDecodeError` handler in `YAMLFileParser._load_config` with a clear actionable
  message ("Save as UTF-8 and try again") instead of a generic crash
- Al2O3 ceramic YAML material configuration
- `__getattr__` on `Material` class to suppress Pylint E1101 false positives for
  dynamically assigned properties
- `docs/_static/.gitkeep` to ensure the directory is tracked in git and Sphinx
  `html_static_path` resolves correctly in CI

### Changed
- **Breaking**: renamed `extrapolate` boundary keyword to `linear` in all YAML files
  to accurately reflect linear extrapolation behaviour; update any custom YAML files
- Enforced strict top-level YAML schema: only `name` and `properties` keys are allowed;
  unknown keys raise `MaterialConfigError` immediately
- Removed `material_type`, `composition`, and other non-essential top-level YAML keys;
  schema is now intentionally minimal
- Replaced fixed property list with dynamic tracker across all processors;
  removed `get_supported_properties` from public API
- All properties on `Material` are now stored in an internal dict via `__setattr__`;
  `evaluate()` returns a new `Material` copy with all expressions substituted at the
  given value - not a dict
- Dependency is now always a `sp.Symbol`; numeric mode removed - use
  `material.evaluate(symbol, value)` for constant properties instead
- `T` restored as the fixed YAML placeholder symbol; equations in YAML files must
  always use `T`, which is substituted at runtime with the caller-provided symbol
- `_validate_config` now raises `MaterialConfigError` (was `ValueError`)
- `_analyze_and_categorize_properties` now raises `PropertyConfigError` (was `ValueError`)
- `create_material` in `material_yaml_parser.py` now wraps failures in
  `MaterialConfigError` (was plain `ValueError`)
- `validate_yaml_file` and `get_material_info` in `api.py` now re-raise typed
  exceptions directly instead of wrapping everything in a new `ValueError`
- `get_material_property_names` and `evaluate_material_properties` now raise
  `TypeError` for wrong argument types (was `ValueError`)
- `open()` call in `YAMLFileParser` now explicitly specifies `encoding='utf-8'`
  for platform-independent behaviour (fixes Pylint W1514)
- `DependencyError.__init__`: `available_props: List[str] = None` →
  `Optional[List[str]] = None` (fixes Pylint type hint warning)
- Simplified materials directory structure

### Fixed
- Removed unused exception variable bindings `as e` in bare `raise` handlers
  (Ruff F841); affected `validate_yaml_file` and `get_material_info` in `api.py`
- Eliminated duplicate error log entries that occurred when exceptions bubbled
  through multiple layers
- Removed unused imports across `api.py`, `material_yaml_parser.py`, and app files

### Testing
- Updated `test_api.py` to match refactored exception types:
  - `test_get_material_property_names_rejects_non_material`: `ValueError` → `TypeError`
  - `test_validate_incomplete_yaml_raises_value_error`: `ValueError` → `MaterialConfigError`
- Moved error class imports to module level in test files to prevent `pytest.raises`
  class identity mismatches caused by local imports resolving to different class objects
- All 295 unit tests passing

### Documentation
- Updated all docstring `Raises` sections in `api.py` and `material_yaml_parser.py`
  to reflect actual exception types; removed stale `ValueError` references
- `evaluate()` docstring in `materials.py`: corrected return description from
  "dict of float values" to "new Material instance with substituted properties"
- `reference/api/material.md`: updated error table entries from `ValueError` to
  `MaterialConfigError`/`PropertyConfigError`; corrected `get_material_property_names`
  description (removed incorrect `sorted()` equivalence claim); updated all `Raises`
  sections; fixed `evaluate()` example to use attribute access on returned `Material`
- `explanation/material_properties.md`: fixed `evaluate()` example (dict access →
  attribute access on `Material`); removed duplicate T-placeholder paragraph
- `how-to/define_materials.md`: fixed `evaluate()` example (dict access →
  attribute access on `Material`)
- `how-to/property_inversion.md`: fixed `evaluate()` example; fixed `create_inverse()`
  calls to pass `sp.Symbol` instances instead of string literals
- `tutorials/getting_started.md`: fixed `evaluate()` iteration - `results.items()`
  replaced with `property_names()` + `getattr()` since `evaluate()` returns `Material`
- `tutorials/first_simulation.md`: fixed `create_material()` to use a plain
  `sp.Symbol` - `u.center()` is not a valid `sp.Symbol` and raises `TypeError`;
  added `.subs(T, u.center())` at the assignment level as the single pystencils
  coupling point; fixed `create_inverse()` string args → `sp.Symbol` instances;
  removed redundant second `create_material()` call in complete example

## [0.6.4] - 2026-02-22

### Added
- User example for a hypothetical alloy demonstrating custom material definition workflow
- Ruff linter validation added to CI pipeline for code quality enforcement

### Changed
- Replaced hardcoded property lists with metadata-tagged fields for cleaner, more maintainable property registration
- Replaced hardcoded physical constants with values sourced from `scipy.constants` (CODATA 2022) for accuracy and traceability
- Migrated `setuptools_scm` `write_to` directive to `version_file` (new setuptools_scm >=8 API)
- Consolidated all dependencies into `pyproject.toml`; removed standalone `docs/requirements.txt`
- Moved shell scripts from root into `/scripts` folder for cleaner project structure
- Removed `[all]` optional dependency group from `pyproject.toml`
- Clean PyPI distribution: excluded `materforge_plots` plot assets and updated license metadata to SPDX string format

### Fixed
- Piecewise inverter final piece detection: replaced `== True` identity check with `is sp.true` for correct SymPy `Piecewise` condition evaluation
- Fixed internal script paths after moving `.sh` scripts to `/scripts` folder
- Fixed inconsistent package naming across modules
- Removed unused imports and non-f-string `f` prefixes (flagged by Ruff)
- Removed unused YAML key from material configuration

### Documentation
- Updated `getting_started.md` with missing `material_type` field and walberla integration notes
- Updated `index.rst` to reflect current module structure
- Added hypothetical alloy example to tutorials
- General documentation refresh across how-to guides and reference pages
- Updated `README.md` with current usage patterns

### Testing
- All unit tests passing

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
  - Previously: `T=300K -> E=1.20e+05 -> T=119597.6K (Error: 1.19e+05)`
  - Now: `T=300K -> E=1.20e+05 -> T=300.0K (Error: ~0)`
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
- Package name: `pymatlib` -> `materforge`
- Import statements: `from pymatlib` -> `from materforge`
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
