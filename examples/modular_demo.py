# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Modular demonstration script for material property evaluation."""
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import sympy as sp
from materforge.algorithms.piecewise_inverter import PiecewiseInverter
from materforge.parsing.api import (
    create_material,
    evaluate_material_properties,
    get_material_info,
    get_material_property_names,
    validate_yaml_file,
)


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

@dataclass
class TestConfig:
    """Controls which tests and subtests to run."""
    yaml_validation: bool = True
    material_info_extraction: bool = True
    material_creation: bool = True
    supported_properties_overview: bool = True
    material_property_analysis: bool = True
    property_evaluation_specific_temp: bool = True
    new_api_methods: bool = True
    batch_temperature_evaluation: bool = True
    inverse_function_testing: bool = True
    advanced_usage: bool = True
    temperature_sweep_analysis: bool = True

    # Subtests for property analysis
    manual_evaluation: bool = True
    object_oriented_evaluation: bool = True
    functional_evaluation: bool = True

    # Material selection
    include_aluminum: bool = True
    include_steel: bool = True

    # Test parameters
    single_test_temperature: float = 500.15
    batch_temperatures: Optional[List[float]] = None
    sweep_temp_range: Tuple[int, int, int] = (300, 3001, 100)

    def __post_init__(self) -> None:
        if self.batch_temperatures is None:
            self.batch_temperatures = [300, 400, 500, 600, 700]


class MaterialPropertyDemonstrator:
    """Modular demonstration class for material property evaluation."""

    def __init__(self, config: TestConfig) -> None:
        self.config = config
        self.materials = []
        self.material_symbols: Dict[str, sp.Symbol] = {}
        self._setup_logging()
        self._setup_paths()

    @staticmethod
    def _setup_logging() -> None:
        logging.basicConfig(
            level=logging.WARNING,
            format="%(asctime)s %(levelname)s %(name)s -> %(message)s",
        )
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        logging.getLogger('PIL').setLevel(logging.WARNING)
        logging.getLogger('fontTools').setLevel(logging.WARNING)

    def _setup_paths(self) -> None:
        base = (Path(__file__).parent.parent / "src" / "materforge" / "data" / "materials")
        self.yaml_paths: Dict[str, Path] = {}
        if self.config.include_aluminum:
            self.yaml_paths["Aluminum"] = base / "pure_metals" / "Al" / "Al.yaml"
        if self.config.include_steel:
            self.yaml_paths["Steel_1.4301"] = base / "alloys" / "1.4301" / "1.4301.yaml"

    @staticmethod
    def _scalar_properties(mat) -> dict:
        result = {}
        for name in mat.property_names():
            val = getattr(mat, name)
            is_symbolic = hasattr(val, 'free_symbols') and val.free_symbols
            if not is_symbolic and isinstance(val, (int, float, sp.Float, sp.Integer)):
                result[name] = val
        return result

    # ------------------------------------------------------------------
    # Main runner
    # ------------------------------------------------------------------

    def run_all_tests(self) -> None:
        print(f"\n{'=' * 80}")
        print("MaterForge MATERIAL PROPERTY DEMONSTRATION")
        print(f"{'=' * 80}")
        self._print_config_summary()

        if self.config.yaml_validation:
            self.test_yaml_validation()
        if self.config.material_info_extraction:
            self.test_material_info_extraction()
        if self.config.material_creation:
            self.test_material_creation()
        if self.config.supported_properties_overview:
            self.test_supported_properties_overview()
        if self.config.material_property_analysis:
            self.test_material_property_analysis()
        if self.config.inverse_function_testing:
            self.test_inverse_functions()
        if self.config.advanced_usage:
            self.test_advanced_usage()

        print(f"\n{'=' * 80}")
        print("DEMONSTRATION COMPLETED")
        print(f"{'=' * 80}")

    def _print_config_summary(self) -> None:
        print("\nTest Configuration Summary:")
        print(f"{'-' * 40}")
        flags = [
            self.config.yaml_validation,
            self.config.material_info_extraction,
            self.config.material_creation,
            self.config.supported_properties_overview,
            self.config.material_property_analysis,
            self.config.inverse_function_testing,
            self.config.advanced_usage,
        ]
        print(f"Enabled Tests:    {sum(flags)}/7")
        mat_str = " ".join(filter(None, [
            "Al"    if self.config.include_aluminum else "",
            "Steel" if self.config.include_steel    else "",
        ]))
        print(f"Materials:        {mat_str}")
        print(f"Test Temperature: {self.config.single_test_temperature} K")
        if self.config.batch_temperature_evaluation:
            print(f"Batch Temps:      {self.config.batch_temperatures}")
        if self.config.temperature_sweep_analysis:
            start, stop, step = self.config.sweep_temp_range
            print(f"Sweep Range:      {start}–{stop} K (step: {step} K)")

    # ------------------------------------------------------------------
    # 1. YAML validation
    # ------------------------------------------------------------------

    def test_yaml_validation(self) -> None:
        print(f"\n{'1. YAML FILE VALIDATION':<50}")
        print(f"{'-' * 50}")
        for name, path in self.yaml_paths.items():
            if path.exists():
                try:
                    is_valid = validate_yaml_file(path)
                    print(f"{name} YAML validation: {'PASSED' if is_valid else 'FAILED'}")
                except Exception as e:
                    print(f"{name} YAML validation: FAILED - {e}")
            else:
                print(f"{name} YAML file not found: {path}")

    # ------------------------------------------------------------------
    # 2. Material info extraction
    # ------------------------------------------------------------------

    def test_material_info_extraction(self) -> None:
        print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
        print(f"{'-' * 50}")
        for name, path in self.yaml_paths.items():
            if not path.exists():
                print(f"{name} YAML file not found: {path}")
                continue
            try:
                info = get_material_info(path)
                print(f"\n{name} Information:")
                print(f"  Name:             {info['name']}")
                print(f"  Total Properties: {info['total_properties']}")
                print(f"  Properties:       {info['properties']}")
                if 'property_types' in info:
                    print("  Property Types:")
                    for pt, count in info['property_types'].items():
                        print(f"    {pt:<30}: {count}")
            except Exception as e:
                print(f"Failed to get {name} info: {e}")

    # ------------------------------------------------------------------
    # 3. Material creation
    # ------------------------------------------------------------------

    def test_material_creation(self) -> None:
        print(f"\n{'3. MATERIAL CREATION':<50}")
        print(f"{'-' * 50}")
        symbol_map = {
            "Aluminum": sp.Symbol('T_Al'),
            "Steel_1.4301": sp.Symbol('T_SS'),
        }
        for name, path in self.yaml_paths.items():
            if not path.exists():
                print(f"{name}: YAML not found")
                continue
            try:
                T = symbol_map.get(name, sp.Symbol(f'T_{name}'))
                mat = create_material(yaml_path=path, dependency=T, enable_plotting=True)
                self.materials.append(mat)
                self.material_symbols[mat.name] = T
                print(f"Successfully created: {mat.name} (symbol: {T})")
            except Exception as e:
                print(f"Failed to create {name}: {e}")

    # ------------------------------------------------------------------
    # 4. Property types overview
    # ------------------------------------------------------------------

    def test_supported_properties_overview(self) -> None:
        print(f"\n{'4. PROPERTY TYPES OVERVIEW':<50}")
        print(f"{'-' * 50}")
        print("MaterForge supports 6 property definition types (any property name is valid):")
        rows = [
            ("CONSTANT_VALUE",     "Single numeric value, e.g. density: 7000.0"),
            ("STEP_FUNCTION",      "Two values split at a scalar-property reference"),
            ("TABULAR_DATA",       "Paired dependency / value lists"),
            ("FILE_IMPORT",        "Load data from .csv / .xlsx / .txt"),
            ("PIECEWISE_EQUATION", "Symbolic equations over temperature breakpoints"),
            ("COMPUTED_PROPERTY",  "Derived from other properties via expression"),
        ]
        for i, (ptype, desc) in enumerate(rows, 1):
            print(f"  {i}. {ptype:<30} - {desc}")

    # ------------------------------------------------------------------
    # 5. Material property analysis
    # ------------------------------------------------------------------

    def test_material_property_analysis(self) -> None:
        print(f"\n{'5. MATERIAL PROPERTY ANALYSIS':<50}")
        print(f"{'-' * 50}")
        for mat in self.materials:
            T_mat = self.material_symbols[mat.name]
            print(f"\n{'=' * 80}")
            print(f"MATERIAL: {mat.name}  (temperature symbol: {T_mat})")
            print(f"{'=' * 80}")
            self._analyze_basic_material_info(mat)
            self._analyze_available_properties(mat)
            if self.config.property_evaluation_specific_temp:
                self._test_property_evaluation_at_specific_temp(mat, T_mat)
            if self.config.batch_temperature_evaluation:
                self._test_batch_temperature_evaluation(mat, T_mat)  # ← pass T_mat

    def _analyze_basic_material_info(self, mat) -> None:
        print(f"Name: {mat.name}")
        scalars = self._scalar_properties(mat)
        if scalars:
            print("Scalar constants:")
            for prop_name, val in sorted(scalars.items()):
                print(f"  {prop_name:<35}: {val}")

    def _analyze_available_properties(self, mat) -> None:
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")
        available_props = get_material_property_names(mat)
        print(f"Material has {len(available_props)} processed properties:")
        for prop_name in sorted(available_props):
            prop_value = getattr(mat, prop_name)
            print(f"  {prop_name:<30}: {prop_value}")
            print(f"  {'':30}  Type: {type(prop_value).__name__}")
            allowed = (sp.Piecewise, sp.Float, sp.Expr, float, int, type(None))
            marker = "Valid type" if isinstance(prop_value, allowed) else "! WARNING: Unexpected type"
            print(f"  {'':30}  {marker}")

    def _test_property_evaluation_at_specific_temp(self, mat, T_mat) -> None:
        test_temp = self.config.single_test_temperature
        available_props = get_material_property_names(mat)

        print(f"\n{'PROPERTY VALUES AT ' + str(test_temp) + 'K:':<50}")
        print(f"{'-' * 50}")

        if self.config.manual_evaluation:
            print("Method 1: Manual Property Evaluation")
            for prop_name in sorted(available_props):
                try:
                    prop_value = getattr(mat, prop_name)
                    if isinstance(prop_value, sp.Expr) and prop_value.free_symbols:
                        numerical = prop_value.subs(T_mat, test_temp).evalf()
                        print(f"  {prop_name:<30}: {numerical} (symbolic)")
                    else:
                        print(f"  {prop_name:<30}: {prop_value} (constant)")
                except Exception as e:
                    print(f"  {prop_name:<30}: Error - {e}")

        if self.config.new_api_methods:
            self._test_new_api_methods(mat, T_mat, test_temp)  # ← pass T_mat

    def _test_new_api_methods(self, mat, T_mat: sp.Symbol, test_temp: float) -> None:
        print(f"\n{'NEW API METHODS:':<50}")
        print(f"{'-' * 50}")
        all_values = None

        if self.config.object_oriented_evaluation:
            print("Method 2: material.evaluate(symbol, value)")
            try:
                all_values = mat.evaluate(T_mat, test_temp)
                print(f"All properties at {test_temp} K:")
                for prop, value in sorted(all_values.items()):
                    print(f"  {prop:<30}: {value:.6e}")
            except Exception as e:
                print(f"  Error: {e}")

        if self.config.functional_evaluation:
            print("\nMethod 3: evaluate_material_properties(material, symbol, value)")
            try:
                all_values_func = evaluate_material_properties(mat, T_mat, test_temp)
                if self.config.object_oriented_evaluation and all_values is not None:
                    print(f"  Results match Method 2: {all_values == all_values_func}")
                else:
                    print("  Functional API evaluation completed successfully")
            except Exception as e:
                print(f"  Error: {e}")

    def _test_batch_temperature_evaluation(self, mat, T_mat: sp.Symbol) -> None:
        print(f"\n{'BATCH TEMPERATURE EVALUATION:':<50}")
        print(f"{'-' * 50}")

        batch_results = {}
        for temp in self.config.batch_temperatures:
            try:
                batch_results[temp] = mat.evaluate(T_mat, temp)
            except Exception as e:
                print(f"  Error at {temp}K: {e}")
                batch_results[temp] = {}

        if not batch_results:
            return

        try:
            df = pd.DataFrame(batch_results).T
            available_props = get_material_property_names(mat)
            display_cols = [c for c in sorted(available_props) if c in df.columns][:3]
            if not display_cols:
                print("  No numeric columns available for batch display")
                return
            header = f"{'Temperature':<12}" + "".join(f"{c:<18}" for c in display_cols)
            print(header)
            print("-" * (12 + 18 * len(display_cols)))
            for temp, row in df.iterrows():
                line = f"{temp:<12.0f}"
                for col in display_cols:
                    val = row.get(col)
                    line += f"{val:.4e}  " if val is not None else f"{'N/A':<18}"
                print(line)
        except Exception as e:
            print(f"  DataFrame creation failed: {e}")

    # ------------------------------------------------------------------
    # 6. Inverse function testing
    # ------------------------------------------------------------------

    def test_inverse_functions(self) -> None:
        print(f"\n{'6. INVERSE FUNCTION TESTING':<50}")
        print(f"{'-' * 50}")
        print(f"\n{'=' * 80}")
        print("INVERSE FUNCTION TESTING")
        print(f"{'=' * 80}")

        for mat in self.materials:
            T_mat = self.material_symbols[mat.name]
            print(f"\n--- {mat.name} (symbol: {T_mat}) ---")

            if 'energy_density' not in mat.property_names():
                print(f"  {mat.name} has no energy_density - skipping")
                continue

            energy_density = getattr(mat, 'energy_density')
            try:
                syms = energy_density.free_symbols
                if len(syms) == 1:
                    temp_sym = next(iter(syms))
                    E = sp.Symbol('E')
                    inv = PiecewiseInverter.create_inverse(energy_density, temp_sym, E)
                    print("  Inverse created successfully")
                    self._test_round_trip_accuracy(energy_density, inv, temp_sym, E)
                else:
                    print(f"  Unexpected free symbols: {syms}")
            except Exception as e:
                print(f"  Inverse testing failed for {mat.name}: {e}")

    def _test_round_trip_accuracy(self, energy_density_expr, inverse_func, temp_symbol, energy_symbol) -> None:
        print("  Round-trip accuracy test:")
        test_temps = [300, 350, 400, 450, 500, 600, 700, 800, 900, 1000,
                      1200, 1500, 1800, 2000, 2500, 3000]
        max_error = 0.0
        successful = 0

        for temp in test_temps:
            try:
                e_val = float(energy_density_expr.subs(temp_symbol, temp))
                recovered = float(inverse_func.subs(energy_symbol, e_val))
                error = abs(temp - recovered)
                max_error = max(max_error, error)
                successful += 1
                print(f"    T={temp:4.0f}K -> E={e_val:.2e} -> T={recovered:6.1f}K  (error: {error:.2e})")
            except Exception as e:
                print(f"    T={temp:4.0f}K -> Error: {e}")

        print(f"  Completed: {successful}/{len(test_temps)}")
        print(f"  Max error: {max_error:.2e} K")
        if max_error < 1e-8:
            print("  Excellent accuracy")
        elif max_error < 1e-4:
            print("  Good accuracy")
        else:
            print("  ! Consider reviewing inverse function accuracy")

    # ------------------------------------------------------------------
    # 7. Advanced usage / temperature sweep
    # ------------------------------------------------------------------

    def test_advanced_usage(self) -> None:
        print(f"\n{'=' * 80}")
        print("ADVANCED USAGE PATTERNS")
        print(f"{'=' * 80}")
        if self.config.temperature_sweep_analysis:
            self._test_temperature_sweep_analysis()

    def _test_temperature_sweep_analysis(self) -> None:
        print("\n7. TEMPERATURE SWEEP ANALYSIS")
        print("-" * 50)

        start, stop, step = self.config.sweep_temp_range
        temp_range = range(start, stop, step)

        for mat in self.materials:
            T_mat = self.material_symbols[mat.name]
            print(f"\n--- Temperature Sweep for {mat.name} ---")

            # Pick up to 3 temperature-dependent (symbolic) properties dynamically
            sweep_props = [
                p for p in sorted(mat.property_names())
                if hasattr(getattr(mat, p, None), 'free_symbols')
                and getattr(mat, p).free_symbols
            ][:3]

            if not sweep_props:
                print(f"  No temperature-dependent properties for {mat.name}")
                continue

            sweep_data = []
            for temp in temp_range:
                try:
                    # Evaluate all, then filter to sweep_props - no kwargs needed
                    all_vals = mat.evaluate(T_mat, temp)
                    row = {p: all_vals[p] for p in sweep_props if p in all_vals}
                    row['temperature'] = temp
                    sweep_data.append(row)
                except Exception as e:
                    print(f"  Error at {temp}K: {e}")

            if not sweep_data:
                continue
            df = pd.DataFrame(sweep_data).set_index('temperature')
            print(f"Sweep ({sweep_props}):")
            print(df.head(5).to_string())

            for col in sweep_props:
                if col in df.columns:
                    series = df[col].dropna()
                    if len(series) >= 2 and series.iloc[0] != 0:
                        change = ((series.iloc[-1] - series.iloc[0]) / series.iloc[0]) * 100
                        print(f"  {col} change {start}K -> {stop - step}K: {change:.2f}%")


# ------------------------------------------------------------------
# Config factory functions
# ------------------------------------------------------------------

def create_quick_test_config() -> TestConfig:
    return TestConfig(
        yaml_validation=True,
        material_creation=True,
        material_info_extraction=False,
        supported_properties_overview=False,
        material_property_analysis=True,
        property_evaluation_specific_temp=True,
        new_api_methods=False,
        batch_temperature_evaluation=False,
        inverse_function_testing=False,
        advanced_usage=False,
        manual_evaluation=True,
        object_oriented_evaluation=True,
        functional_evaluation=False,
        include_aluminum=True,
        include_steel=False,
        single_test_temperature=500.0,
        batch_temperatures=[400, 500, 600],
    )

def create_comprehensive_test_config() -> TestConfig:
    return TestConfig()

def create_everything_config() -> TestConfig:
    return TestConfig(
        single_test_temperature=500.15,
        batch_temperatures=[200, 273, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000],
        sweep_temp_range=(200, 3001, 50),
    )

def create_inverse_only_config() -> TestConfig:
    return TestConfig(
        yaml_validation=False,
        material_info_extraction=False,
        material_creation=True,
        supported_properties_overview=False,
        material_property_analysis=False,
        inverse_function_testing=True,
        advanced_usage=False,
    )

def create_api_comparison_config() -> TestConfig:
    return TestConfig(
        yaml_validation=False,
        material_info_extraction=False,
        material_creation=True,
        supported_properties_overview=False,
        material_property_analysis=True,
        property_evaluation_specific_temp=True,
        new_api_methods=True,
        batch_temperature_evaluation=False,
        inverse_function_testing=False,
        advanced_usage=False,
        manual_evaluation=True,
        object_oriented_evaluation=True,
        functional_evaluation=True,
    )

def run_everything() -> None:
    print(f"\n{'#' * 80}")
    print("RUNNING EVERYTHING - MAXIMUM COMPREHENSIVE TEST SUITE")
    print(f"{'#' * 80}")
    MaterialPropertyDemonstrator(create_everything_config()).run_all_tests()
    print(f"\n{'#' * 80}")
    print("EVERYTHING COMPLETED")
    print(f"{'#' * 80}")


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

if __name__ == "__main__":
    print("MaterForge Material Property Demonstration")
    print("=" * 50)
    print("Available test configurations:")
    print("1. Quick test (minimal, fast)")
    print("2. Comprehensive test (all standard features)")
    print("3. Inverse functions only")
    print("4. API comparison")
    print("5. Custom configuration")
    print("6. RUN EVERYTHING (maximum comprehensive test)")
    print("0. Exit")

    choice = input("\nSelect configuration (0-6, default=2): ").strip() or "2"

    configs = {
        "1": create_quick_test_config,
        "2": create_comprehensive_test_config,
        "3": create_inverse_only_config,
        "4": create_api_comparison_config,
    }

    if choice == "0":
        print("Exiting...")
    elif choice in configs:
        MaterialPropertyDemonstrator(configs[choice]()).run_all_tests()
    elif choice == "5":
        print("\nUsing custom configuration...")
        MaterialPropertyDemonstrator(TestConfig(
            yaml_validation=True,
            material_creation=True,
            material_property_analysis=True,
            inverse_function_testing=True,
            include_aluminum=True,
            include_steel=True,
            single_test_temperature=600.0,
            batch_temperatures=[300, 600, 900, 1200],
        )).run_all_tests()
    elif choice == "6":
        run_everything()
    else:
        print(f"Invalid choice '{choice}' - using default comprehensive config.")
        MaterialPropertyDemonstrator(create_comprehensive_test_config()).run_all_tests()
