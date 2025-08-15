"""Demonstration script for material property evaluation."""
import logging
from pathlib import Path
import sympy as sp
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from materforge.parsing.api import (
    create_material,
    get_supported_properties,
    evaluate_material_properties,
    get_material_property_names,
    validate_yaml_file,
    get_material_info
)
from materforge.algorithms.piecewise_inverter import PiecewiseInverter


@dataclass
class TestConfig:
    """Configuration for controlling which tests to run."""
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

    # Detailed subtests for property analysis
    manual_evaluation: bool = True
    object_oriented_evaluation: bool = True
    functional_evaluation: bool = True
    specific_properties_only: bool = True
    temperature_dependent_only: bool = True

    # Material selection
    include_aluminum: bool = True
    include_steel: bool = True

    # Test parameters
    single_test_temperature: float = 500.15  # Kelvin
    batch_temperatures: List[float] = None # type: ignore
    sweep_temp_range: Tuple[int, int, int] = (300, 3001, 100)  # start, stop, step

    def __post_init__(self):
        if self.batch_temperatures is None:
            self.batch_temperatures = [300, 400, 500, 600, 700]


class MaterialPropertyDemonstrator:
    """Modular demonstration class for material property evaluation."""

    def __init__(self, config: TestConfig):
        self.config = config
        self.materials = []
        self.material_symbols = {}
        self.setup_logging()
        self.setup_paths()

    def setup_logging(self):
        """Setup logging configuration."""
        logging.basicConfig(
            level=logging.WARNING,  # DEBUG/INFO/WARNING/ERROR/CRITICAL
            format="%(asctime)s %(levelname)s %(name)s -> %(message)s"
        )
        # Silence noisy libraries
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        logging.getLogger('PIL').setLevel(logging.WARNING)
        logging.getLogger('fontTools').setLevel(logging.WARNING)

    def setup_paths(self):
        """Setup file paths for materials."""
        current_file = Path(__file__)
        base_path = current_file.parent.parent / "src" / "materforge" / "data" / "materials"

        self.yaml_paths = {}
        if self.config.include_aluminum:
            self.yaml_paths["Aluminum"] = base_path / "pure_metals" / "Al" / "Al.yaml"
        if self.config.include_steel:
            self.yaml_paths["Steel_1.4301"] = base_path / "alloys" / "1.4301" / "1.4301.yaml"

    def run_all_tests(self):
        """Run all enabled tests."""
        print(f"\n{'=' * 80}")
        print("PYMATLIB MATERIAL PROPERTY DEMONSTRATION")
        print(f"{'=' * 80}")

        # Print configuration summary
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

    def _print_config_summary(self):
        """Print a summary of the current configuration."""
        print(f"\nTest Configuration Summary:")
        print(f"{'-' * 40}")

        # Count enabled tests
        test_flags = [
            self.config.yaml_validation,
            self.config.material_info_extraction,
            self.config.material_creation,
            self.config.supported_properties_overview,
            self.config.material_property_analysis,
            self.config.inverse_function_testing,
            self.config.advanced_usage
        ]
        enabled_tests = sum(test_flags)

        print(f"Enabled Tests: {enabled_tests}/7")
        print(f"Materials: {'Al' if self.config.include_aluminum else ''} {'Steel' if self.config.include_steel else ''}")
        print(f"Test Temperature: {self.config.single_test_temperature}K")

        if self.config.batch_temperature_evaluation:
            print(f"Batch Temperatures: {self.config.batch_temperatures}")

        if self.config.temperature_sweep_analysis:
            start, stop, step = self.config.sweep_temp_range
            print(f"Sweep Range: {start}-{stop}K (step: {step}K)")

    def test_yaml_validation(self):
        """Test YAML file validation."""
        print(f"\n{'1. YAML FILE VALIDATION':<50}")
        print(f"{'-' * 50}")

        for name, yaml_path in self.yaml_paths.items():
            if yaml_path.exists():
                try:
                    is_valid = validate_yaml_file(yaml_path)
                    print(f"✓ {name} YAML validation: {'PASSED' if is_valid else 'FAILED'}")
                except Exception as e:
                    print(f"✗ {name} YAML validation: FAILED - {e}")
            else:
                print(f"✗ {name} YAML file not found: {yaml_path}")

    def test_material_info_extraction(self):
        """Test material info extraction."""
        print(f"\n{'2. MATERIAL INFO EXTRACTION':<50}")
        print(f"{'-' * 50}")

        for name, yaml_path in self.yaml_paths.items():
            if yaml_path.exists():
                try:
                    info = get_material_info(yaml_path)
                    print(f"\n{name} Information:")
                    print(f"  Name: {info['name']}")
                    print(f"  Type: {info['material_type']}")
                    print(f"  Total Properties: {info['total_properties']}")
                    print(f"  Composition: {info['composition']}")

                    # Temperature properties
                    if info['material_type'] == 'pure_metal':
                        print(f"  Melting Temperature: {info.get('melting_temperature', 'N/A')} K")
                        print(f"  Boiling Temperature: {info.get('boiling_temperature', 'N/A')} K")
                    else:
                        print(f"  Solidus Temperature: {info.get('solidus_temperature', 'N/A')} K")
                        print(f"  Liquidus Temperature: {info.get('liquidus_temperature', 'N/A')} K")

                except Exception as e:
                    print(f"✗ Failed to get {name} info: {e}")

    def test_material_creation(self):
        """Test material creation."""
        print(f"\n{'3. MATERIAL CREATION':<50}")
        print(f"{'-' * 50}")

        symbol_mapping = {
            "Aluminum": sp.Symbol('T_Al'),
            "Steel_1.4301": sp.Symbol('T_SS')
        }

        for name, yaml_path in self.yaml_paths.items():
            if yaml_path.exists():
                try:
                    T_symbol = symbol_mapping.get(name, sp.Symbol(f'T_{name}'))
                    mat = create_material(yaml_path=yaml_path, dependency=T_symbol, enable_plotting=True)
                    self.materials.append(mat)
                    self.material_symbols[mat.name] = T_symbol
                    print(f"✓ Successfully created: {mat.name} (using symbol {T_symbol})")
                except Exception as e:
                    print(f"✗ Failed to create {name}: {e}")

    def test_supported_properties_overview(self):
        """Test supported properties overview."""
        print(f"\n{'4. SUPPORTED PROPERTIES OVERVIEW':<50}")
        print(f"{'-' * 50}")

        supported_props = get_supported_properties()
        print(f"PyMatLib supports {len(supported_props)} property types:")
        for i, prop in enumerate(sorted(supported_props), 1):
            print(f"  {i:2d}. {prop}")

    def test_material_property_analysis(self):
        """Test material property analysis."""
        print(f"\n{'5. MATERIAL PROPERTY ANALYSIS':<50}")
        print(f"{'-' * 50}")

        for mat in self.materials:
            T_mat = self.material_symbols[mat.name]

            print(f"\n{'=' * 80}")
            print(f"MATERIAL: {mat.name} (Temperature symbol: {T_mat})")
            print(f"{'=' * 80}")

            self._analyze_basic_material_info(mat)
            self._analyze_available_properties(mat)

            if self.config.property_evaluation_specific_temp:
                self._test_property_evaluation_at_specific_temp(mat, T_mat)

            if self.config.batch_temperature_evaluation:
                self._test_batch_temperature_evaluation(mat)

    def _analyze_basic_material_info(self, mat):
        """Analyze basic material information."""
        print(f"Name: {mat.name}")
        print(f"Type: {mat.material_type}")
        print(f"Elements: {[elem.name for elem in mat.elements]}")
        print(f"Composition: {mat.composition}")
        for i in range(len(mat.composition)):
            print(f"  {mat.elements[i].name}: {mat.composition[i]}")

        # Temperature properties
        for attr in ['solidus_temperature', 'liquidus_temperature', 'melting_temperature', 'boiling_temperature']:
            if hasattr(mat, attr):
                print(f"{attr.replace('_', ' ').title()}: {getattr(mat, attr)}")

    def _analyze_available_properties(self, mat):
        """Analyze available properties on material."""
        print(f"\n{'AVAILABLE PROPERTIES:':<50}")
        print(f"{'-' * 50}")

        available_props = get_material_property_names(mat)
        print(f"Material has {len(available_props)} properties:")

        for prop_name in sorted(available_props):
            prop_value = getattr(mat, prop_name)
            print(f"{prop_name:<30}: {prop_value}")
            print(f"{' ' * 24} Type : {type(prop_value)}")

            # Validate property types
            allowed_types = (sp.Piecewise, sp.Float, type(None))
            if not isinstance(prop_value, allowed_types):
                print(f"{' ' * 30}  !  WARNING: Unexpected type")
            else:
                print(f"{' ' * 30}  ✓ Valid type")

    def _test_property_evaluation_at_specific_temp(self, mat, T_mat):
        """Test property evaluation at specific temperature."""
        test_temp = self.config.single_test_temperature
        available_props = get_material_property_names(mat)

        print(f"\n{'PROPERTY VALUES AT ' + str(test_temp) + 'K:':<50}")
        print(f"{'-' * 50}")

        # Method 1: Manual evaluation
        if self.config.manual_evaluation:
            print("Method 1: Manual Property Evaluation")
            for prop_name in sorted(available_props):
                try:
                    prop_value = getattr(mat, prop_name)
                    if isinstance(prop_value, sp.Expr):
                        numerical_value = prop_value.subs(T_mat, test_temp).evalf() # type: ignore
                        print(f"{prop_name:<30}: {numerical_value} (symbolic)")
                    else:
                        print(f"{prop_name:<30}: {prop_value} (constant)")
                except Exception as e:
                    print(f"{prop_name:<30}: Error - {str(e)}")

        if self.config.new_api_methods:
            self._test_new_api_methods(mat, test_temp, available_props)

    def _test_new_api_methods(self, mat, test_temp, available_props):
        """Test new API methods for property evaluation."""
        print(f"\n{'NEW API METHODS:':<50}")
        print(f"{'-' * 50}")

        # Method 2: Object-oriented approach
        if self.config.object_oriented_evaluation:
            print("Method 2: Material.evaluate_properties_at_temperature()")
            try:
                all_values = mat.evaluate_properties_at_temperature(test_temp)
                print(f"All properties at {test_temp}K:")
                for prop, value in sorted(all_values.items()):
                    print(f"  {prop:<28}: {value:.6e}")
            except Exception as e:
                print(f"Error: {e}")
                return

        # Method 3: Functional approach
        if self.config.functional_evaluation:
            print(f"\nMethod 3: evaluate_material_properties()")
            try:
                all_values_func = evaluate_material_properties(mat, test_temp)
                if self.config.object_oriented_evaluation:
                    print(f"Functional API results match: {all_values == all_values_func}") # type: ignore
                else:
                    print("Functional API evaluation completed successfully")
            except Exception as e:
                print(f"Error: {e}")

        # Method 4: Specific properties only
        if self.config.specific_properties_only:
            print(f"\nMethod 4: Specific Properties Only")
            try:
                specific_props = ['density', 'heat_capacity', 'energy_density']
                available_specific = [p for p in specific_props if p in available_props]

                if available_specific:
                    specific_values = mat.evaluate_properties_at_temperature(
                        test_temp,
                        properties=available_specific
                    )
                    print(f"Specific properties at {test_temp}K:")
                    for prop, value in sorted(specific_values.items()):
                        print(f"  {prop:<28}: {value:.6e}")
                else:
                    print("None of the requested specific properties available")
            except Exception as e:
                print(f"Error: {e}")

        # Method 5: Temperature-dependent properties only
        if self.config.temperature_dependent_only:
            print(f"\nMethod 5: Temperature-Dependent Properties Only")
            try:
                temp_dependent = mat.evaluate_properties_at_temperature(
                    test_temp,
                    include_constants=False
                )
                print(f"Temperature-dependent properties at {test_temp}K:")
                for prop, value in sorted(temp_dependent.items()):
                    print(f"  {prop:<28}: {value:.6e}")
            except Exception as e:
                print(f"Error: {e}")

    def _test_batch_temperature_evaluation(self, mat):
        """Test batch temperature evaluation."""
        print(f"\n{'BATCH TEMPERATURE EVALUATION:':<50}")
        print(f"{'-' * 50}")

        batch_results = {}
        for temp in self.config.batch_temperatures:
            try:
                batch_results[temp] = mat.evaluate_properties_at_temperature(temp)
            except Exception as e:
                print(f"Error at {temp}K: {e}")
                batch_results[temp] = {}

        # Create DataFrame for analysis
        if batch_results:
            try:
                df = pd.DataFrame(batch_results).T  # Transpose so temperatures are rows
                print("Temperature variation analysis:")
                print(f"{'Temperature':<12} {'density':<12} {'heat_capacity':<15} {'energy_density':<15}")
                print("-" * 55)

                for temp, row in df.iterrows():
                    density_val = row.get('density', 'N/A')
                    heat_cap_val = row.get('heat_capacity', 'N/A')
                    energy_val = row.get('energy_density', 'N/A')

                    density_str = f"{density_val:.2e}" if density_val != 'N/A' else 'N/A'
                    heat_str = f"{heat_cap_val:.2e}" if heat_cap_val != 'N/A' else 'N/A'
                    energy_str = f"{energy_val:.2e}" if energy_val != 'N/A' else 'N/A'

                    print(f"{temp:<12.0f} {density_str:<12} {heat_str:<15} {energy_str:<15}")

            except Exception as e:
                print(f"DataFrame creation failed: {e}")

    def test_inverse_functions(self):
        """Test inverse function creation and accuracy."""
        print(f"\n{'6. INVERSE FUNCTION TESTING':<50}")
        print(f"{'-' * 50}")

        print(f"\n{'=' * 80}")
        print("INVERSE FUNCTION TESTING")
        print(f"{'=' * 80}")

        for mat in self.materials:
            T_mat = self.material_symbols[mat.name]

            print(f"\n--- Testing Inverse Functions for {mat.name} (symbol: {T_mat}) ---")

            if not hasattr(mat, 'energy_density') or mat.energy_density is None:
                print(f"!  {mat.name} has no energy_density property - skipping inverse test")
                continue

            try:
                # Method 1: Try convenience function
                print("Method 1: Convenience function...")
                try:
                    E = sp.Symbol('E')
                    inverse_func1 = PiecewiseInverter.create_energy_density_inverse(mat, 'E')
                    print("✓ Method 1 succeeded!")
                    self._test_round_trip_accuracy(mat, inverse_func1, T_mat, E, method="Method 1")
                except Exception as e:
                    print(f"✗ Method 1 failed: {e}")

                # Method 2: Direct approach
                print("Method 2: Direct approach...")
                try:
                    energy_symbols = mat.energy_density.free_symbols
                    if len(energy_symbols) == 1:
                        temp_symbol = list(energy_symbols)[0]
                        E_symbol = sp.Symbol('E')
                        inverse_func2 = PiecewiseInverter.create_inverse(
                            mat.energy_density, temp_symbol, E_symbol
                        )
                        print(f"✓ Method 2 succeeded!")
                        self._test_round_trip_accuracy(mat, inverse_func2, temp_symbol, E_symbol, method="Method 2")
                    else:
                        print(f"✗ Unexpected symbols in energy density: {energy_symbols}")
                except Exception as e:
                    print(f"✗ Method 2 failed: {e}")

            except Exception as e:
                print(f"✗ Inverse testing failed for {mat.name}: {e}")

    def _test_round_trip_accuracy(self, material, inverse_func, temp_symbol, energy_symbol, method="Unknown"):
        """Test round-trip accuracy: T -> E -> T."""
        print(f"  Round-trip accuracy test ({method}):")

        test_temperatures = [300, 350, 400, 450, 500, 600, 700, 800, 900, 1000,
                             1200, 1500, 1800, 2000, 2500, 3000]

        max_error = 0.0
        successful_tests = 0

        for temp in test_temperatures:
            try:
                # Forward: T -> E
                energy_val = float(material.energy_density.subs(temp_symbol, temp))

                # Backward: E -> T
                recovered_temp = float(inverse_func.subs(energy_symbol, energy_val))

                # Calculate error
                error = abs(temp - recovered_temp)
                max_error = max(max_error, error)
                successful_tests += 1

                print(f"    T={temp:4.0f}K -> E={energy_val:.2e} -> T={recovered_temp:6.1f}K (Error: {error:.2e})")

            except Exception as e:
                print(f"    T={temp:4.0f}K -> Error: {e}")

        print(f"  Tests completed: {successful_tests}/{len(test_temperatures)}")
        print(f"  Maximum error: {max_error:.2e} K")

        if max_error < 1e-8:
            print("  ✓ Excellent accuracy!")
        elif max_error < 1e-4:
            print("  ✓ Good accuracy")
        else:
            print("  !  Consider reviewing inverse function accuracy")

    def test_advanced_usage(self):
        """Test advanced usage patterns."""
        print(f"\n{'=' * 80}")
        print("ADVANCED USAGE PATTERNS")
        print(f"{'=' * 80}")

        if self.config.temperature_sweep_analysis:
            self._test_temperature_sweep_analysis()

    def _test_temperature_sweep_analysis(self):
        """Test temperature sweep analysis."""
        print("\n7. TEMPERATURE SWEEP ANALYSIS")
        print("-" * 50)

        start, stop, step = self.config.sweep_temp_range
        temp_range = range(start, stop, step)

        for mat in self.materials:
            T_symbol = self.material_symbols[mat.name]
            mat_name = mat.name

            print(f"\n--- Temperature Sweep for {mat_name} ---")

            sweep_data = []
            for temp in temp_range:
                try:
                    props = mat.evaluate_properties_at_temperature(
                        temp,
                        properties=['density', 'heat_capacity', 'thermal_diffusivity'],
                        include_constants=False
                    )
                    props['temperature'] = temp
                    sweep_data.append(props)
                except Exception as e:
                    print(f"Error at {temp}K: {e}")

            if sweep_data:
                df = pd.DataFrame(sweep_data)
                print(f"Temperature sweep results for {mat_name}:")
                print(df.head(5))

                # Simple analysis
                if 'density' in df.columns:
                    density_change = ((df['density'].iloc[-1] - df['density'].iloc[0]) / df['density'].iloc[0]) * 100
                    print(f"Density change from {temp_range[0]}K to {temp_range[-1]}K: {density_change:.2f}%")


def create_quick_test_config() -> TestConfig:
    """Create a configuration for quick testing (minimal tests)."""
    return TestConfig(
        yaml_validation=True,
        material_creation=True,
        material_property_analysis=True,
        property_evaluation_specific_temp=True,
        new_api_methods=False,  # Skip detailed API tests
        batch_temperature_evaluation=False,
        inverse_function_testing=False,
        advanced_usage=False,
        # Subtests
        manual_evaluation=True,
        object_oriented_evaluation=True,
        functional_evaluation=False,
        specific_properties_only=False,
        temperature_dependent_only=False,
        # Materials
        include_aluminum=True,
        include_steel=False,  # Only test aluminum for speed
        # Parameters
        single_test_temperature=500.0,
        batch_temperatures=[400, 500, 600]
    )


def create_comprehensive_test_config() -> TestConfig:
    """Create a configuration for comprehensive testing (all tests)."""
    return TestConfig()  # All defaults are True


def create_everything_config() -> TestConfig:
    """Create a configuration that runs absolutely everything with maximum detail."""
    return TestConfig(
        # Enable all main tests
        yaml_validation=True,
        material_info_extraction=True,
        material_creation=True,
        supported_properties_overview=True,
        material_property_analysis=True,
        property_evaluation_specific_temp=True,
        new_api_methods=True,
        batch_temperature_evaluation=True,
        inverse_function_testing=True,
        advanced_usage=True,
        temperature_sweep_analysis=True,

        # Enable all sub-tests
        manual_evaluation=True,
        object_oriented_evaluation=True,
        functional_evaluation=True,
        specific_properties_only=True,
        temperature_dependent_only=True,

        # Include all materials
        include_aluminum=True,
        include_steel=True,

        # Extended test parameters for maximum coverage
        single_test_temperature=500.15,
        batch_temperatures=[200, 273, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000],
        sweep_temp_range=(200, 3001, 50)  # More detailed sweep: 200K to 3000K in 50K steps
    )


def create_inverse_only_config() -> TestConfig:
    """Create a configuration for testing only inverse functions."""
    return TestConfig(
        yaml_validation=False,
        material_info_extraction=False,
        material_creation=True,  # Need this for inverse tests
        supported_properties_overview=False,
        material_property_analysis=False,
        inverse_function_testing=True,
        advanced_usage=False
    )


def create_api_comparison_config() -> TestConfig:
    """Create a configuration for comparing different API methods."""
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
        # All API comparison methods
        manual_evaluation=True,
        object_oriented_evaluation=True,
        functional_evaluation=True,
        specific_properties_only=True,
        temperature_dependent_only=True
    )


def run_with_config(config_func, description: str):
    """Helper function to run tests with a specific configuration."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")

    config = config_func()
    demonstrator = MaterialPropertyDemonstrator(config)
    demonstrator.run_all_tests()


def run_everything():
    """Run absolutely everything - the most comprehensive test possible."""
    print(f"\n{'#'*80}")
    print("RUNNING EVERYTHING - MAXIMUM COMPREHENSIVE TEST SUITE")
    print(f"{'#'*80}")
    print("\nThis will run all available tests with maximum detail and coverage.")
    print("This may take several minutes to complete...")

    config = create_everything_config()
    demonstrator = MaterialPropertyDemonstrator(config)
    demonstrator.run_all_tests()

    print(f"\n{'#'*80}")
    print("EVERYTHING COMPLETED - ALL TESTS FINISHED")
    print(f"{'#'*80}")


if __name__ == "__main__":
    # Enhanced menu with "Run Everything" option
    print("PyMatLib Material Property Demonstration")
    print("="*50)
    print("Available test configurations:")
    print("1. Quick test (minimal, fast execution)")
    print("2. Comprehensive test (all standard features)")
    print("3. Inverse functions only")
    print("4. API comparison")
    print("5. Custom configuration")
    print("6. RUN EVERYTHING (maximum comprehensive test)")
    print("0. Exit")

    choice = input("\nSelect configuration (0-6, default=2): ").strip() or "2"

    if choice == "0":
        print("Exiting...")
        exit()
    elif choice == "1":
        config = create_quick_test_config()
        demonstrator = MaterialPropertyDemonstrator(config)
        demonstrator.run_all_tests()
    elif choice == "2":
        config = create_comprehensive_test_config()
        demonstrator = MaterialPropertyDemonstrator(config)
        demonstrator.run_all_tests()
    elif choice == "3":
        config = create_inverse_only_config()
        demonstrator = MaterialPropertyDemonstrator(config)
        demonstrator.run_all_tests()
    elif choice == "4":
        config = create_api_comparison_config()
        demonstrator = MaterialPropertyDemonstrator(config)
        demonstrator.run_all_tests()
    elif choice == "5":
        # Custom configuration example
        print("\nUsing custom configuration...")
        config = TestConfig(
            yaml_validation=True,
            material_creation=True,
            material_property_analysis=True,
            inverse_function_testing=True,
            include_aluminum=True,
            include_steel=True,
            single_test_temperature=600.0,
            batch_temperatures=[300, 600, 900, 1200]
        )
        demonstrator = MaterialPropertyDemonstrator(config)
        demonstrator.run_all_tests()
    elif choice == "6":
        # RUN EVERYTHING option
        """confirm = input("\nThis will run the most comprehensive test suite possible.\n"
                        "It may take several minutes. Continue? (y/N): ").strip().lower()
        if confirm in ['y', 'yes']:
            run_everything()
        else:
            print("Cancelled.")"""
        run_everything()
    else:
        print(f"Invalid choice: {choice}")
        print("Using default comprehensive configuration...")
        config = create_comprehensive_test_config()
        demonstrator = MaterialPropertyDemonstrator(config)
        demonstrator.run_all_tests()
