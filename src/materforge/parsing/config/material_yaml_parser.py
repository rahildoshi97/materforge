import logging
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import sympy as sp
from ruamel.yaml import YAML, constructor, scanner

from materforge.core.materials import Material
from materforge.data import ProcessingConstants
from materforge.parsing.processors.property_processor import PropertyProcessor
from materforge.parsing.validation.property_type_detector import PropertyType, PropertyTypeDetector
from materforge.visualization.plotters import PropertyVisualizer
from materforge.parsing.config.yaml_keys import PROPERTIES_KEY, MATERIAL_TYPE_KEY, \
    COMPOSITION_KEY, PURE_METAL_KEY, MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY, SOLIDUS_TEMPERATURE_KEY, \
    LIQUIDUS_TEMPERATURE_KEY, INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY, ALLOY_KEY, NAME_KEY, \
    INDEPENDENT_VARIABLES_KEY

logger = logging.getLogger(__name__)


class BaseFileParser:
    """Base class for parsing configuration files."""

    def __init__(self, config_path: Union[str, Path]) -> None:
        self.config_path = Path(config_path)
        self.base_dir = self.config_path.parent
        self.config = self._load_config()
        logger.info("Successfully loaded configuration from: %s", self.config_path)

    def _load_config(self) -> Dict[str, Any]:
        raise NotImplementedError("Subclasses must implement _load_config method")


class YAMLFileParser(BaseFileParser):
    """Parser for YAML configuration files."""

    def _load_config(self) -> Dict[str, Any]:
        yaml = YAML(typ='safe')
        yaml.allow_duplicate_keys = False
        try:
            logger.debug("Loading YAML file: %s", self.config_path)
            with open(self.config_path, 'r') as f:
                config = yaml.load(f)
            logger.debug("YAML file loaded successfully, found %d top-level keys", len(config) if config else 0)
            return config
        except FileNotFoundError as e:
            logger.error("YAML file not found: %s", self.config_path)
            raise FileNotFoundError(f"YAML file not found: {self.config_path}") from e
        except constructor.DuplicateKeyError as e:
            logger.error("Duplicate key found in YAML file %s: %s", self.config_path, e)
            raise constructor.DuplicateKeyError(f"Duplicate key in {self.config_path}: {str(e)}") from e
        except scanner.ScannerError as e:
            logger.error("YAML syntax error in file %s: %s", self.config_path, e)
            raise scanner.ScannerError(f"YAML syntax error in {self.config_path}: {str(e)}") from e
        except Exception as e:
            logger.error("Unexpected error parsing YAML file %s: %s", self.config_path, e, exc_info=True)
            raise ValueError(f"Error parsing {self.config_path}: {str(e)}") from e


class MaterialYAMLParser(YAMLFileParser):
    """Parser for material configuration files in YAML format."""

    VALID_YAML_PROPERTIES = {
        "bulk_modulus",
        "density",
        "dynamic_viscosity",
        "elastic_modulus",
        "electrical_conductivity",
        "electrical_resistivity",
        "energy_density",
        "fracture_toughness",
        "hardness",
        "heat_capacity",
        "heat_conductivity",
        "kinematic_viscosity",
        "latent_heat_of_fusion",
        "latent_heat_of_vaporization",
        "magnetic_permeability",
        "poisson_ratio",
        "shear_modulus",
        "specific_enthalpy",
        "surface_tension",
        "thermal_diffusivity",
        "thermal_expansion_coefficient",
        "ultimate_tensile_strength",
        "viscosity",
        "yield_strength",
        # Extend to include other material properties as needed
    }

    # --- Constructor ---
    def __init__(self, yaml_path: Union[str, Path]) -> None:
        super().__init__(yaml_path)
        logger.info("Initializing MaterialYAMLParser for: %s", yaml_path)
        self._validate_config()
        self.categorized_properties = self._analyze_and_categorize_properties(self.config[PROPERTIES_KEY])
        self.property_processor = PropertyProcessor()
        self.visualizer = PropertyVisualizer(self)
        logger.info("MaterialYAMLParser initialized successfully with %d property categories",
                    len(self.categorized_properties))

    # --- Public API ---
    def create_material(self, dependency: Union[float, sp.Symbol], enable_plotting: bool = True) -> Material:
        """Create a Material instance from the parsed configuration and temperature."""
        logger.info("Creating material from configuration: %s", self.config_path)
        try:
            name = self.config.get(NAME_KEY, "Unnamed Material")
            material_type = self.config[MATERIAL_TYPE_KEY]
            logger.info("Creating material: %s (type: %s)", name, material_type)
            elements = self._get_elements()
            composition = [val for val in self.config[COMPOSITION_KEY].values()]
            # Create material with different parameters based on material_type
            if material_type == PURE_METAL_KEY:
                material = Material(
                    name=name,
                    elements=elements,
                    composition=composition,
                    material_type=material_type,
                    melting_temperature=sp.Float(self.config[MELTING_TEMPERATURE_KEY]),
                    boiling_temperature=sp.Float(self.config[BOILING_TEMPERATURE_KEY]),
                )
                logger.debug("Created pure metal with melting temp: %s K, boiling temp: %s K",
                             self.config[MELTING_TEMPERATURE_KEY], self.config[BOILING_TEMPERATURE_KEY])
            else:  # alloy
                material = Material(
                    name=name,
                    elements=elements,
                    composition=composition,
                    material_type=material_type,
                    solidus_temperature=sp.Float(self.config[SOLIDUS_TEMPERATURE_KEY]),
                    liquidus_temperature=sp.Float(self.config[LIQUIDUS_TEMPERATURE_KEY]),
                    initial_boiling_temperature=sp.Float(self.config[INITIAL_BOILING_TEMPERATURE_KEY]),
                    final_boiling_temperature=sp.Float(self.config[FINAL_BOILING_TEMPERATURE_KEY]),
                )
                logger.debug("Created alloy with solidus: %s K, liquidus: %s K",
                             self.config[SOLIDUS_TEMPERATURE_KEY], self.config[LIQUIDUS_TEMPERATURE_KEY])
            # Initialize visualizer only if plotting is enabled AND dependency is symbolic
            visualizer = None
            should_visualize = enable_plotting and isinstance(dependency, sp.Symbol)
            if should_visualize:
                self.visualizer.initialize_plots()
                self.visualizer.reset_visualization_tracking()
                visualizer = self.visualizer
                logger.info("Visualization enabled for symbolic temperature")
            else:
                if not isinstance(dependency, sp.Symbol):
                    logger.debug("Visualization disabled - numeric temperature provided")
                else:
                    logger.debug("Visualization disabled - plotting not enabled")
            # Process properties
            logger.info("Starting property processing for material: %s", name)
            self.property_processor.process_properties(
                material=material,
                dependency=dependency,
                properties=self.config[PROPERTIES_KEY],
                categorized_properties=self.categorized_properties,
                base_dir=self.base_dir,
                visualizer=visualizer
            )
            # Save plots only if visualization was actually enabled
            if should_visualize and visualizer is not None:
                logger.info("Saving property plots for material: %s", name)
                self.visualizer.save_property_plots()
                logger.info(f"Property plots saved for {material}")
            logger.info(f"Successfully created material: {name}")
            return material
        except KeyError as e:
            logger.error("Configuration error for material creation - missing key: %s", e, exc_info=True)
            raise ValueError(f"Configuration error: Missing {str(e)}") from e
        except Exception as e:
            logger.error("Failed to create material from %s: %s", self.config_path, e, exc_info=True)
            raise ValueError(f"Failed to create material \n -> {str(e)}") from e

    # --- Validation Methods ---
    def _validate_config(self) -> None:
        """Validate the configuration structure and content."""
        logger.debug("Starting configuration validation")
        if not isinstance(self.config, dict):
            logger.error("Invalid YAML structure - expected dictionary at root level")
            raise ValueError("The YAML file must start with a dictionary/object structure with key-value pairs,"
                             "not a list or scalar value")
        self._validate_required_fields()
        properties = self.config.get(PROPERTIES_KEY, {})
        if not isinstance(properties, dict):
            logger.error("Properties section is not a dictionary: %s", type(properties))
            raise ValueError("The 'properties' section in your YAML file must be a dictionary with key-value pairs")
        self._validate_property_names(properties)
        logger.info("Configuration validation completed successfully")

    def _validate_required_fields(self) -> None:
        """Validate that all required fields are present."""
        logger.debug("Validating required fields")
        # Check for material_type first
        if MATERIAL_TYPE_KEY not in self.config:
            logger.error("Missing required field: %s", MATERIAL_TYPE_KEY)
            raise ValueError("Missing required field: material_type")
        material_type = self.config[MATERIAL_TYPE_KEY]
        logger.debug("Material type: %s", material_type)
        if material_type not in [ALLOY_KEY, PURE_METAL_KEY]:
            logger.error("Invalid material_type: %s", material_type)
            raise ValueError(f"Invalid material_type: {material_type}. Must be {PURE_METAL_KEY} or {ALLOY_KEY}")
        # Common required fields
        common_fields = {NAME_KEY, MATERIAL_TYPE_KEY, COMPOSITION_KEY, PROPERTIES_KEY}
        # Material-type specific fields
        if material_type == PURE_METAL_KEY:
            required_fields = common_fields | {MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY}
        elif material_type == ALLOY_KEY:
            required_fields = common_fields | {SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
                                               INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY}
        else:
            logger.error("Unsupported material_type: %s", material_type)
            raise ValueError(f"Unsupported material_type: {material_type}. "
                             f"Supported types are: {PURE_METAL_KEY}, {ALLOY_KEY}.")
        missing_fields = required_fields - set(self.config.keys())
        if missing_fields:
            logger.error("Missing required fields for %s: %s", material_type, missing_fields)
            raise ValueError(f"Missing required fields for {material_type}: {', '.join(missing_fields)}")
        extra_fields = set(self.config.keys()) - required_fields
        if extra_fields:
            logger.warning("Extra fields found in configuration: %s", extra_fields)
            suggestions = {
                field: get_close_matches(field, required_fields, n=1, cutoff=0.6)
                for field in extra_fields
            }
            error_msg = "Extra fields found in configuration: \n ->"
            for field, matches in suggestions.items():
                suggestion = f" (did you mean '{matches[0]}'?)" if matches else ""
                error_msg += f" - '{field}'{suggestion}\n"
            raise ValueError(error_msg)
        self._validate_name()
        self._validate_composition()
        logger.debug("Required fields validation completed")

    def _validate_name(self) -> None:
        """Validate the material name."""
        logger.debug("Validating material name")
        name = self.config[NAME_KEY]
        if name is None:
            logger.error("Material name is None")
            raise ValueError("Material name cannot be None")
        if not isinstance(name, str):
            logger.error("Material name is not a string: %s (type: %s)", name, type(name))
            raise ValueError(f"Material name must be a string, got {type(name).__name__} (value: {name})")
        if len(name) > 100:
            logger.warning("Material name '%s' exceeds 100 characters", name)
        logger.debug("Material name validation completed successfully")

    def _validate_composition(self) -> None:
        """Validate composition for both pure metals and alloys."""
        logger.debug("Validating composition")
        composition = self.config.get(COMPOSITION_KEY, {})
        material_type = self.config[MATERIAL_TYPE_KEY]
        if not isinstance(composition, dict):
            logger.error("Composition is not a dictionary: %s", type(composition))
            raise ValueError("Composition must be a dictionary")
        if not composition:
            logger.error("Composition is empty")
            raise ValueError("Composition cannot be empty")
        logger.debug("Composition contains %d elements: %s", len(composition), list(composition.keys()))
        # Check that all fractions are valid numbers
        for element, fraction in composition.items():
            if not isinstance(fraction, (int, float)):
                logger.error("Invalid composition fraction for '%s': %s (type: %s)",
                             element, fraction, type(fraction))
                raise ValueError(
                    f"Composition fraction for '{element}' must be a number, got {type(fraction).__name__}")
            if fraction < 0:
                logger.error("Negative composition fraction for '%s': %s", element, fraction)
                raise ValueError(f"Composition fraction for '{element}' cannot be negative, got {fraction}")
            if fraction > 1:
                logger.error("Composition fraction exceeds 1.0 for '%s': %s", element, fraction)
                raise ValueError(f"Composition fraction for '{element}' cannot exceed 1.0, got {fraction}")
        # Check that fractions sum to 1.0
        total = sum(composition.values())
        if not abs(total - 1.0) < ProcessingConstants.COMPOSITION_THRESHOLD:
            logger.error("Composition fractions sum to %s, expected 1.0", total)
            raise ValueError(f"Composition fractions must sum to 1.0, got {total}")
        # Material-type specific validation
        if material_type == PURE_METAL_KEY:
            self._validate_pure_metal_composition_rules(composition)
        else:  # alloy
            self._validate_alloy_composition_rules(composition)
        logger.debug("Composition validation completed successfully")

    @staticmethod
    def _validate_pure_metal_composition_rules(composition: dict) -> None:
        """Validate composition rules specific to pure metals."""
        logger.debug("Validating pure metal composition rules")
        # Count non-zero elements
        non_zero_elements = {element: fraction for element, fraction in composition.items()
                             if fraction > ProcessingConstants.COMPOSITION_THRESHOLD}
        if len(non_zero_elements) == 0:
            logger.error("Pure metal has no elements with non-zero composition")
            raise ValueError("Pure metals must have at least one element with non-zero composition")

        if len(non_zero_elements) > 1:
            element_list = ", ".join(f"{elem}: {frac}" for elem, frac in non_zero_elements.items())
            logger.error("Pure metal has multiple non-zero elements: %s", element_list)
            raise ValueError(
                f"Pure metals must contain exactly one element with composition 1.0. "
                f"Found multiple non-zero elements: {element_list}. "
                f"Use material_type: 'alloy' for multi-element materials."
            )
        # Check that the single element has composition 1.0
        single_element, single_fraction = list(non_zero_elements.items())[0]
        if not abs(single_fraction - 1.0) < ProcessingConstants.COMPOSITION_THRESHOLD:
            logger.error("Pure metal element '%s' has composition %s, expected 1.0",
                         single_element, single_fraction)
            raise ValueError(
                f"Pure metal element '{single_element}' must have composition 1.0, "
                f"got {single_fraction}. Use material_type: 'alloy' for fractional compositions."
            )
        # ERROR for zero-valued elements in pure metals
        zero_elements = [element for element, fraction in composition.items() if fraction == 0.0]
        if zero_elements:
            logger.error("Pure metal composition includes zero-valued elements: %s", zero_elements)
            raise ValueError(
                f"Pure metal composition should not include zero-valued elements: {zero_elements}. "
                f"Remove these elements from the composition dictionary."
            )
        logger.debug("Pure metal composition rules validated successfully")

    @staticmethod
    def _validate_alloy_composition_rules(composition: dict) -> None:
        """Validate composition rules specific to alloys."""
        logger.debug("Validating alloy composition rules")
        non_zero_elements = {element: fraction for element, fraction in composition.items()
                             if fraction > 1e-10}
        if len(non_zero_elements) < 2:
            if len(non_zero_elements) == 1:
                single_element = list(non_zero_elements.keys())[0]
                logger.error("Alloy has only one non-zero element: %s", single_element)
                raise ValueError(
                    f"Alloys must have at least 2 elements with non-zero composition. "
                    f"Found only '{single_element}'. Use material_type: 'pure_metal' for single elements."
                )
            else:
                logger.error("Alloy has no elements with non-zero composition")
                raise ValueError("Alloys must have at least 2 elements with non-zero composition")
        # Warning for zero-valued elements in alloys - they might be intentional
        zero_elements = [element for element, fraction in composition.items() if fraction == 0.0]
        if zero_elements:
            logger.warning(
                "Alloy composition includes zero-valued elements: %s. Consider removing if not needed",
                zero_elements
            )
        logger.debug("Alloy composition rules validated successfully")

    def _validate_property_names(self, properties: Dict[str, Any]) -> None:
        """Validate that all property names are supported."""
        logger.debug("Validating property names for %d properties", len(properties))
        invalid_props = set(properties.keys()) - self.VALID_YAML_PROPERTIES
        if invalid_props:
            logger.error("Invalid properties found: %s", invalid_props)
            suggestions = {
                prop: get_close_matches(prop, self.VALID_YAML_PROPERTIES, n=1, cutoff=0.6)
                for prop in invalid_props
            }
            error_msg = "Invalid properties found: \n ->"
            for prop, matches in suggestions.items():
                suggestion = f" (did you mean '{matches[0]}'?)" if matches else ""
                error_msg += f" - '{prop}'{suggestion}\n"
            raise ValueError(error_msg)
    logger.debug("Property names validation completed successfully")

    # --- Processing Methods ---
    def _get_elements(self) -> List:
        """Get element objects from composition keys."""
        from materforge.data.elements.element_data import element_map
        element_symbols = list(self.config[COMPOSITION_KEY].keys())
        logger.debug("Looking up elements: %s", element_symbols)
        try:
            elements = [element_map[sym] for sym in element_symbols]
            logger.debug("Successfully found all %d elements", len(elements))
            return elements
        except KeyError as e:
            logger.error("Invalid element symbol: %s", e)
            raise ValueError(f"Invalid element symbol: {str(e)}") from e

    @staticmethod
    def _analyze_and_categorize_properties(properties: Dict[str, Any]) -> Dict[PropertyType, List[Tuple[str, Any]]]:
        """Categorizes properties after detecting and validating their types."""
        logger.debug("Analyzing and categorizing %d properties", len(properties))
        categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]] = {
            prop_type: [] for prop_type in PropertyType
        }
        for prop_name, config in properties.items():
            try:
                logger.debug("Processing property: %s", prop_name)
                prop_type = PropertyTypeDetector.determine_property_type(prop_name, config)
                PropertyTypeDetector.validate_property_config(prop_name, config, prop_type)
                categorized_properties[prop_type].append((prop_name, config))
                logger.debug("Property '%s' categorized as: %s", prop_name, prop_type.name)
            except ValueError as e:
                logger.error("Configuration error for property '%s': %s", prop_name, e)
                raise ValueError(f"Configuration error for property '{prop_name}': {str(e)}") from e
        # Log summary
        for prop_type, prop_list in categorized_properties.items():
            if prop_list:
                logger.info("Found %d properties of type %s: %s",
                            len(prop_list), prop_type.name, [p[0] for p in prop_list])
        logger.debug(f"Categorized properties: {categorized_properties}")
        return categorized_properties
