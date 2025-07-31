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
from materforge.parsing.processors.dependency_resolver import DependencyResolver
from materforge.visualization.plotters import PropertyVisualizer
from materforge.parsing.config.yaml_keys import (
    PROPERTIES_KEY, MATERIAL_TYPE_KEY, COMPOSITION_KEY, PURE_METAL_KEY,
    MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY, SOLIDUS_TEMPERATURE_KEY,
    LIQUIDUS_TEMPERATURE_KEY, INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY,
    ALLOY_KEY, NAME_KEY, INDEPENDENT_VARIABLES_KEY
)

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
    """Enhanced parser for material configuration files with multi-dependency support."""

    VALID_YAML_PROPERTIES = {
        "bulk_modulus", "density", "dynamic_viscosity", "elastic_modulus",
        "electrical_conductivity", "electrical_resistivity", "energy_density",
        "fracture_toughness", "hardness", "heat_capacity", "heat_conductivity",
        "kinematic_viscosity", "latent_heat_of_fusion", "latent_heat_of_vaporization",
        "magnetic_permeability", "poisson_ratio", "shear_modulus", "specific_enthalpy",
        "surface_tension", "thermal_diffusivity", "thermal_expansion_coefficient",
        "ultimate_tensile_strength", "viscosity", "yield_strength",
    }

    def __init__(self, yaml_path: Union[str, Path]) -> None:
        super().__init__(yaml_path)
        logger.info("Initializing MaterialYAMLParser for: %s", yaml_path)

        # Initialize dependency resolver
        self.dependency_resolver = DependencyResolver(self.config)

        # Validate configuration
        self._validate_config()

        # Analyze and categorize properties
        self.categorized_properties = self._analyze_and_categorize_properties(self.config[PROPERTIES_KEY])

        # Initialize processors
        self.property_processor = PropertyProcessor()
        self.visualizer = PropertyVisualizer(self)

        logger.info("MaterialYAMLParser initialized successfully with %d property categories",
                    len(self.categorized_properties))

    def create_material(self, symbol_mapping: Dict[str, sp.Symbol], enable_plotting: bool = True) -> Material:
        """
        Create a Material instance from the parsed configuration.

        Args:
            symbol_mapping: Dictionary mapping dependency names to SymPy symbols
                          e.g., {'temperature': sp.Symbol('T'), 'pressure': sp.Symbol('P')}
            enable_plotting: Whether to generate visualization plots
        """
        logger.info("Creating material from configuration: %s", self.config_path)

        try:
            # Validate symbol mapping
            self._validate_symbol_mapping(symbol_mapping)

            # Extract basic material info
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

            # Initialize visualization if enabled
            visualizer = None
            if enable_plotting:
                self.visualizer.initialize_plots()
                self.visualizer.reset_visualization_tracking()
                visualizer = self.visualizer
                logger.info("Visualization enabled")
            else:
                logger.debug("Visualization disabled")

            # Process properties with multi-dependency support
            logger.info("Starting property processing for material: %s", name)
            self.property_processor.process_properties(
                material=material,
                symbol_mapping=symbol_mapping,
                dependency_resolver=self.dependency_resolver,
                properties=self.config[PROPERTIES_KEY],
                categorized_properties=self.categorized_properties,
                base_dir=self.base_dir,
                visualizer=visualizer
            )

            # Save plots if visualization was enabled
            if enable_plotting and visualizer is not None:
                logger.info("Saving property plots for material: %s", name)
                self.visualizer.save_property_plots()

            logger.info("Successfully created material: %s", name)
            return material

        except KeyError as e:
            logger.error("Configuration error for material creation - missing key: %s", e)
            raise ValueError(f"Configuration error: Missing {str(e)}") from e
        except Exception as e:
            logger.error("Failed to create material from %s: %s", self.config_path, e)
            raise ValueError(f"Failed to create material -> {str(e)}") from e

    def _validate_symbol_mapping(self, symbol_mapping: Dict[str, sp.Symbol]) -> None:
        """Validate that symbol mapping matches the configuration."""
        independent_vars = self.dependency_resolver.independent_vars

        # Check that all defined independent variables have corresponding symbols
        for dep_name in independent_vars.keys():
            if dep_name not in symbol_mapping:
                raise ValueError(f"Missing symbol for dependency '{dep_name}'. "
                                 f"Required symbols: {list(independent_vars.keys())}")

        # Check that all provided symbols are actually used
        unused_symbols = set(symbol_mapping.keys()) - set(independent_vars.keys())
        if unused_symbols:
            logger.warning("Unused symbols provided: %s", unused_symbols)

        # Validate that all values are SymPy symbols
        for dep_name, symbol in symbol_mapping.items():
            if not isinstance(symbol, sp.Symbol):
                raise TypeError(f"Symbol for dependency '{dep_name}' must be a SymPy Symbol, got {type(symbol)}")

    def _validate_config(self) -> None:
        """Validate the configuration structure and content."""
        logger.debug("Starting configuration validation")

        if not isinstance(self.config, dict):
            raise ValueError("The YAML file must start with a dictionary/object structure")

        self._validate_required_fields()

        properties = self.config.get(PROPERTIES_KEY, {})
        if not isinstance(properties, dict):
            raise ValueError("The 'properties' section must be a dictionary")

        self._validate_property_names(properties)

        # Validate multi-dependency structure
        for prop_name, prop_config in properties.items():
            self.dependency_resolver.validate_dependency_structure(prop_name, prop_config)

        logger.info("Configuration validation completed successfully")

    def _validate_required_fields(self) -> None:
        """Validate that all required fields are present."""
        logger.debug("Validating required fields")

        if MATERIAL_TYPE_KEY not in self.config:
            raise ValueError("Missing required field: material_type")

        material_type = self.config[MATERIAL_TYPE_KEY]
        if material_type not in [ALLOY_KEY, PURE_METAL_KEY]:
            raise ValueError(f"Invalid material_type: {material_type}. Must be {PURE_METAL_KEY} or {ALLOY_KEY}")

        # Common required fields
        common_fields = {NAME_KEY, MATERIAL_TYPE_KEY, COMPOSITION_KEY, PROPERTIES_KEY}

        # Type-specific required fields
        if material_type == PURE_METAL_KEY:
            type_specific_fields = {MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY}
        else:  # alloy
            type_specific_fields = {
                SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
                INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY
            }

        required_fields = common_fields | type_specific_fields
        missing_fields = required_fields - set(self.config.keys())

        if missing_fields:
            raise ValueError(f"Missing required fields: {sorted(list(missing_fields))}")

    def _validate_property_names(self, properties: Dict[str, Any]) -> None:
        """Validate property names against known properties."""
        for prop_name in properties.keys():
            if prop_name not in self.VALID_YAML_PROPERTIES:
                suggestions = get_close_matches(prop_name, self.VALID_YAML_PROPERTIES, n=3, cutoff=0.6)
                error_msg = f"Unknown property: '{prop_name}'"
                if suggestions:
                    error_msg += f". Did you mean: {', '.join(suggestions)}?"
                raise ValueError(error_msg)

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
