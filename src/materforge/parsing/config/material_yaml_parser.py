# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
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
    LIQUIDUS_TEMPERATURE_KEY, INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY, ALLOY_KEY, NAME_KEY

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
    def create_material(self, dependency: sp.Symbol, enable_plotting: bool = True) -> Material:
        """Create a Material instance from the parsed configuration."""
        logger.info("Creating material from configuration: %s", self.config_path)
        try:
            name = self.config.get(NAME_KEY, "Unnamed Material")
            material_type = self.config[MATERIAL_TYPE_KEY]
            logger.info("Creating material: %s (type: %s)", name, material_type)
            elements = self._get_elements()
            composition = list(self.config[COMPOSITION_KEY].values())
            # Build temperature kwargs dynamically — any top-level key that is not
            # one of the four structural keys is treated as a temperature/scalar param.
            # This allows future material types to define their own temperature fields
            # without changing this method.
            reserved_keys = {NAME_KEY, MATERIAL_TYPE_KEY, COMPOSITION_KEY, PROPERTIES_KEY}
            temperature_params = {
                k: sp.Float(v)
                for k, v in self.config.items()
                if k not in reserved_keys
            }
            material = Material(
                name=name,
                elements=elements,
                composition=composition,
                material_type=material_type,
                **temperature_params,
            )
            logger.debug("Created material '%s' with params: %s", name, list(temperature_params.keys()))
            # Visualizer setup
            visualizer = None
            should_visualize = enable_plotting and isinstance(dependency, sp.Symbol)
            if should_visualize:
                self.visualizer.initialize_plots()
                self.visualizer.reset_visualization_tracking()
                visualizer = self.visualizer
                logger.info("Visualization enabled for symbolic dependency")
            else:
                logger.debug("Visualization disabled - %s",
                             "numeric dependency" if not isinstance(dependency, sp.Symbol)
                             else "plotting not enabled")
            # Process properties
            logger.info("Starting property processing for material: %s", name)
            self.property_processor.process_properties(
                material=material,
                dependency=dependency,
                properties=self.config[PROPERTIES_KEY],
                categorized_properties=self.categorized_properties,
                base_dir=self.base_dir,
                visualizer=visualizer,
            )
            if should_visualize and visualizer is not None:
                logger.info("Saving property plots for material: %s", name)
                self.visualizer.save_property_plots()
            logger.info(f"Successfully created material: {name}")
            return material
        except KeyError as e:
            logger.error("Missing configuration key: %s", e, exc_info=True)
            raise ValueError(f"Configuration error: Missing {str(e)}") from e
        except Exception as e:
            logger.error("Failed to create material from %s: %s", self.config_path, e, exc_info=True)
            raise ValueError(f"Failed to create material \n -> {str(e)}") from e

    # --- Validation ---
    def _validate_config(self) -> None:
        """Validate top-level structure only. Property names are not restricted."""
        logger.debug("Starting configuration validation")
        if not isinstance(self.config, dict):
            raise ValueError("The YAML file must start with a dictionary/object structure, not a list or scalar value")
        self._validate_required_fields()
        properties = self.config.get(PROPERTIES_KEY, {})
        if not isinstance(properties, dict):
            raise ValueError("The 'properties' section in your YAML file must be a dictionary")
        logger.info("Configuration validation completed successfully")

    def _validate_required_fields(self) -> None:
        """Validate structural required fields. Unknown extra keys are permitted."""
        logger.debug("Validating required fields")
        if MATERIAL_TYPE_KEY not in self.config:
            raise ValueError("Missing required field: material_type")
        material_type = self.config[MATERIAL_TYPE_KEY]
        logger.debug("Material type: %s", material_type)
        if material_type not in [ALLOY_KEY, PURE_METAL_KEY]:
            raise ValueError(f"Invalid material_type: '{material_type}'. Must be '{PURE_METAL_KEY}' or '{ALLOY_KEY}'")
        common_fields = {NAME_KEY, MATERIAL_TYPE_KEY, COMPOSITION_KEY, PROPERTIES_KEY}
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
            raise ValueError(f"Missing required fields for {material_type}: {', '.join(sorted(missing_fields))}")
        # NOTE: extra top-level keys are intentionally allowed.
        # Users may add arbitrary scalar fields (e.g. reference_pressure, custom_param).
        self._validate_name()
        self._validate_composition()
        logger.debug("Required fields validation completed")

    def _validate_name(self) -> None:
        logger.debug("Validating material name")
        name = self.config[NAME_KEY]
        if name is None:
            raise ValueError("Material name cannot be None")
        if not isinstance(name, str):
            raise ValueError(f"Material name must be a string, got {type(name).__name__} (value: {name})")
        if len(name) > 100:
            logger.warning("Material name '%s' exceeds 100 characters", name)

    def _validate_composition(self) -> None:
        logger.debug("Validating composition")
        composition = self.config.get(COMPOSITION_KEY, {})
        material_type = self.config[MATERIAL_TYPE_KEY]
        if not isinstance(composition, dict):
            raise ValueError("Composition must be a dictionary")
        if not composition:
            raise ValueError("Composition cannot be empty")
        for element, fraction in composition.items():
            if not isinstance(fraction, (int, float)):
                raise ValueError(
                    f"Composition fraction for '{element}' must be a number, got {type(fraction).__name__}")
            if fraction < 0:
                raise ValueError(f"Composition fraction for '{element}' cannot be negative, got {fraction}")
            if fraction > 1:
                raise ValueError(f"Composition fraction for '{element}' cannot exceed 1.0, got {fraction}")
        total = sum(composition.values())
        if not abs(total - 1.0) < ProcessingConstants.COMPOSITION_THRESHOLD:
            raise ValueError(f"Composition fractions must sum to 1.0, got {total}")
        if material_type == PURE_METAL_KEY:
            self._validate_pure_metal_composition_rules(composition)
        else:  # alloy
            self._validate_alloy_composition_rules(composition)

    @staticmethod
    def _validate_pure_metal_composition_rules(composition: dict) -> None:
        logger.debug("Validating pure metal composition rules")
        non_zero = {e: f for e, f in composition.items() if f > ProcessingConstants.COMPOSITION_THRESHOLD}
        if len(non_zero) == 0:
            raise ValueError("Pure metals must have at least one element with non-zero composition")
        if len(non_zero) > 1:
            element_list = ", ".join(f"{e}: {f}" for e, f in non_zero.items())
            raise ValueError(
                f"Pure metals must contain exactly one element with composition 1.0. "
                f"Found: {element_list}. Use material_type: 'alloy' for multi-element materials.")
        single_element, single_fraction = list(non_zero.items())[0]
        if not abs(single_fraction - 1.0) < ProcessingConstants.COMPOSITION_THRESHOLD:
            raise ValueError(
                f"Pure metal element '{single_element}' must have composition 1.0, "
                f"got {single_fraction}.")
        zero_elements = [e for e, f in composition.items() if f == 0.0]
        if zero_elements:
            raise ValueError(f"Pure metal composition should not include zero-valued elements: {zero_elements}.")

    @staticmethod
    def _validate_alloy_composition_rules(composition: dict) -> None:
        logger.debug("Validating alloy composition rules")
        non_zero = {e: f for e, f in composition.items() if f > 1e-10}
        if len(non_zero) < 2:
            if len(non_zero) == 1:
                raise ValueError(
                    f"Alloys must have at least 2 elements. "
                    f"Found only '{list(non_zero.keys())[0]}'. Use material_type: 'pure_metal' for single elements.")
            else:
                raise ValueError("Alloys must have at least 2 elements with non-zero composition")
        zero_elements = [e for e, f in composition.items() if f == 0.0]
        if zero_elements:
            logger.warning("Alloy has zero-valued elements: %s. Consider removing.", zero_elements)

    # --- Processing helpers ---
    def _get_elements(self) -> List:
        from materforge.data.elements.element_data import element_map
        element_symbols = list(self.config[COMPOSITION_KEY].keys())
        logger.debug("Looking up elements: %s", element_symbols)
        try:
            elements = [element_map[sym] for sym in element_symbols]
            return elements
        except KeyError as e:
            raise ValueError(f"Invalid element symbol: {str(e)}") from e

    @staticmethod
    def _analyze_and_categorize_properties(properties: Dict[str, Any]) -> Dict[PropertyType, List[Tuple[str, Any]]]:
        logger.debug("Analyzing and categorizing %d properties", len(properties))
        categorized: Dict[PropertyType, List[Tuple[str, Any]]] = {
            pt: [] for pt in PropertyType
        }
        for prop_name, config in properties.items():
            try:
                prop_type = PropertyTypeDetector.determine_property_type(prop_name, config)
                PropertyTypeDetector.validate_property_config(prop_name, config, prop_type)
                categorized[prop_type].append((prop_name, config))
                logger.debug("Property '%s' -> %s", prop_name, prop_type.name)
            except ValueError as e:
                raise ValueError(f"Configuration error for property '{prop_name}': {str(e)}") from e
        for prop_type, prop_list in categorized.items():
            if prop_list:
                logger.info("Found %d %s properties: %s",
                            len(prop_list), prop_type.name, [p[0] for p in prop_list])
        return categorized
