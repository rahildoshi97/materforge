# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import sympy as sp
from ruamel.yaml import YAML, constructor, scanner
from materforge.core.materials import Material
from materforge.parsing.processors.property_processor import PropertyProcessor
from materforge.parsing.validation.errors import MaterialConfigError, PropertyConfigError
from materforge.parsing.validation.property_type_detector import PropertyType, PropertyTypeDetector
from materforge.visualization.plotters import PropertyVisualizer
from materforge.parsing.config.yaml_keys import NAME_KEY, PROPERTIES_KEY

logger = logging.getLogger(__name__)

_REQUIRED_TOP_LEVEL_KEYS = frozenset({NAME_KEY, PROPERTIES_KEY})


class BaseFileParser:
    """Base class for parsing configuration files."""

    def __init__(self, config_path: Union[str, Path]) -> None:
        self.config_path = Path(config_path)
        self.base_dir = self.config_path.parent
        self.config = self._load_config()
        logger.info("Successfully loaded configuration from: %s", self.config_path)

    def _load_config(self) -> Dict[str, Any]:
        raise NotImplementedError("Subclasses must implement _load_config")


class YAMLFileParser(BaseFileParser):
    """Parser for YAML configuration files."""

    def _load_config(self) -> Dict[str, Any]:
        yaml = YAML(typ='safe')
        yaml.allow_duplicate_keys = False
        try:
            logger.debug("Loading YAML file: %s", self.config_path)
            with open(self.config_path, 'r', encoding='utf-8') as f:
                config = yaml.load(f)
            logger.debug("YAML loaded - %d top-level keys", len(config) if config else 0)
            return config
        except FileNotFoundError as e:
            logger.error("YAML file not found: %s", self.config_path)
            raise FileNotFoundError(f"YAML file not found: {self.config_path}") from e
        except constructor.DuplicateKeyError as e:
            logger.error("Duplicate key in %s: %s", self.config_path, e)
            raise constructor.DuplicateKeyError(f"Duplicate key in {self.config_path}: {str(e)}") from e
        except scanner.ScannerError as e:
            logger.error("YAML syntax error in %s: %s", self.config_path, e)
            raise scanner.ScannerError(f"YAML syntax error in {self.config_path}: {str(e)}") from e
        except Exception as e:
            logger.error("Unexpected error parsing %s: %s", self.config_path, e, exc_info=True)
            raise ValueError(f"Error parsing {self.config_path}: {str(e)}") from e


class MaterialYAMLParser(YAMLFileParser):
    """Parser for material YAML configuration files.

    Validates and processes a schema-agnostic YAML file containing a name
    and a properties block. All thermophysical properties live in the properties
    block and are assigned dynamically to the Material instance.
    """

    def __init__(self, yaml_path: Union[str, Path]) -> None:
        super().__init__(yaml_path)
        logger.info("Initializing MaterialYAMLParser for: %s", yaml_path)
        self._validate_config()
        self.categorized_properties = self._analyze_and_categorize_properties(self.config[PROPERTIES_KEY])
        self.property_processor = PropertyProcessor()
        self.visualizer = PropertyVisualizer(self)
        logger.info("MaterialYAMLParser ready - %d property categories",
                    len(self.categorized_properties))

    # --- Public API ---
    def create_material(self, dependency: sp.Symbol, enable_plotting: bool = True) -> Material:
        """Creates a Material instance from the parsed configuration.

        Args:
            dependency:      SymPy symbol used as the independent variable in
                             property expressions (e.g. sp.Symbol('T')).
            enable_plotting: Whether to generate and save property plots.
        Returns:
            Fully initialised Material with all properties assigned.
        Raises:
            MaterialConfigError: If configuration is invalid or property processing fails.
        """
        name = self.config.get(NAME_KEY, "Unnamed Material")
        logger.info("Creating material: '%s'", name)
        try:
            material = Material(name=name)
            should_visualize = enable_plotting and isinstance(dependency, sp.Symbol)
            visualizer = None
            if should_visualize:
                self.visualizer.initialize_plots()
                self.visualizer.reset_visualization_tracking()
                visualizer = self.visualizer
                logger.info("Visualization enabled")
            else:
                logger.debug("Visualization disabled - %s",
                    "numeric dependency" if not isinstance(dependency, sp.Symbol)
                    else "plotting not requested")
            logger.info("Processing properties for '%s'", name)
            self.property_processor.process_properties(
                material=material,
                dependency=dependency,
                properties=self.config[PROPERTIES_KEY],
                categorized_properties=self.categorized_properties,
                base_dir=self.base_dir,
                visualizer=visualizer,
            )
            if should_visualize and visualizer is not None:
                logger.info("Saving property plots for '%s'", name)
                self.visualizer.save_property_plots()
            logger.info("Successfully created material: '%s'", name)
            return material
        except Exception as e:
            # Intermediate layer - do not log here, bubble up to api.py
            raise MaterialConfigError(f"Failed to create material \n -> {str(e)}") from e

    # --- Validation ---
    def _validate_config(self) -> None:
        """Validates the top-level YAML structure.

        Raises:
            MaterialConfigError: If the config is not a dict, contains unknown
                top-level keys, or if name/properties fields are missing or invalid.
        """
        logger.debug("Validating configuration structure")
        if not isinstance(self.config, dict):
            raise MaterialConfigError(f"YAML file must contain a top-level mapping, got {type(self.config).__name__}")
        unknown_keys = set(self.config.keys()) - _REQUIRED_TOP_LEVEL_KEYS
        if unknown_keys:
            raise MaterialConfigError(f"Unknown top-level keys: {sorted(unknown_keys)}. "
                f"Only {sorted(_REQUIRED_TOP_LEVEL_KEYS)} are allowed.")
        if NAME_KEY not in self.config:
            raise MaterialConfigError(f"Missing required field: '{NAME_KEY}'")
        name = self.config[NAME_KEY]
        if name is None or not isinstance(name, str) or not name.strip():
            raise MaterialConfigError(f"'{NAME_KEY}' must be a non-empty string, "
                f"got {type(name).__name__!r}: {name!r}")
        if len(name) > 100:
            logger.warning("Material name '%s' exceeds 100 characters", name)
        if PROPERTIES_KEY not in self.config:
            raise MaterialConfigError(f"Missing required field: '{PROPERTIES_KEY}'")
        props = self.config[PROPERTIES_KEY]
        if not isinstance(props, dict):
            raise MaterialConfigError(f"'{PROPERTIES_KEY}' must be a mapping, got {type(props).__name__}")
        if not props:
            raise MaterialConfigError(f"'{PROPERTIES_KEY}' block cannot be empty")
        logger.info("Configuration validation passed")

    # --- Property categorisation ---
    @staticmethod
    def _analyze_and_categorize_properties(properties: Dict[str, Any]) -> Dict[PropertyType, List[Tuple[str, Any]]]:
        """Detects and categorises all properties by their type.

        Args:
            properties: Raw properties mapping from the YAML config.
        Returns:
            Dict keyed by PropertyType, each entry a list of
            (property_name, config) tuples.
        Raises:
            PropertyConfigError: If any property config is structurally invalid.
        """
        logger.debug("Categorising %d properties", len(properties))
        categorized: Dict[PropertyType, List[Tuple[str, Any]]] = {
            pt: [] for pt in PropertyType
        }
        for prop_name, config in properties.items():
            try:
                prop_type = PropertyTypeDetector.determine_property_type(prop_name, config)
                PropertyTypeDetector.validate_property_config(prop_name, config, prop_type)
                categorized[prop_type].append((prop_name, config))
                logger.debug("'%s' -> %s", prop_name, prop_type.name)
            except ValueError as e:
                raise PropertyConfigError(f"Configuration error for property '{prop_name}': {str(e)}") from e
        for prop_type, prop_list in categorized.items():
            if prop_list:
                logger.info("%d %s: %s", len(prop_list), prop_type.name, [p[0] for p in prop_list])
        return categorized
