import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.validation.property_type_detector import PropertyType
from materforge.parsing.processors.property_processor_base import PropertyProcessorBase
from materforge.parsing.processors.property_handlers import (
    ConstantValuePropertyHandler,
    StepFunctionPropertyHandler,
    FileImportPropertyHandler,
    TabularDataPropertyHandler,
    PiecewiseEquationPropertyHandler,
    ComputedPropertyHandler
)
from materforge.parsing.processors.post_processor import PropertyPostProcessor

logger = logging.getLogger(__name__)


class PropertyProcessor(PropertyProcessorBase):
    """
    Main orchestrator for processing different property types for material objects.

    This class coordinates the processing of various property types by delegating
    to specialized handlers for each property type.
    """

    def __init__(self) -> None:
        """Initialize processor with specialized handlers."""
        super().__init__()
        # Initialize property handlers
        self.handlers = {
            PropertyType.CONSTANT_VALUE: ConstantValuePropertyHandler(),
            PropertyType.STEP_FUNCTION: StepFunctionPropertyHandler(),
            PropertyType.FILE_IMPORT: FileImportPropertyHandler(),
            PropertyType.TABULAR_DATA: TabularDataPropertyHandler(),
            PropertyType.PIECEWISE_EQUATION: PiecewiseEquationPropertyHandler(),
            PropertyType.COMPUTED_PROPERTY: ComputedPropertyHandler()
        }
        # Initialize post-processor
        self.post_processor = PropertyPostProcessor()
        # Processing state
        self.properties: Optional[Dict[str, Any]] = None
        self.categorized_properties: Optional[Dict[PropertyType, List[Tuple[str, Any]]]] = None
        self.base_dir: Optional[Path] = None
        logger.debug("PropertyProcessor initialized with %d handlers", len(self.handlers))

    def process_properties(self, material: Material, dependency: Union[float, sp.Symbol],
                           properties: Dict[str, Any],
                           categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                           base_dir: Path, visualizer) -> None:
        """Process all properties for the material."""
        logger.info("Starting property processing for material: %s", material.name)
        # Set processing context
        self._initialize_processing_context(material, dependency, properties, categorized_properties, base_dir, visualizer)
        try:
            # Process properties by type
            self._process_by_category(material, dependency)
            # Post-process properties (regression, etc.)
            logger.info("Starting post-processing for material: %s", material.name)
            self.post_processor.post_process_properties(
                material, dependency, self.properties, self.categorized_properties, self.processed_properties
            )
            logger.info("Successfully processed all properties for material: %s", material.name)
        except Exception as e:
            logger.error("Property processing failed for material '%s': %s", material.name, e, exc_info=True)
            raise ValueError(f"Failed to process properties \n -> {str(e)}") from e

    def _initialize_processing_context(self, material: Material, dependency: Union[float, sp.Symbol],
                                       properties: Dict[str, Any],
                                       categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                       base_dir: Path, visualizer) -> None:
        """Initialize processing context and handler dependencies."""
        logger.debug("Initializing processing context for material: %s", material.name)
        self.properties = properties
        self.categorized_properties = categorized_properties
        self.base_dir = base_dir
        self.visualizer = visualizer
        self.processed_properties = set()
        # Set context for all handlers
        for handler_type, handler in self.handlers.items():
            handler.set_processing_context(self.base_dir, visualizer, self.processed_properties)
            logger.debug("Set processing context for handler: %s", handler_type.name)
        # Initialize dependency processor for computed properties
        computed_handler = self.handlers.get(PropertyType.COMPUTED_PROPERTY)
        if computed_handler:
            computed_handler.set_computed_property_processor(properties)
            logger.debug("Dependency processor initialized for computed properties")

    def _process_by_category(self, material: Material, dependency: Union[float, sp.Symbol]) -> None:
        """Process properties grouped by category."""
        total_properties = sum(len(prop_list) for prop_list in self.categorized_properties.values())
        active_categories = len([cat for cat, props in self.categorized_properties.items() if props])
        logger.info("Processing %d properties across %d categories", total_properties, active_categories)
        for prop_type, prop_list in self.categorized_properties.items():
            if not prop_list:
                continue
            logger.info("Processing %d properties of type: %s", len(prop_list), prop_type.name)
            handler = self.handlers.get(prop_type)
            if handler is None:
                logger.error("No handler available for property type: %s", prop_type.name)
                raise ValueError(f"No handler available for property type: {prop_type.name}")
            # Sort properties to prioritize temperature references
            sorted_props = self._sort_properties_by_priority(prop_list)
            for prop_name, config in sorted_props:
                logger.debug("Processing property: %s", prop_name)
                try:
                    handler.process_property(material, prop_name, config, dependency)
                    logger.debug("Successfully processed property: %s", prop_name)
                except Exception as e:
                    logger.error("Failed to process property '%s': %s", prop_name, e, exc_info=True)
                    raise
            logger.info("Completed processing %s properties", prop_type.name)

    @staticmethod
    def _sort_properties_by_priority(prop_list: List[Tuple[str, Any]]) -> List[Tuple[str, Any]]:
        """Sort properties to process temperature references first."""
        from materforge.parsing.config.yaml_keys import (
            MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
            SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
            INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY
        )
        priority_props = {
            MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
            SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
            INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY
        }
        sorted_props = sorted(prop_list, key=lambda x: 0 if x[0] in priority_props else 1)
        priority_count = len([p for p in prop_list if p[0] in priority_props])
        if priority_count > 0:
            logger.debug("Prioritized %d temperature reference properties", priority_count)
        return sorted_props
