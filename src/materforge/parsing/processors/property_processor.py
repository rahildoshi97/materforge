import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import sympy as sp

from materforge.core.materials import Material
from materforge.parsing.validation.property_type_detector import PropertyType
from materforge.parsing.processors.property_processor_base import PropertyProcessorBase
from materforge.parsing.processors.dependency_resolver import DependencyResolver
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
    Enhanced to support multi-dependency properties through symbol mapping.
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
        self.dependency_resolver: Optional[DependencyResolver] = None
        self.symbol_mapping: Optional[Dict[str, sp.Symbol]] = None
        self.base_dir: Optional[Path] = None

        logger.debug("PropertyProcessor initialized with %d handlers", len(self.handlers))

    def process_properties(self, material: Material,
                           symbol_mapping: Dict[str, sp.Symbol],
                           dependency_resolver: DependencyResolver,
                           properties: Dict[str, Any],
                           categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                           base_dir: Path,
                           visualizer) -> None:
        """
        Process all properties for the material with multi-dependency support.

        Args:
            material: Material instance to process properties for
            symbol_mapping: Dictionary mapping dependency names to SymPy symbols
            dependency_resolver: Resolver for handling multi-dependency configurations
            properties: Raw property configurations from YAML
            categorized_properties: Properties organized by type
            base_dir: Base directory for file operations
            visualizer: Visualizer instance for plotting
        """
        logger.info("Starting multi-dependency property processing for material: %s", material.name)
        logger.debug("Symbol mapping: %s", {k: str(v) for k, v in symbol_mapping.items()})

        # Set processing context
        self._initialize_processing_context(
            material, symbol_mapping, dependency_resolver, properties,
            categorized_properties, base_dir, visualizer
        )

        try:
            # Process properties by type with multi-dependency support
            self._process_by_category(material)

            # Post-process properties (regression, etc.)
            logger.info("Starting post-processing for material: %s", material.name)
            self.post_processor.post_process_multi_dependency_properties(
                material, self.symbol_mapping, self.dependency_resolver,
                self.properties, self.categorized_properties, self.processed_properties
            )

            logger.info("Successfully processed all properties for material: %s", material.name)

        except Exception as e:
            logger.error("Property processing failed for material '%s': %s", material.name, e, exc_info=True)
            raise ValueError(f"Failed to process properties -> {str(e)}") from e

    def _initialize_processing_context(self, material: Material,
                                       symbol_mapping: Dict[str, sp.Symbol],
                                       dependency_resolver: DependencyResolver,
                                       properties: Dict[str, Any],
                                       categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                       base_dir: Path,
                                       visualizer) -> None:
        """Initialize processing context with multi-dependency support."""
        logger.debug("Initializing multi-dependency processing context for material: %s", material.name)

        self.properties = properties
        self.categorized_properties = categorized_properties
        self.dependency_resolver = dependency_resolver
        self.symbol_mapping = symbol_mapping
        self.base_dir = base_dir
        self.visualizer = visualizer
        self.processed_properties = set()

        # Set context for all handlers with multi-dependency support
        for handler_type, handler in self.handlers.items():
            handler.set_multi_dependency_context(
                base_dir=self.base_dir,
                visualizer=visualizer,
                processed_properties=self.processed_properties,
                dependency_resolver=self.dependency_resolver,
                symbol_mapping=self.symbol_mapping
            )
            logger.debug("Set multi-dependency context for handler: %s", handler_type.name)

        # Initialize dependency processor for computed properties
        computed_handler = self.handlers.get(PropertyType.COMPUTED_PROPERTY)
        if computed_handler:
            computed_handler.set_dependency_processor(properties, self.dependency_resolver)
            logger.debug("Multi-dependency processor initialized for computed properties")

    def _process_by_category(self, material: Material) -> None:
        """Process properties grouped by category with multi-dependency support."""
        total_properties = sum(len(prop_list) for prop_list in self.categorized_properties.values())
        active_categories = len([cat for cat, props in self.categorized_properties.items() if props])

        logger.info("Processing %d properties across %d categories with multi-dependency support",
                    total_properties, active_categories)

        for prop_type, prop_list in self.categorized_properties.items():
            if not prop_list:
                continue

            logger.info("Processing %d properties of type: %s", len(prop_list), prop_type.name)

            handler = self.handlers.get(prop_type)
            if handler is None:
                logger.error("No handler available for property type: %s", prop_type.name)
                raise ValueError(f"No handler available for property type: {prop_type.name}")

            # Sort properties to prioritize dependencies
            sorted_props = self._sort_properties_by_dependency_priority(prop_list)

            for prop_name, config in sorted_props:
                logger.debug("Processing multi-dependency property: %s", prop_name)

                try:
                    # Process property with multi-dependency support
                    handler.process_multi_dependency_property(material, prop_name, config)
                    logger.debug("Successfully processed multi-dependency property: %s", prop_name)

                except Exception as e:
                    logger.error("Failed to process multi-dependency property '%s': %s",
                                 prop_name, e, exc_info=True)
                    raise

            logger.info("Completed processing %s properties", prop_type.name)

    @staticmethod
    def _sort_properties_by_dependency_priority(prop_list: List[Tuple[str, Any]]) -> List[Tuple[str, Any]]:
        """Sort properties by dependency complexity and reference priority."""
        from materforge.parsing.config.yaml_keys import (
            MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
            SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
            INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY,
            DEPENDENCIES_KEY
        )

        # High priority: temperature references
        temperature_refs = {
            MELTING_TEMPERATURE_KEY, BOILING_TEMPERATURE_KEY,
            SOLIDUS_TEMPERATURE_KEY, LIQUIDUS_TEMPERATURE_KEY,
            INITIAL_BOILING_TEMPERATURE_KEY, FINAL_BOILING_TEMPERATURE_KEY
        }

        def get_priority(prop_tuple):
            prop_name, config = prop_tuple

            # Highest priority: temperature references
            if prop_name in temperature_refs:
                return 0

            # High priority: constants (no dependencies)
            if not isinstance(config, dict):
                return 1

            # Medium priority: single dependency
            if isinstance(config, dict) and DEPENDENCIES_KEY in config:
                dependencies = config[DEPENDENCIES_KEY]
                if isinstance(dependencies, list):
                    return 2 + len(dependencies)  # More dependencies = lower priority

            # Low priority: legacy format or complex dependencies
            return 10

        sorted_props = sorted(prop_list, key=get_priority)

        # Log priority information
        priority_counts = {}
        for prop_name, config in sorted_props:
            priority = get_priority((prop_name, config))
            priority_counts[priority] = priority_counts.get(priority, 0) + 1

        logger.debug("Property processing priority distribution: %s", priority_counts)

        return sorted_props

    def get_processing_statistics(self) -> Dict[str, Any]:
        """Get statistics about the processing session."""
        if not self.categorized_properties:
            return {}

        total_properties = sum(len(prop_list) for prop_list in self.categorized_properties.values())
        processed_count = len(self.processed_properties)

        stats = {
            'total_properties': total_properties,
            'processed_properties': processed_count,
            'processing_rate': processed_count / total_properties if total_properties > 0 else 0,
            'property_types': {
                prop_type.name: len(props)
                for prop_type, props in self.categorized_properties.items()
                if props
            },
            'dependencies_used': list(self.symbol_mapping.keys()) if self.symbol_mapping else [],
            'multi_dependency_enabled': self.dependency_resolver.is_multi_dependency if self.dependency_resolver else False
        }

        return stats
