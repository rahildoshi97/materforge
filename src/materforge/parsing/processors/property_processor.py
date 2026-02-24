# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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
    ComputedPropertyHandler,
)
from materforge.parsing.processors.post_processor import PropertyPostProcessor

logger = logging.getLogger(__name__)


class PropertyProcessor(PropertyProcessorBase):
    """Orchestrates property processing by delegating to specialised handlers.
    Properties are processed in PropertyType enum definition order, which
    guarantees CONSTANT_VALUE is resolved before any type that may reference
    those scalars (STEP_FUNCTION, TABULAR_DATA, etc.).
    """

    def __init__(self) -> None:
        super().__init__()
        self.handlers = {
            PropertyType.CONSTANT_VALUE: ConstantValuePropertyHandler(),
            PropertyType.STEP_FUNCTION: StepFunctionPropertyHandler(),
            PropertyType.FILE_IMPORT: FileImportPropertyHandler(),
            PropertyType.TABULAR_DATA: TabularDataPropertyHandler(),
            PropertyType.PIECEWISE_EQUATION: PiecewiseEquationPropertyHandler(),
            PropertyType.COMPUTED_PROPERTY: ComputedPropertyHandler(),
        }
        self.post_processor = PropertyPostProcessor()
        self.properties: Optional[Dict[str, Any]] = None
        self.categorized_properties: Optional[Dict[PropertyType, List[Tuple[str, Any]]]] = None
        self.base_dir: Optional[Path] = None
        logger.debug("PropertyProcessor initialised with %d handlers", len(self.handlers))

    # --- Public API ---
    def process_properties(self, material: Material, dependency: sp.Symbol,
                           properties: Dict[str, Any],
                           categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                           base_dir: Path, visualizer) -> None:
        """Processes all properties for the material.
        Args:
            material: Target Material instance.
            dependency: SymPy symbol used as the independent variable.
            properties: Raw properties dict from the YAML config.
            categorized_properties: Properties pre-grouped by PropertyType.
            base_dir: Base directory for resolving relative file paths.
            visualizer: PropertyVisualizer instance, or None.
        """
        logger.info("Starting property processing for '%s'", material.name)
        self._initialize_processing_context(material, dependency, properties, categorized_properties, base_dir, visualizer)
        try:
            self._process_by_category(material, dependency)
            logger.info("Starting post-processing for '%s'", material.name)
            self.post_processor.post_process_properties(
                material, dependency, self.properties, self.categorized_properties, self.processed_properties)
            logger.info("Finished processing all properties for '%s'", material.name)
        except Exception as e:
            logger.error("Property processing failed for '%s': %s", material.name, e, exc_info=True)
            raise ValueError(f"Failed to process properties\n -> {str(e)}") from e

    # --- Private helpers ---
    def _initialize_processing_context(self, material: Material, dependency: sp.Symbol,
                                       properties: Dict[str, Any],
                                       categorized_properties: Dict[PropertyType, List[Tuple[str, Any]]],
                                       base_dir: Path, visualizer) -> None:
        """Sets up shared state and injects context into all handlers."""
        logger.debug("Initialising processing context for '%s'", material.name)
        self.properties = properties
        self.categorized_properties = categorized_properties
        self.base_dir = base_dir
        self.visualizer = visualizer
        self.processed_properties = set()
        for handler_type, handler in self.handlers.items():
            handler.set_processing_context(self.base_dir, visualizer, self.processed_properties)
            logger.debug("Context set for handler: %s", handler_type.name)
        computed_handler = self.handlers.get(PropertyType.COMPUTED_PROPERTY)
        if computed_handler:
            computed_handler.set_computed_property_processor(properties)

    def _process_by_category(self, material: Material, dependency: sp.Symbol) -> None:
        """Iterates PropertyType enum order and processes each category.
        CONSTANT_VALUE is the first enum value, so all scalar constants are
        assigned to the material before any other type runs - enabling safe
        forward references from STEP_FUNCTION, TABULAR_DATA, etc.
        """
        total = sum(len(v) for v in self.categorized_properties.values())
        active = sum(1 for v in self.categorized_properties.values() if v)
        logger.info("Processing %d properties across %d categories", total, active)
        for prop_type, prop_list in self.categorized_properties.items():
            if not prop_list:
                continue
            handler = self.handlers.get(prop_type)
            if handler is None:
                raise ValueError(f"No handler for property type: {prop_type.name}")
            logger.info("Processing %d %s properties", len(prop_list), prop_type.name)
            for prop_name, config in prop_list:
                logger.debug("Processing '%s'", prop_name)
                try:
                    handler.process_property(material, prop_name, config, dependency)
                except Exception as e:
                    logger.error("Failed to process '%s': %s", prop_name, e, exc_info=True)
                    raise
            logger.info("Completed %s", prop_type.name)
