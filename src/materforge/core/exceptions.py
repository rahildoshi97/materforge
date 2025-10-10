# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Custom exceptions for materforge core functionality."""
import logging

logger = logging.getLogger(__name__)


class MaterialError(Exception):
    """Base exception for all material-related errors."""

    def __init__(self, message):
        super().__init__(message)
        logger.error("MaterialError raised: %s", message)


class MaterialCompositionError(MaterialError):
    """Exception raised when material composition validation fails."""

    def __init__(self, message):
        super().__init__(message)
        logger.error("MaterialCompositionError raised: %s", message)


class MaterialTemperatureError(MaterialError):
    """Exception raised when material temperature validation fails."""

    def __init__(self, message):
        super().__init__(message)
        logger.error("MaterialTemperatureError raised: %s", message)
