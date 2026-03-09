# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

from typing import List, Optional


class MaterialConfigError(ValueError):
    """Raised when the top-level YAML structure is invalid.

    Subclasses ValueError so existing callers catching ValueError are unaffected.
    """
    pass


class PropertyConfigError(MaterialConfigError):
    """Raised when an individual property block is structurally invalid.

    Narrower than MaterialConfigError - callers that need to distinguish a
    bad property config from a bad top-level structure can catch this directly.
    """
    pass


class PropertyError(Exception):
    """Base exception for property runtime errors."""
    pass


class DependencyError(PropertyError):
    """Exception for dependency-related errors."""

    def __init__(self, expression: str, missing_deps: List[str], available_props: Optional[List[str]] = None):
        self.expression = expression
        self.missing_deps = missing_deps
        self.available_props = available_props or []
        message = f"Missing dependencies in expression '{expression}': {', '.join(missing_deps)}"
        if available_props:
            message += f"\nAvailable properties: {', '.join(available_props)}"
            message += "\nPlease check for typos or add the missing properties to your configuration."
        super().__init__(message)


class CircularDependencyError(PropertyError):
    """Exception for circular dependency errors."""

    def __init__(self, dependency_path: List[str]):
        self.dependency_path = dependency_path
        cycle_str = " -> ".join(dependency_path)
        message = (f"Circular dependency detected: {cycle_str}\n"
            "Please resolve this cycle in your configuration.")
        super().__init__(message)
