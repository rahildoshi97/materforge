from typing import List


class PropertyError(Exception):
    """Base exception for property-related errors."""
    pass


class DependencyError(PropertyError):
    """Exception for dependency-related errors."""

    def __init__(self, expression: str, missing_deps: List[str], available_props: List[str] = None):
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
        message = f"Circular dependency detected: {cycle_str}\nPlease resolve this cycle in your configuration."
        super().__init__(message)
