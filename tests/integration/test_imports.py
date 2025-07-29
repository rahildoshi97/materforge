"""Test imports work correctly."""

from pathlib import Path
import pytest
import sympy as sp
from pymatlib import create_material

def test_all_imports():
    """Test that all modules can be imported without circular dependencies."""
    try:
        # Test main package imports
        import pymatlib
        import pymatlib.core
        import pymatlib.algorithms
        import pymatlib.parsing
        import pymatlib.visualization
        # Test core imports
        from pymatlib.core.materials import Material
        from pymatlib.core.elements import ChemicalElement
        # Test parsing imports
        from pymatlib.parsing.validation.property_validator import validate_monotonic_energy_density
        # Test algorithm imports
        from pymatlib.algorithms.interpolation import interpolate_value
        from pymatlib.algorithms.piecewise_builder import PiecewiseBuilder
        # Test visualization imports
        from pymatlib.visualization.plotters import PropertyVisualizer
    except ImportError as e:
        pytest.fail(f"Import failed: {e}")

def test_basic_material_creation():
    """Test basic material creation functionality."""
    T = sp.Symbol('T')  # Use standard temperature symbol
    # Construct paths more reliably
    current_file = Path(__file__)
    project_root = current_file.parent.parent.parent
    yaml_path_Al = project_root / "src" / "pymatlib" / "data" / "materials" / "pure_metals" / "Al" / "Al.yaml"
    yaml_path_SS304L = project_root / "src" / "pymatlib" / "data" / "materials" / "alloys" / "1.4301" / "1.4301.yaml"
    # Test aluminum material creation if file exists
    if yaml_path_Al.exists():
        try:
            mat_Al = create_material(yaml_path=yaml_path_Al, dependency=T, enable_plotting=False)
            assert mat_Al is not None
            assert mat_Al.name == "Aluminum"
            assert mat_Al.material_type == "pure_metal"
        except Exception as e:
            pytest.fail(f"Failed to create aluminum material: {e}")
    else:
        pytest.fail(f"Aluminum YAML file not found: {yaml_path_Al}")
    # Test steel material creation if file exists
    if yaml_path_SS304L.exists():
        try:
            mat_SS304L = create_material(yaml_path=yaml_path_SS304L, dependency=T, enable_plotting=False)
            assert mat_SS304L is not None
            assert "Steel" in mat_SS304L.name or "1.4301" in mat_SS304L.name
            assert mat_SS304L.material_type == "alloy"
        except Exception as e:
            pytest.fail(f"Failed to create steel material: {e}")
    else:
        pytest.fail(f"Steel YAML file not found: {yaml_path_SS304L}")

def test_circular_dependencies():
    """Test specifically for circular import dependencies."""
    import sys
    # Clear any previously imported pymatlib modules
    modules_to_remove = [name for name in sys.modules.keys() if name.startswith('pymatlib')]
    for module_name in modules_to_remove:
        if module_name in sys.modules:
            del sys.modules[module_name]
    try:
        # Import in different orders to catch circular dependencies
        import pymatlib.core.materials
        import pymatlib.parsing.api
        import pymatlib.algorithms.interpolation
        import pymatlib.visualization.plotters
        # Import main package last
        import pymatlib
    except ImportError as e:
        pytest.fail(f"Circular dependency detected: {e}")

if __name__ == "__main__":
    test_all_imports()
    test_basic_material_creation()
    test_circular_dependencies()
    print("✅ All import tests passed - no circular dependencies detected")
