# test_validation_functions.py
import numpy as np
from materforge.parsing.validation.property_validator import validate_monotonic_energy_density
from materforge.parsing.validation.array_validator import is_monotonic

def test_moved_functions():
    """Test that moved validation functions work correctly."""
    temp_array = np.array([300, 400, 500])
    prop_array = np.array([1000, 1100, 1200])
    # Test monotonicity
    assert is_monotonic(prop_array, mode="strictly_increasing")
    # Test energy density validation
    validate_monotonic_energy_density("energy_density", temp_array, prop_array)
    print("✅ All validation functions working correctly")

if __name__ == "__main__":
    test_moved_functions()
