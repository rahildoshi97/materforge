"""Unit tests for Material class."""
import pytest
import numpy as np
import sympy as sp

from materforge.core.materials import Material
from materforge.core.exceptions import MaterialCompositionError, MaterialTemperatureError
from materforge.core.elements import ChemicalElement

class TestMaterial:
    """Test cases for Material class."""
    def test_pure_metal_creation_valid(self, sample_aluminum_element):
        """Test valid pure metal creation."""
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        assert material.name == "Test Aluminum"
        assert material.material_type == "pure_metal"
        assert len(material.elements) == 1
        assert material.composition[0] == 1.0
        assert float(material.melting_temperature) == 933.47

    def test_alloy_creation_valid(self, sample_steel_elements):
        """Test valid alloy creation."""
        material = Material(
            name="Test Steel",
            material_type="alloy",
            elements=sample_steel_elements,
            composition=[0.68, 0.02, 0.20, 0.10],
            solidus_temperature=sp.Float(1400.0),
            liquidus_temperature=sp.Float(1450.0),
            initial_boiling_temperature=sp.Float(2800.0),
            final_boiling_temperature=sp.Float(2900.0)
        )
        assert material.name == "Test Steel"
        assert material.material_type == "alloy"
        assert len(material.elements) == 4
        assert np.isclose(sum(material.composition), 1.0)

    def test_composition_validation_sum_not_one(self, sample_aluminum_element):
        """Test composition validation when sum is not 1.0."""
        with pytest.raises(MaterialCompositionError, match="sum of the composition array must be 1.0"):
            Material(
                name="Invalid Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[0.8],  # Sum is not 1.0
                melting_temperature=sp.Float(933.47),
                boiling_temperature=sp.Float(2792.0)
            )

    def test_pure_metal_multiple_elements(self, sample_steel_elements):
        """Test pure metal with multiple elements raises error."""
        with pytest.raises(MaterialCompositionError, match="Pure metals must have exactly 1 element"):
            Material(
                name="Invalid Pure Metal",
                material_type="pure_metal",
                elements=sample_steel_elements,
                composition=[0.25, 0.25, 0.25, 0.25],
                melting_temperature=sp.Float(933.47),
                boiling_temperature=sp.Float(2792.0)
            )

    def test_alloy_single_element(self, sample_aluminum_element):
        """Test alloy with single element raises error."""
        with pytest.raises(MaterialCompositionError, match="Alloys must have at least 2 elements"):
            Material(
                name="Invalid Alloy",
                material_type="alloy",
                elements=[sample_aluminum_element],
                composition=[1.0],
                solidus_temperature=sp.Float(1400.0),
                liquidus_temperature=sp.Float(1450.0),
                initial_boiling_temperature=sp.Float(2800.0),
                final_boiling_temperature=sp.Float(2900.0)
            )

    def test_temperature_validation_valid(self, sample_aluminum_element):
        """Test valid temperature combinations."""
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        assert material is not None
        assert float(material.melting_temperature) == 933.47
        assert float(material.boiling_temperature) == 2792.0

    def test_temperature_validation_melting_greater_than_boiling(self, sample_aluminum_element):
        """Test that melting temperature greater than boiling temperature raises error."""
        exception_caught = False
        try:
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(2792.0),
                boiling_temperature=sp.Float(933.47)
            )
        except Exception as e:
            # Check exception by name instead of isinstance
            if e.__class__.__name__ == "MaterialTemperatureError":
                exception_caught = True
                assert "melting_temperature" in str(e)
                assert "must be less than" in str(e)
            else:
                pytest.fail(f"Wrong exception type raised: {type(e).__name__}: {e}")
        if not exception_caught:
            pytest.fail("Expected MaterialTemperatureError was not raised")

    def test_temperature_validation_melting_greater_than_boiling_1(self, sample_aluminum_element):
        """Test that melting temperature greater than boiling temperature raises error."""
        try:
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(2792.0),
                boiling_temperature=sp.Float(933.47)
            )
            pytest.fail("Expected MaterialTemperatureError was not raised")
        except Exception as e:
            assert "MaterialTemperatureError" in str(type(e))
            assert "melting_temperature" in str(e)
            assert "must be less than" in str(e)

    def test_temperature_validation_melting_greater_than_boiling_2(self, sample_aluminum_element):
        """Test that melting temperature greater than boiling temperature raises error."""
        with pytest.raises(MaterialTemperatureError, match=r"melting_temperature.*must be less than.*boiling_temperature"):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(2792.0),
                boiling_temperature=sp.Float(933.47)
            )

    def test_temperature_validation_negative_temperature(self, sample_aluminum_element):
        """Test that negative temperatures raise error."""
        exception_caught = False
        try:
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(-10.0),
                boiling_temperature=sp.Float(2792.0)
            )
        except Exception as e:
            if e.__class__.__name__ == "MaterialTemperatureError":
                exception_caught = True
                assert "must be above absolute zero" in str(e)
                assert "-10.0K" in str(e)
            else:
                pytest.fail(f"Wrong exception type raised: {type(e).__name__}: {e}")
        if not exception_caught:
            pytest.fail("Expected MaterialTemperatureError was not raised")

    def test_temperature_validation_negative_temperature_1(self, sample_aluminum_element):
        """Test that negative temperatures raise error."""
        with pytest.raises(MaterialTemperatureError, match=r"must be above absolute zero, got -10\.0K"):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(-10.0),
                boiling_temperature=sp.Float(2792.0)
            )

    def test_temperature_validation_zero_kelvin(self, sample_aluminum_element):
        """Test that zero Kelvin is invalid."""
        exception_caught = False
        try:
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(0.0),
                boiling_temperature=sp.Float(2792.0)
            )
        except Exception as e:
            if e.__class__.__name__ == "MaterialTemperatureError":
                exception_caught = True
                assert "must be above absolute zero" in str(e)
                assert "0.0K" in str(e)
            else:
                pytest.fail(f"Wrong exception type raised: {type(e).__name__}: {e}")
        if not exception_caught:
            pytest.fail("Expected MaterialTemperatureError was not raised")

    def test_temperature_validation_very_small_positive(self, sample_aluminum_element):
        """Test that temperatures below minimum allowed value raise error."""
        exception_caught = False
        try:
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(0.1),
                boiling_temperature=sp.Float(2792.0)
            )
        except Exception as e:
            if e.__class__.__name__ == "MaterialTemperatureError":
                exception_caught = True
                assert "is below minimum allowed value" in str(e)
                assert "0.1K" in str(e)
            else:
                pytest.fail(f"Wrong exception type raised: {type(e).__name__}: {e}")
        if not exception_caught:
            pytest.fail("Expected MaterialTemperatureError was not raised")

    def test_exception_handling_works_properly(self, sample_aluminum_element):
        """Verify that proper exception handling works now."""
        with pytest.raises(MaterialTemperatureError):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(-10.0),
                boiling_temperature=sp.Float(2792.0)
            )

    def test_temperature_validation_above_minimum_allowed(self, sample_aluminum_element):
        """Test that temperatures above minimum allowed value are valid."""
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(350.0),  # Above minimum allowed (302.0K)
            boiling_temperature=sp.Float(2792.0)
        )
        assert material is not None
        assert float(material.melting_temperature) == 350.0

    def test_material_calculated_properties(self, sample_aluminum_element):
        """Test that calculated properties are set correctly."""
        material = Material(
            name="Test Aluminum",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        # Test that calculated properties exist
        assert hasattr(material, 'atomic_number')
        assert hasattr(material, 'atomic_mass')
        assert material.atomic_number > 0
        assert material.atomic_mass > 0

    def test_alloy_calculated_properties(self, sample_steel_elements):
        """Test calculated properties for alloys."""
        material = Material(
            name="Test Steel",
            material_type="alloy",
            elements=sample_steel_elements,
            composition=[0.68, 0.02, 0.20, 0.10],
            solidus_temperature=sp.Float(1400.0),
            liquidus_temperature=sp.Float(1450.0),
            initial_boiling_temperature=sp.Float(2800.0),
            final_boiling_temperature=sp.Float(2900.0)
        )
        # Test interpolated properties
        assert hasattr(material, 'atomic_number')
        assert hasattr(material, 'atomic_mass')
        assert material.atomic_number > 0
        assert material.atomic_mass > 0
        # For steel, atomic number should be weighted average (mostly Fe = 26)
        assert 20 < material.atomic_number < 27  # Should be close to iron's atomic number

    def test_solidification_interval_alloy(self, sample_steel_elements):
        """Test solidification interval calculation for alloys."""
        material = Material(
            name="Test Steel",
            material_type="alloy",
            elements=sample_steel_elements,
            composition=[0.68, 0.02, 0.20, 0.10],
            solidus_temperature=sp.Float(1400.0),
            liquidus_temperature=sp.Float(1450.0),
            initial_boiling_temperature=sp.Float(2800.0),
            final_boiling_temperature=sp.Float(2900.0)
        )
        solidus, liquidus = material.solidification_interval()
        assert float(solidus) == 1400.0
        assert float(liquidus) == 1450.0
        assert float(liquidus) > float(solidus)

    def test_temperature_validation_minimum_boundary(self, sample_aluminum_element):
        """Test temperature validation at the minimum boundary."""
        # Test exactly at the minimum allowed temperature
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(302.0),  # Exactly at minimum allowed
            boiling_temperature=sp.Float(2792.0)
        )
        assert material is not None
        assert float(material.melting_temperature) == 302.0

    def test_temperature_validation_alloy_solidus_liquidus_order(self, sample_steel_elements):
        """Test that solidus temperature must be less than liquidus temperature."""
        with pytest.raises(MaterialTemperatureError, match=r"solidus_temperature.*must be less than.*liquidus_temperature"):
            Material(
                name="Test Alloy",
                material_type="alloy",
                elements=sample_steel_elements,
                composition=[0.68, 0.02, 0.20, 0.10],
                solidus_temperature=sp.Float(1500.0),  # Higher than liquidus
                liquidus_temperature=sp.Float(1400.0),  # Lower than solidus
                initial_boiling_temperature=sp.Float(2800.0),
                final_boiling_temperature=sp.Float(2900.0)
            )

    def test_temperature_validation_alloy_boiling_order(self, sample_steel_elements):
        """Test that initial boiling must be less than final boiling temperature."""
        with pytest.raises(MaterialTemperatureError, match=r"initial_boiling_temperature.*must be less than.*final_boiling_temperature"):
            Material(
                name="Test Alloy",
                material_type="alloy",
                elements=sample_steel_elements,
                composition=[0.68, 0.02, 0.20, 0.10],
                solidus_temperature=sp.Float(1400.0),
                liquidus_temperature=sp.Float(1450.0),
                initial_boiling_temperature=sp.Float(2900.0),  # Higher than final
                final_boiling_temperature=sp.Float(2800.0)    # Lower than initial
            )

    def test_temperature_validation_liquidus_boiling_order(self, sample_steel_elements):
        """Test that liquidus temperature must be less than initial boiling temperature."""
        with pytest.raises(MaterialTemperatureError, match=r"liquidus_temperature.*must be less than.*initial_boiling_temperature"):
            Material(
                name="Test Alloy",
                material_type="alloy",
                elements=sample_steel_elements,
                composition=[0.68, 0.02, 0.20, 0.10],
                solidus_temperature=sp.Float(1400.0),
                liquidus_temperature=sp.Float(2100.0),  # Higher than initial boiling
                initial_boiling_temperature=sp.Float(2000.0),
                final_boiling_temperature=sp.Float(2900.0)
            )

    def test_temperature_validation_pure_metal_missing_melting(self, sample_aluminum_element):
        """Test that pure metals must specify melting temperature."""
        with pytest.raises(MaterialTemperatureError, match="Pure metals must specify melting_temperature"):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=None,  # Missing
                boiling_temperature=sp.Float(2792.0)
            )

    def test_temperature_validation_pure_metal_missing_boiling(self, sample_aluminum_element):
        """Test that pure metals must specify boiling temperature."""
        with pytest.raises(MaterialTemperatureError, match="Pure metals must specify boiling_temperature"):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(933.47),
                boiling_temperature=None  # Missing
            )

    def test_temperature_validation_alloy_missing_temperatures(self, sample_steel_elements):
        """Test that alloys must specify all temperature properties."""
        with pytest.raises(MaterialTemperatureError, match="Alloys must specify all temperature properties"):
            Material(
                name="Test Alloy",
                material_type="alloy",
                elements=sample_steel_elements,
                composition=[0.68, 0.02, 0.20, 0.10],
                solidus_temperature=sp.Float(1400.0),
                liquidus_temperature=sp.Float(1450.0),
                initial_boiling_temperature=None,  # Missing
                final_boiling_temperature=sp.Float(2900.0)
            )

    def test_temperature_validation_wrong_type(self, sample_aluminum_element):
        """Test that temperatures must be SymPy Float type."""
        with pytest.raises(MaterialTemperatureError, match="must be a SymPy Float"):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=933.47,  # Regular float instead of sp.Float
                boiling_temperature=sp.Float(2792.0)
            )

    def test_temperature_validation_above_maximum(self, sample_aluminum_element):
        """Test that temperatures cannot exceed maximum allowed values."""
        with pytest.raises(MaterialTemperatureError, match="is above maximum allowed value"):
            Material(
                name="Test Material",
                material_type="pure_metal",
                elements=[sample_aluminum_element],
                composition=[1.0],
                melting_temperature=sp.Float(5000.0),  # Above MAX_MELTING_TEMP (3695.0)
                boiling_temperature=sp.Float(7000.0)
            )
