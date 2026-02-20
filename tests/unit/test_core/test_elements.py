"""Unit tests for ChemicalElement class."""

import pytest
from materforge.core.elements import ChemicalElement, interpolate, interpolate_atomic_number
from materforge.data import element_map, get_element

class TestChemicalElement:
    """Test cases for ChemicalElement class."""
    def test_element_creation_valid(self):
        """Test valid element creation."""
        element = ChemicalElement(
            name="Aluminum",
            atomic_number=13.0,
            atomic_mass=26.9815385,
            melting_temperature=933.47,
            boiling_temperature=2792.0,
            latent_heat_of_fusion=397000.0,
            latent_heat_of_vaporization=10900000.0
        )
        assert element.name == "Aluminum"
        assert element.atomic_number == 13.0
        assert element.atomic_mass == 26.9815385
        assert element.melting_temperature == 933.47
        assert element.boiling_temperature == 2792.0

    def test_element_from_data_map(self):
        """Test element creation from data map."""
        al_element = get_element('Al')
        assert al_element.atomic_number == 13.0
        assert al_element.name in ["Aluminum", "Aluminium"]
        assert al_element.atomic_mass > 0
        assert al_element.melting_temperature is None or al_element.melting_temperature > 0
        if al_element.boiling_temperature is not None and al_element.melting_temperature is not None:
            assert al_element.boiling_temperature > al_element.melting_temperature

    def test_element_interpolation_functions(self):
        """Test element interpolation utility functions."""
        # Test basic interpolation
        values = [10.0, 20.0, 30.0]
        composition = [0.5, 0.3, 0.2]
        result = interpolate(values, composition)
        expected = 0.5*10 + 0.3*20 + 0.2*30  # 17.0
        assert result == expected
        # Test atomic number interpolation
        fe = get_element('Fe')
        c = get_element('C')
        elements = [fe, c]
        composition = [0.98, 0.02]
        result = interpolate_atomic_number(elements, composition)
        expected = 0.98*26.0 + 0.02*6.0  # Steel approximation
        assert abs(result - expected) < 1e-10

    def test_element_properties_validation(self):
        """Test that element properties are valid where present."""
        # Force-load a reasonable subset
        symbols = ['C', 'N', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']
        for s in symbols:
            _ = get_element(s)

        for symbol, element in element_map.items():
            assert element.atomic_number > 0, f"Element {symbol} has invalid atomic number"
            assert element.atomic_mass > 0, f"Element {symbol} has invalid atomic mass"
            if element.melting_temperature is not None:
                assert element.melting_temperature > 0, f"Element {symbol} has invalid melting temperature"
            if element.boiling_temperature is not None and element.melting_temperature is not None:
                assert element.boiling_temperature > element.melting_temperature, \
                    f"Element {symbol} has boiling temp <= melting temp"
            if element.latent_heat_of_fusion is not None:
                assert element.latent_heat_of_fusion > 0, \
                    f"Element {symbol} has invalid latent heat of fusion"
            if element.latent_heat_of_vaporization is not None:
                assert element.latent_heat_of_vaporization > 0, \
                    f"Element {symbol} has invalid latent heat of vaporization"

    def test_common_elements_exist(self):
        """Test that common elements can be loaded and have valid core properties."""
        common_elements = ['C', 'N', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']
        for symbol in common_elements:
            element = get_element(symbol)
            assert symbol in element_map, f"Element {symbol} not found in element_map"
            assert isinstance(element, ChemicalElement)
            assert element.atomic_number > 0
            assert element.atomic_mass > 0
            assert element.melting_temperature is None or element.melting_temperature > 0
            if element.boiling_temperature is not None and element.melting_temperature is not None:
                assert element.boiling_temperature > element.melting_temperature
                
    def test_validate_element_completeness(self):
        """Test that incomplete elements are rejected for pure_metal context."""
        from materforge.core.elements import validate_element_completeness
        # Carbon has no melting point at 1 atm — must fail for pure_metal
        carbon = get_element('C')
        with pytest.raises(ValueError, match="missing required properties"):
            validate_element_completeness(carbon, context="pure_metal")
        # Iron is complete - must pass
        iron = get_element('Fe')
        validate_element_completeness(iron, context="pure_metal")  # no exception
        # Carbon is fine for alloy context (only atomic props needed)
        validate_element_completeness(carbon, context="alloy")  # no exception
