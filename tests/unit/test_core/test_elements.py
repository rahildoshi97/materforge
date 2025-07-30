"""Unit tests for ChemicalElement class."""

import pytest
from materforge.core.elements import ChemicalElement, interpolate, interpolate_atomic_number
from materforge.data.elements.element_data import element_map

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
        al_element = element_map['Al']
        assert al_element.atomic_number == 13.0
        assert al_element.name in ["Aluminum", "Aluminium"]
        assert al_element.atomic_mass > 0
        assert al_element.melting_temperature > 0
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
        elements = [element_map['Fe'], element_map['C']]
        composition = [0.98, 0.02]
        result = interpolate_atomic_number(elements, composition)
        expected = 0.98*26.0 + 0.02*6.0  # Steel approximation
        assert abs(result - expected) < 1e-10

    def test_element_properties_validation(self):
        """Test that element properties are valid."""
        for symbol, element in element_map.items():
            assert element.atomic_number > 0, f"Element {symbol} has invalid atomic number"
            assert element.atomic_mass > 0, f"Element {symbol} has invalid atomic mass"
            assert element.melting_temperature > 0, f"Element {symbol} has invalid melting temperature"
            assert element.boiling_temperature > element.melting_temperature, \
                f"Element {symbol} has boiling temp <= melting temp"
            assert element.latent_heat_of_fusion > 0, f"Element {symbol} has invalid latent heat of fusion"
            assert element.latent_heat_of_vaporization > 0, f"Element {symbol} has invalid latent heat of vaporization"

    def test_common_elements_exist(self):
        """Test that common elements exist in the map."""
        common_elements = ['C', 'N', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']
        for symbol in common_elements:
            assert symbol in element_map, f"Element {symbol} not found in element_map"
            element = element_map[symbol]
            assert isinstance(element, ChemicalElement)
            assert element.atomic_number > 0
            assert element.atomic_mass > 0
