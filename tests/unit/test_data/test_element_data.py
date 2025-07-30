"""Unit tests for element data."""

import pytest
from materforge.data.elements.element_data import element_map, get_element
from materforge.core.elements import ChemicalElement

class TestElementData:
    """Test cases for element data."""
    def test_element_map_exists(self):
        """Test that element map is properly defined."""
        assert isinstance(element_map, dict)
        assert len(element_map) > 0

    def test_common_elements_exist(self):
        """Test that common elements exist in the map."""
        common_elements = ['C', 'N', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']
        for symbol in common_elements:
            assert symbol in element_map, f"Element {symbol} not found in element_map"
            assert isinstance(element_map[symbol], ChemicalElement)

    def test_aluminum_element_properties(self):
        """Test aluminum element properties."""
        al = element_map['Al']
        assert al.atomic_number == 13.0
        assert al.name in ['Aluminum', 'Aluminium']  # Both spellings accepted
        assert al.atomic_mass > 0
        assert al.melting_temperature > 0
        assert al.boiling_temperature > al.melting_temperature

    def test_iron_element_properties(self):
        """Test iron element properties."""
        fe = element_map['Fe']
        assert fe.atomic_number == 26.0
        assert fe.name == 'Iron'
        assert fe.atomic_mass > 0
        assert fe.melting_temperature > 0
        assert fe.boiling_temperature > fe.melting_temperature

    def test_carbon_element_properties(self):
        """Test carbon element properties."""
        c = element_map['C']
        assert c.atomic_number == 6.0
        assert c.name == 'Carbon'
        assert c.atomic_mass > 0
        assert c.melting_temperature > 0
        assert c.boiling_temperature > c.melting_temperature

    def test_element_access_by_symbol(self):
        """Test element access by symbol."""
        # Test valid symbols
        al = element_map['Al']
        assert isinstance(al, ChemicalElement)
        assert al.atomic_number == 13.0
        fe = element_map['Fe']
        assert isinstance(fe, ChemicalElement)
        assert fe.atomic_number == 26.0

    def test_element_access_invalid_symbol(self):
        """Test element access with invalid symbols."""
        with pytest.raises(KeyError):
            element_map['Xx']  # Non-existent element
        with pytest.raises(KeyError):
            element_map['invalid']

    def test_element_access_case_sensitivity(self):
        """Test element access case sensitivity."""
        # Should work with correct case
        al = element_map['Al']
        assert al.atomic_number == 13.0
        # Should fail with incorrect case
        with pytest.raises(KeyError):
            element_map['al']  # lowercase
        with pytest.raises(KeyError):
            element_map['AL']  # uppercase

    def test_element_atomic_numbers_unique(self):
        """Test that all elements have unique atomic numbers."""
        atomic_numbers = [element.atomic_number for element in element_map.values()]
        assert len(atomic_numbers) == len(set(atomic_numbers)), "Duplicate atomic numbers found"

    def test_element_symbols_unique(self):
        """Test that all element symbols are unique."""
        symbols = list(element_map.keys())
        assert len(symbols) == len(set(symbols)), "Duplicate symbols found"

    def test_element_atomic_masses_positive(self):
        """Test that all atomic masses are positive."""
        for symbol, element in element_map.items():
            assert element.atomic_mass > 0, f"Element {symbol} has non-positive atomic mass"

    def test_element_atomic_numbers_positive(self):
        """Test that all atomic numbers are positive."""
        for symbol, element in element_map.items():
            assert element.atomic_number > 0, f"Element {symbol} has non-positive atomic number"

    def test_element_names_non_empty(self):
        """Test that all element names are non-empty strings."""
        for symbol, element in element_map.items():
            assert isinstance(element.name, str), f"Element {symbol} has non-string name"
            assert len(element.name) > 0, f"Element {symbol} has empty name"

    def test_element_symbols_valid_format(self):
        """Test that element symbols follow valid format."""
        for symbol in element_map.keys():
            assert isinstance(symbol, str), f"Symbol {symbol} is not a string"
            assert 1 <= len(symbol) <= 2, f"Symbol {symbol} has invalid length"
            assert symbol[0].isupper(), f"Symbol {symbol} doesn't start with uppercase"
            if len(symbol) == 2:
                assert symbol[1].islower(), f"Symbol {symbol} second character not lowercase"

    def test_steel_alloy_elements_available(self):
        """Test that common steel alloy elements are available."""
        steel_elements = ['Fe', 'Cr', 'Ni', 'Mn', 'Mo']
        for symbol in steel_elements:
            assert symbol in element_map, f"Steel element {symbol} not available"
            element = element_map[symbol]
            assert isinstance(element, ChemicalElement)
            assert element.atomic_number > 0

    def test_aluminum_alloy_elements_available(self):
        """Test that common aluminum alloy elements are available."""
        al_alloy_elements = ['Al', 'Cu', 'Si', 'Mn']
        for symbol in al_alloy_elements:
            assert symbol in element_map, f"Aluminum alloy element {symbol} not available"
            element = element_map[symbol]
            assert isinstance(element, ChemicalElement)
            assert element.atomic_number > 0

    def test_element_temperature_properties(self):
        """Test that elements have valid temperature properties."""
        for symbol, element in element_map.items():
            assert element.melting_temperature > 0, f"Element {symbol} has invalid melting temperature"
            assert element.boiling_temperature > 0, f"Element {symbol} has invalid boiling temperature"
            assert element.boiling_temperature > element.melting_temperature, \
                f"Element {symbol} has boiling temperature <= melting temperature"

    def test_element_latent_heat_properties(self):
        """Test that elements have valid latent heat properties."""
        for symbol, element in element_map.items():
            assert element.latent_heat_of_fusion > 0, f"Element {symbol} has invalid latent heat of fusion"
            assert element.latent_heat_of_vaporization > 0, f"Element {symbol} has invalid latent heat of vaporization"
            assert element.latent_heat_of_vaporization > element.latent_heat_of_fusion, \
                f"Element {symbol} has latent heat of vaporization <= latent heat of fusion"

    def test_get_element_function_valid_symbol(self):
        """Test get_element function with valid symbols."""
        al = get_element('Al')
        assert isinstance(al, ChemicalElement)
        assert al.atomic_number == 13.0
        fe = get_element('Fe')
        assert isinstance(fe, ChemicalElement)
        assert fe.atomic_number == 26.0

    def test_get_element_function_invalid_symbol(self):
        """Test get_element function with invalid symbols."""
        with pytest.raises(KeyError):
            get_element('Xx')  # Non-existent element
        with pytest.raises(KeyError):
            get_element('invalid')

    def test_get_element_function_case_sensitivity(self):
        """Test get_element function case sensitivity."""
        # Should work with correct case
        al = get_element('Al')
        assert al.atomic_number == 13.0
        # Should fail with incorrect case
        with pytest.raises(KeyError):
            get_element('al')  # lowercase
        with pytest.raises(KeyError):
            get_element('AL')  # uppercase
