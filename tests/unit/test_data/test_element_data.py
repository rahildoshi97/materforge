"""Unit tests for element data."""

import pytest
from materforge.data import element_map, get_element
from materforge.core.elements import ChemicalElement

class TestElementData:
    """Test cases for element data."""

    def test_element_map_exists(self):
        """Test that element map is properly defined and can be populated."""
        assert isinstance(element_map, dict)
        # Force-load a couple of elements to ensure population works
        get_element('Al')
        get_element('Fe')
        assert len(element_map) >= 2

    def test_common_elements_exist(self):
        """Test that common elements can be loaded."""
        common_elements = ['C', 'N', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']
        for symbol in common_elements:
            el = get_element(symbol)
            assert symbol in element_map, f"Element {symbol} not found in element_map"
            assert isinstance(el, ChemicalElement)

    def test_aluminum_element_properties(self):
        """Test aluminum element properties."""
        al = get_element('Al')
        assert al.atomic_number == 13
        assert al.name in ['Aluminum', 'Aluminium']
        assert al.atomic_mass > 0
        if al.melting_temperature is not None:
            assert al.melting_temperature > 0
        if al.boiling_temperature is not None and al.melting_temperature is not None:
            assert al.boiling_temperature > al.melting_temperature

    def test_iron_element_properties(self):
        """Test iron element properties."""
        fe = get_element('Fe')
        assert fe.atomic_number == 26
        assert fe.name == 'Iron'
        assert fe.atomic_mass > 0
        if fe.melting_temperature is not None:
            assert fe.melting_temperature > 0
        if fe.boiling_temperature is not None and fe.melting_temperature is not None:
            assert fe.boiling_temperature > fe.melting_temperature

    def test_carbon_element_properties(self):
        """Test carbon element properties."""
        c = get_element('C')
        assert c.atomic_number == 6
        assert c.name == 'Carbon'
        assert c.atomic_mass > 0
        # Carbon often has no melting point at 1 atm; allow None
        if c.melting_temperature is not None:
            assert c.melting_temperature > 0
        if c.boiling_temperature is not None and c.melting_temperature is not None:
            assert c.boiling_temperature > c.melting_temperature

    def test_element_access_by_symbol(self):
        """Test element access by symbol via map after loading."""
        al = get_element('Al')
        fe = get_element('Fe')
        assert element_map['Al'] is al
        assert element_map['Fe'] is fe
        assert isinstance(al, ChemicalElement)
        assert isinstance(fe, ChemicalElement)
        assert al.atomic_number == 13
        assert fe.atomic_number == 26

    def test_element_access_invalid_symbol(self):
        """Test element access with invalid symbols via get_element."""
        with pytest.raises(KeyError):
            get_element('Xx')
        with pytest.raises(KeyError):
            get_element('invalid')

    def test_element_access_case_sensitivity(self):
        """Test element access case sensitivity."""
        al = get_element('Al')
        assert al.atomic_number == 13
        with pytest.raises(KeyError):
            get_element('al')
        with pytest.raises(KeyError):
            get_element('AL')

    def test_element_atomic_numbers_unique(self):
        """Test that all cached elements have unique atomic numbers."""
        # Load a representative set
        for s in ['H', 'He', 'C', 'N', 'O', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']:
            get_element(s)
        atomic_numbers = [el.atomic_number for el in element_map.values()]
        assert len(atomic_numbers) == len(set(atomic_numbers)), "Duplicate atomic numbers found"

    def test_element_symbols_unique(self):
        """Test that all element symbols are unique."""
        symbols = list(element_map.keys())
        assert len(symbols) == len(set(symbols)), "Duplicate symbols found"

    def test_element_atomic_masses_positive(self):
        """Test that all cached atomic masses are positive."""
        for symbol, element in element_map.items():
            assert element.atomic_mass > 0, f"Element {symbol} has non-positive atomic mass"

    def test_element_atomic_numbers_positive(self):
        """Test that all cached atomic numbers are positive."""
        for symbol, element in element_map.items():
            assert element.atomic_number > 0, f"Element {symbol} has non-positive atomic number"

    def test_element_names_non_empty(self):
        """Test that all cached element names are non-empty strings."""
        for symbol, element in element_map.items():
            assert isinstance(element.name, str), f"Element {symbol} has non-string name"
            assert element.name, f"Element {symbol} has empty name"

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
            el = get_element(symbol)
            assert symbol in element_map, f"Steel element {symbol} not available"
            assert isinstance(el, ChemicalElement)
            assert el.atomic_number > 0

    def test_aluminum_alloy_elements_available(self):
        """Test that common aluminum alloy elements are available."""
        al_alloy_elements = ['Al', 'Cu', 'Si', 'Mn']
        for symbol in al_alloy_elements:
            el = get_element(symbol)
            assert symbol in element_map, f"Aluminum alloy element {symbol} not available"
            assert isinstance(el, ChemicalElement)
            assert el.atomic_number > 0

    def test_element_temperature_properties(self):
        """Test that cached elements have valid temperature properties where present."""
        # Load a representative set
        for s in ['C', 'N', 'Al', 'Fe', 'Cu', 'Ni', 'Cr']:
            get_element(s)
        for symbol, element in element_map.items():
            if element.melting_temperature is not None:
                assert element.melting_temperature > 0, f"{symbol}: invalid melting temperature"
            if element.boiling_temperature is not None:
                assert element.boiling_temperature > 0, f"{symbol}: invalid boiling temperature"
            if (
                element.boiling_temperature is not None
                and element.melting_temperature is not None
            ):
                assert element.boiling_temperature > element.melting_temperature, \
                    f"{symbol}: boiling temperature <= melting temperature"

    def test_element_latent_heat_properties(self):
        """Test that cached elements have valid latent heat properties where present."""
        for symbol in ['Al', 'Fe', 'Cu', 'Ni', 'Cr']:
            get_element(symbol)
        for symbol, element in element_map.items():
            if element.latent_heat_of_fusion is not None:
                assert element.latent_heat_of_fusion > 0, f"{symbol}: invalid latent heat of fusion"
            if element.latent_heat_of_vaporization is not None:
                assert element.latent_heat_of_vaporization > 0, \
                    f"{symbol}: invalid latent heat of vaporization"
            if (
                element.latent_heat_of_fusion is not None
                and element.latent_heat_of_vaporization is not None
            ):
                assert element.latent_heat_of_vaporization > element.latent_heat_of_fusion, \
                    f"{symbol}: latent heat of vaporization <= latent heat of fusion"

    def test_get_element_function_valid_symbol(self):
        """Test get_element function with valid symbols."""
        al = get_element('Al')
        fe = get_element('Fe')
        assert isinstance(al, ChemicalElement)
        assert isinstance(fe, ChemicalElement)
        assert al.atomic_number == 13
        assert fe.atomic_number == 26

    def test_get_element_function_invalid_symbol(self):
        """Test get_element function with invalid symbols."""
        with pytest.raises(KeyError):
            get_element('Xx')  # Non-existent element
        with pytest.raises(KeyError):
            get_element('invalid')

    def test_get_element_function_case_sensitivity(self):
        """Test get_element function case sensitivity."""
        al = get_element('Al')
        assert al.atomic_number == 13
        with pytest.raises(KeyError):
            get_element('al')
        with pytest.raises(KeyError):
            get_element('AL')
