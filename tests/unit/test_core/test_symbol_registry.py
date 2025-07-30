"""Unit tests for SymbolRegistry class."""

import pytest
import sympy as sp
from materforge.core.symbol_registry import SymbolRegistry

class TestSymbolRegistry:
    """Test cases for SymbolRegistry class."""
    def test_get_symbol(self):
        """Test symbol retrieval and creation."""
        SymbolRegistry.clear()
        symbol = SymbolRegistry.get("test_property")
        assert isinstance(symbol, sp.Symbol)
        assert str(symbol) == "test_property"
        symbol2 = SymbolRegistry.get("test_property")
        assert symbol == symbol2
        assert id(symbol) == id(symbol2)

    def test_get_all_symbols(self):
        """Test getting all registered symbols."""
        SymbolRegistry.clear()
        SymbolRegistry.get("prop1")
        SymbolRegistry.get("prop2")
        all_symbols = SymbolRegistry.get_all()
        assert len(all_symbols) == 2
        assert "prop1" in all_symbols
        assert "prop2" in all_symbols
        assert isinstance(all_symbols["prop1"], sp.Symbol)
        assert isinstance(all_symbols["prop2"], sp.Symbol)

    def test_clear_symbols(self):
        """Test clearing all symbols."""
        SymbolRegistry.clear()
        SymbolRegistry.get("temp_symbol")
        assert len(SymbolRegistry._symbols) == 1
        SymbolRegistry.clear()
        assert len(SymbolRegistry._symbols) == 0

    def test_symbol_persistence(self):
        """Test that symbols persist across calls."""
        SymbolRegistry.clear()
        symbol1 = SymbolRegistry.get("persistent_symbol")
        symbol2 = SymbolRegistry.get("persistent_symbol")
        assert symbol1 is symbol2
        assert id(symbol1) == id(symbol2)
