import logging
import sympy as sp
from typing import Dict

logger = logging.getLogger(__name__)


class SymbolRegistry:
    """Registry for SymPy symbols to ensure uniqueness."""
    _symbols = {}

    @classmethod
    def get(cls, name: str) -> sp.Symbol:
        """Get or create a symbol with the given name."""
        if name not in cls._symbols:
            cls._symbols[name] = sp.Symbol(name)
            logger.debug("Created new symbol: %s", name)
        else:
            logger.debug("Retrieved existing symbol: %s", name)
        return cls._symbols[name]

    @classmethod
    def get_all(cls) -> Dict[str, sp.Symbol]:
        """Get all registered symbols."""
        logger.debug("Retrieved all symbols, count: %d", len(cls._symbols))
        return cls._symbols.copy()

    @classmethod
    def clear(cls) -> None:
        """Clear all registered symbols."""
        count = len(cls._symbols)
        cls._symbols.clear()
        logger.info("Cleared %d symbols from registry", count)
