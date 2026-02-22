# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import scipy.constants as _sc


class PhysicalConstants:
    """
    Fundamental physical constants sourced from scipy.constants (CODATA 2022).

    References:
        - CODATA 2022: https://physics.nist.gov/cuu/Constants/
        - SciPy Constants: https://docs.scipy.org/doc/scipy/reference/constants.html

    Attributes:
        ABSOLUTE_ZERO (float): Absolute zero temperature, K.
        ROOM_TEMPERATURE (float): Standard room temperature, K (298.15 K).
        TRIPLE_POINT_WATER (float): Triple point of water, K (273.16 K).
        AMU (float): Atomic mass unit, kg.
        N_A (float): Avogadro's number, mol⁻¹.
        AVOGADRO_NUMBER (float): Alias for N_A.
        BOLTZMANN_CONSTANT (float): Boltzmann constant, J/K.
        GAS_CONSTANT (float): Universal gas constant, J/(mol·K).
        STEFAN_BOLTZMANN_CONSTANT (float): Stefan-Boltzmann constant, W/(m²·K⁴).
        GRAVITY (float): Acceleration due to gravity, m/s².
    """

    # -------------------------------------------------------------------------
    # Temperature references (conventional values, not in scipy)
    # -------------------------------------------------------------------------
    ABSOLUTE_ZERO: float = 0.0        # Absolute zero, K
    ROOM_TEMPERATURE: float = 298.15  # Standard room temperature, K
    TRIPLE_POINT_WATER: float = 273.16  # Triple point of water, K

    # -------------------------------------------------------------------------
    # Core constants used in MaterForge (sourced from scipy.constants)
    # -------------------------------------------------------------------------
    AMU: float = _sc.atomic_mass               # Atomic mass unit, kg
    N_A: float = _sc.Avogadro                  # Avogadro's number, mol⁻¹
    AVOGADRO_NUMBER: float = _sc.Avogadro      # Alias for N_A
    BOLTZMANN_CONSTANT: float = _sc.Boltzmann  # Boltzmann constant, J/K
    GAS_CONSTANT: float = _sc.R                # Universal gas constant, J/(mol·K)
    STEFAN_BOLTZMANN_CONSTANT: float = _sc.Stefan_Boltzmann  # W/(m²·K⁴)
    GRAVITY: float = _sc.g                     # m/s²

    @classmethod
    def get_all_constants(cls) -> dict:
        """Return a dictionary of all constants with their values."""
        return {name: getattr(cls, name) for name in dir(cls)
                if not name.startswith('_') and not callable(getattr(cls, name))}

    @classmethod
    def get_constant(cls, name: str) -> float:
        """Get a specific constant by name."""
        if hasattr(cls, name):
            return getattr(cls, name)
        available = [
            name for name in dir(cls)
            if not name.startswith('_') and not callable(getattr(cls, name))
        ]
        raise AttributeError(
            f"Constant '{name}' not found. Available constants: {available}"
        )
