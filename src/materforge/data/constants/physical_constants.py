from dataclasses import dataclass


@dataclass
class PhysicalConstants:
    """
    A dataclass to store fundamental constants.

    Attributes:
        ROOM_TEMPERATURE (float): Room temperature in Kelvin.
        N_A (float): Avogadro's number,
            the number of constituent particles (usually atoms or molecules) in one mole of a given substance.
        AMU (float): Atomic mass unit in kilograms.
        E (float): Elementary charge,
            the electric charge carried by a single proton
            or the magnitude of the electric charge carried by a single electron.
        SPEED_OF_LIGHT (float): Speed of light in vacuum in meters per second.
        BOLTZMANN_CONSTANT (float): Boltzmann constant in J/K.
        GAS_CONSTANT (float): Universal gas constant in J/(mol·K).
        PLANCK_CONSTANT (float): Planck constant in J·s.
        STEFAN_BOLTZMANN_CONSTANT (float): Stefan-Boltzmann constant in W/(m²·K⁴).
        GRAVITATIONAL_CONSTANT (float): Gravitational constant in m³/(kg·s²).
        PERMITTIVITY_VACUUM (float): Permittivity of free space in F/m.
        PERMEABILITY_VACUUM (float): Permeability of free space in H/m.
        ELECTRON_MASS (float): Electron rest mass in kg.
        PROTON_MASS (float): Proton rest mass in kg.
        NEUTRON_MASS (float): Neutron rest mass in kg.
        FINE_STRUCTURE_CONSTANT (float): Fine structure constant (dimensionless).
        RYDBERG_CONSTANT (float): Rydberg constant in m⁻¹.
        BOHR_RADIUS (float): Bohr radius in m.
        ELECTRON_VOLT (float): Electron volt in J.
        STANDARD_ATMOSPHERE (float): Standard atmospheric pressure in Pa.
        STANDARD_GRAVITY (float): Standard acceleration due to gravity in m/s².
        ABSOLUTE_ZERO (float): Absolute zero temperature in Kelvin.
        TRIPLE_POINT_WATER (float): Triple point of water in Kelvin.
        MOLAR_VOLUME_STP (float): Molar volume of ideal gas at STP in m³/mol.
        FARADAY_CONSTANT (float): Faraday constant in C/mol.
        WIEN_DISPLACEMENT_CONSTANT (float): Wien displacement constant in m·K.
    """

    # Temperature constants
    ROOM_TEMPERATURE: float = 298.15  # Room temperature in Kelvin
    ABSOLUTE_ZERO: float = 0.0  # Absolute zero in Kelvin
    TRIPLE_POINT_WATER: float = 273.16  # Triple point of water in Kelvin

    # Fundamental physical constants (CODATA 2018 values)
    SPEED_OF_LIGHT: float = 299792458.0  # Speed of light in vacuum, m/s (exact)
    PLANCK_CONSTANT: float = 6.62607015e-34  # Planck constant, J·s (exact)
    BOLTZMANN_CONSTANT: float = 1.380649e-23  # Boltzmann constant, J/K (exact)
    AVOGADRO_NUMBER: float = 6.02214076e23  # Avogadro's number, /mol (exact)
    N_A: float = 6.02214076e23  # Alias for Avogadro's number
    ELEMENTARY_CHARGE: float = 1.602176634e-19  # Elementary charge, C (exact)
    E: float = 1.602176634e-19  # Alias for elementary charge

    # Derived constants
    GAS_CONSTANT: float = 8.314462618  # Universal gas constant, J/(mol·K)
    FARADAY_CONSTANT: float = 96485.33212  # Faraday constant, C/mol

    # Electromagnetic constants
    PERMITTIVITY_VACUUM: float = 8.8541878128e-12  # Permittivity of free space, F/m
    PERMEABILITY_VACUUM: float = 1.25663706212e-6  # Permeability of free space, H/m
    FINE_STRUCTURE_CONSTANT: float = 7.2973525693e-3  # Fine structure constant (dimensionless)

    # Particle masses
    ELECTRON_MASS: float = 9.1093837015e-31  # Electron rest mass, kg
    PROTON_MASS: float = 1.67262192369e-27  # Proton rest mass, kg
    NEUTRON_MASS: float = 1.67492749804e-27  # Neutron rest mass, kg
    AMU: float = 1.66053906660e-27  # Atomic mass unit, kg

    # Atomic and quantum constants
    BOHR_RADIUS: float = 5.29177210903e-11  # Bohr radius, m
    RYDBERG_CONSTANT: float = 1.0973731568160e7  # Rydberg constant, m⁻¹

    # Energy conversion
    ELECTRON_VOLT: float = 1.602176634e-19  # Electron volt, J

    # Gravitational and mechanical constants
    GRAVITATIONAL_CONSTANT: float = 6.67430e-11  # Gravitational constant, m³/(kg·s²)
    STANDARD_GRAVITY: float = 9.80665  # Standard acceleration due to gravity, m/s² (exact)
    STANDARD_ATMOSPHERE: float = 101325.0  # Standard atmospheric pressure, Pa (exact)

    # Thermodynamic constants
    STEFAN_BOLTZMANN_CONSTANT: float = 5.670374419e-8  # Stefan-Boltzmann constant, W/(m²·K⁴)
    WIEN_DISPLACEMENT_CONSTANT: float = 2.897771955e-3  # Wien displacement constant, m·K
    MOLAR_VOLUME_STP: float = 0.02241396954  # Molar volume of ideal gas at STP, m³/mol

    # Mathematical constants (for convenience)
    PI: float = 3.141592653589793
    E_EULER: float = 2.718281828459045

    # Material science specific constants
    VACUUM_IMPEDANCE: float = 376.730313668  # Impedance of free space, Ω
    CONDUCTANCE_QUANTUM: float = 7.748091729e-5  # Conductance quantum, S
    MAGNETIC_FLUX_QUANTUM: float = 2.067833848e-15  # Magnetic flux quantum, Wb

    # Additional useful constants for materials modeling
    CLASSICAL_ELECTRON_RADIUS: float = 2.8179403262e-15  # Classical electron radius, m
    COMPTON_WAVELENGTH: float = 2.42631023867e-12  # Compton wavelength, m
    THOMSON_CROSS_SECTION: float = 6.6524587321e-29  # Thomson scattering cross-section, m²

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
        else:
            available = [name for name in dir(cls)
                         if not name.startswith('_') and not callable(getattr(cls, name))]
            raise AttributeError(f"Constant '{name}' not found. Available constants: {available}")
