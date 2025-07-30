import logging
from materforge.data.constants import PhysicalConstants
from materforge.core.elements import ChemicalElement

logger = logging.getLogger(__name__)

# NIST: National Institute of Standards and Technology
# RSC: Royal Society of Chemistry
# CRC: CRC Handbook of Chemistry and Physics

CARBON = ChemicalElement(
    name="Carbon",
    atomic_number=6,
    atomic_mass=12.0107 * PhysicalConstants.AMU,
    melting_temperature=3915,  # Melting temperature = 3915 K
    boiling_temperature=4300,  # Boiling temperature = 4300 K
    latent_heat_of_fusion=117000,  # Latent heat of fusion = 117 kJ/mol
    latent_heat_of_vaporization=355000  # Latent heat of vaporization = 355 kJ/mol
)

NITROGEN = ChemicalElement(
    name="Nitrogen",
    atomic_number=7,
    atomic_mass=14.0067 * PhysicalConstants.AMU,
    melting_temperature=63.15,  # Melting temperature = 63.15 K
    boiling_temperature=77.36,  # Boiling temperature = 77.36 K
    latent_heat_of_fusion=720,  # Latent heat of fusion = 0.72 kJ/mol
    latent_heat_of_vaporization=5570  # Latent heat of vaporization = 5.57 kJ/mol
)

ALUMINIUM = ChemicalElement(
    name="Aluminium",
    atomic_number=13,  # Atomic number = 13 / Source: Periodic Table
    atomic_mass=26.9815384 * PhysicalConstants.AMU,  # Atomic mass = 26.9815384 amu / Source: NIST
    melting_temperature=933.35,  # Melting temperature = 933.35 K / Source: RSC
    boiling_temperature=2743,  # Boiling temperature = 2743 K / Source: RSC
    latent_heat_of_fusion=10700,  # Latent heat of fusion = 10700 J/kg / Source: CRC
    latent_heat_of_vaporization=284000  # Latent heat of vaporization = 284000 J/kg / Source: CRC
)

SILICON = ChemicalElement(
    name="Silicon",
    atomic_number=14,
    atomic_mass=28.0855 * PhysicalConstants.AMU,
    melting_temperature=1687,  # Melting temperature = 1687 K
    boiling_temperature=3538,  # Boiling temperature = 3538 K
    latent_heat_of_fusion=50200,  # Latent heat of fusion = 50.2 kJ/mol
    latent_heat_of_vaporization=359000  # Latent heat of vaporization = 359 kJ/mol
)

PHOSPHORUS = ChemicalElement(
    name="Phosphorus",
    atomic_number=15,
    atomic_mass=30.973762 * PhysicalConstants.AMU,
    melting_temperature=317.3,  # Melting temperature = 317.3 K
    boiling_temperature=553.7,  # Boiling temperature = 553.7 K
    latent_heat_of_fusion=2510,  # Latent heat of fusion = 2.51 kJ/mol
    latent_heat_of_vaporization=12400  # Latent heat of vaporization = 12.4 kJ/mol
)

SULFUR = ChemicalElement(
    name="Sulfur",
    atomic_number=16,
    atomic_mass=32.065 * PhysicalConstants.AMU,
    melting_temperature=388.36,  # Melting temperature = 388.36 K
    boiling_temperature=717.8,  # Boiling temperature = 717.8 K
    latent_heat_of_fusion=1730,  # Latent heat of fusion = 1.73 kJ/mol
    latent_heat_of_vaporization=9800  # Latent heat of vaporization = 9.8 kJ/mol
)

TITANIUM = ChemicalElement(
    name="Titanium",
    atomic_number=22,  # Atomic number = 22 / Source: Periodic Table
    atomic_mass=47.867 * PhysicalConstants.AMU,  # Atomic mass = 47.867 amu / Source: NIST
    melting_temperature=1941,  # Melting temperature = 1941 K / Source: RSC
    boiling_temperature=3560,  # Boiling temperature = 3560 K / Source: RSC
    latent_heat_of_fusion=18700,  # Latent heat of fusion = 18700 J/kg / Source: CRC
    latent_heat_of_vaporization=427000  # Latent heat of vaporization = 427000 J/kg / Source: CRC
)

VANADIUM = ChemicalElement(
    name="Vanadium",
    atomic_number=23,  # Atomic number = 23 / Source: Periodic Table
    atomic_mass=50.9415 * PhysicalConstants.AMU,  # Atomic mass = 50.9415 amu / Source: NIST
    melting_temperature=2183,  # Melting temperature = 2183 K / Source: RSC
    boiling_temperature=3680,  # Boiling temperature = 3680 K / Source: RSC
    latent_heat_of_fusion=21500,  # Latent heat of fusion = 21500 J/kg / Source: CRC
    latent_heat_of_vaporization=444000  # Latent heat of vaporization = 444000 J/kg / Source: CRC
)

CHROMIUM = ChemicalElement(
    name="Chromium",
    atomic_number=24,  # Atomic number = 24 / Source: Periodic Table
    atomic_mass=51.9961 * PhysicalConstants.AMU,  # Atomic mass = 51.9961 amu / Source: NIST
    melting_temperature=2180,  # Melting temperature = 2180 K / Source: RSC
    boiling_temperature=2944,  # Boiling temperature = 2944 K / Source: RSC
    latent_heat_of_fusion=16500,  # Latent heat of fusion = 16500 J/kg / Source: CRC
    latent_heat_of_vaporization=344000  # Latent heat of vaporization = 344000 J/kg / Source: CRC
)

MANGANESE = ChemicalElement(
    name="Manganese",
    atomic_number=25,  # Atomic number = 25 / Source: Periodic Table
    atomic_mass=54.938045 * PhysicalConstants.AMU,  # Atomic mass = 54.938045 amu / Source: NIST
    melting_temperature=1519,  # Melting temperature = 1519 K / Source: RSC
    boiling_temperature=2334,  # Boiling temperature = 2334 K / Source: RSC
    latent_heat_of_fusion=12500,  # Latent heat of fusion = 12500 J/kg / Source: CRC
    latent_heat_of_vaporization=220000  # Latent heat of vaporization = 220000 J/kg / Source: CRC
)

IRON = ChemicalElement(
    name="Iron",
    atomic_number=26,  # Atomic number = 26 / Source: Periodic Table
    atomic_mass=55.845 * PhysicalConstants.AMU,  # Atomic mass = 55.845 amu / Source: NIST
    melting_temperature=1809,  # Melting temperature = 1809 K / Source: RSC
    boiling_temperature=3134,  # Boiling temperature = 3134 K / Source: RSC
    latent_heat_of_fusion=13800,  # Latent heat of fusion = 13800 J/kg / Source: CRC
    latent_heat_of_vaporization=340000  # Latent heat of vaporization = 340000 J/kg / Source: CRC
)

NICKEL = ChemicalElement(
    name="Nickel",
    atomic_number=28,  # Atomic number = 28 / Source: Periodic Table
    atomic_mass=58.6934 * PhysicalConstants.AMU,  # Atomic mass = 58.6934 amu / Source: NIST
    melting_temperature=1728,  # Melting temperature = 1728 K / Source: RSC
    boiling_temperature=3186,  # Boiling temperature = 3186 K / Source: RSC
    latent_heat_of_fusion=17200,  # Latent heat of fusion = 17200 J/kg / Source: CRC
    latent_heat_of_vaporization=377000  # Latent heat of vaporization = 377000 J/kg / Source: CRC
)

COPPER = ChemicalElement(
    name="Copper",
    atomic_number=29,  # Atomic number = 29 / Source: Periodic Table
    atomic_mass=63.546 * PhysicalConstants.AMU,  # Atomic mass = 63.546 amu / Source: NIST
    melting_temperature=1357.77,  # Melting temperature = 1357.77 K / Source: RSC
    boiling_temperature=2835,  # Boiling temperature = 2835 K / Source: RSC
    latent_heat_of_fusion=209000,  # Latent heat of fusion = 20500 J/kg / Source: CRC
    latent_heat_of_vaporization=4730000.0  # Latent heat of vaporization = 453000 J/kg / Source: CRC
)

MOLYBDENUM = ChemicalElement(
    name="Molybdenum",
    atomic_number=42,
    atomic_mass=95.96 * PhysicalConstants.AMU,
    melting_temperature=2896,  # Melting temperature = 2896K (2623°C)
    boiling_temperature=4912,  # Boiling temperature = 4912K (4639°C)
    latent_heat_of_fusion=37480,  # Latent heat of fusion = 37.48 kJ/mol
    latent_heat_of_vaporization=598000  # Latent heat of vaporization = 598 kJ/mol
)

# This dictionary maps chemical symbols (strings) to their corresponding ChemicalElement instances,
# allowing the parser to convert composition keys from the YAML file (like 'Fe': 0.675)
# to actual ChemicalElement objects needed by the Alloy class.
# TODO: Update dictionary as more elements are added!
element_map = {
    'C': CARBON,
    'N': NITROGEN,
    'Al': ALUMINIUM,
    'Si': SILICON,
    'P': PHOSPHORUS,
    'S': SULFUR,
    'Ti': TITANIUM,
    'V': VANADIUM,
    'Cr': CHROMIUM,
    'Mn': MANGANESE,
    'Fe': IRON,
    'Ni': NICKEL,
    'Cu': COPPER,
    'Mo': MOLYBDENUM,
}


def get_element(symbol: str) -> ChemicalElement:
    """Get element by symbol with error handling."""
    logger.debug("Looking up element: %s", symbol)
    if symbol not in element_map:
        available_symbols = list(element_map.keys())
        logger.error("Element symbol '%s' not found. Available: %s", symbol, available_symbols)
        raise KeyError(f"Element with symbol '{symbol}' not found. Available elements: {available_symbols}")
    element = element_map[symbol]
    logger.debug("Found element: %s (atomic number: %d)", element.name, element.atomic_number)
    return element


# Log available elements on module import
logger.info("Element database loaded with %d elements: %s", len(element_map), list(element_map.keys()))
