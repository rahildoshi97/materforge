Welcome to materforge's documentation!
=======================================

materforge is a Python library for materials science and thermophysical property calculations.

Installation
============

Install materforge using pip:

.. code-block:: bash

   pip install materforge

Or install from source:

.. code-block:: bash

   git clone https://github.com/rahildoshi97/materforge.git
   cd materforge
   pip install -e .

Quick Start
===========

Here's a simple example to get you started:

.. code-block:: python

   import sympy as sp
   from materforge import create_material

   # Create a material with symbolic temperature
   T = sp.Symbol('T')
   material = create_material('steel.yaml', T)

   # Evaluate properties at a specific temperature
   properties = material.evaluate_properties_at_temperature(500.0)
   print(properties)

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   tutorials/getting_started
   tutorials/first_simulation
   how-to/define_materials
   how-to/energy_temperature_conversion
   explanation/design_philosophy
   explanation/material_properties
   reference/yaml_schema
   reference/api
   reference/api/material

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
