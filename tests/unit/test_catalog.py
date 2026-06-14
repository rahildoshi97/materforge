# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for the bundled-material catalog helpers."""

import pytest
import sympy as sp

from materforge import list_materials, load_material, get_material_path
from materforge.core.materials import Material


def test_list_materials_returns_sorted_names():
    names = list_materials()
    assert names == sorted(names)
    assert all(isinstance(name, str) for name in names)


def test_list_materials_includes_shipped_examples():
    # The reference YAMLs that ship with the package.
    assert {"Al", "Al2O3", "1.4301"} <= set(list_materials())


def test_get_material_path_points_at_existing_yaml():
    path = get_material_path("Al")
    assert path.is_file()
    assert path.name == "Al.yaml"


def test_get_material_path_is_case_insensitive():
    assert get_material_path("al") == get_material_path("Al")


def test_get_material_path_unknown_name_lists_options():
    with pytest.raises(ValueError) as exc:
        get_material_path("unobtainium")
    message = str(exc.value)
    assert "unobtainium" in message
    assert "Al" in message  # the error should point at what is available


def test_load_material_returns_populated_material(temp_symbol):
    material = load_material("Al", temp_symbol)
    assert isinstance(material, Material)
    assert material.property_names()


def test_load_material_resolves_companion_file(temp_symbol):
    # 1.4301 pulls heat_capacity and density from the sibling .xlsx via
    # FILE_IMPORT, so a successful load proves companion-file resolution works.
    material = load_material("1.4301", temp_symbol)
    assert "heat_capacity" in material.property_names()


def test_load_material_case_insensitive(temp_symbol):
    material = load_material("al", temp_symbol)
    assert isinstance(material, Material)


def test_load_material_rejects_non_symbol_dependency():
    with pytest.raises(TypeError):
        load_material("Al", "T")


def test_load_material_unknown_name_raises():
    with pytest.raises(ValueError):
        load_material("not_a_material", sp.Symbol("T"))
