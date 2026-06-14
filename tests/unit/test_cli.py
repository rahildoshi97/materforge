# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for the materforge command-line interface."""

import logging
import shutil

import pytest

from materforge import __version__
from materforge.cli import main, validate_entry


@pytest.fixture(autouse=True)
def _restore_materforge_logger():
    """Snapshots and restores the materforge logger.

    The CLI reconfigures the package logger (handlers, level, propagate) on every
    run; without this, a CLI test would leave it muted and break other tests that
    rely on caplog capturing materforge warnings.
    """
    package_logger = logging.getLogger("materforge")
    saved = (package_logger.level, list(package_logger.handlers), package_logger.propagate)
    try:
        yield
    finally:
        package_logger.setLevel(saved[0])
        package_logger.handlers = saved[1]
        package_logger.propagate = saved[2]


# --- list ---------------------------------------------------------------

def test_list_shows_bundled_materials(capsys):
    assert main(["list"]) == 0
    out = capsys.readouterr().out
    assert "Al" in out
    assert "1.4301" in out


def test_list_paths_shows_file_locations(capsys):
    assert main(["list", "--paths"]) == 0
    out = capsys.readouterr().out
    assert "Al.yaml" in out
    assert ".yaml" in out


# --- validate -----------------------------------------------------------

def test_validate_valid_file(capsys, aluminum_yaml_path):
    assert main(["validate", str(aluminum_yaml_path)]) == 0
    assert "OK" in capsys.readouterr().out


def test_validate_missing_file_is_usage_error(aluminum_yaml_path):
    # A missing path is rejected by the argparse type, which exits with code 2.
    with pytest.raises(SystemExit) as exc:
        main(["validate", "does_not_exist.yaml"])
    assert exc.value.code == 2


def test_validate_invalid_yaml_returns_error(tmp_path, capsys):
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: Broken\nproperties:\n  density: {}\n")
    assert main(["validate", str(bad)]) == 1
    assert "Error:" in capsys.readouterr().err


def test_validate_entry_shim_forwards_to_validate(capsys, aluminum_yaml_path):
    assert validate_entry([str(aluminum_yaml_path)]) == 0
    assert "OK" in capsys.readouterr().out


# --- info ---------------------------------------------------------------

def test_info_reports_name_and_property_types(capsys, aluminum_yaml_path):
    assert main(["info", str(aluminum_yaml_path)]) == 0
    out = capsys.readouterr().out
    assert "Aluminum" in out
    assert "Property types:" in out
    assert "density" in out


# --- evaluate -----------------------------------------------------------

def test_evaluate_prints_numeric_values(capsys, aluminum_yaml_path):
    assert main(["evaluate", str(aluminum_yaml_path), "500"]) == 0
    out = capsys.readouterr().out
    assert "@ T = 500.0" in out
    assert "density" in out


def test_evaluate_respects_custom_symbol(capsys, aluminum_yaml_path):
    assert main(["evaluate", str(aluminum_yaml_path), "500", "-s", "u_C"]) == 0
    assert "@ u_C = 500.0" in capsys.readouterr().out


def test_evaluate_rejects_non_numeric_value(aluminum_yaml_path):
    with pytest.raises(SystemExit) as exc:
        main(["evaluate", str(aluminum_yaml_path), "not_a_number"])
    assert exc.value.code == 2


# --- plot ---------------------------------------------------------------

def test_plot_writes_png(tmp_path, capsys, aluminum_yaml_path):
    local_yaml = tmp_path / "Al.yaml"
    shutil.copy(aluminum_yaml_path, local_yaml)
    assert main(["plot", str(local_yaml)]) == 0
    out = capsys.readouterr().out
    assert "Built material" in out
    plots = list((tmp_path / "materforge_plots").glob("*.png"))
    assert plots, "expected a PNG to be written to the materforge_plots directory"


# --- top level ----------------------------------------------------------

def test_no_command_prints_help_and_fails(capsys):
    assert main([]) == 1
    assert "usage:" in capsys.readouterr().out


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out
