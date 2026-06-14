# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

"""Command-line interface for MaterForge.

Exposes the most common library operations as a ``materforge`` command so YAML
material files can be inspected, validated, plotted, and evaluated without
writing any Python. Every subcommand is a thin wrapper over the public API in
:mod:`materforge`; see ``materforge --help`` for the full list.
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional, Sequence

import sympy as sp

from materforge import __version__
from materforge.catalog import get_material_path, list_materials
from materforge.parsing.api import (
    create_material,
    get_material_info,
    validate_yaml_file,
)

logger = logging.getLogger(__name__)

# Exit codes: 0 success, 1 a handled error (bad file, invalid YAML, ...),
# 2 is reserved by argparse for usage errors.
_EXIT_OK = 0
_EXIT_ERROR = 1


# ====================================================================
# SUBCOMMANDS
# ====================================================================

def _cmd_list(args: argparse.Namespace) -> int:
    """Lists the example material files bundled with the package."""
    names = list_materials()
    if not names:
        print("No bundled materials found.", file=sys.stderr)
        return _EXIT_ERROR
    if args.paths:
        width = max(len(name) for name in names)
        for name in names:
            print(f"{name:<{width}}  {get_material_path(name)}")
    else:
        print(f"Bundled example materials ({len(names)}):")
        for name in names:
            print(f"  {name}")
    return _EXIT_OK


def _cmd_validate(args: argparse.Namespace) -> int:
    """Validates a YAML material file without building the full material."""
    validate_yaml_file(args.yaml_path)
    print(f"OK: {args.yaml_path} is valid.")
    return _EXIT_OK


def _cmd_info(args: argparse.Namespace) -> int:
    """Prints a material's metadata without fully processing its properties."""
    info = get_material_info(args.yaml_path)
    print(f"Material: {info.get('name', 'Unknown')}")
    print(f"Source:   {args.yaml_path}")

    properties = info.get("properties", [])
    print(f"Properties ({len(properties)}):")
    for name in properties:
        print(f"  - {name}")

    property_types = info.get("property_types")
    if property_types:
        print("Property types:")
        for type_name, count in property_types.items():
            print(f"  {type_name}: {count}")

    reserved = {"name", "properties", "total_properties", "property_types"}
    extras = {key: value for key, value in info.items() if key not in reserved}
    if extras:
        print("Other fields:")
        for key, value in extras.items():
            print(f"  {key}: {value}")
    return _EXIT_OK


def _cmd_plot(args: argparse.Namespace) -> int:
    """Builds a material with plotting enabled and reports where plots landed."""
    symbol = sp.Symbol(args.symbol)
    material = create_material(args.yaml_path, symbol, enable_plotting=True)
    plot_dir = args.yaml_path.resolve().parent / "materforge_plots"
    print(f"Built material '{material.name}' ({len(material.property_names())} properties).")
    print(f"Plots saved under: {plot_dir}")
    return _EXIT_OK


def _cmd_evaluate(args: argparse.Namespace) -> int:
    """Evaluates every property at a single dependency value."""
    symbol = sp.Symbol(args.symbol)
    material = create_material(args.yaml_path, symbol, enable_plotting=False)
    evaluated = material.evaluate(symbol, args.value)

    names = sorted(evaluated.property_names())
    if not names:
        print(f"No properties could be evaluated at {args.symbol}={args.value}.",
              file=sys.stderr)
        return _EXIT_ERROR

    print(f"{material.name} @ {args.symbol} = {args.value}")
    width = max(len(name) for name in names)
    for name in names:
        print(f"  {name:<{width}}  {_format_value(evaluated.properties[name])}")
    return _EXIT_OK


def _format_value(value: sp.Basic) -> str:
    """Formats an evaluated property value for display, preferring a short float."""
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return str(value)


# ====================================================================
# PARSER
# ====================================================================

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="materforge",
        description="Inspect, validate, plot, and evaluate MaterForge YAML materials.",
    )
    parser.add_argument("--version", action="version",
                        version=f"materforge {__version__}")

    # Shared options that apply to every subcommand.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-v", "--verbose", action="count", default=0,
                        help="increase log verbosity (-v for info, -vv for debug)")

    # Options shared by subcommands that build a material from a symbol.
    symbol_opt = argparse.ArgumentParser(add_help=False)
    symbol_opt.add_argument("-s", "--symbol", default="T",
                            help="dependency symbol to substitute for the YAML "
                                 "placeholder 'T' (default: T)")

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    p_list = subparsers.add_parser(
        "list", parents=[common], help="list bundled example materials")
    p_list.add_argument("--paths", action="store_true",
                        help="also show the YAML file path for each material")
    p_list.set_defaults(func=_cmd_list)

    p_validate = subparsers.add_parser(
        "validate", parents=[common], help="validate a YAML material file")
    p_validate.add_argument("yaml_path", type=_existing_path,
                            help="path to the YAML material file")
    p_validate.set_defaults(func=_cmd_validate)

    p_info = subparsers.add_parser(
        "info", parents=[common], help="show a material's metadata")
    p_info.add_argument("yaml_path", type=_existing_path,
                        help="path to the YAML material file")
    p_info.set_defaults(func=_cmd_info)

    p_plot = subparsers.add_parser(
        "plot", parents=[common, symbol_opt],
        help="build a material and save its property plots")
    p_plot.add_argument("yaml_path", type=_existing_path,
                        help="path to the YAML material file")
    p_plot.set_defaults(func=_cmd_plot)

    p_evaluate = subparsers.add_parser(
        "evaluate", parents=[common, symbol_opt],
        help="evaluate all properties at a dependency value")
    p_evaluate.add_argument("yaml_path", type=_existing_path,
                            help="path to the YAML material file")
    p_evaluate.add_argument("value", type=float,
                            help="numeric value of the dependency to evaluate at")
    p_evaluate.set_defaults(func=_cmd_evaluate)

    return parser


def _existing_path(raw: str):
    """argparse type that resolves a path and fails early if it is missing."""
    from pathlib import Path

    path = Path(raw)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"file not found: {raw}")
    return path


def _configure_logging(verbosity: int) -> None:
    """Routes library logs to stderr, quiet by default and louder with -v/-vv."""
    level = {0: logging.CRITICAL, 1: logging.INFO}.get(verbosity, logging.DEBUG)
    package_logger = logging.getLogger("materforge")
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    package_logger.handlers = [handler]
    package_logger.setLevel(level)
    package_logger.propagate = False


# ====================================================================
# ENTRY POINTS
# ====================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    """Runs the MaterForge CLI and returns a process exit code.

    Args:
        argv: Argument list to parse (defaults to ``sys.argv[1:]``).
    Returns:
        0 on success, 1 on a handled error.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging(getattr(args, "verbose", 0))

    if not getattr(args, "func", None):
        parser.print_help()
        return _EXIT_ERROR

    try:
        return args.func(args)
    except Exception as exc:  # noqa: BLE001 - surface any failure as a clean CLI error
        logger.debug("Command '%s' failed", args.command, exc_info=True)
        print(f"Error: {exc}", file=sys.stderr)
        return _EXIT_ERROR


def validate_entry(argv: Optional[Sequence[str]] = None) -> int:
    """Console-script shim for the legacy ``materforge-validate`` command.

    Forwards its arguments to ``materforge validate`` so the standalone entry
    point keeps working alongside the unified CLI.
    """
    raw = list(sys.argv[1:] if argv is None else argv)
    return main(["validate", *raw])


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
