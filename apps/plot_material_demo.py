#!/usr/bin/env python3
"""Standalone figure: temperature-dependent properties of stainless steel 1.4301
as produced by MaterForge from the shipped data/materials/1.4301.yaml file.

This figure illustrates MaterForge's capability to convert a single declarative
YAML file (mixing tabular data, file imports, piecewise equations, and computed
properties) into smooth symbolic expressions usable downstream.

The Couette LBM benchmark uses a synthetic LBM-scaled viscosity expression for
numerical stability; this figure makes the framework's behaviour on a real
engineering alloy independently visible.
"""
from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from materforge import create_material

# Silence the noisy pwlf RuntimeWarnings - they don't affect the fitted output.
warnings.filterwarnings("ignore", category=RuntimeWarning)

_HERE       = Path(__file__).parent
REPO_ROOT   = _HERE.parent
MATERIAL    = REPO_ROOT / "src/materforge/data/materials/1.4301.yaml"
OUT_PATH    = _HERE / "material_1.4301_properties.png"
T_MIN, T_MAX, NT = 300.0, 3000.0, 271   # 10 K resolution

PANELS = [
    # (property name, y-axis label, title, optional unit-conversion factor + label)
    ("density",                       r"$\rho$ (kg/m$^3$)",              "Density",                  1.0,   ""),
    ("heat_capacity",                 r"$c_p$ (J/(kg$\cdot$K))",         "Specific heat capacity",  1.0,   ""),
    ("heat_conductivity",             r"$k$ (W/(m$\cdot$K))",            "Thermal conductivity",    1.0,   ""),
    ("thermal_expansion_coefficient", r"$\alpha$ ($10^{-5}$ 1/K)",       "Thermal expansion coef.", 1e5,   ""),
    ("thermal_diffusivity",           r"$\alpha_T$ ($10^{-6}$ m$^2$/s)", "Thermal diffusivity",     1e6,   ""),
    ("specific_enthalpy",             r"$h$ (MJ/kg)",                    "Specific enthalpy",       1e-6,  ""),
    ("energy_density",                r"$\rho h$ (GJ/m$^3$)",            "Energy density",          1e-9,  ""),
    ("latent_heat_of_fusion",         r"$L_f$ (kJ/kg)",                  "Latent heat of fusion",   1e-3,  ""),
]


def main() -> None:
    if not MATERIAL.exists():
        raise SystemExit(f"Missing material YAML: {MATERIAL}")

    T = sp.Symbol("T")
    print(f"Loading {MATERIAL.name}...")
    mat = create_material(MATERIAL, dependency=T, enable_plotting=False)
    print("  done.")

    Ts = np.linspace(T_MIN, T_MAX, NT)
    samples: dict[str, np.ndarray] = {}
    for name, *_ in PANELS:
        expr = getattr(mat, name)
        try:
            f = sp.lambdify(T, expr, "numpy")
            y = np.asarray(f(Ts), dtype=float)
            if y.ndim == 0:  # constant property
                y = np.full_like(Ts, float(y))
        except Exception as e:
            print(f"  skip {name}: {e}")
            continue
        samples[name] = y

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5))
    fig.suptitle("MaterForge processing of stainless steel 1.4301 (AISI 304)\n"
                 "Properties derived from data/materials/1.4301.yaml",
                 fontsize=13, y=1.00)

    for ax, (name, ylabel, title, scale, _) in zip(axes.flat, PANELS):
        if name not in samples:
            ax.text(0.5, 0.5, f"({name} unavailable)",
                    transform=ax.transAxes, ha="center", va="center", color="gray")
            ax.set_title(title, fontsize=11)
            continue
        ax.plot(Ts, samples[name] * scale, color="#1f3c88", lw=1.8)
        ax.set_xlim(T_MIN, T_MAX)
        ax.set_xlabel("T (K)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.35)
        ax.tick_params(direction="in")

        # Mark solidus/liquidus from the YAML
        try:
            ts = float(mat.solidus_temperature)
            tl = float(mat.liquidus_temperature)
            ax.axvspan(ts, tl, color="#ffcc66", alpha=0.30, zorder=0)
        except (AttributeError, TypeError):
            pass

    # Discrete legend strip across the bottom
    fig.text(0.5, 0.02,
             "Shaded band: solidus-liquidus interval (1605 - 1735 K)",
             ha="center", va="bottom", fontsize=9, color="dimgray", style="italic")

    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {OUT_PATH}")


if __name__ == "__main__":
    main()
