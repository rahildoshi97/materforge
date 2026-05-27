#!/usr/bin/env python3
"""Strong-scaling plots for the Couette benchmark.

Reads apps/scaling_data.csv produced by parse_scaling.py and produces:
  perf_scaling_mlups.png       Total MLUPS vs MPI ranks (log-log)
  perf_scaling_speedup.png     Speedup vs ideal
  perf_scaling_efficiency.png  Parallel efficiency
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

_HERE = Path(__file__).parent
CSV   = _HERE / "scaling_data.csv"

LABEL_CONST   = "const_0.08"
LABEL_TEMPDEP = "tempdep"
DISP_CONST    = r"Constant $\nu = 0.08$"
DISP_TEMPDEP  = "Temperature-dependent (MaterForge)"
COLOR_CONST   = "#1f77b4"
COLOR_TEMPDEP = "#ff7f0e"
CI_LEVEL      = 0.95
FIG_DPI       = 300


def _apply_style() -> None:
    plt.rcParams.update({
        "font.size":         11,
        "axes.titlesize":    12,
        "axes.labelsize":    11,
        "axes.linewidth":    1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.direction":   "in",
        "ytick.direction":   "in",
        "legend.fontsize":   9,
        "savefig.dpi":       FIG_DPI,
        "savefig.bbox":      "tight",
        "axes.grid":         True,
        "grid.alpha":        0.35,
    })


def _ci(values: np.ndarray, level: float = CI_LEVEL) -> tuple[float, float]:
    n = len(values)
    if n < 2:
        return float(values.mean()), 0.0
    m = values.mean()
    s = values.std(ddof=1)
    t = stats.t.ppf((1 + level) / 2, df=n - 1)
    return m, t * s / np.sqrt(n)


def _aggregate(df: pd.DataFrame, value_col: str) -> dict[str, pd.DataFrame]:
    """Return per-case dataframe with columns ranks, mean, halfwidth."""
    out: dict[str, pd.DataFrame] = {}
    for case, sub in df.groupby("case"):
        rows = []
        for ranks, g in sub.groupby("ranks"):
            m, h = _ci(g[value_col].to_numpy())
            rows.append({"ranks": ranks, "mean": m, "halfwidth": h})
        out[case] = pd.DataFrame(rows).sort_values("ranks").reset_index(drop=True)
    return out


def plot_mlups(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    agg = _aggregate(df, "mlups_total")
    ranks_all = sorted(df["ranks"].unique())

    for case, color, disp in [
        (LABEL_CONST,   COLOR_CONST,   DISP_CONST),
        (LABEL_TEMPDEP, COLOR_TEMPDEP, DISP_TEMPDEP),
    ]:
        if case not in agg:
            continue
        a = agg[case]
        ax.errorbar(a["ranks"], a["mean"], yerr=a["halfwidth"],
                    marker="o", lw=1.6, ms=7, capsize=4,
                    color=color, label=disp, zorder=3)
        for _, row in a.iterrows():
            ax.annotate(f"{row['mean']:.1f}",
                        xy=(row["ranks"], row["mean"]),
                        xytext=(6, 6), textcoords="offset points",
                        fontsize=8, color=color)

    # Ideal-scaling reference: extend the lowest-rank const measurement linearly
    if LABEL_CONST in agg and len(agg[LABEL_CONST]) > 0:
        base = agg[LABEL_CONST].iloc[0]
        ideal_ranks = np.array(ranks_all, dtype=float)
        ideal_mlups = base["mean"] * ideal_ranks / base["ranks"]
        ax.plot(ideal_ranks, ideal_mlups, "--", color="gray", lw=1.2,
                label="Ideal (linear)", zorder=2)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(ranks_all)
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("MPI ranks (one Xeon Gold 6326 node)")
    ax.set_ylabel("Total MLUPS")
    ax.set_title("Strong Scaling - Total MLUPS\n"
                 r"Fixed 128$\times$64$\times$64 domain, 10 000 timesteps")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(_HERE / "perf_scaling_mlups.png")
    plt.close(fig)
    print("  Saved: perf_scaling_mlups.png")


def plot_speedup(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    agg = _aggregate(df, "mlups_total")
    ranks_all = sorted(df["ranks"].unique())

    for case, color, disp in [
        (LABEL_CONST,   COLOR_CONST,   DISP_CONST),
        (LABEL_TEMPDEP, COLOR_TEMPDEP, DISP_TEMPDEP),
    ]:
        if case not in agg:
            continue
        a = agg[case]
        base = a.iloc[0]
        speedup = a["mean"].to_numpy() / base["mean"] * base["ranks"]
        # Half-width propagation: linear in mean
        half = a["halfwidth"].to_numpy() / base["mean"] * base["ranks"]
        ax.errorbar(a["ranks"], speedup, yerr=half,
                    marker="o", lw=1.6, ms=7, capsize=4,
                    color=color, label=disp, zorder=3)

    ax.plot(ranks_all, ranks_all, "--", color="gray", lw=1.2, label="Ideal", zorder=2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(ranks_all)
    ax.set_yticks(ranks_all)
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Speedup vs. 1 rank")
    ax.set_title("Strong Scaling - Speedup")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(_HERE / "perf_scaling_speedup.png")
    plt.close(fig)
    print("  Saved: perf_scaling_speedup.png")


def plot_efficiency(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    agg = _aggregate(df, "mlups_total")
    ranks_all = sorted(df["ranks"].unique())

    for case, color, disp in [
        (LABEL_CONST,   COLOR_CONST,   DISP_CONST),
        (LABEL_TEMPDEP, COLOR_TEMPDEP, DISP_TEMPDEP),
    ]:
        if case not in agg:
            continue
        a = agg[case]
        base = a.iloc[0]
        efficiency = (a["mean"].to_numpy() / a["ranks"].to_numpy()) / (base["mean"] / base["ranks"]) * 100
        half = (a["halfwidth"].to_numpy() / a["ranks"].to_numpy()) / (base["mean"] / base["ranks"]) * 100
        ax.errorbar(a["ranks"], efficiency, yerr=half,
                    marker="o", lw=1.6, ms=7, capsize=4,
                    color=color, label=disp, zorder=3)
        for r, e in zip(a["ranks"], efficiency):
            ax.annotate(f"{e:.0f}%",
                        xy=(r, e), xytext=(6, -12), textcoords="offset points",
                        fontsize=8, color=color)

    ax.axhline(100, color="gray", ls="--", lw=1.2, label="Ideal (100 %)")
    ax.set_xscale("log", base=2)
    ax.set_xticks(ranks_all)
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Parallel efficiency (%)")
    ax.set_title("Strong Scaling - Parallel Efficiency")
    ax.set_ylim(0, 115)
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(_HERE / "perf_scaling_efficiency.png")
    plt.close(fig)
    print("  Saved: perf_scaling_efficiency.png")


def main() -> None:
    if not CSV.exists():
        raise SystemExit(f"Missing {CSV}; run parse_scaling.py first.")
    _apply_style()
    df = pd.read_csv(CSV)
    print(f"Loaded {len(df)} rows from {CSV}")
    plot_mlups(df)
    plot_speedup(df)
    plot_efficiency(df)


if __name__ == "__main__":
    main()
