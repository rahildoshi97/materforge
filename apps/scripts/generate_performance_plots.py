#!/usr/bin/env python3
"""
generate_performance_plots.py

Publication-ready performance comparison plots for the Couette flow benchmark
(constant viscosity nu=0.08  vs  temperature-dependent viscosity via MaterForge).

Data sources
------------
Log files are parsed automatically when present.  If a log file is missing or
a value cannot be parsed, the hardcoded FALLBACK_DATA block is used instead.
Edit FALLBACK_DATA to update results from a new benchmark campaign without
rerunning the parser.

Configuration
-------------
All tunable parameters are in the CONFIGURATION block below.  No other part
of the script needs editing for routine use.

Outputs (all written to OUTPUT_DIR)
------------------------------------
perf_mlups.png                   MLUPS bars + individual trial scatter + 95 % CI
perf_wall_time.png               Wall-time bars (same style)
perf_timer_breakdown.png         Stacked absolute time: StreamCollide / Comms / BCs
perf_timer_breakdown_pct.png     Same, normalised to 100 %
perf_overhead_decomposition.png  Predicted vs measured overhead bar chart
perf_op_counts.png               Symbolic operation counts (adds / muls / divs)
perf_bandwidth.png               Algorithmic bandwidth vs hardware peak
perf_summary.png                 Multi-panel summary figure (6 sub-plots)
perf_data.csv                    Parsed per-trial data table
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import stats

matplotlib.use("Agg")


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION  — edit here; nowhere else needs changing
# ══════════════════════════════════════════════════════════════════════════════

# ── I/O paths ─────────────────────────────────────────────────────────────────
_HERE             = Path(__file__).parent              # apps/scripts/
APPS_DIR          = _HERE.parent                       # apps/
LOG_DIR           = APPS_DIR / "logs" / "performance"
LOG_CONST_FILE    = LOG_DIR / "run_perf_const.log"
LOG_TEMPDEP_FILE  = LOG_DIR / "run_perf_tempdep.log"
OUTPUT_DIR        = APPS_DIR / "output" / "plots" / "performance"
DATA_DIR          = APPS_DIR / "output" / "data"
FIGURE_DPI        = 300    # 300 for submission, 150 for draft review

# ── Benchmark labels and colours ──────────────────────────────────────────────
LABEL_CONST  = r"Constant $\nu = 0.08$"
LABEL_TEMPDEP = "Temperature-dependent\n(MaterForge)"
COLOR_CONST  = "#1f77b4"   # matplotlib tab:blue
COLOR_TEMPDEP = "#ff7f0e"  # matplotlib tab:orange
CI_LEVEL      = 0.95       # confidence interval coverage

# ── Hardware specification (Xeon Gold 6326, Ice Lake) ─────────────────────────
HW_PEAK_BW_GB_S   = 204.8  # theoretical peak per socket (DDR4-3200, 8 ch)
HW_CLOCK_GHZ      = 2.90
HW_CORES_USED     = 4      # MPI ranks used in the benchmark

# ── Algorithmic bandwidth model ───────────────────────────────────────────────
BW_CONST_BYTES_CELL   = 336   # 19 PDFs r+w (304) + density w (8) + velocity w (24)
BW_TEMPDEP_BYTES_CELL = 352   # + temperature r (8) + viscosity w (8)

# ── Symbolic operation counts (from cmake configure log) ──────────────────────
OPS = {
    # { label: {adds, muls, divs, total} }
    LABEL_CONST: {
        "adds": 167, "muls": 212, "divs": 1,
        "total": 380,
        "description": "omega baked as compile-time literal",
    },
    LABEL_TEMPDEP: {
        "adds": 238, "muls": 295, "divs": 4,
        "total": 537,
        "description": "runtime piecewise-polynomial nu(T)",
    },
}

# ── Overhead decomposition estimates ──────────────────────────────────────────
# Sources: bandwidth model, perf-stat hardware-counter analysis,
#          waLBerla timer reductions (see couette_performance_analysis.md)
OVERHEAD_COMPONENTS = {
    "Algorithmic\nbandwidth\n(+4.8 % bytes/cell)": 5.0,   # % wall-time contribution
    "Extra L1 traffic\n(+20.9 % loads)":            3.5,
    "Extra instructions\n(+21.3 % retired,\npartly hidden by IPC)": 2.0,
    "Communication\nscheduling\n(net negative)":    -2.0,
}
OVERHEAD_MEASURED_PCT = 8.52   # measured total wall-clock overhead (%)

# ── Fallback data (used when log files are missing or unparseable) ─────────────
# Update this block after a new benchmark campaign.
FALLBACK_DATA = {
    LABEL_CONST: {
        "mlups_total":   [67.1908, 66.5277, 66.7587, 66.7707, 66.5366],
        "wall_time_s":   [468.178, 472.845, 471.209, 471.124, 472.726],
        "timer_total_s": {            # sum across 4 ranks (REDUCE_TOTAL)
            "StreamCollide":    [1648.18, 1665.32, 1660.05, 1659.73, 1665.12],
            "LBM Communication":[181.63,  183.06,  181.99,  181.88,  183.31],
            "UBB":              [31.62,   31.57,   31.61,   31.58,   31.59],
            "NoSlip":           [10.30,   10.41,   10.24,   10.22,   10.38],
        },
    },
    LABEL_TEMPDEP: {
        "mlups_total":   [62.1174, 61.3147, 61.5239, 61.3884, 61.2328],
        "wall_time_s":   [506.417, 513.047, 511.302, 512.431, 513.733],
        "timer_total_s": {
            "StreamCollide":    [1829.14, 1841.26, 1834.71, 1838.52, 1840.85],
            "LBM Communication":[182.71,  183.42,  183.51,  183.79,  183.08],
            "UBB":              [31.64,   31.58,   31.60,   31.61,   31.62],
            "NoSlip":           [10.35,   10.42,   10.38,   10.40,   10.37],
        },
    },
}

# ══════════════════════════════════════════════════════════════════════════════
# MATPLOTLIB STYLE
# ══════════════════════════════════════════════════════════════════════════════

def _apply_style():
    plt.rcParams.update({
        "font.family":         "sans-serif",
        "font.size":           11,
        "axes.titlesize":      12,
        "axes.labelsize":      11,
        "axes.linewidth":      1.2,
        "xtick.major.width":   1.2,
        "ytick.major.width":   1.2,
        "xtick.minor.width":   0.8,
        "ytick.minor.width":   0.8,
        "xtick.major.size":    5,
        "ytick.major.size":    5,
        "xtick.direction":     "in",
        "ytick.direction":     "in",
        "legend.fontsize":     9,
        "legend.framealpha":   0.85,
        "figure.dpi":          100,
        "savefig.dpi":         FIGURE_DPI,
        "savefig.bbox":        "tight",
        "axes.grid":           True,
        "grid.alpha":          0.35,
        "grid.linewidth":      0.7,
    })


# ══════════════════════════════════════════════════════════════════════════════
# LOG PARSER
# ══════════════════════════════════════════════════════════════════════════════

_RE_MLUPS_TOTAL    = re.compile(r'\[0\]\[RESULT\s*\].*?Total MLUPS:\s*([\d.]+)')
_RE_WALL_TIME      = re.compile(r'\[0\]\[RESULT\s*\].*?Total simulation time:\s*([\d.]+)\s+seconds')
# Timer table row: name | pct | total_time | avg | count | min | max | var
# "TOTAL reduction" rows have count >= 100 000; AVG/MAX reduction rows have count == n_ranks (4)
_RE_TIMER_ROW      = re.compile(
    r'\[0\]\s+([\w ]+?)\s*\|\s*[\d.]+%\s*\|\s*([\d.]+)\s*\|\s*[\d.]+\s*\|\s*(\d+)'
)
_TOTAL_COUNT_MIN   = 100_000   # distinguishes REDUCE_TOTAL from REDUCE_AVG/MAX rows


def _parse_log(path: Path) -> dict | None:
    """Parse a SLURM run log, returning per-trial MLUPS, wall time, and timer totals."""
    if not path.exists():
        return None
    text = path.read_text(errors="replace")

    mlups  = [float(x) for x in _RE_MLUPS_TOTAL.findall(text)]
    wall   = [float(x) for x in _RE_WALL_TIME.findall(text)]

    # Collect all TOTAL-reduction timer rows grouped by trial.
    # Strategy: scan lines; when we hit a "Total MLUPS" result line we know a trial
    # ended.  Timer rows between two successive "Total MLUPS" lines belong to the
    # preceding trial.
    timers_per_trial: list[dict[str, float]] = []
    current: dict[str, float] = {}
    for line in text.splitlines():
        m_timer = _RE_TIMER_ROW.search(line)
        if m_timer:
            name, total, count = m_timer.group(1).strip(), float(m_timer.group(2)), int(m_timer.group(3))
            if count >= _TOTAL_COUNT_MIN:
                current[name] = total
        if _RE_MLUPS_TOTAL.search(line):
            if current:
                timers_per_trial.append(current)
                current = {}

    # Invert: { timer_name -> [val_trial_1, ...] }
    timer_totals: dict[str, list[float]] = {}
    for trial_d in timers_per_trial:
        for name, val in trial_d.items():
            timer_totals.setdefault(name, []).append(val)

    if not mlups:
        return None
    return {"mlups_total": mlups, "wall_time_s": wall, "timer_total_s": timer_totals}


def load_data() -> dict[str, dict]:
    """Load data from logs; fall back to FALLBACK_DATA per case if needed."""
    data: dict[str, dict] = {}
    for label, log_path in [
        (LABEL_CONST,   LOG_CONST_FILE),
        (LABEL_TEMPDEP, LOG_TEMPDEP_FILE),
    ]:
        parsed = _parse_log(log_path)
        if parsed and len(parsed["mlups_total"]) >= 2:
            data[label] = parsed
            print(f"  Loaded from log: {log_path.name} "
                  f"({len(parsed['mlups_total'])} trials)")
        else:
            data[label] = FALLBACK_DATA[label]
            status = "log missing" if not log_path.exists() else "parse failed or < 2 trials"
            print(f"  Using fallback data for '{label}' ({status})")
    return data


# ══════════════════════════════════════════════════════════════════════════════
# STATISTICS HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _ci(values: list[float], level: float = CI_LEVEL) -> tuple[float, float, float, float]:
    """Return (mean, std, ci_lo, ci_hi) using t-distribution."""
    arr = np.asarray(values, dtype=float)
    n   = len(arr)
    m   = arr.mean()
    s   = arr.std(ddof=1) if n > 1 else 0.0
    if n > 1:
        t = stats.t.ppf((1 + level) / 2, df=n - 1)
        h = t * s / np.sqrt(n)
    else:
        h = 0.0
    return m, s, m - h, m + h


# ══════════════════════════════════════════════════════════════════════════════
# PLOT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _save(fig: plt.Figure, name: str) -> None:
    path = OUTPUT_DIR / name
    fig.savefig(path)
    plt.close(fig)
    print(f"  Saved: {name}")


def _add_value_label(ax, x, y, text, *, va="bottom", offset=0.01):
    ax.text(x, y + offset, text, ha="center", va=va, fontsize=9, fontweight="bold")


def _bar_pair_with_trials(ax, labels, values_list, colors,
                           ylabel, title, unit="", fmt=".2f"):
    """Two bars with individual trial scatter dots and 95 % CI error bar."""
    x = np.arange(len(labels))
    for i, (label, vals, color) in enumerate(zip(labels, values_list, colors)):
        m, s, lo, hi = _ci(vals)
        err_lo = m - lo
        err_hi = hi - m
        ax.bar(x[i], m, color=color, alpha=0.85, zorder=2, label=label,
               error_kw={"elinewidth": 1.5, "capsize": 5},
               yerr=[[err_lo], [err_hi]])
        # Individual trial dots
        jitter = (np.random.RandomState(42 + i).rand(len(vals)) - 0.5) * 0.12
        ax.scatter(x[i] + jitter, vals, color="black", s=20, zorder=5, alpha=0.7)
        # Annotate mean
        _add_value_label(ax, x[i], m + err_hi, f"{m:{fmt}}{unit}",
                         offset=err_hi * 0.05 + 0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


# ══════════════════════════════════════════════════════════════════════════════
# INDIVIDUAL PLOTS
# ══════════════════════════════════════════════════════════════════════════════

def plot_mlups(data: dict) -> None:
    fig, ax = plt.subplots(figsize=(6, 5.5))
    labels = [LABEL_CONST, LABEL_TEMPDEP]
    colors = [COLOR_CONST, COLOR_TEMPDEP]
    values = [data[l]["mlups_total"] for l in labels]
    _bar_pair_with_trials(ax, labels, values, colors,
                           "Total MLUPS", "Throughput — Total MLUPS", unit=" MLUPS")

    m_c  = np.mean(data[LABEL_CONST]["mlups_total"])
    m_td = np.mean(data[LABEL_TEMPDEP]["mlups_total"])
    overhead_pct = (m_c - m_td) / m_c * 100
    ax.annotate(
        f"Overhead: {overhead_pct:.2f} %\n"
        f"(ratio {m_c/m_td:.4f}×)",
        xy=(0.97, 0.97), xycoords="axes fraction",
        ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
    )
    ci_label = f"{int(CI_LEVEL*100)} % CI, n={len(data[LABEL_CONST]['mlups_total'])} trials"
    ax.text(0.50, 0.03, ci_label, transform=ax.transAxes,
            fontsize=8, va="bottom", ha="center", color="gray")
    # Extra headroom above value labels
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi * 1.12)
    fig.tight_layout()
    _save(fig, "perf_mlups.png")


def plot_wall_time(data: dict) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    labels = [LABEL_CONST, LABEL_TEMPDEP]
    colors = [COLOR_CONST, COLOR_TEMPDEP]
    values = [data[l]["wall_time_s"] for l in labels]
    _bar_pair_with_trials(ax, labels, values, colors,
                           "Wall-clock time (s)",
                           f"Wall-clock Time — {data[LABEL_CONST]['wall_time_s'].__len__()} Trials",
                           unit=" s")

    m_c  = np.mean(data[LABEL_CONST]["wall_time_s"])
    m_td = np.mean(data[LABEL_TEMPDEP]["wall_time_s"])
    overhead_pct = (m_td - m_c) / m_c * 100
    ax.annotate(
        f"Overhead: +{overhead_pct:.2f} %\n"
        f"(ratio {m_td/m_c:.4f}×)",
        xy=(0.97, 0.07), xycoords="axes fraction",
        ha="right", va="bottom", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
    )
    fig.tight_layout()
    _save(fig, "perf_wall_time.png")


def plot_timer_breakdown(data: dict) -> None:
    timer_names = ["StreamCollide", "LBM Communication", "UBB", "NoSlip"]
    timer_labels = ["StreamCollide", "Communication", "UBB wall", "NoSlip wall"]
    timer_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    labels = [LABEL_CONST, LABEL_TEMPDEP]
    x = np.arange(len(labels))
    width = 0.55

    # Mean timer totals (sum across 4 ranks, averaged over trials)
    means: dict[str, list[float]] = {t: [] for t in timer_names}
    errs:  dict[str, list[float]] = {t: [] for t in timer_names}
    for label in labels:
        td = data[label]["timer_total_s"]
        for t in timer_names:
            vals = td.get(t, [0.0])
            m, _, lo, hi = _ci(vals)
            means[t].append(m)
            errs[t].append(hi - m)

    # ── absolute time ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    bottoms = np.zeros(len(labels))
    for t, label_t, color in zip(timer_names, timer_labels, timer_colors):
        vals = np.array(means[t])
        ax.bar(x, vals, width, bottom=bottoms, label=label_t,
               color=color, alpha=0.85, zorder=2)
        # Annotate segment centre
        for j, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 15:
                ax.text(x[j], b + v / 2, f"{v:.0f} s",
                        ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("CPU-seconds (sum across 4 ranks)")
    ax.set_title("Timer Breakdown — Total CPU Time per Component")
    ax.legend(loc="upper right")
    fig.tight_layout()
    _save(fig, "perf_timer_breakdown.png")

    # ── percentage ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    totals = [sum(means[t][j] for t in timer_names) for j in range(len(labels))]
    bottoms = np.zeros(len(labels))
    for t, label_t, color in zip(timer_names, timer_labels, timer_colors):
        pcts = np.array([means[t][j] / totals[j] * 100 for j in range(len(labels))])
        ax.bar(x, pcts, width, bottom=bottoms, label=label_t, color=color, alpha=0.85)
        for j, (p, b) in enumerate(zip(pcts, bottoms)):
            if p > 3:
                ax.text(x[j], b + p / 2, f"{p:.1f} %",
                        ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        bottoms += pcts

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Fraction of total CPU time (%)")
    ax.set_title("Timer Breakdown — Percentage of Total")
    ax.set_ylim(0, 110)
    ax.legend(loc="upper right")
    fig.tight_layout()
    _save(fig, "perf_timer_breakdown_pct.png")


def plot_overhead_decomposition() -> None:
    """Estimated overhead components vs measured total."""
    fig, ax = plt.subplots(figsize=(8, 5.5))

    components = list(OVERHEAD_COMPONENTS.keys())
    values     = list(OVERHEAD_COMPONENTS.values())

    # Short x-tick labels (details already shown in bar annotations and CONFIGURATION)
    short_labels = [
        "Algorithmic\nbandwidth",
        "Extra L1\ntraffic",
        "Extra\ninstructions",
        "Comm.\nscheduling",
    ]

    # Separate positive and negative
    pos_colors = ["#1f77b4", "#aec7e8", "#ffbb78", "#98df8a"]
    neg_color  = "#d62728"
    bar_colors = [neg_color if v < 0 else pos_colors[i % len(pos_colors)]
                  for i, v in enumerate(values)]

    bars = ax.bar(range(len(components)), values, color=bar_colors, alpha=0.85, zorder=2)
    for bar, val in zip(bars, values):
        sign = "+" if val >= 0 else ""
        va   = "bottom" if val >= 0 else "top"
        y_off = 0.18 if val >= 0 else -0.18
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + y_off,
                f"{sign}{val:.1f} %",
                ha="center", va=va,
                fontsize=9, fontweight="bold")

    # Measured total - dashed line with inline label; no legend box
    ax.axhline(OVERHEAD_MEASURED_PCT, color="black", lw=2, ls="--")
    # x in axes fraction, y in data coords: label hugs the left edge of the line
    ax.text(0.02, OVERHEAD_MEASURED_PCT,
            f"Measured: +{OVERHEAD_MEASURED_PCT:.2f} %",
            transform=ax.get_yaxis_transform(),
            ha="left", va="bottom", fontsize=9, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.axhline(0, color="gray", lw=0.8, ls="-")

    ax.set_xticks(range(len(components)))
    ax.set_xticklabels(short_labels, fontsize=9)
    ax.set_ylabel("Estimated wall-clock overhead (%)")
    ax.set_title("Overhead Decomposition: Estimated vs Measured")

    # Symbolic FLOP note outside axes at figure bottom (avoids clashing with bars)
    symbolic_pct = (OPS[LABEL_TEMPDEP]["total"] / OPS[LABEL_CONST]["total"] - 1) * 100
    fig.text(0.50, 0.01,
             f"Symbolic FLOP prediction: +{symbolic_pct:.1f} % (count_operations ratio)",
             ha="center", va="bottom", fontsize=8, color="gray", style="italic")

    # Extra headroom so bar value labels don't clip the top
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo * 1.3, yhi * 1.25)
    # rect reserves 7 % at the bottom of the figure for the fig.text note
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    _save(fig, "perf_overhead_decomposition.png")


def plot_op_counts() -> None:
    """Grouped bar chart of symbolic operation counts."""
    op_types   = ["adds", "muls", "divs", "total"]
    op_labels  = ["Additions", "Multiplications", "Divisions", "Total ops"]
    x = np.arange(len(op_types))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for i, (label, color) in enumerate(
        [(LABEL_CONST, COLOR_CONST), (LABEL_TEMPDEP, COLOR_TEMPDEP)]
    ):
        vals = [OPS[label][k] for k in op_types]
        offset = (i - 0.5) * width
        bars = ax.bar(x + offset, vals, width, label=label, color=color, alpha=0.85)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                    str(val), ha="center", va="bottom", fontsize=8)

    # Delta annotations — placed just above the taller of each pair
    for j, key in enumerate(op_types):
        c = OPS[LABEL_CONST][key]
        t = OPS[LABEL_TEMPDEP][key]
        delta_pct = (t - c) / c * 100
        ax.text(x[j], max(c, t) + 30,
                f"+{delta_pct:.0f} %", ha="center", va="bottom",
                fontsize=8, color="dimgray", style="italic")

    # Explicit y-limit: 28 % headroom above the tallest bar + delta label
    y_top = max(OPS[LABEL_TEMPDEP][k] for k in op_types)
    ax.set_ylim(0, y_top * 1.28)

    ax.set_xticks(x)
    ax.set_xticklabels(op_labels)
    ax.set_ylabel("Symbolic operation count\n(count_operations, pre-compilation)")
    ax.set_title("Symbolic Operation Counts — Code-Generation Metric")
    # Legend above the axes so it never overlaps bars or delta labels
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    # Short note as fig.text outside the axes at the figure bottom
    note = ("count_operations: AST-level diagnostic. "
            "Measured instruction overhead: +21.3 % "
            "(compiler fuses ops via AVX-512 FMA).")
    fig.text(0.50, 0.01, note, ha="center", va="bottom",
             fontsize=7.5, color="gray", style="italic")
    # rect reserves 7 % at bottom for note, 10 % at top for legend
    fig.tight_layout(rect=[0, 0.07, 1, 0.90])
    _save(fig, "perf_op_counts.png")


def plot_bandwidth(data: dict) -> None:
    """Algorithmic bandwidth per case vs hardware peak."""
    m_mlups_c  = np.mean(data[LABEL_CONST]["mlups_total"])
    m_mlups_td = np.mean(data[LABEL_TEMPDEP]["mlups_total"])

    bw_const_gs  = m_mlups_c  * 1e6 * BW_CONST_BYTES_CELL  / 1e9
    bw_tempdep_gs = m_mlups_td * 1e6 * BW_TEMPDEP_BYTES_CELL / 1e9

    fig, ax = plt.subplots(figsize=(6, 5))
    labels = [LABEL_CONST, LABEL_TEMPDEP]
    bws    = [bw_const_gs, bw_tempdep_gs]
    colors = [COLOR_CONST, COLOR_TEMPDEP]

    bars = ax.bar(range(2), bws, color=colors, alpha=0.85, zorder=2)
    for bar, bw, label, b_cell in zip(bars, bws, labels,
                                       [BW_CONST_BYTES_CELL, BW_TEMPDEP_BYTES_CELL]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                f"{bw:.2f} GB/s\n({b_cell} B/cell)",
                ha="center", va="bottom", fontsize=9)

    # Hardware peak reference
    ax.axhline(HW_PEAK_BW_GB_S, color="red", lw=2, ls="--",
               label=f"Hardware peak: {HW_PEAK_BW_GB_S:.0f} GB/s/socket")

    # Utilisation labels
    for i, (bw, label) in enumerate(zip(bws, labels)):
        util = bw / HW_PEAK_BW_GB_S * 100
        ax.text(i, HW_PEAK_BW_GB_S * 0.55,
                f"{util:.1f} %\nof peak",
                ha="center", va="center", fontsize=9, color="dimgray",
                bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", ec="gray", alpha=0.7))

    ax.set_xticks(range(2))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Effective DRAM bandwidth (GB/s)")
    ax.set_title("Algorithmic Bandwidth vs Hardware Peak\n"
                 f"(Xeon Gold 6326, {HW_CORES_USED} cores, DDR4-3200)")
    ax.set_ylim(0, HW_PEAK_BW_GB_S * 1.15)
    ax.legend(loc="upper right")
    fig.tight_layout()
    _save(fig, "perf_bandwidth.png")


def plot_summary(data: dict) -> None:
    """6-panel summary figure combining key results."""
    fig = plt.figure(figsize=(15, 9))
    fig.suptitle(
        "Couette Flow Benchmark — MaterForge Overhead on Xeon Gold 6326 (Ice Lake, woody NHR@FAU)\n"
        "D3Q19 TRT LBM, 128×64×64 cells, 4 MPI ranks, 60 000 timesteps",
        fontsize=12, y=0.99,
    )

    # 2×3 grid
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.38)
    ax_mlups   = fig.add_subplot(gs[0, 0])
    ax_wt      = fig.add_subplot(gs[0, 1])
    ax_bw      = fig.add_subplot(gs[0, 2])
    ax_timer   = fig.add_subplot(gs[1, 0])
    ax_ops     = fig.add_subplot(gs[1, 1])
    ax_overhead = fig.add_subplot(gs[1, 2])

    labels  = [LABEL_CONST, LABEL_TEMPDEP]
    colors  = [COLOR_CONST, COLOR_TEMPDEP]

    # ── Panel 1: MLUPS ────────────────────────────────────────────────────────
    _bar_pair_with_trials(ax_mlups, labels,
                           [data[l]["mlups_total"] for l in labels],
                           colors, "Total MLUPS", "Throughput (MLUPS)")

    # ── Panel 2: Wall time ────────────────────────────────────────────────────
    _bar_pair_with_trials(ax_wt, labels,
                           [data[l]["wall_time_s"] for l in labels],
                           colors, "Wall time (s)", "Wall-clock Time (s)")

    # ── Panel 3: Bandwidth ────────────────────────────────────────────────────
    m_c  = np.mean(data[LABEL_CONST]["mlups_total"])
    m_td = np.mean(data[LABEL_TEMPDEP]["mlups_total"])
    bws = [m_c * 1e6 * BW_CONST_BYTES_CELL / 1e9,
           m_td * 1e6 * BW_TEMPDEP_BYTES_CELL / 1e9]
    bars = ax_bw.bar(range(2), bws, color=colors, alpha=0.85, zorder=2)
    for bar, bw, b_cell in zip(bars, bws, [BW_CONST_BYTES_CELL, BW_TEMPDEP_BYTES_CELL]):
        ax_bw.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                   f"{bw:.1f} GB/s\n{b_cell} B/cell", ha="center", va="bottom", fontsize=8)
    ax_bw.axhline(HW_PEAK_BW_GB_S, color="red", lw=1.5, ls="--",
                  label=f"Peak {HW_PEAK_BW_GB_S:.0f} GB/s")
    ax_bw.set_xticks(range(2))
    ax_bw.set_xticklabels(["Const", "Temp-dep"], fontsize=9)
    ax_bw.set_ylabel("Bandwidth (GB/s)")
    ax_bw.set_title("Algorithmic Bandwidth")
    ax_bw.set_ylim(0, HW_PEAK_BW_GB_S * 1.15)
    ax_bw.legend(fontsize=8)

    # ── Panel 4: Timer breakdown (%) ──────────────────────────────────────────
    timer_names  = ["StreamCollide", "LBM Communication", "UBB", "NoSlip"]
    timer_lbl    = ["StreamCollide", "Comms", "UBB", "NoSlip"]
    timer_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    x = np.arange(len(labels))
    means = {t: [np.mean(data[l]["timer_total_s"].get(t, [0])) for l in labels]
             for t in timer_names}
    totals = [sum(means[t][j] for t in timer_names) for j in range(len(labels))]
    bottoms = np.zeros(len(labels))
    for t, tl, tc in zip(timer_names, timer_lbl, timer_colors):
        pcts = np.array([means[t][j] / totals[j] * 100 for j in range(len(labels))])
        ax_timer.bar(x, pcts, 0.55, bottom=bottoms, label=tl, color=tc, alpha=0.85)
        for j, (p, b) in enumerate(zip(pcts, bottoms)):
            if p > 4:
                ax_timer.text(x[j], b + p / 2, f"{p:.1f}%",
                              ha="center", va="center", fontsize=7.5,
                              color="white", fontweight="bold")
        bottoms += pcts
    ax_timer.set_xticks(x)
    ax_timer.set_xticklabels(["Const", "Temp-dep"], fontsize=9)
    ax_timer.set_ylabel("Fraction of CPU time (%)")
    ax_timer.set_title("Timer Breakdown")
    ax_timer.legend(fontsize=7.5, loc="upper right")

    # ── Panel 5: Symbolic op counts (total only) ──────────────────────────────
    op_vals = {"Const": OPS[LABEL_CONST]["total"], "Temp-dep": OPS[LABEL_TEMPDEP]["total"]}
    bars = ax_ops.bar(range(2), list(op_vals.values()), color=colors, alpha=0.85, zorder=2)
    for bar, (lbl, val) in zip(bars, op_vals.items()):
        ax_ops.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                    str(val), ha="center", va="bottom", fontsize=9, fontweight="bold")
    delta_pct = (OPS[LABEL_TEMPDEP]["total"] / OPS[LABEL_CONST]["total"] - 1) * 100
    ax_ops.text(0.97, 0.97, f"Symbolic +{delta_pct:.1f} %",
                transform=ax_ops.transAxes, ha="right", va="top",
                fontsize=9, style="italic", color="gray")
    ax_ops.set_xticks(range(2))
    ax_ops.set_xticklabels(["Const", "Temp-dep"], fontsize=9)
    ax_ops.set_ylabel("Total symbolic ops")
    ax_ops.set_title("Symbolic Op Count\n(count_operations)")

    # ── Panel 6: Overhead summary bars ───────────────────────────────────────
    items   = list(OVERHEAD_COMPONENTS.items()) + [("Measured\ntotal", OVERHEAD_MEASURED_PCT)]
    names   = [k for k, _ in items]
    vals    = [v for _, v in items]
    bcolors = ["#d62728" if v < 0 else ("#2ca02c" if k == "Measured\ntotal" else "#aec7e8")
               for k, v in items]
    ax_overhead.bar(range(len(names)), vals, color=bcolors, alpha=0.85, zorder=2)
    ax_overhead.axhline(0, color="gray", lw=0.8)
    ax_overhead.axhline(OVERHEAD_MEASURED_PCT, color="#2ca02c", lw=1.5, ls="--", alpha=0.7)
    for j, (name, val) in enumerate(zip(names, vals)):
        sign = "+" if val >= 0 else ""
        ax_overhead.text(j, val + (0.2 if val >= 0 else -0.5),
                         f"{sign}{val:.1f}%",
                         ha="center", va="bottom" if val >= 0 else "top",
                         fontsize=7.5)
    ax_overhead.set_xticks(range(len(names)))
    ax_overhead.set_xticklabels(names, fontsize=7)
    ax_overhead.set_ylabel("Wall-clock overhead (%)")
    ax_overhead.set_title("Overhead Decomposition")

    fig.savefig(OUTPUT_DIR / "perf_summary.png")
    plt.close(fig)
    print("  Saved: perf_summary.png")


# ══════════════════════════════════════════════════════════════════════════════
# DATA EXPORT
# ══════════════════════════════════════════════════════════════════════════════

def export_csv(data: dict) -> None:
    rows = []
    for label, d in data.items():
        for i, mlups in enumerate(d["mlups_total"]):
            row = {
                "Case":       label,
                "Trial":      i + 1,
                "MLUPS_total": mlups,
                "Wall_time_s": d["wall_time_s"][i] if i < len(d["wall_time_s"]) else float("nan"),
            }
            for timer in ["StreamCollide", "LBM Communication", "UBB", "NoSlip"]:
                vals = d["timer_total_s"].get(timer, [])
                row[f"Timer_{timer.replace(' ','_')}"] = vals[i] if i < len(vals) else float("nan")
            rows.append(row)
    df = pd.DataFrame(rows)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / "perf_data.csv"
    df.to_csv(path, index=False, float_format="%.6f")
    print(f"  Saved: {path}  ({len(df)} rows)")

    # Print summary statistics
    print("\n  Summary statistics:")
    for label in [LABEL_CONST, LABEL_TEMPDEP]:
        sub = df[df["Case"] == label]
        m, _, lo, hi = _ci(sub["MLUPS_total"].tolist())
        print(f"    {label[:30]:30s}  MLUPS: {m:.4f}  95%CI [{lo:.4f}, {hi:.4f}]")

    # Welch's t-test
    g1 = df[df["Case"] == LABEL_CONST]["MLUPS_total"].values
    g2 = df[df["Case"] == LABEL_TEMPDEP]["MLUPS_total"].values
    t_stat, p_val = stats.ttest_ind(g1, g2, equal_var=False)
    overhead = (g1.mean() - g2.mean()) / g1.mean() * 100
    print(f"\n  Welch t-test:  t = {t_stat:.4f},  p = {p_val:.3e}")
    print(f"  MLUPS overhead:  {overhead:.4f} %  "
          f"(const {g1.mean():.4f} vs tempdep {g2.mean():.4f})")
    print(f"  Wall-time overhead:  "
          f"{(np.mean(df[df['Case']==LABEL_TEMPDEP]['Wall_time_s'].values) / np.mean(df[df['Case']==LABEL_CONST]['Wall_time_s'].values) - 1)*100:.4f} %")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    _apply_style()
    np.random.seed(42)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("COUETTE FLOW — PERFORMANCE PLOTS")
    print("=" * 70)

    print("\nLoading data...")
    data = load_data()

    print("\nGenerating plots:")
    plot_mlups(data)
    plot_wall_time(data)
    plot_timer_breakdown(data)
    plot_overhead_decomposition()
    plot_op_counts()
    plot_bandwidth(data)
    plot_summary(data)

    print("\nExporting data table:")
    export_csv(data)

    print(f"\nAll outputs written to: {OUTPUT_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()
