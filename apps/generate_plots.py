#!/usr/bin/env python3
"""
generate_plots.py

Self-contained publication workflow for thermal Couette flow validation.

Expects CSV files produced by extract_vtk_profiles.py, named:
    cf_cpu_mf{const_NNN|tempdep}_dat_<last_index>_<last_step>.csv

Missing CSV files are skipped gracefully.

Notes on the CSV format
-----------------------
ParaView's "Plot Over Line" filter probes at z = 0, 1, …, N (N+1 points) for
an N-cell domain.  The last sample point at z = N falls on the top-wall boundary
of cell N-1, so ParaView clamps it to that cell and returns a duplicate of the
last real cell's values.  That spurious row is dropped during loading.

ParaView exports VTK *cell* data (not interpolated point data): probing at
integer z = i returns cell i's stored value, whose physical centre is at
z = i + 0.5.  The code accounts for this with y_norm[i] = (i + 0.5) / N.

The simulation initialises temperature at cell centres:
    T(cell i) = T_bottom + (T_top - T_bottom) * (i + 0.5) / N
so the bottom cell reads ~321 K (not 300 K) and the top cell ~2979 K (not 3000 K).
This is correct cell-centred behaviour - the walls impose 300 K and 3000 K.

Alternative ParaView export (cell-centre format)
-------------------------------------------------
To probe exactly at cell centres without the spurious-point workaround, set
Plot Over Line as follows (for a 64-cell domain):
    Point1 z = 0.5,  Point2 z = 63.5,  Resolution = 63  →  64 points
Arc-length then runs 0, 1, …, 63 and actual z = arc_length + 0.5.
In load_csv_and_extract_steady_state, replace the drop + y_norm lines with:
    # (no drop needed - line ends at last cell centre, no spurious point)
    # y_norm = (ss["arc_length"].values + 0.5) / (ss["arc_length"].max() + 1)
Numerically equivalent to the current approach; included here for reference.
"""

import re
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d


# ── Configurable ───────────────────────────────────────────────────────────────
# Plot velocity profiles at every TIME_EVERY_N-th CSV entry.
# The CSV stores one entry per vtkWriteFrequency simulation steps, so the
# actual simulation-step interval is TIME_EVERY_N * vtkWriteFrequency.
# Example: TIME_EVERY_N = 10, vtkWriteFrequency = 100  →  every 1000 sim steps.
TIME_EVERY_N: int = 20

# Maximum CSV entry index (inclusive) to use across all plots and steady-state
# selection.  Set to None to use all available timesteps.
# Example: TIME_MAX_INDEX = 99 restricts a 200-entry CSV to the first 100 entries
# (indices 0–99), equivalent to the first half of the simulation.
TIME_MAX_INDEX: int | None = None


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _apply_log_yticks(ax):
    """Apply fine log-scale ticks (2–9 × 10^n minor, 10^n major) to ax.yaxis."""
    ax.yaxis.set_major_locator(plt.LogLocator(base=10, numticks=15))
    ax.yaxis.set_minor_locator(plt.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterSciNotation(base=10, labelOnlyBase=True))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.grid(True, which="major", alpha=0.4)
    ax.grid(True, which="minor", alpha=0.15)


# ──────────────────────────────────────────────────────────────────────────────
# Plotter
# ──────────────────────────────────────────────────────────────────────────────

class CouetteFlowPlotter:
    """Accumulates per-case data and generates publication-quality plots."""

    _COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    def __init__(self, dpi: int = 300):
        self.dpi = dpi
        self.cases: dict = {}

    def add_case(self, name, y_norm, u_sim, u_analytical, T_profile, nu_profile):
        u_wall = float(np.asarray(u_analytical).max())
        error_abs = np.asarray(u_sim) - np.asarray(u_analytical)
        self.cases[name] = {
            "y_norm":       np.asarray(y_norm,       dtype=float),
            "u_sim":        np.asarray(u_sim,         dtype=float),
            "u_analytical": np.asarray(u_analytical,  dtype=float),
            "T_profile":    np.asarray(T_profile,     dtype=float),
            "nu_profile":   np.asarray(nu_profile,    dtype=float),
            "error_abs":    error_abs,
            # Relative error normalised by u_wall so it stays bounded near the
            # no-slip wall where u_analytical → 0.
            "error_rel":    np.abs(error_abs) / u_wall * 100.0,
            "u_wall":       u_wall,
        }

    # ── internal helpers ──────────────────────────────────────────────────────

    def _save(self, fig, fname: str):
        fig.savefig(fname, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {fname}")

    # ── per-case plots ────────────────────────────────────────────────────────

    def plot_single_case_detailed(self, case_name: str, prefix: str):
        """4-panel: temperature · viscosity · velocity profile · absolute error."""
        d = self.cases[case_name]
        y = d["y_norm"]

        fig, axes = plt.subplots(1, 4, figsize=(17, 4.8))
        fig.suptitle(f"Couette Flow - {case_name}", fontsize=12)

        axes[0].plot(d["T_profile"],    y, "k-", lw=1.5)
        axes[0].set(xlabel="Temperature (K)",          ylabel="y / H", title="Temperature", ylim=(0, 1))

        axes[1].plot(d["nu_profile"],   y, "b-", lw=1.5)
        axes[1].set(xlabel=r"Kinematic viscosity $\nu$", ylabel="y / H", title="Viscosity",    ylim=(0, 1))

        axes[2].plot(d["u_analytical"], y, "k-",  lw=1.5,       label="Analytical")
        axes[2].plot(d["u_sim"],        y, "o", color="#d62728", markersize=3, alpha=0.7, label="LBM")
        axes[2].set(xlabel=r"Velocity $u_x$",          ylabel="y / H", title="Velocity Profile", ylim=(0, 1))
        axes[2].legend(fontsize=9)

        err = np.abs(d["error_abs"])
        axes[3].plot(err, y, color="#2ca02c", lw=1.5)
        axes[3].set(xlabel=r"$|u_\mathrm{LBM}-u_\mathrm{ana}|$", ylabel="y / H", title="Absolute Error", ylim=(0, 1))

        fig.tight_layout()
        self._save(fig, f"{prefix}_detailed.png")

    def plot_single_case_errors(self, case_name: str, prefix: str):
        """2-panel: absolute and relative velocity errors."""
        d = self.cases[case_name]
        y = d["y_norm"]

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        fig.suptitle(f"Error Analysis - {case_name}", fontsize=12)

        err_abs = np.abs(d["error_abs"])
        axes[0].plot(err_abs, y, "b-", lw=1.5)
        axes[0].set(xlabel=r"$|u_\mathrm{LBM}-u_\mathrm{ana}|$", ylabel="y / H",
                    title="Absolute Velocity Error", ylim=(0, 1))

        axes[1].plot(d["error_rel"], y, "r-", lw=1.5)
        axes[1].set(xlabel=r"$|u_\mathrm{LBM}-u_\mathrm{ana}|\,/\,u_\mathrm{wall}$ (%)",
                    ylabel="y / H", title=r"Relative Error (% of $u_\mathrm{wall}$)", ylim=(0, 1))

        fig.tight_layout()
        self._save(fig, f"{prefix}_errors.png")

    # ── multi-case plots ──────────────────────────────────────────────────────

    def plot_all_cases_comparison(self, filename: str):
        """Velocity profiles + absolute errors for all cases on one figure.

        The two constant-ν cases share the same analytical solution (linear
        Couette profile, beta=0).  Their analytical lines are identical and
        would overlap.  This is handled by suppressing duplicate legend entries:
        only the first occurrence of each unique analytical curve is labelled.
        """
        if not self.cases:
            return
        fig, (ax_v, ax_e) = plt.subplots(1, 2, figsize=(13, 5))

        seen_analytical: list[np.ndarray] = []
        for i, (name, d) in enumerate(self.cases.items()):
            c = self._COLORS[i % len(self._COLORS)]
            y = d["y_norm"]

            # Check if this analytical curve is a duplicate of one already plotted.
            is_dup = any(np.allclose(d["u_analytical"], prev, rtol=1e-6) for prev in seen_analytical)
            if is_dup:
                # Plot the line (so its colour/style is visible) but hide from legend.
                ax_v.plot(d["u_analytical"], y, "-", color=c, lw=1.5, label="_nolegend_")
                lbm_label = f"{name} (LBM) *"
            else:
                ax_v.plot(d["u_analytical"], y, "-", color=c, lw=1.5, label=f"{name} (ana)")
                lbm_label = f"{name} (LBM)"
                seen_analytical.append(d["u_analytical"])

            ax_v.plot(d["u_sim"], y, "--", color=c, lw=1.0, label=lbm_label)
            ax_e.plot(np.abs(d["error_abs"]), y, "-", color=c, lw=1.5, label=name)

        has_duplicate = any("*" in lbl for lbl in [h.get_label() for h in ax_v.get_lines()])

        ax_v.set(xlabel=r"Velocity $u_x$", ylabel="y / H",
                 title="Velocity Profiles - All Cases", ylim=(0, 1))
        ax_v.legend(fontsize=8)
        ax_e.set(xlabel=r"$|u_\mathrm{LBM}-u_\mathrm{ana}|$", ylabel="y / H",
                 title="Absolute Errors - All Cases", ylim=(0, 1))
        ax_e.legend(fontsize=8)

        # Footnote below both subplots, clear of the x-axis labels.
        if has_duplicate:
            fig.tight_layout()
            fig.subplots_adjust(bottom=0.14)
            fig.text(0.01, 0.02,
                     "* Analytical solutions for all constant-ν cases are identical "
                     "(linear Couette profile - viscosity cancels in steady state)",
                     fontsize=7, color="grey", style="italic")
        else:
            fig.tight_layout()
        self._save(fig, filename)

    def plot_all_cases_direct_comparison(self, filename: str):
        """One sub-panel per case: analytical vs LBM velocity."""
        if not self.cases:
            return
        n = len(self.cases)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.8), sharey=True)
        if n == 1:
            axes = [axes]
        for ax, (name, d) in zip(axes, self.cases.items()):
            y = d["y_norm"]
            ax.plot(d["u_analytical"], y, "k-",  lw=1.5,       label="Analytical")
            ax.plot(d["u_sim"],        y, "o", color="#d62728", markersize=3, alpha=0.7, label="LBM")
            ax.set(xlabel=r"Velocity $u_x$", title=name, ylim=(0, 1))
            ax.legend(fontsize=9)
        axes[0].set_ylabel("y / H")
        fig.suptitle("LBM vs Analytical - Direct Comparison", fontsize=12)
        fig.tight_layout()
        self._save(fig, filename)

    def plot_error_metrics_summary(self, filename: str):
        """Bar chart of L² / L∞ errors and max relative error per case."""
        if not self.cases:
            return
        names = list(self.cases.keys())
        l2   = [float(np.sqrt(np.mean(d["error_abs"] ** 2))) for d in self.cases.values()]
        linf = [float(np.max(np.abs(d["error_abs"])))         for d in self.cases.values()]
        mrel = [float(d["error_rel"].max())                    for d in self.cases.values()]
        x, w = np.arange(len(names)), 0.35

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
        ax1.bar(x - w/2, l2,   w, label="L² error",  color="#1f77b4")
        ax1.bar(x + w/2, linf, w, label="L∞ error", color="#ff7f0e")
        ax1.set_yscale("log")
        ax1.set_xticks(x); ax1.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax1.set(ylabel="Velocity error", title="L² and L∞ Velocity Errors")
        ax1.legend()

        ax2.bar(x, mrel, color="#2ca02c")
        ax2.set_xticks(x); ax2.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax2.set(ylabel=r"Max relative error (% of $u_\mathrm{wall}$)",
                title="Maximum Relative Error")

        fig.tight_layout()
        self._save(fig, filename)

    def export_summary_table(self, filename: str) -> pd.DataFrame:
        """Write per-case error metrics to CSV and return the DataFrame."""
        rows = [
            {
                "Case":             name,
                "N_cells":          int(len(d["y_norm"])),
                "L2_error":         float(np.sqrt(np.mean(d["error_abs"] ** 2))),
                "Linf_error":       float(np.max(np.abs(d["error_abs"]))),
                "Max_rel_error_%":  float(d["error_rel"].max()),
                "Mean_rel_error_%": float(d["error_rel"].mean()),
            }
            for name, d in self.cases.items()
        ]
        df = pd.DataFrame(rows)
        df.to_csv(filename, index=False, float_format="%.6e")
        print(f"  Saved: {filename}")
        return df

    # ── time-evolution plots ──────────────────────────────────────────────────

    def plot_time_evolution(self, case_name: str, time_series: dict,
                            analytical_func, prefix: str):
        """2-panel: velocity profiles coloured by time + L² convergence curve.

        time_series : output of load_all_timesteps()
        analytical_func : callable y_norm → u_analytical (from make_analytical_func)
        """
        if not time_series:
            return

        sorted_keys = sorted(time_series)
        sim_times   = [time_series[k]["sim_time"] for k in sorted_keys]
        t_min, t_max = sim_times[0], sim_times[-1]

        cmap = plt.cm.plasma  # dark-purple (early) → bright-yellow (late)
        norm = plt.Normalize(vmin=t_min, vmax=t_max)

        fig, (ax_vel, ax_err) = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f"Time Evolution - {case_name}", fontsize=12)

        l2_errors = []
        for k in sorted_keys:
            d       = time_series[k]
            y       = d["y_norm"]
            u_sim   = d["u_sim"]
            u_ana   = analytical_func(y)
            color   = cmap(norm(d["sim_time"]))
            alpha   = 0.4 + 0.6 * (d["sim_time"] - t_min) / max(t_max - t_min, 1)
            ax_vel.plot(u_sim, y, color=color, lw=0.9, alpha=alpha)
            l2_errors.append(float(np.sqrt(np.mean((u_sim - u_ana) ** 2))))

        # Overlay the analytical solution as a reference
        y_ref = time_series[sorted_keys[-1]]["y_norm"]
        ax_vel.plot(analytical_func(y_ref), y_ref, "k--", lw=2.0, label="Analytical")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax_vel, label="Simulation timestep")

        ax_vel.set(xlabel=r"Velocity $u_x$", ylabel="y / H",
                   title="Velocity profiles (dark=early, bright=late)", ylim=(0, 1))
        ax_vel.legend(fontsize=9)

        ax_err.semilogy(sim_times, l2_errors, "o-", color="#1f77b4",
                        markersize=4, lw=1.5)
        ax_err.set(xlabel="Simulation timestep", ylabel=r"$L^2$ error",
                   title=r"Convergence: $L^2(u_\mathrm{LBM} - u_\mathrm{ana})$")
        _apply_log_yticks(ax_err)

        fig.tight_layout()
        self._save(fig, f"{prefix}_time_evolution.png")

    def plot_convergence_comparison(self, all_time_series: dict,
                                    all_analytical_funcs: dict, filename: str):
        """L² error vs simulation time for all cases on a single log-scale plot.

        all_time_series     : {case_name: output of load_all_timesteps()}
        all_analytical_funcs: {case_name: callable from make_analytical_func()}
        """
        if not all_time_series:
            return

        fig, ax = plt.subplots(figsize=(8, 5))

        for i, (name, ts) in enumerate(all_time_series.items()):
            if not ts:
                continue
            ana_func    = all_analytical_funcs[name]
            sorted_keys = sorted(ts)
            sim_times   = [ts[k]["sim_time"] for k in sorted_keys]
            l2_errors   = [float(np.sqrt(np.mean(
                               (ts[k]["u_sim"] - ana_func(ts[k]["y_norm"])) ** 2)))
                           for k in sorted_keys]
            c = self._COLORS[i % len(self._COLORS)]
            ax.semilogy(sim_times, l2_errors, "o-", color=c,
                        markersize=4, lw=1.5, label=name)

        ax.set(xlabel="Simulation timestep", ylabel=r"$L^2$ error",
               title=r"Convergence comparison - $L^2(u_\mathrm{LBM} - u_\mathrm{ana})$")
        ax.legend(fontsize=9)
        _apply_log_yticks(ax)
        fig.tight_layout()
        self._save(fig, filename)


# ──────────────────────────────────────────────────────────────────────────────
# Analytical solver
# ──────────────────────────────────────────────────────────────────────────────

class CouetteAnalyticalSolver:
    """Analytical solution for Couette flow with exponential viscosity ν(T)."""

    @staticmethod
    def velocity_profile_exponential(y_norm, T_bottom, T_top, u_wall, nu_0, beta):
        """
        u(y) = u_wall * F(y) / F(1),  F(y) = ∫₀ʸ dy' / ν(T(y'))

        y_norm : cell-centre positions in (0, 1)
        beta=0 reduces to the linear Couette profile.
        """
        y_norm = np.asarray(y_norm)
        y_hi = np.linspace(0.0, 1.0, 2000)
        T_hi = T_bottom + (T_top - T_bottom) * y_hi
        nu_hi = nu_0 * np.exp(beta * (T_hi - 300.0))
        integral = cumulative_trapezoid(1.0 / nu_hi, y_hi, initial=0.0)
        integral_norm = integral / integral[-1]
        f = interp1d(y_hi, integral_norm, kind="cubic",
                     bounds_error=False, fill_value=(0.0, 1.0))
        return u_wall * f(np.clip(y_norm, 0.0, 1.0))

    @staticmethod
    def temperature_profile_linear(y_norm, T_bottom, T_top):
        return T_bottom + (T_top - T_bottom) * np.asarray(y_norm)

    @staticmethod
    def viscosity_profile_exponential(T, nu_0, beta):
        return nu_0 * np.exp(beta * (np.asarray(T) - 300.0))

    @staticmethod
    def velocity_profile_from_nu(y_norm, nu_sim, u_wall):
        """Analytical Couette velocity by integrating the actual 1/nu(y) from the simulation.

        Avoids the wavy error that arises when the simulation uses a piecewise
        polynomial ν(T) (e.g. from MaterForge) that deviates from the idealised
        exponential by up to ±2e-4 — integrating that deviation produces visible
        oscillations in u_analytical when the exponential formula is used instead.
        """
        y_norm = np.asarray(y_norm)
        nu_f = interp1d(y_norm, nu_sim, kind="cubic", fill_value="extrapolate")
        y_hi = np.linspace(0.0, 1.0, 2000)
        nu_hi = np.maximum(nu_f(y_hi), 1e-12)
        integral = cumulative_trapezoid(1.0 / nu_hi, y_hi, initial=0.0)
        f = interp1d(y_hi, integral / integral[-1], kind="cubic",
                     bounds_error=False, fill_value=(0.0, 1.0))
        return u_wall * f(np.clip(y_norm, 0.0, 1.0))


def make_analytical_func(T_bottom, T_top, u_wall, nu_0, beta, nu_sim=None):
    """Return a callable y_norm → u_analytical for use in time-evolution plots.

    When beta != 0 and nu_sim is provided, integrates the actual simulation
    viscosity profile instead of the idealised exponential formula.
    """
    if beta != 0.0 and nu_sim is not None:
        def _f(y_norm):
            return CouetteAnalyticalSolver.velocity_profile_from_nu(y_norm, nu_sim, u_wall)
    else:
        def _f(y_norm):
            return CouetteAnalyticalSolver.velocity_profile_exponential(
                y_norm, T_bottom, T_top, u_wall, nu_0, beta)
    return _f


# ──────────────────────────────────────────────────────────────────────────────
# CSV loading
# ──────────────────────────────────────────────────────────────────────────────

def load_csv_and_extract_steady_state(csv_file, max_index: int | None = TIME_MAX_INDEX):
    """
    Return (y_norm, u_sim, T_sim) at steady state from a ParaView CSV export.

    ParaView's Plot Over Line samples at z = 0, 1, …, N for an N-cell domain,
    producing N+1 rows.  The last row (z = N) is outside the domain and
    duplicates cell N-1's values - it is dropped here.

    ParaView returns VTK cell data (nearest cell, not interpolated), so probing
    at integer z = i gives cell i's stored value.  Cell i's physical centre is
    at z = i + 0.5, so y_norm[i] = (i + 0.5) / N correctly places each sample.

    max_index : if set, only timesteps with index <= max_index are considered,
                so the "steady state" is the last entry up to that index.
    """
    df = pd.read_csv(csv_file)
    all_ts = sorted(df["TimeStep"].unique())
    if max_index is not None:
        all_ts = [t for t in all_ts if t <= max_index]
    last_timestep = all_ts[-1]
    ss = df[df["TimeStep"] == last_timestep].sort_values("arc_length")

    # Drop the spurious duplicate at arc_length == N (outside domain)
    ss = ss[ss["arc_length"] < ss["arc_length"].max()]

    # --- Alternative for cell-centre ParaView export (start=0.5, end=N-0.5, res=N-1) ---
    # No drop needed; arc_length runs 0, 1, …, N-1 and actual z = arc_length + 0.5.
    # y_norm = (ss["arc_length"].values + 0.5) / (ss["arc_length"].max() + 1)
    # ---------------------------------------------------------------------------------

    u_sim  = ss["velocity:0"].values
    T_sim  = ss["temperature"].values
    nu_sim = ss["viscosity"].values
    n = len(u_sim)
    y_norm = (np.arange(n, dtype=float) + 0.5) / n
    return y_norm, u_sim, T_sim, nu_sim


def load_all_timesteps(csv_file: str, every_n: int = TIME_EVERY_N,
                       max_index: int | None = TIME_MAX_INDEX) -> dict:
    """
    Return time-series data at every every_n-th CSV entry.

    Returns a dict  {csv_timestep_index: {"y_norm", "u_sim", "T_sim", "sim_time"}}
    Always includes the last kept timestep so the terminal profile is present.

    every_n   : sampling stride in CSV entries (not simulation steps).
    max_index : if set, discard all CSV entries with index > max_index before
                sampling, so the range [0, max_index] is treated as the full data.
    """
    df = pd.read_csv(csv_file)
    all_ts = sorted(df["TimeStep"].unique())
    if max_index is not None:
        all_ts = [t for t in all_ts if t <= max_index]

    # Select every every_n-th entry, always keeping the last one.
    selected = all_ts[::every_n]
    if all_ts[-1] not in selected:
        selected = list(selected) + [all_ts[-1]]

    result = {}
    for ts in selected:
        ss = df[df["TimeStep"] == ts].sort_values("arc_length")
        ss = ss[ss["arc_length"] < ss["arc_length"].max()]
        n = len(ss)
        y_norm = (np.arange(n, dtype=float) + 0.5) / n
        result[int(ts)] = {
            "y_norm":   y_norm,
            "u_sim":    ss["velocity:0"].values,
            "T_sim":    ss["temperature"].values,
            "sim_time": int(ss["Time"].iloc[0]),
        }
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Per-case processing
# ──────────────────────────────────────────────────────────────────────────────

def process_case(csv_file, case_name, T_bottom, T_top, u_wall, nu_0, beta):
    print(f"\nProcessing: {case_name}")
    print(f"  CSV: {csv_file}")

    y_norm, u_sim, T_sim, nu_sim = load_csv_and_extract_steady_state(csv_file)
    print(f"  {len(y_norm)} cells  |  y ∈ ({y_norm[0]:.4f}, {y_norm[-1]:.4f})")

    if beta != 0.0:
        # Use actual simulation viscosity to avoid wavy errors from piecewise-polynomial
        # ν(T) deviating from the idealised exponential.
        u_analytical = CouetteAnalyticalSolver.velocity_profile_from_nu(y_norm, nu_sim, u_wall)
    else:
        u_analytical = CouetteAnalyticalSolver.velocity_profile_exponential(
            y_norm, T_bottom, T_top, u_wall, nu_0, beta)
    T_profile  = CouetteAnalyticalSolver.temperature_profile_linear(y_norm, T_bottom, T_top)
    nu_profile = nu_sim  # actual piecewise-polynomial viscosity used by the simulation

    error_abs = u_sim - u_analytical
    error_rel = np.abs(error_abs) / u_wall * 100.0
    print(f"  L∞ error:          {np.max(np.abs(error_abs)):.3e}")
    print(f"  L2 error:          {np.sqrt(np.mean(error_abs**2)):.3e}")
    print(f"  Max relative error:{error_rel.max():.3f}%  (% of u_wall={u_wall})")

    return dict(y_norm=y_norm, u_sim=u_sim, u_analytical=u_analytical,
                T_profile=T_profile, nu_profile=nu_profile, T_sim=T_sim, nu_sim=nu_sim)


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def _resolve_csv(base: str, search_dir: Path) -> Path | None:
    """
    Find the CSV file for a given base name, ignoring the _dat_NNN_TTT suffix.

    Globs for  {base}_dat_*.csv  and returns the file with the highest step
    count (the number after the first underscore in the dat suffix), so the
    longest / most-converged run is always picked automatically.

    Falls back to {base}.csv if no dat-suffixed file exists.
    """
    candidates = sorted(search_dir.glob(f"{base}_dat_*.csv"))
    if candidates:
        # Pick the file whose last-step index is highest
        def _step(p: Path) -> int:
            m = re.search(r"_dat_(\d+)_", p.name)
            return int(m.group(1)) if m else 0
        return max(candidates, key=_step)
    direct = search_dir / f"{base}.csv"
    return direct if direct.exists() else None


def main():
    print("=" * 70)
    print("COUETTE FLOW - VALIDATION PLOTS (PUBLICATION VERSION)")
    print("=" * 70)
    print(f"Time-evolution sampling: every {TIME_EVERY_N} CSV entries  "
          f"(change TIME_EVERY_N at the top of this file)")

    script_dir = Path(__file__).parent

    # Wall temperatures and u_max come from CouetteFlowScaling.prm:
    #   T_bottom = 300 K,  T_top = 3000 K,  u_max = 0.025
    # The exponential viscosity model is from CouetteFlowMaterial.yaml:
    #   ν(T) = 0.1667 * exp(-0.0005 * (T - 300))
    # csv_base: stem without the _dat_NNN_TTT suffix — resolved automatically.
    cases_config = [
        dict(csv_base="cf_cpu_mfconst_0.04",   name="Constant ν=0.04",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.04,   beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.06",   name="Constant ν=0.06",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.06,   beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.08",   name="Constant ν=0.08",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.08,   beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.1",    name="Constant ν=0.1",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.1,    beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.1667", name="Constant ν=0.1667",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.1667, beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.2",    name="Constant ν=0.2",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.2,    beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.4",    name="Constant ν=0.4",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.4,    beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.6",    name="Constant ν=0.6",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.6,    beta=0.0),
        dict(csv_base="cf_cpu_mfconst_0.8",    name="Constant ν=0.8",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.8,    beta=0.0),
        dict(csv_base="cf_cpu_mfconst_1.0",    name="Constant ν=1.0",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=1.0,    beta=0.0),
        dict(csv_base="cf_cpu_mftempdep",      name="Temperature-Dependent",
             T_bottom=300.0, T_top=3000.0, u_wall=0.025, nu_0=0.1667, beta=-0.0005),
    ]

    # Resolve each base name to the actual CSV path
    for cfg in cases_config:
        cfg["csv_file"] = _resolve_csv(cfg["csv_base"], script_dir)

    plotter = CouetteFlowPlotter(dpi=300)
    all_cases       = {}
    all_time_series = {}
    all_ana_funcs   = {}

    # ── steady-state processing ───────────────────────────────────────────────
    print("\n--- Loading and processing cases (steady state) ---")
    for cfg in cases_config:
        p = cfg["csv_file"]
        if p is None:
            print(f"  WARNING: no CSV found for {cfg['csv_base']}, skipping")
            continue
        try:
            data = process_case(str(p), cfg["name"], cfg["T_bottom"], cfg["T_top"],
                                cfg["u_wall"], cfg["nu_0"], cfg["beta"])
            all_cases[cfg["name"]] = data
            plotter.add_case(cfg["name"], data["y_norm"], data["u_sim"],
                             data["u_analytical"], data["T_profile"], data["nu_profile"])
        except Exception:
            import traceback; traceback.print_exc()

    if not all_cases:
        print("No cases processed. Exiting.")
        return
    print(f"\n{len(all_cases)} case(s) loaded successfully.")

    # ── time-series loading ───────────────────────────────────────────────────
    print(f"\n--- Loading time-series data (every {TIME_EVERY_N} CSV entries) ---")
    for cfg in cases_config:
        p = cfg["csv_file"]
        if p is None:
            continue
        try:
            ts = load_all_timesteps(str(p), every_n=TIME_EVERY_N)
            all_time_series[cfg["name"]] = ts
            nu_sim_ss = all_cases[cfg["name"]].get("nu_sim") if cfg["name"] in all_cases else None
            all_ana_funcs[cfg["name"]]   = make_analytical_func(
                cfg["T_bottom"], cfg["T_top"], cfg["u_wall"], cfg["nu_0"], cfg["beta"],
                nu_sim=nu_sim_ss)
            print(f"  {cfg['name']}: {len(ts)} snapshots loaded")
        except Exception:
            import traceback; traceback.print_exc()

    # ── steady-state plots ────────────────────────────────────────────────────
    print("\n--- Generating per-case steady-state plots ---")
    for name in all_cases:
        prefix = name.replace(" ", "_").replace("=", "")
        plotter.plot_single_case_detailed(name, prefix)
        plotter.plot_single_case_errors(name, prefix)

    print("\n--- Generating comparison plots ---")
    plotter.plot_all_cases_comparison("couette_all_cases_comparison.png")
    plotter.plot_all_cases_direct_comparison("couette_all_cases_direct_comparison.png")
    plotter.plot_error_metrics_summary("couette_error_summary.png")

    # ── time-evolution plots ──────────────────────────────────────────────────
    print("\n--- Generating time-evolution plots ---")
    for name, ts in all_time_series.items():
        prefix = name.replace(" ", "_").replace("=", "")
        plotter.plot_time_evolution(name, ts, all_ana_funcs[name], prefix)

    if len(all_time_series) > 1:
        plotter.plot_convergence_comparison(
            all_time_series, all_ana_funcs,
            "couette_convergence_comparison.png")

    # ── summary table ─────────────────────────────────────────────────────────
    print("\n--- Summary table ---")
    df = plotter.export_summary_table("couette_validation_summary.csv")
    print(df.to_string(index=False))
    print("\nDone.")


if __name__ == "__main__":
    main()
