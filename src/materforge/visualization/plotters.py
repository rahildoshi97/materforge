# SPDX-FileCopyrightText: 2025 - 2026 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-FileCopyrightText: 2026 Matthias Markl, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause

import logging
from datetime import datetime
from typing import Optional, Union
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from matplotlib.gridspec import GridSpec
from materforge.core.materials import Material
from materforge.algorithms.regression_processor import RegressionProcessor
from materforge.parsing.config.yaml_keys import CONSTANT_KEY, NAME_KEY, POST_KEY, PRE_KEY
from materforge.data.constants import ProcessingConstants

logger = logging.getLogger(__name__)


class PropertyVisualizer:
    """Handles visualization of material properties."""

    # --- Constructor ---
    def __init__(self, parser) -> None:
        self.parser = parser
        self.fig = None
        self.gs = None
        self.current_subplot = 0
        self.plot_directory = self.parser.base_dir / "materforge_plots"
        self.visualized_properties = set()
        self.is_enabled = True
        self.setup_style()
        logger.debug("PropertyVisualizer initialized for: %s", parser.config_path)

    @staticmethod
    def setup_style() -> None:
        """Applies global matplotlib style settings."""
        plt.rcParams.update({
            'font.size': 10,
            'font.family': 'sans-serif',
            'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica', 'Liberation Sans'],
            'axes.titlesize': 12,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9,
            'figure.titlesize': 14,
            'axes.grid': True,
            'grid.alpha': 0.3,
            'grid.linestyle': '--',
            'axes.axisbelow': True,
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'savefig.facecolor': 'white',
            'savefig.edgecolor': 'none',
            'savefig.dpi': 300,
            'figure.autolayout': True,
        })

    # --- Public API ---
    def is_visualization_enabled(self) -> bool:
        """Returns True if visualization is active and a figure exists."""
        return self.is_enabled and self.fig is not None

    def initialize_plots(self) -> None:
        """Initialises the figure and grid layout for all properties."""
        if not self.is_enabled:
            logger.debug("Visualization disabled - skipping plot initialization")
            return
        if self.parser.categorized_properties is None:
            raise ValueError("No properties to plot.")
        property_count = sum(len(props) for props in self.parser.categorized_properties.values())
        logger.info("Initializing visualization for %d properties", property_count)
        fig_width = 12
        fig_height = max(4 * property_count, 4)  # Minimum height for readability
        self.fig = plt.figure(figsize=(fig_width, fig_height))
        self.gs = GridSpec(property_count, 1, figure=self.fig, )
        self.current_subplot = 0
        self.plot_directory.mkdir(exist_ok=True)
        logger.debug("Plot directory ready: %s", self.plot_directory)

    def reset_visualization_tracking(self) -> None:
        """Clears the set of already-visualized properties."""
        logger.debug("Resetting visualization tracking - clearing %d tracked properties",
                     len(self.visualized_properties))
        self.visualized_properties = set()

    def visualize_property(
        self,
        material: Material,
        prop_name: str,
        dependency: Union[float, sp.Symbol],
        prop_type: str,
        x_data: Optional[np.ndarray] = None,
        y_data: Optional[np.ndarray] = None,
        has_regression: bool = False,
        simplify_type: Optional[str] = None,
        degree: int = 1,
        segments: int = 1,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        lower_bound_type: str = CONSTANT_KEY,
        upper_bound_type: str = CONSTANT_KEY,) -> None:
        """Visualizes a single material property on the next available subplot."""
        if prop_name in self.visualized_properties:
            logger.debug("Property %r already visualized - skipping", prop_name)
            return
        if self.fig is None:
            logger.warning("No figure available for %r - visualization skipped", prop_name)
            return
        if not isinstance(dependency, sp.Symbol):
            logger.debug("Numeric dependency for %r - visualization skipped", prop_name)
            return
        logger.info("Visualizing %r (%s) for material %r", prop_name, prop_type, material.name)
        try:
            ax = self.fig.add_subplot(self.gs[self.current_subplot])
            self.current_subplot += 1
            ax.set_aspect('auto')
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.set_axisbelow(True)
            for spine in ax.spines.values():
                spine.set_color('#CCCCCC')
                spine.set_linewidth(1.2)
            current_prop = getattr(material, prop_name)
            # Build dependency evaluation range
            if x_data is not None and len(x_data) > 0:
                data_lower, data_upper = np.min(x_data), np.max(x_data)
                step = (data_upper - data_lower) / 1000
            else:
                data_lower = ProcessingConstants.DEFAULT_DEPENDENCY_LOWER
                data_upper = ProcessingConstants.DEFAULT_DEPENDENCY_UPPER
                step = (data_upper - data_lower) / 1000
                logger.debug("Using default dependency range: %.1f-%.1f", data_lower, data_upper)
            if lower_bound is None:
                lower_bound = data_lower
            if upper_bound is None:
                upper_bound = data_upper
            padding = (upper_bound - lower_bound) * ProcessingConstants.DEPENDENCY_PADDING_FACTOR
            padded_lower = lower_bound - padding
            padded_upper = upper_bound + padding
            num_points = int(np.ceil((padded_upper - padded_lower) / step)) + 1
            extended_dep = np.linspace(padded_lower, padded_upper, num_points)
            ax.set_title(f"{prop_name} ({prop_type})", fontsize=14, fontweight="bold", pad=15)
            ax.set_xlabel(str(dependency), fontsize=12, fontweight="bold")
            ax.set_ylabel(prop_name, fontsize=12, fontweight="bold")
            colors = {
                'constant': '#1f77b4',  # blue
                'raw': '#ff7f0e',  # orange
                'regression_pre': '#2ca02c',  # green
                'regression_post': '#d62728',  # red
                'bounds': '#9467bd',  # purple
                'extended': '#8c564b',  # brown
            }
            _y_value = 0.0
            if prop_type == 'CONSTANT_VALUE':
                value = float(current_prop)
                ax.plot(extended_dep, np.full_like(extended_dep, value), color=colors['constant'], linestyle='-',
                           linewidth=2.5, label='constant', alpha=0.8)
                ax.text(0.5, 0.9, f"Value: {value:.3e}", transform=ax.transAxes,
                        horizontalalignment="center", fontweight="bold",
                        bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5',
                                  edgecolor=colors['constant']))
                ax.set_ylim(value * 0.9, value * 1.1)
                y_range = ax.get_ylim()
                offset = (y_range[1] - y_range[0]) * 0.1
                _y_value = value + offset
                logger.debug("Plotted constant %r: %g", prop_name, value)
            elif prop_type == 'STEP_FUNCTION':
                try:
                    f_current = sp.lambdify(dependency, current_prop, 'numpy')
                    y_extended = f_current(extended_dep)
                    ax.plot(extended_dep, y_extended, color=colors['extended'],
                            linestyle='-', linewidth=2.5, label='extended behavior',
                            zorder=1, alpha=0.6)
                    if x_data is not None and y_data is not None:
                        ax.plot(x_data, y_data, color=colors['raw'], linestyle='-',
                                linewidth=2.5, marker="o", markersize=6,
                                label='step function', zorder=3, alpha=0.8)
                        transition_idx = len(x_data) // 2
                        transition_point = x_data[transition_idx]
                        ax.axvline(x=transition_point, color='red', linestyle='--',
                                   alpha=0.7, linewidth=2, label='transition point')
                        ax.text(transition_point, y_data[0], f" Before: {y_data[0]:.2e}",
                                verticalalignment="bottom", horizontalalignment="left",
                                fontweight="bold",
                                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
                        ax.text(transition_point, y_data[-1], f" After: {y_data[-1]:.2e}",
                                verticalalignment="top", horizontalalignment="left",
                                fontweight="bold",
                                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
                        _y_value = float(np.mean(y_data))
                        logger.debug("Step %r - transition at %.1f", prop_name, transition_point)
                    else:
                        _y_value = float(f_current(lower_bound))
                except Exception as e:
                    logger.warning("Could not evaluate step function %r: %s", prop_name, e)
            else:
                try:
                    f_current = sp.lambdify(dependency, current_prop, 'numpy')
                    if has_regression and simplify_type == PRE_KEY:
                        main_color = colors['regression_pre']
                        main_label = 'regression (pre)'
                    else:
                        main_color = colors['extended']
                        main_label = 'raw (extended)'
                    try:
                        y_extended = f_current(extended_dep)
                        ax.plot(extended_dep, y_extended, color=main_color,
                                linestyle='-', linewidth=2.5, label=main_label,
                                zorder=2, alpha=0.8)
                    except Exception as e:
                        logger.warning("Could not evaluate extended range for %r: %s", prop_name, e)
                        if x_data is not None and y_data is not None:
                            ax.plot(x_data, y_data, color=colors['raw'],
                                    linestyle='-', linewidth=2, label='data points', zorder=2)
                    if y_data is not None and len(y_data) > 0:
                        _y_value = float(np.percentile(y_data, 25))
                    else:
                        try:
                            midpoint = (lower_bound + upper_bound) / 2
                            _y_value = float(f_current(midpoint))
                        except (ValueError, TypeError, AttributeError) as e:
                            logger.error("Could not evaluate midpoint for %r: %s", prop_name, e)
                            _y_value = 0.0
                    if has_regression and simplify_type == POST_KEY and x_data is not None and y_data is not None:
                        try:
                            preview_pw = RegressionProcessor.process_regression(
                                dep_array=x_data, prop_array=y_data, dependency=dependency,
                                lower_bound_type=lower_bound_type, upper_bound_type=upper_bound_type,
                                degree=degree, segments=segments,
                                seed=ProcessingConstants.DEFAULT_REGRESSION_SEED
                            )
                            f_preview = sp.lambdify(dependency, preview_pw, 'numpy')
                            ax.plot(extended_dep, f_preview(extended_dep), color=colors['regression_post'],
                                    linestyle='--', linewidth=2.5, label='regression (post)',
                                    zorder=4, alpha=0.8)
                        except Exception as e:
                            logger.warning("Post-regression preview failed for %r: %s", prop_name, e)
                except Exception as e:
                    logger.error("Error creating function for %r: %s", prop_name, e)
                    ax.text(0.5, 0.5, f"Error: {str(e)}", transform=ax.transAxes,
                            horizontalalignment="center", fontweight="bold",
                            bbox=dict(facecolor='red', alpha=0.2))
            # --- Boundary annotations ---
            ax.axvline(x=lower_bound, color=colors['bounds'], linestyle='--',
                       alpha=0.6, linewidth=1.5, label="_nolegend_")
            ax.axvline(x=upper_bound, color=colors['bounds'], linestyle='--',
                       alpha=0.6, linewidth=1.5, label="_nolegend_")
            if _y_value is None or not np.isfinite(_y_value):
                try:
                    if hasattr(current_prop, 'subs') and hasattr(current_prop, 'evalf'):
                        _y_value = float(current_prop.subs(dependency, lower_bound).evalf())
                    elif hasattr(current_prop, '__float__'):
                        _y_value = float(current_prop)
                    else:
                        _y_value = 0.0
                except (ValueError, TypeError, AttributeError):
                    _y_value = 0.0
                    logger.warning("Could not determine y_value for annotations of %r", prop_name)
            ax.text(lower_bound, _y_value, f" {lower_bound_type}",
                    verticalalignment="top", horizontalalignment="right", fontweight="bold",
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3',
                              edgecolor=colors['bounds']))
            ax.text(upper_bound, _y_value, f" {upper_bound_type}",
                    verticalalignment="top", horizontalalignment="left", fontweight="bold",
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3',
                              edgecolor=colors['bounds']))
            if has_regression and degree is not None:
                ax.text(0.5, 0.98, f"Simplify: {simplify_type} | Degree: {degree} | Segments: {segments}",
                        transform=ax.transAxes, horizontalalignment="center", fontweight="bold",
                        bbox=dict(facecolor='lightblue', alpha=0.8, boxstyle='round,pad=0.3'))
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                legend = ax.legend(handles, labels, loc='best', framealpha=0.9,
                                   fancybox=True, shadow=True, edgecolor="gray")
                legend.get_frame().set_linewidth(1.2)
            self.visualized_properties.add(prop_name)
            logger.info("Successfully visualized %r", prop_name)
        except Exception as e:
            logger.error("Unexpected error visualizing %r: %s", prop_name, e, exc_info=True)
            raise ValueError(f"Unexpected error visualizing {prop_name!r}: {e}") from e

    def save_property_plots(self) -> None:
        """Saves the composed property figure to disk and closes it."""
        if not self.is_enabled or self.fig is None:
            logger.debug("No plots to save - visualization disabled or no figure")
            return
        try:
            material_name = self.parser.config[NAME_KEY]
            self.fig.suptitle(f"Material Properties: {material_name}", fontsize=16, fontweight="bold", y=0.98)
            try:
                plt.tight_layout(rect=[0, 0.01, 1, 0.98], pad=1.0)
            except Exception as e:
                logger.warning("tight_layout failed: %s - using subplots_adjust", e)
                plt.subplots_adjust(left=0.08, bottom=0.08, right=0.92, top=0.88, hspace=0.8)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{material_name.replace(' ', '_')}_properties_{timestamp}.png"
            filepath = self.plot_directory / filename
            self.fig.savefig(str(filepath), dpi=300, bbox_inches='tight',
                             facecolor='white', edgecolor='none', pad_inches=0.4)
            total = sum(len(p) for p in self.parser.categorized_properties.values())
            visualized = len(self.visualized_properties)
            if visualized != total:
                logger.warning("Not all properties visualized - %d/%d", visualized, total)
            else:
                logger.info("All %d properties visualized successfully", total)
            logger.info("Property plots saved: %s", filepath)
        finally:
            if self.fig is not None:
                plt.close(self.fig)
                self.fig = None
                logger.debug("Figure closed and memory freed")
