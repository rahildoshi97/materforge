import logging
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from matplotlib.gridspec import GridSpec
from datetime import datetime

from pymatlib.core.materials import Material
from pymatlib.algorithms.regression_processor import RegressionProcessor
from pymatlib.parsing.config.yaml_keys import CONSTANT_KEY, PRE_KEY, POST_KEY, NAME_KEY, MATERIAL_TYPE_KEY
from pymatlib.data.constants import PhysicalConstants, ProcessingConstants

logger = logging.getLogger(__name__)


class PropertyVisualizer:
    """Handles visualization of material properties."""

    # --- Constructor ---
    def __init__(self, parser) -> None:
        self.parser = parser
        self.fig = None
        self.gs = None
        self.current_subplot = 0
        yaml_dir = self.parser.base_dir
        self.plot_directory = yaml_dir / "pymatlib_plots"
        self.visualized_properties = set()
        self.is_enabled = True
        self.setup_style()
        logger.debug("PropertyVisualizer initialized for parser: %s", parser.config_path)

    @staticmethod
    def setup_style() -> None:
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
            'figure.autolayout': True
        })

    def is_visualization_enabled(self) -> bool:
        """Check if visualization is currently enabled."""
        return self.is_enabled and self.fig is not None

    # --- Public API Methods ---
    def initialize_plots(self) -> None:
        """Initialize plots only if visualization is enabled."""
        if not self.is_enabled:
            logger.debug("Visualization disabled, skipping plot initialization")
            return
        if self.parser.categorized_properties is None:
            logger.error("categorized_properties is None - cannot initialize plots")
            raise ValueError("No properties to plot.")
        property_count = sum(len(props) for props in self.parser.categorized_properties.values())
        logger.info("Initializing visualization for %d properties", property_count)
        fig_width = 10
        fig_height = max(6 * property_count, 6)  # Minimum height for readability (5.6 for thermal_diffusivity)
        self.fig = plt.figure(figsize=(fig_width, fig_height))
        self.gs = GridSpec(property_count, 1, figure=self.fig, )
        # self.gs = GridSpec(property_count, ncols=2, figure=self.fig, width_ratios=[1, 1], wspace=0.2)
        self.current_subplot = 0
        self.plot_directory.mkdir(exist_ok=True)
        logger.debug("Plot directory created: %s", self.plot_directory)

    def reset_visualization_tracking(self) -> None:
        logger.debug("Resetting visualization tracking - clearing %d tracked properties",
                     len(self.visualized_properties))
        self.visualized_properties = set()

    def generate_yaml_snippet(self, prop_name: str, prop_type: str,
                              has_regression: bool = False,
                              simplify_type: Optional[str] = None,
                              degree: int = 1, segments: int = 1) -> str:
        """Generate YAML snippet based on property configuration."""

        # Method-specific YAML templates
        yaml_templates = {
            'CONSTANT': f'''{prop_name.lower()}: 6950.0''',

            'STEP_FUNCTION': f'''{prop_name.lower()}:
      temperature: melting_temperature
      value: [7500.0, 6500.0]
      bounds: [extrapolate, extrapolate]''',

            'KEY_VAL': f'''{prop_name.lower()}:
      temperature: [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1500, 2000, 2500, 3000]
      value: [7747, 7716, 7685, 7652, 7617, 7582, 7545, 7508, 7469, 7305, 6825, 6388, 5925]
      bounds: [extrapolate, extrapolate]''',

            'FILE': f'''{prop_name.lower()}:
      file_path: ./SS304L.xlsx
      temperature_header: T (K)
      value_header: Density (kg/(m)^3)
      bounds: [extrapolate, extrapolate]''',

            'PIECEWISE_EQUATION': f'''{prop_name.lower()}:
      temperature: [300, 1660, 1736, 3000] 
      equation: [7877.39163826692 - 0.377781577789007*T, 11816.6337569868 - 2.74041499241854*T, 8596.40178865677 - 0.885849240373116*T]
      bounds: [extrapolate, extrapolate]''',

            'COMPUTE': f'''{prop_name.lower()}:
      temperature: (300, 3000, 5.0)
      equation: 2700 * (1 - 3*thermal_expansion_coefficient * (T - 293))
      bounds: [extrapolate, extrapolate]'''
        }

        # Get base YAML
        yaml_snippet = yaml_templates.get(prop_type, f'''{prop_name.lower()}:
      # Configuration for {prop_type} property
      bounds: [constant, constant]''')

        # Add regression configuration if applicable
        if has_regression and prop_type != 'CONSTANT':
            regression_config = f'''
      regression:
        simplify: {simplify_type}
        degree: {degree}
        segments: {segments}'''
            yaml_snippet += regression_config

        return yaml_snippet

    def visualize_property(
            self,
            material: Material,
            prop_name: str,
            T: Union[float, sp.Symbol],
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
            upper_bound_type: str = CONSTANT_KEY) -> None:
        """Visualize a single property."""
        if prop_name == 'heat_conductivity':
            return
        if prop_name == 'density':
            return
        if prop_name == 'heat_capacity':
            return
        if prop_name in self.visualized_properties:
            logger.debug("Property '%s' already visualized, skipping", prop_name)
            return
        if not hasattr(self, 'fig') or self.fig is None:
            logger.warning("No figure available for property '%s' - visualization skipped", prop_name)
            return
        if not isinstance(T, sp.Symbol):
            logger.debug("Temperature is not symbolic for property '%s' - visualization skipped", prop_name)
            return
        logger.info("Visualizing property: %s (type: %s) for material: %s",
                    prop_name, prop_type, material.name)
        try:
            """# Method-specific colors
            method_colors = {
                'CONSTANT': '#FF6B6B',
                'STEP_FUNCTION': '#4ECDC4',
                'KEY_VAL': '#45B7D1',
                'FILE': '#96CEB4',
                'PIECEWISE_EQUATION': '#FFEAA7',
                'COMPUTE': '#DDA0DD'
            }
            plot_color = method_colors.get(prop_type, '#8c564b')

            # LEFT PANEL: YAML Configuration
            ax_yaml = self.fig.add_subplot(self.gs[self.current_subplot, 0])
            ax_yaml.axis('off')

            # Generate YAML snippet
            yaml_text = self.generate_yaml_snippet(
                prop_name=prop_name,
                prop_type=prop_type,
                has_regression=has_regression,
                simplify_type=simplify_type,
                degree=degree,
                segments=segments
            )

            # Add YAML text
            ax_yaml.text(0.05, 0.95, yaml_text,
                         fontsize=10,
                         family='monospace',
                         verticalalignment='top',
                         horizontalalignment='left',
                         bbox=dict(boxstyle="round,pad=0.4",
                                   facecolor="lightgray",
                                   alpha=0.2,
                                   edgecolor=plot_color,
                                   linewidth=1.5))

            # Add method type label
            ax_yaml.text(0.05, 0.05, f"Method: {prop_type}",
                         fontsize=10,
                         fontweight='bold',
                         color=plot_color,
                         verticalalignment='bottom',
                         horizontalalignment='left')

            # Add title for YAML section
            # ax_yaml.set_title(f"YAML Configuration\n{prop_name}", fontsize=12, fontweight='bold', pad=20)
            ax_yaml.set_title(f"YAML Configuration", fontsize=12, fontweight='bold', pad=20)"""

            ax = self.fig.add_subplot(self.gs[self.current_subplot, 0])
            # RIGHT PANEL: Property Plot
            # ax = self.fig.add_subplot(self.gs[self.current_subplot, 0])
            ax.tick_params(axis='both', which='major', labelsize=16)
            self.current_subplot += 1
            ax.set_aspect('auto')
            # Grid and border styling
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.set_axisbelow(True)
            # Border styling
            for spine in ax.spines.values():
                spine.set_color('#CCCCCC')
                spine.set_linewidth(1.2)
            # Get property and prepare temperature array
            current_prop = getattr(material, prop_name)
            if x_data is not None and len(x_data) > 0:
                # Use property's own temperature range
                data_lower, data_upper = np.min(x_data), np.max(x_data)
                temp_range = data_upper - data_lower
                step = temp_range / 1000  # Create 1000 points for smooth visualization
            else:
                # Fallback for properties without explicit temperature data
                data_lower, data_upper = (ProcessingConstants.DEFAULT_TEMP_LOWER,
                                          ProcessingConstants.DEFAULT_TEMP_UPPER)
                step = (data_upper - data_lower) / 1000
                logger.debug("Using data temperature range: %.1f - %.1f K", data_lower, data_upper)
            # Set bounds with property-specific defaults
            if lower_bound is None:
                lower_bound = data_lower
            if upper_bound is None:
                upper_bound = data_upper
                logger.debug("Using default temperature range: %.1f - %.1f K", data_lower, data_upper)
            # Create extended temperature range for visualization
            padding = (upper_bound - lower_bound) * ProcessingConstants.TEMPERATURE_PADDING_FACTOR
            padded_lower = max(lower_bound - padding, PhysicalConstants.ABSOLUTE_ZERO)
            padded_upper = upper_bound + padding
            num_points = int(np.ceil((padded_upper - padded_lower) / step)) + 1
            extended_temp = np.linspace(padded_lower, padded_upper, num_points)
            # Title and labels
            ax.set_title(f"{prop_name} ({prop_type})", fontsize=24, fontweight='bold', pad=10)
            ax.set_xlabel("Temperature", fontsize=22, fontweight='bold')
            ax.set_ylabel(f"{prop_name}", fontsize=22,fontweight='bold')
            # Color scheme
            colors = {
                'constant': '#1f77b4',  # blue
                'raw': '#ff7f0e',  # orange
                'regression_pre': '#2ca02c',  # green
                'regression_post': '#e62728',  # red
                'bounds': '#9467bd',  # purple
                'extended': '#8c564b',  # brown
            }
            # Initialize y_value for annotations
            _y_value = 0.0
            if prop_type == 'CONSTANT_VALUE':
                value = float(current_prop)
                ax.axhline(y=value, color='#bcbd22', linestyle='-',
                           linewidth=3.5, label='constant', alpha=0.8)
                # Annotation
                ax.text(0.5, 0.9, f"Value: {value:.3e}", transform=ax.transAxes,
                        horizontalalignment='center', fontweight='bold',
                        bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5',
                                  edgecolor=colors['constant']))
                ax.set_ylim(value * 0.9, value * 1.1)
                # Add small offset to avoid overlap with horizontal line
                y_range = ax.get_ylim()
                offset = (y_range[1] - y_range[0]) * 0.1
                _y_value = value + offset
                logger.debug("Plotted constant property '%s' with value: %g", prop_name, value)
            elif prop_type == 'STEP_FUNCTION':
                try:
                    f_current = sp.lambdify(T, current_prop, 'numpy')
                    # Always plot the extended behavior first (background)
                    y_extended = f_current(extended_temp)
                    """ax.plot(extended_temp, y_extended, color=colors['extended'],
                            linestyle='-', linewidth=2.5, label='extended behavior',
                            zorder=1, alpha=0.6)"""
                    # Overlay data points if available (foreground)
                    if x_data is not None and y_data is not None:
                        ax.plot(x_data, y_data, color='#bcbd22', linestyle='-',
                                linewidth=3.5, marker='o', markersize=3.5,
                                label='step function', zorder=3, alpha=0.8)
                        # Add vertical line at transition point
                        transition_idx = len(x_data) // 2
                        transition_temp = x_data[transition_idx]
                        """ax.axvline(x=transition_temp, color='red', linestyle='--',
                                   alpha=0.7, linewidth=2, label='transition point')
                        # Annotations
                        ax.text(transition_temp, y_data[0], f' Before: {y_data[0]:.2e}',
                                verticalalignment='bottom', horizontalalignment='left',
                                fontweight='bold',
                                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
                        ax.text(transition_temp, y_data[-1], f' After: {y_data[-1]:.2e}',
                                verticalalignment='top', horizontalalignment='left',
                                fontweight='bold',
                                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))"""
                        _y_value = np.mean(y_data)
                        logger.debug("Plotted step function '%s' with transition at %.1f K",
                                     prop_name, transition_temp)
                    else:
                        # No data points available, use function evaluation
                        _y_value = f_current(lower_bound)
                except Exception as e:
                    logger.warning("Could not evaluate step function '%s': %s", prop_name, e)
                    _y_value = 0.0
            else:  # Handle all other property types (FILE_IMPORT, TABULAR_DATA, PIECEWISE_EQUATION, COMPUTED_PROPERTY)
                try:
                    f_current = sp.lambdify(T, current_prop, 'numpy')
                    # Determine the appropriate label and color based on regression status
                    if has_regression and simplify_type == PRE_KEY:
                        main_color = colors['regression_pre']
                        main_label = 'regression (pre)'
                    else:
                        main_color = colors['extended']
                        main_label = 'raw data'  #
                    try:  # Plot the main function over extended range
                        y_extended = f_current(extended_temp)
                        ax.plot(extended_temp, y_extended, color='#bcbd22',  # Yellow-Green
                                linestyle='-', linewidth=3.5, label=main_label,
                                marker='o', markersize=3.5, zorder=2, alpha=0.8)  # markersize=3.5 for boundary_behavior.png
                        logger.debug("Plotted extended range for property '%s'", prop_name)
                    except Exception as e:
                        logger.warning("Could not evaluate function over extended range for '%s': %s",
                                       prop_name, e)
                        # Fallback to data range if available
                        if x_data is not None and y_data is not None:
                            ax.plot(x_data, y_data, color=colors['regression_pre'],  #
                                    linestyle='None', linewidth=2, label='data points', marker='o', markersize=4, zorder=2)
                    # Plot data points if available (for FILE, KEY_VAL properties)
                    if x_data is not None and y_data is not None and prop_type in ['FILE', 'KEY_VAL']:
                        # marker_size = 6 if prop_type == 'KEY_VAL' else 3
                        # ax.scatter(x_data, y_data, color=colors['raw'], marker='o', s=marker_size**2,
                        #            alpha=0.7, label='data points', zorder=3)
                        pass
                    # Set y_value for boundary annotations to avoid overlap
                    if y_data is not None and len(y_data) > 0:
                        # Use 25th percentile instead of max to avoid high regions
                        _y_value = np.percentile(y_data, 25)
                    else:
                        try:
                            # Use midpoint instead of upper_bound
                            midpoint = (lower_bound + upper_bound) / 2
                            _y_value = f_current(midpoint)
                        except (ValueError, TypeError, AttributeError) as e:
                            logger.error(f"Could not evaluate function at midpoint for '%s': %s", prop_name, e)
                            _y_value = 0.0
                    # Post-regression overlay
                    if has_regression and simplify_type == POST_KEY and x_data is not None and y_data is not None:
                        try:
                            preview_pw = RegressionProcessor.process_regression(
                                temp_array=x_data, prop_array=y_data, T=T,
                                lower_bound_type=lower_bound_type, upper_bound_type=upper_bound_type,
                                degree=degree, segments=segments, seed=ProcessingConstants.DEFAULT_REGRESSION_SEED
                            )
                            f_preview = sp.lambdify(T, preview_pw, 'numpy')
                            extended_temp = np.linspace(lower_bound - 270.0, upper_bound + 270.0, num_points)  # For boundary_behavior.png
                            y_preview = f_preview(extended_temp)
                            print(f"x_data: {x_data.tolist()}")
                            print(f"y_data: {y_data.tolist()}")
                            print(f"extended_temp: {extended_temp.tolist()}")
                            print(f"y_preview: {y_preview.tolist()}")
                            ax.plot(extended_temp, y_preview, color=colors['regression_post'],  #
                                    linestyle='--', linewidth=2.5, label='Degree: 1 | Segments = 1',
                                    zorder=4, alpha=0.8)
                            """temp5 = np.asarray([30.0, 33.24, 36.480000000000004, 39.72, 42.96, 46.2, 49.44, 52.68, 55.92, 59.160000000000004, 62.400000000000006, 65.64, 68.88, 72.12, 75.36, 78.6, 81.84, 85.08000000000001, 88.32000000000001, 91.56, 94.80000000000001, 98.04, 101.28, 104.52000000000001, 107.76, 111.0, 114.24000000000001, 117.48, 120.72, 123.96000000000001, 127.2, 130.44, 133.68, 136.92000000000002, 140.16000000000003, 143.4, 146.64000000000001, 149.88, 153.12, 156.36, 159.60000000000002, 162.84, 166.08, 169.32000000000002, 172.56, 175.8, 179.04000000000002, 182.28, 185.52, 188.76000000000002, 192.0, 195.24, 198.48000000000002, 201.72, 204.96, 208.20000000000002, 211.44, 214.68, 217.92000000000002, 221.16000000000003, 224.4, 227.64000000000001, 230.88000000000002, 234.12, 237.36, 240.60000000000002, 243.84, 247.08, 250.32000000000002, 253.56, 256.8, 260.04, 263.28000000000003, 266.52, 269.76, 273.0, 276.24, 279.48, 282.72, 285.96000000000004, 289.20000000000005, 292.44, 295.68, 298.92, 302.16, 305.40000000000003, 308.64000000000004, 311.88, 315.12, 318.36, 321.6, 324.84000000000003, 328.08000000000004, 331.32, 334.56, 337.8, 341.04, 344.28000000000003, 347.52000000000004, 350.76000000000005, 354.0, 357.24, 360.48, 363.72, 366.96000000000004, 370.20000000000005, 373.44, 376.68, 379.92, 383.16, 386.40000000000003, 389.64000000000004, 392.88, 396.12, 399.36, 402.6, 405.84000000000003, 409.08000000000004, 412.32000000000005, 415.56, 418.8, 422.04, 425.28000000000003, 428.52000000000004, 431.76000000000005, 435.0, 438.24, 441.48, 444.72, 447.96000000000004, 451.20000000000005, 454.44000000000005, 457.68, 460.92, 464.16, 467.40000000000003, 470.64000000000004, 473.88000000000005, 477.12, 480.36, 483.6, 486.84000000000003, 490.08000000000004, 493.32000000000005, 496.56000000000006, 499.8, 503.04, 506.28000000000003, 509.52000000000004, 512.76, 516.0, 519.24, 522.48, 525.72, 528.96, 532.2, 535.44, 538.6800000000001, 541.9200000000001, 545.1600000000001, 548.4000000000001, 551.64, 554.88, 558.12, 561.36, 564.6, 567.84, 571.08, 574.32, 577.5600000000001, 580.8000000000001, 584.0400000000001, 587.2800000000001, 590.52, 593.76, 597.0, 600.24, 603.48, 606.72, 609.96, 613.2, 616.44, 619.6800000000001, 622.9200000000001, 626.1600000000001, 629.4000000000001, 632.64, 635.88, 639.12, 642.36, 645.6, 648.84, 652.08, 655.32, 658.5600000000001, 661.8000000000001, 665.0400000000001, 668.2800000000001, 671.5200000000001, 674.76, 678.0, 681.24, 684.48, 687.72, 690.96, 694.2, 697.44, 700.6800000000001, 703.9200000000001, 707.1600000000001, 710.4000000000001, 713.6400000000001, 716.88, 720.12, 723.36, 726.6, 729.84, 733.08, 736.32, 739.5600000000001, 742.8000000000001, 746.0400000000001, 749.2800000000001, 752.5200000000001, 755.76, 759.0, 762.24, 765.48, 768.72, 771.96, 775.2, 778.44, 781.6800000000001, 784.9200000000001, 788.1600000000001, 791.4000000000001, 794.6400000000001, 797.88, 801.12, 804.36, 807.6, 810.84, 814.08, 817.32, 820.5600000000001, 823.8000000000001, 827.0400000000001, 830.2800000000001, 833.5200000000001, 836.7600000000001, 840.0, 843.24, 846.48, 849.72, 852.96, 856.2, 859.44, 862.6800000000001, 865.9200000000001, 869.1600000000001, 872.4000000000001, 875.6400000000001, 878.8800000000001, 882.12, 885.36, 888.6, 891.84, 895.08, 898.32, 901.5600000000001, 904.8000000000001, 908.0400000000001, 911.2800000000001, 914.5200000000001, 917.7600000000001, 921.0000000000001, 924.24, 927.48, 930.72, 933.96, 937.2, 940.44, 943.6800000000001, 946.9200000000001, 950.1600000000001, 953.4000000000001, 956.6400000000001, 959.8800000000001, 963.1200000000001, 966.36, 969.6, 972.84, 976.08, 979.32, 982.5600000000001, 985.8000000000001, 989.0400000000001, 992.2800000000001, 995.5200000000001, 998.7600000000001, 1002.0000000000001, 1005.24, 1008.48, 1011.72, 1014.96, 1018.2, 1021.44, 1024.68, 1027.92, 1031.16, 1034.4, 1037.64, 1040.88, 1044.1200000000001, 1047.3600000000001, 1050.6, 1053.8400000000001, 1057.0800000000002, 1060.3200000000002, 1063.5600000000002, 1066.8000000000002, 1070.04, 1073.28, 1076.52, 1079.76, 1083.0, 1086.24, 1089.48, 1092.72, 1095.96, 1099.2, 1102.44, 1105.68, 1108.92, 1112.16, 1115.4, 1118.64, 1121.88, 1125.1200000000001, 1128.3600000000001, 1131.6000000000001, 1134.8400000000001, 1138.0800000000002, 1141.3200000000002, 1144.5600000000002, 1147.8000000000002, 1151.04, 1154.28, 1157.52, 1160.76, 1164.0, 1167.24, 1170.48, 1173.72, 1176.96, 1180.2, 1183.44, 1186.68, 1189.92, 1193.16, 1196.4, 1199.64, 1202.88, 1206.1200000000001, 1209.3600000000001, 1212.6000000000001, 1215.8400000000001, 1219.0800000000002, 1222.3200000000002, 1225.5600000000002, 1228.8000000000002, 1232.0400000000002, 1235.28, 1238.52, 1241.76, 1245.0, 1248.24, 1251.48, 1254.72, 1257.96, 1261.2, 1264.44, 1267.68, 1270.92, 1274.16, 1277.4, 1280.64, 1283.88, 1287.1200000000001, 1290.3600000000001, 1293.6000000000001, 1296.8400000000001, 1300.0800000000002, 1303.3200000000002, 1306.5600000000002, 1309.8000000000002, 1313.0400000000002, 1316.28, 1319.52, 1322.76, 1326.0, 1329.24, 1332.48, 1335.72, 1338.96, 1342.2, 1345.44, 1348.68, 1351.92, 1355.16, 1358.4, 1361.64, 1364.88, 1368.1200000000001, 1371.3600000000001, 1374.6000000000001, 1377.8400000000001, 1381.0800000000002, 1384.3200000000002, 1387.5600000000002, 1390.8000000000002, 1394.0400000000002, 1397.2800000000002, 1400.52, 1403.76, 1407.0, 1410.24, 1413.48, 1416.72, 1419.96, 1423.2, 1426.44, 1429.68, 1432.92, 1436.16, 1439.4, 1442.64, 1445.88, 1449.1200000000001, 1452.3600000000001, 1455.6000000000001, 1458.8400000000001, 1462.0800000000002, 1465.3200000000002, 1468.5600000000002, 1471.8000000000002, 1475.0400000000002, 1478.2800000000002, 1481.52, 1484.76, 1488.0, 1491.24, 1494.48, 1497.72, 1500.96, 1504.2, 1507.44, 1510.68, 1513.92, 1517.16, 1520.4, 1523.64, 1526.88, 1530.1200000000001, 1533.3600000000001, 1536.6000000000001, 1539.8400000000001, 1543.0800000000002, 1546.3200000000002, 1549.5600000000002, 1552.8000000000002, 1556.0400000000002, 1559.2800000000002, 1562.5200000000002, 1565.76, 1569.0, 1572.24, 1575.48, 1578.72, 1581.96, 1585.2, 1588.44, 1591.68, 1594.92, 1598.16, 1601.4, 1604.64, 1607.88, 1611.1200000000001, 1614.3600000000001, 1617.6000000000001, 1620.8400000000001, 1624.0800000000002, 1627.3200000000002, 1630.5600000000002, 1633.8000000000002, 1637.0400000000002, 1640.2800000000002, 1643.5200000000002, 1646.7600000000002, 1650.0, 1653.24, 1656.48, 1659.72, 1662.96, 1666.2, 1669.44, 1672.68, 1675.92, 1679.16, 1682.4, 1685.64, 1688.88, 1692.1200000000001, 1695.3600000000001, 1698.6000000000001, 1701.8400000000001, 1705.0800000000002, 1708.3200000000002, 1711.5600000000002, 1714.8000000000002, 1718.0400000000002, 1721.2800000000002, 1724.5200000000002, 1727.7600000000002, 1731.0, 1734.24, 1737.48, 1740.72, 1743.96, 1747.2, 1750.44, 1753.68, 1756.92, 1760.16, 1763.4, 1766.64, 1769.88, 1773.1200000000001, 1776.3600000000001, 1779.6000000000001, 1782.8400000000001, 1786.0800000000002, 1789.3200000000002, 1792.5600000000002, 1795.8000000000002, 1799.0400000000002, 1802.2800000000002, 1805.5200000000002, 1808.7600000000002, 1812.0000000000002, 1815.24, 1818.48, 1821.72, 1824.96, 1828.2, 1831.44, 1834.68, 1837.92, 1841.16, 1844.4, 1847.64, 1850.88, 1854.1200000000001, 1857.3600000000001, 1860.6000000000001, 1863.8400000000001, 1867.0800000000002, 1870.3200000000002, 1873.5600000000002, 1876.8000000000002, 1880.0400000000002, 1883.2800000000002, 1886.5200000000002, 1889.7600000000002, 1893.0000000000002, 1896.2400000000002, 1899.48, 1902.72, 1905.96, 1909.2, 1912.44, 1915.68, 1918.92, 1922.16, 1925.4, 1928.64, 1931.88, 1935.1200000000001, 1938.3600000000001, 1941.6000000000001, 1944.8400000000001, 1948.0800000000002, 1951.3200000000002, 1954.5600000000002, 1957.8000000000002, 1961.0400000000002, 1964.2800000000002, 1967.5200000000002, 1970.7600000000002, 1974.0000000000002, 1977.2400000000002, 1980.48, 1983.72, 1986.96, 1990.2, 1993.44, 1996.68, 1999.92, 2003.16, 2006.4, 2009.64, 2012.88, 2016.1200000000001, 2019.3600000000001, 2022.6000000000001, 2025.8400000000001, 2029.0800000000002, 2032.3200000000002, 2035.5600000000002, 2038.8000000000002, 2042.0400000000002, 2045.2800000000002, 2048.5200000000004, 2051.76, 2055.0, 2058.2400000000002, 2061.4800000000005, 2064.7200000000003, 2067.96, 2071.2, 2074.44, 2077.6800000000003, 2080.92, 2084.1600000000003, 2087.4, 2090.6400000000003, 2093.88, 2097.1200000000003, 2100.36, 2103.6000000000004, 2106.84, 2110.08, 2113.32, 2116.56, 2119.8, 2123.04, 2126.28, 2129.52, 2132.76, 2136.0, 2139.2400000000002, 2142.48, 2145.7200000000003, 2148.96, 2152.2000000000003, 2155.44, 2158.6800000000003, 2161.92, 2165.1600000000003, 2168.4, 2171.6400000000003, 2174.88, 2178.1200000000003, 2181.36, 2184.6000000000004, 2187.84, 2191.08, 2194.32, 2197.56, 2200.8, 2204.04, 2207.28, 2210.52, 2213.76, 2217.0, 2220.2400000000002, 2223.48, 2226.7200000000003, 2229.96, 2233.2000000000003, 2236.44, 2239.6800000000003, 2242.92, 2246.1600000000003, 2249.4, 2252.6400000000003, 2255.88, 2259.1200000000003, 2262.36, 2265.6000000000004, 2268.84, 2272.08, 2275.32, 2278.56, 2281.8, 2285.04, 2288.28, 2291.52, 2294.76, 2298.0, 2301.2400000000002, 2304.48, 2307.7200000000003, 2310.96, 2314.2000000000003, 2317.44, 2320.6800000000003, 2323.92, 2327.1600000000003, 2330.4, 2333.6400000000003, 2336.88, 2340.1200000000003, 2343.36, 2346.6000000000004, 2349.84, 2353.08, 2356.32, 2359.56, 2362.8, 2366.04, 2369.28, 2372.52, 2375.76, 2379.0, 2382.2400000000002, 2385.48, 2388.7200000000003, 2391.96, 2395.2000000000003, 2398.44, 2401.6800000000003, 2404.92, 2408.1600000000003, 2411.4, 2414.6400000000003, 2417.88, 2421.1200000000003, 2424.36, 2427.6000000000004, 2430.84, 2434.0800000000004, 2437.32, 2440.56, 2443.8, 2447.04, 2450.28, 2453.52, 2456.76, 2460.0, 2463.2400000000002, 2466.48, 2469.7200000000003, 2472.96, 2476.2000000000003, 2479.44, 2482.6800000000003, 2485.92, 2489.1600000000003, 2492.4, 2495.6400000000003, 2498.88, 2502.1200000000003, 2505.36, 2508.6000000000004, 2511.84, 2515.0800000000004, 2518.32, 2521.56, 2524.8, 2528.04, 2531.28, 2534.52, 2537.76, 2541.0, 2544.2400000000002, 2547.48, 2550.7200000000003, 2553.96, 2557.2000000000003, 2560.44, 2563.6800000000003, 2566.92, 2570.1600000000003, 2573.4, 2576.6400000000003, 2579.88, 2583.1200000000003, 2586.36, 2589.6000000000004, 2592.84, 2596.0800000000004, 2599.32, 2602.56, 2605.8, 2609.04, 2612.28, 2615.52, 2618.76, 2622.0, 2625.2400000000002, 2628.48, 2631.7200000000003, 2634.96, 2638.2000000000003, 2641.44, 2644.6800000000003, 2647.92, 2651.1600000000003, 2654.4, 2657.6400000000003, 2660.88, 2664.1200000000003, 2667.36, 2670.6000000000004, 2673.84, 2677.0800000000004, 2680.32, 2683.5600000000004, 2686.8, 2690.04, 2693.28, 2696.52, 2699.76, 2703.0, 2706.2400000000002, 2709.48, 2712.7200000000003, 2715.96, 2719.2000000000003, 2722.44, 2725.6800000000003, 2728.92, 2732.1600000000003, 2735.4, 2738.6400000000003, 2741.88, 2745.1200000000003, 2748.36, 2751.6000000000004, 2754.84, 2758.0800000000004, 2761.32, 2764.5600000000004, 2767.8, 2771.04, 2774.28, 2777.52, 2780.76, 2784.0, 2787.2400000000002, 2790.48, 2793.7200000000003, 2796.96, 2800.2000000000003, 2803.44, 2806.6800000000003, 2809.92, 2813.1600000000003, 2816.4, 2819.6400000000003, 2822.88, 2826.1200000000003, 2829.36, 2832.6000000000004, 2835.84, 2839.0800000000004, 2842.32, 2845.5600000000004, 2848.8, 2852.04, 2855.28, 2858.52, 2861.76, 2865.0, 2868.2400000000002, 2871.48, 2874.7200000000003, 2877.96, 2881.2000000000003, 2884.44, 2887.6800000000003, 2890.92, 2894.1600000000003, 2897.4, 2900.6400000000003, 2903.88, 2907.1200000000003, 2910.36, 2913.6000000000004, 2916.84, 2920.0800000000004, 2923.32, 2926.5600000000004, 2929.8, 2933.04, 2936.28, 2939.52, 2942.76, 2946.0, 2949.2400000000002, 2952.48, 2955.7200000000003, 2958.96, 2962.2000000000003, 2965.44, 2968.6800000000003, 2971.92, 2975.1600000000003, 2978.4, 2981.6400000000003, 2984.88, 2988.1200000000003, 2991.36, 2994.6000000000004, 2997.84, 3001.0800000000004, 3004.32, 3007.5600000000004, 3010.8, 3014.0400000000004, 3017.28, 3020.52, 3023.76, 3027.0, 3030.2400000000002, 3033.48, 3036.7200000000003, 3039.96, 3043.2000000000003, 3046.44, 3049.6800000000003, 3052.92, 3056.1600000000003, 3059.4, 3062.6400000000003, 3065.88, 3069.1200000000003, 3072.36, 3075.6000000000004, 3078.84, 3082.0800000000004, 3085.32, 3088.5600000000004, 3091.8, 3095.0400000000004, 3098.28, 3101.52, 3104.76, 3108.0, 3111.2400000000002, 3114.48, 3117.7200000000003, 3120.96, 3124.2000000000003, 3127.44, 3130.6800000000003, 3133.92, 3137.1600000000003, 3140.4, 3143.6400000000003, 3146.88, 3150.1200000000003, 3153.36, 3156.6000000000004, 3159.84, 3163.0800000000004, 3166.32, 3169.5600000000004, 3172.8, 3176.0400000000004, 3179.28, 3182.52, 3185.76, 3189.0, 3192.2400000000002, 3195.48, 3198.7200000000003, 3201.96, 3205.2000000000003, 3208.44, 3211.6800000000003, 3214.92, 3218.1600000000003, 3221.4, 3224.6400000000003, 3227.88, 3231.1200000000003, 3234.36, 3237.6000000000004, 3240.84, 3244.0800000000004, 3247.32, 3250.5600000000004, 3253.8, 3257.0400000000004, 3260.28, 3263.5200000000004, 3266.76, 3270.0]
                                               )  #
                            y5 = np.asarray([7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7747.3508470962, 7746.693989509354, 7745.707755991116, 7744.720385907332, 7743.731879258001, 7742.742236043124, 7741.751456262699, 7740.759539916729, 7739.766487005213, 7738.7722975281495, 7737.77697148554, 7736.780508877384, 7735.782909703682, 7734.784173964433, 7733.784301659638, 7732.783292789296, 7731.7811473534075, 7730.777865351974, 7729.773446784992, 7728.767891652465, 7727.761199954391, 7726.75337169077, 7725.744406861603, 7724.734305466891, 7723.72306750663, 7722.710692980824, 7721.697181889472, 7720.682534232573, 7719.6667500101275, 7718.649829222136, 7717.631771868597, 7716.612577949512, 7715.592247464881, 7714.570780414703, 7713.5481767989795, 7712.5244366177085, 7711.499559870892, 7710.473546558528, 7709.446396680619, 7708.418110237162, 7707.388687228159, 7706.35812765361, 7705.326431513515, 7704.2935988078725, 7703.259629536684, 7702.224523699949, 7701.188281297667, 7700.1509023298395, 7699.112386796465, 7698.072734697545, 7697.0319460330775, 7695.990020803063, 7694.946959007503, 7693.902760646397, 7692.857425719744, 7691.810954227544, 7690.763346169799, 7689.714601546506, 7688.664720357668, 7687.613702603282, 7686.561548283351, 7685.5082573978725, 7684.453829946849, 7683.398265930277, 7682.34156534816, 7681.283728200496, 7680.224754487286, 7679.164644208529, 7678.103397364226, 7677.0410139543765, 7675.97749397898, 7674.912837438038, 7673.847044331549, 7672.780114659513, 7671.712048421931, 7670.642845618803, 7669.572506250129, 7668.501030315907, 7667.42841781614, 7666.354668750826, 7665.279783119965, 7664.203760923559, 7663.126602161605, 7662.0483068341055, 7660.968874941059, 7659.888306482466, 7658.806601458327, 7657.7237598686415, 7656.639781713409, 7655.554666992631, 7654.468415706306, 7653.381027854435, 7652.292503437017, 7651.202842454053, 7650.1120449055425, 7649.0201107914845, 7647.927040111881, 7646.832832866731, 7645.737489056035, 7644.641008679792, 7643.543391738002, 7642.444638230667, 7641.344748157785, 7640.243721519356, 7639.141558315381, 7638.038258545859, 7636.933822210792, 7635.828249310177, 7634.721539844016, 7633.613693812309, 7632.504711215055, 7631.394592052255, 7630.283336323908, 7629.170944030015, 7628.057415170576, 7626.94274974559, 7625.8269477550575, 7624.710009198979, 7623.5919340773535, 7622.472722390181, 7621.352374137463, 7620.2308893191985, 7619.108267935388, 7617.98450998603, 7616.859615471126, 7615.733584390676, 7614.606416744678, 7613.478112533136, 7612.348671756045, 7611.21809441341, 7610.086380505227, 7608.953530031498, 7607.8195429922225, 7606.6844193874, 7605.548159217032, 7604.410762481118, 7603.272229179656, 7602.1325593126485, 7600.991752880094, 7599.849809881994, 7598.706730318347, 7597.562514189153, 7596.417161494413, 7595.270672234127, 7594.123046408294, 7592.974284016915, 7591.824385059989, 7590.673349537517, 7589.521177449498, 7588.367868795934, 7587.213423576822, 7586.057841792164, 7584.9011234419595, 7583.743268526209, 7582.584277044912, 7581.424148998068, 7580.262884385678, 7579.1004832077415, 7577.936945464258, 7576.772271155229, 7575.606460280653, 7574.43951284053, 7573.271428834862, 7572.102208263646, 7570.931851126885, 7569.760357424577, 7568.587727156722, 7567.413960323322, 7566.239056924374, 7565.06301695988, 7563.88584042984, 7562.707527334253, 7561.52807767312, 7560.3474914464405, 7559.165768654214, 7557.982909296442, 7556.798913373123, 7555.613780884258, 7554.427511829846, 7553.240106209887, 7552.051564024382, 7550.861885273332, 7549.671069956734, 7548.47911807459, 7547.286029626899, 7546.091804613662, 7544.896443034879, 7543.699944890549, 7542.502310180673, 7541.30353890525, 7540.10363106428, 7538.902586657765, 7537.700405685703, 7536.497088148095, 7535.292634044939, 7534.087043376238, 7532.88031614199, 7531.672452342196, 7530.463451976855, 7529.2533150459685, 7528.042041549535, 7526.829631487555, 7525.616084860028, 7524.401401666955, 7523.185581908336, 7521.96862558417, 7520.750532694458, 7519.5313032392, 7518.310937218394, 7517.089434632043, 7515.866795480145, 7514.6430197627005, 7513.418107479709, 7512.192058631172, 7510.964873217088, 7509.736551237458, 7508.507092692282, 7507.276497581558, 7506.044765905289, 7504.811897663472, 7503.57789285611, 7502.342751483201, 7501.106473544746, 7499.869059040744, 7498.6305079711965, 7497.390820336102, 7496.14999613546, 7494.9080353692725, 7493.664938037538, 7492.420704140259, 7491.175333677432, 7489.9288266490585, 7488.681183055139, 7487.432402895673, 7486.182486170659, 7484.9314328801, 7483.679243023995, 7482.425916602343, 7481.171453615145, 7479.9158540624, 7478.659117944108, 7477.401245260271, 7476.142236010886, 7474.882090195956, 7473.620807815479, 7472.358388869456, 7471.094833357885, 7469.83014128077, 7468.5643126381065, 7467.297347429897, 7466.029245656141, 7464.76000731684, 7463.489632411991, 7462.218120941596, 7460.945472905654, 7459.671688304166, 7458.396767137132, 7457.120709404551, 7455.843515106424, 7454.56518424275, 7453.28571681353, 7452.005112818763, 7450.723372258451, 7449.4404951325905, 7448.156481441185, 7446.871331184233, 7445.585044361734, 7444.297620973689, 7443.009061020097, 7441.719364500959, 7440.428531416274, 7439.136561766043, 7437.843455550266, 7436.549212768942, 7435.253833422072, 7433.957317509655, 7432.6596650316915, 7431.360875988182, 7430.060950379126, 7428.759888204523, 7427.457689464374, 7426.154354158679, 7424.849882287437, 7423.544273850649, 7422.237528848314, 7420.929647280433, 7419.620629147005, 7418.310474448031, 7416.9991831835105, 7415.686755353444, 7414.37319095783, 7413.058489996671, 7411.742652469964, 7410.425678377711, 7409.107567719912, 7407.7883204965665, 7406.467936707674, 7405.1464163532355, 7403.823759433251, 7402.49996594772, 7401.175035896642, 7399.848969280018, 7398.521766097847, 7397.19342635013, 7395.863950036866, 7394.5333371580555, 7393.201587713699, 7391.868701703796, 7390.534679128347, 7389.199519987351, 7387.863224280809, 7386.525792008721, 7385.187223171085, 7383.847517767904, 7382.506675799175, 7381.164697264901, 7379.82158216508, 7378.4773304997125, 7377.131942268799, 7375.785417472339, 7374.437756110332, 7373.088958182779, 7371.739023689679, 7370.387952631033, 7369.03574500684, 7367.682400817102, 7366.327920061816, 7364.972302740985, 7363.6155488546065, 7362.257658402682, 7360.898631385211, 7359.538467802193, 7358.177167653629, 7356.814730939519, 7355.451157659862, 7354.086447814659, 7352.720601403908, 7351.353618427613, 7349.98549888577, 7348.616242778381, 7347.245850105445, 7345.874320866963, 7344.501655062935, 7343.12785269336, 7341.752913758239, 7340.376838257571, 7338.999626191357, 7337.621277559596, 7336.241792362289, 7334.861170599435, 7333.4794122710355, 7332.096517377089, 7330.7124859175965, 7329.327317892557, 7327.941013301971, 7326.553572145839, 7325.164994424161, 7323.775280136935, 7322.384429284164, 7320.992441865847, 7319.599317881982, 7318.2050573325705, 7316.809660217614, 7315.41312653711, 7314.015456291059, 7312.616649479463, 7311.21670610232, 7309.8156261596305, 7308.413409651394, 7307.010056577612, 7305.605566938283, 7304.199940733408, 7302.793177962986, 7301.385278627018, 7299.976242725504, 7298.566070258443, 7297.1547612258355, 7295.742315627682, 7294.328733463981, 7292.914014734734, 7291.498159439941, 7290.081167579601, 7288.6630391537155, 7287.2437741622825, 7285.823372605304, 7284.401834482778, 7282.979159794706, 7281.555348541087, 7280.130400721923, 7278.704316337212, 7277.277095386954, 7275.848737871151, 7274.4192437898, 7272.988613142903, 7271.5568459304595, 7270.12394215247, 7268.689901808933, 7267.254724899851, 7265.818411425221, 7264.380961385045, 7262.9423747793235, 7261.502651608054, 7260.06179187124, 7258.6197955688785, 7257.17666270097, 7255.732393267516, 7254.286987268515, 7252.840444703968, 7251.392765573874, 7249.943949878235, 7248.493997617048, 7247.042908790315, 7245.590683398035, 7244.13732144021, 7242.682822916838, 7241.227187827919, 7239.770416173454, 7238.312507953442, 7236.853463167884, 7235.39328181678, 7233.869969268111, 7231.865976551591, 7229.303541288871, 7226.182663479878, 7222.50334312467, 7218.265580223204, 7213.469374775523, 7208.114726781569, 7202.201636241414, 7195.730103154987, 7188.700127522359, 7181.111709343459, 7172.964848618358, 7164.259545346969, 7154.995799529395, 7145.173611165534, 7134.792980255472, 7123.853906799151, 7112.356390796616, 7100.300432247808, 7087.6860311528, 7074.513187511504, 7060.781901324008, 7046.492172590282, 7039.197636125616, 7036.604274794377, 7034.009492000925, 7031.413287745259, 7028.815662027382, 7026.216614847292, 7023.616146204988, 7021.014256100471, 7018.410944533742, 7015.8062115048, 7013.200057013645, 7010.592481060277, 7007.983483644696, 7005.373064766903, 7002.761224426897, 7000.147962624676, 6997.533279360245, 6994.9171746336, 6992.299648444741, 6989.68070079367, 6987.060331680387, 6984.438541104891, 6981.815329067181, 6979.190695567258, 6976.564640605124, 6973.937164180776, 6971.308266294214, 6968.677946945441, 6966.046206134454, 6963.413043861255, 6960.778460125843, 6958.1424549282165, 6955.505028268379, 6952.866180146328, 6950.225910562063, 6947.584219515587, 6944.941107006897, 6942.296573035995, 6939.65061760288, 6937.0032407075505, 6934.35444235001, 6931.704222530256, 6929.052581248289, 6926.399518504109, 6923.745034297716, 6921.089128629111, 6918.431801498293, 6915.7730529052615, 6913.112882850018, 6910.451291332561, 6907.788278352891, 6905.1238439110075, 6902.457988006912, 6899.790710640604, 6897.122011812083, 6894.451891521348, 6891.780349768402, 6889.107386553242, 6886.433001875868, 6883.7571957362825, 6881.079968134484, 6878.401319070473, 6875.721248544249, 6873.039756555811, 6870.356843105161, 6867.672508192298, 6864.986751817222, 6862.2995739799335, 6859.610974680432, 6856.920953918718, 6854.229511694791, 6851.53664800865, 6848.842362860298, 6846.146656249732, 6843.449528176952, 6840.750978641961, 6838.051007644757, 6835.34961518534, 6832.646801263709, 6829.942565879866, 6827.23690903381, 6824.529830725542, 6821.821330955059, 6819.111409722365, 6816.400067027457, 6813.687302870337, 6810.973117251004, 6808.257510169457, 6805.540481625699, 6802.822031619727, 6800.102160151542, 6797.380867221144, 6794.658152828534, 6791.934016973711, 6789.208459656675, 6786.481480877425, 6783.753080635963, 6781.023258932289, 6778.292015766401, 6775.5593511383, 6772.825265047986, 6770.089757495461, 6767.352828480722, 6764.61447800377, 6761.874706064605, 6759.133512663227, 6756.390897799636, 6753.646861473832, 6750.901403685816, 6748.154524435587, 6745.406223723145, 6742.65650154849, 6739.905357911623, 6737.1527928125415, 6734.398806251247, 6731.643398227741, 6728.886568742022, 6726.1283177940895, 6723.368645383945, 6720.607551511586, 6717.845036177016, 6715.0810993802315, 6712.315741121235, 6709.548961400026, 6706.780760216603, 6704.0111375709685, 6701.2400934631205, 6698.467627893059, 6695.6937408607855, 6692.918432366299, 6690.141702409599, 6687.363550990687, 6684.583978109561, 6681.802983766223, 6679.0205679606715, 6676.236730692908, 6673.4514719629315, 6670.664791770741, 6667.876690116339, 6665.087166999724, 6662.296222420896, 6659.503856379854, 6656.7100688766, 6653.914859911134, 6651.118229483454, 6648.320177593561, 6645.520704241456, 6642.719809427137, 6639.917493150606, 6637.113755411861, 6634.308596210905, 6631.502015547735, 6628.694013422352, 6625.884589834757, 6623.073744784948, 6620.261478272927, 6617.447790298693, 6614.632680862245, 6611.816149963586, 6608.998197602712, 6606.178823779627, 6603.358028494328, 6600.535811746817, 6597.712173537093, 6594.887113865156, 6592.060632731005, 6589.232730134643, 6586.403406076067, 6583.572660555278, 6580.740493572277, 6577.906905127062, 6575.071895219635, 6572.235463849995, 6569.397611018142, 6566.558336724076, 6563.717640967797, 6560.875523749305, 6558.031985068601, 6555.1870249256835, 6552.340643320553, 6549.49284025321, 6546.643615723654, 6543.792969731885, 6540.940902277904, 6538.087413361709, 6535.232502983301, 6532.376171142681, 6529.518417839848, 6526.659243074801, 6523.798646847543, 6520.936629158071, 6518.073190006386, 6515.208329392489, 6512.342047316378, 6509.474343778054, 6506.6052187775185, 6503.734672314769, 6500.862704389808, 6497.989315002633, 6495.114504153245, 6492.238271841645, 6489.360618067831, 6486.481542831805, 6483.601046133565, 6480.719127973113, 6477.835788350449, 6474.951027265571, 6472.064844718479, 6469.177240709177, 6466.288215237661, 6463.397768303931, 6460.505899907988, 6457.612610049834, 6454.717898729466, 6451.821765946885, 6448.924211702091, 6446.025235995085, 6443.124838825866, 6440.223020194433, 6437.319780100788, 6434.41511854493, 6431.509035526859, 6428.6015310465755, 6425.692605104079, 6422.782257699369, 6419.870488832447, 6416.957298503312, 6414.042686711964, 6411.126653458403, 6408.209198742629, 6405.290322564642, 6402.370024924443, 6399.44830582203, 6396.525165257405, 6393.600603230567, 6390.674619741516, 6387.747214790252, 6384.818388376775, 6381.888140501085, 6378.956471163183, 6376.023380363067, 6373.088868100739, 6370.152934376199, 6367.215579189444, 6364.276802540478, 6361.336604429298, 6358.394984855905, 6355.451943820299, 6352.507481322481, 6349.56159736245, 6346.6142919402055, 6343.665565055748, 6340.715416709078, 6337.763846900196, 6334.8108556291, 6331.856442895792, 6328.900608700271, 6325.943353042536, 6322.984675922589, 6320.0245773404295, 6317.063057296056, 6314.100115789471, 6311.135752820672, 6308.16996838966, 6305.202762496436, 6302.234135140999, 6299.264086323348, 6296.292616043486, 6293.31972430141, 6290.3454110971215, 6287.36967643062, 6284.392520301905, 6281.4139427109785, 6278.433943657838, 6275.452523142485, 6272.469681164919, 6269.48541772514, 6266.499732823148, 6263.5126264589435, 6260.524098632526, 6257.534149343896, 6254.542778593053, 6251.549986379997, 6248.555772704728, 6245.560137567247, 6242.563080967552, 6239.564602905644, 6236.564703381524, 6233.56338239519, 6230.560639946645, 6227.556476035886, 6224.550890662914, 6221.543883827729, 6218.535455530331, 6215.525605770721, 6212.514334548898, 6209.501641864861, 6206.487527718613, 6203.4719921101505, 6200.455035039476, 6197.436656506588, 6194.416856511487, 6191.395635054174, 6188.372992134648, 6185.348927752909, 6182.323441908957, 6179.296534602791, 6176.2682058344135, 6173.238455603823, 6170.207283911019, 6167.174690756003, 6164.140676138773, 6161.105240059332, 6158.068382517677, 6155.030103513809, 6151.990403047728, 6148.949281119434, 6145.9067377289275, 6142.8627728762085, 6139.817386561276, 6136.770578784131, 6133.722349544773, 6130.672698843202, 6127.621626679418, 6124.569133053422, 6121.515217965212, 6118.45988141479, 6115.403123402155, 6112.344943927306, 6109.285342990246, 6106.224320590972, 6103.161876729485, 6100.098011405786, 6097.032724619873, 6093.966016371747, 6090.89788666141, 6087.828335488858, 6084.757362854094, 6081.6849687571175, 6078.611153197928, 6075.535916176525, 6072.45925769291, 6069.381177747082, 6066.301676339041, 6063.220753468787, 6060.13840913632, 6057.05464334164, 6053.969456084747, 6050.882847365641, 6047.794817184324, 6044.705365540793, 6041.614492435047, 6038.522197867091, 6035.428481836921, 6032.333344344538, 6029.236785389942, 6026.138804973134, 6023.0394030941125, 6019.9385797528785, 6016.836334949431, 6013.732668683771, 6010.627580955898, 6007.521071765813, 6004.413141113514, 6001.303788999003, 5998.193015422278, 5995.080820383341, 5991.967203882191, 5988.852165918828, 5985.735706493253, 5982.617825605464, 5979.498523255463, 5976.377799443248, 5973.255654168821, 5970.13208743218, 5967.007099233328, 5963.880689572261, 5960.752858448983, 5957.62360586349, 5954.4929318157865, 5951.360836305868, 5948.227319333739, 5945.092380899396, 5941.956021002839, 5938.81823964407, 5935.679036823089, 5932.538412539894, 5929.396366794486, 5926.252899586866, 5923.1080109170325, 5919.9617007849865, 5916.813969190727, 5913.6648161342555, 5910.514241615571, 5907.362245634673, 5904.208828191563, 5901.05398928624, 5897.897728918703, 5894.740047088954, 5891.580943796993, 5888.420419042817, 5885.25847282643, 5882.09510514783, 5878.930316007016, 5875.76410540399, 5872.59647333875, 5869.427419811298, 5866.256944821633, 5863.085048369756, 5859.911730455664, 5856.736991079361, 5853.560830240845, 5850.383247940115, 5847.204244177173, 5844.023818952019, 5840.84197226465, 5837.65870411507, 5834.474014503276, 5831.287903429269, 5828.100370893049, 5824.911416894618, 5821.721041433973, 5818.5292445111145, 5815.336026126044, 5812.1413862787595, 5808.945324969263, 5805.747842197554, 5802.548937963631, 5799.348612267496, 5796.146865109148, 5792.943696488587, 5789.739106405814, 5786.533094860826, 5783.325661853627, 5780.116807384215, 5776.906531452589, 5773.694834058751, 5770.4817152027, 5767.267174884436, 5764.0512131039595, 5760.833829861269, 5757.615025156367, 5754.394798989251, 5751.1731513599225, 5747.950082268382, 5744.725591714628, 5741.499679698662, 5738.272346220481, 5735.043591280089, 5731.8134148774825, 5728.581817012665, 5725.348797685633, 5722.114356896389, 5718.878494644932, 5715.641210931262, 5712.402505755379, 5709.162379117284, 5705.920831016975, 5702.677861454454, 5699.43347042972, 5696.187657942773, 5692.940423993612, 5689.69176858224, 5686.441691708654, 5683.190193372855, 5679.937273574844, 5676.682932314619, 5673.427169592182, 5670.169985407532, 5666.911379760669, 5663.651352651594, 5660.389904080304, 5657.127034046804, 5653.862742551089]
                                            )  #
                            ax.plot(temp5, y5, color='#000000',  # Brown: #8c564b  #
                                    linestyle='-.', linewidth=2., label='Degree: 2 | Segments = 3',
                                    zorder=4, alpha=0.8)"""  #
                            logger.debug("Added post-regression preview for property '%s'", prop_name)
                        except Exception as e:
                            logger.warning("Could not generate post-regression preview for '%s': %s",
                                           prop_name, e)
                except Exception as e:
                    logger.error("Error creating function for property '%s': %s", prop_name, e)
                    ax.text(0.5, 0.5, f"Error: {str(e)}", transform=ax.transAxes,
                            horizontalalignment='center', fontweight='bold',
                            bbox=dict(facecolor='red', alpha=0.2))
                    _y_value = 0.0
            # Add boundary lines and annotations
            """ax.axvline(x=lower_bound, color=colors['bounds'], linestyle='--',
                       alpha=0.6, linewidth=1.5, label='_nolegend_')
            ax.axvline(x=upper_bound, color=colors['bounds'], linestyle='--',
                       alpha=0.6, linewidth=1.5, label='_nolegend_')"""
            # Ensure _y_value is valid for annotations
            if _y_value is None or not np.isfinite(_y_value):
                try:
                    if hasattr(current_prop, 'subs') and hasattr(current_prop, 'evalf'):
                        _y_value = float(current_prop.subs(T, lower_bound).evalf())
                    else:
                        _y_value = float(current_prop) if hasattr(current_prop, '__float__') else 0.0
                except (ValueError, TypeError, AttributeError):
                    _y_value = 0.0
                    logger.warning("Could not determine y_value for annotations for property '%s'", prop_name)
            # Add boundary type annotations
            '''ax.text(lower_bound, _y_value, f' {lower_bound_type}',
                    verticalalignment='top', horizontalalignment='right',
                    fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3',
                              edgecolor=colors['bounds']))
            ax.text(upper_bound, _y_value, f' {upper_bound_type}',
                    verticalalignment='top', horizontalalignment='left',
                    fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3',
                              edgecolor=colors['bounds']))'''
            # Add regression info
            """if has_regression and degree is not None:
                ax.text(0.5, 0.98, f"Simplify: {simplify_type} | Degree: {degree} | Segments: {segments}",
                        transform=ax.transAxes, horizontalalignment='center',
                        fontweight='bold',
                        bbox=dict(facecolor='lightblue', alpha=0.8, boxstyle='round,pad=0.3'))"""  #
            # Add legend
            '''legend = ax.legend(loc='best', framealpha=0.9, fancybox=True,
                               shadow=True, edgecolor='gray')
            legend.get_frame().set_linewidth(1.2)'''
            # Add property to visualized set
            self.visualized_properties.add(prop_name)
            logger.info("Successfully visualized property: %s", prop_name)
        except Exception as e:
            logger.error("Unexpected error visualizing property '%s': %s", prop_name, e, exc_info=True)
            raise ValueError(f"Unexpected error in property {prop_name}: {e}")

    def save_property_plots(self) -> None:
        """Save plots only if visualization is enabled and plots exist."""
        if not self.is_enabled or not hasattr(self, 'fig') or self.fig is None:
            logger.debug("No plots to save - visualization disabled or no plots created")
            return
        try:
            if hasattr(self, 'fig') and self.fig is not None:
                material_type = self.parser.config[MATERIAL_TYPE_KEY]
                title = f"Material Properties: {self.parser.config[NAME_KEY]} ({material_type})"
                #self.fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
                try:
                    plt.tight_layout(rect=[0, 0.01, 1, 0.98], pad=1.0)
                except Exception as e:
                    logger.warning("tight_layout failed: %s. Using subplots_adjust as fallback", e)
                    plt.subplots_adjust(
                        left=0.08,  # Left margin
                        bottom=0.08,  # Bottom margin
                        right=0.92,  # Right margin
                        top=0.88,  # Top margin (leave space for subtitle)
                        hspace=0.8  # Height spacing between subplots
                    )
                material_name = self.parser.config[NAME_KEY].replace(' ', '_')
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename = f"{material_name}_properties_{timestamp}.png"
                filepath = self.plot_directory / filename
                # Save settings
                self.fig.savefig(
                    str(filepath),
                    dpi=600,  # High resolution
                    bbox_inches="tight",  # Cropping
                    facecolor='white',  # Clean background
                    edgecolor='none',  # No border
                    pad_inches=0.4  # Padding
                )
                # Also save as vector format for infinite scalability
                # svg_filepath = filepath.with_suffix('.svg')
                # self.fig.savefig(str(svg_filepath), format='svg', bbox_inches="tight",
                #                  facecolor='white', edgecolor='none', pad_inches=0.4)
                total_properties = sum(len(props) for props in self.parser.categorized_properties.values())
                visualized_count = len(self.visualized_properties)
                if visualized_count != total_properties:
                    logger.warning(
                        f"Not all properties visualized! "
                        f"Visualized: {visualized_count}, "
                        f"Total: {total_properties}"
                    )
                else:
                    logger.info(f"All properties ({total_properties}) visualized successfully.")
                logger.info(f"All property plots saved as {filepath}")
        finally:  # Always close the figure to prevent memory leaks
            if hasattr(self, 'fig') and self.fig is not None:
                plt.close(self.fig)
                self.fig = None
                logger.debug("Figure closed and memory cleaned up")
