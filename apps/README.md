# MaterForge Applications

Demonstration applications integrating MaterForge with
[waLBerla](https://i10git.cs.fau.de/walberla/walberla),
[lbmpy](https://pypi.org/project/lbmpy/), and
[pystencils](https://pypi.org/project/pystencils/).

The headline example is a 3D thermal Couette flow benchmark used to
validate the integration and to quantify the wall-clock overhead of a
MaterForge-driven temperature-dependent viscosity vs. a baked-in constant
viscosity.

## License

All files in this directory (except the `walberla/` submodule) are licensed
under **GPL-3.0-or-later** (inherited from waLBerla). See [LICENSE](LICENSE).
MaterForge itself (`../src/materforge/`) is licensed under
**BSD-3-Clause** - see the repository root for details.

The `walberla/` subdirectory is a Git submodule containing waLBerla (GPLv3),
maintained separately at:
<https://i10git.cs.fau.de/walberla/walberla>

## Directory layout

```
apps/
├── CMakeLists.txt, CMakePresets.json, CMakeUserPresets.json
├── CouetteFlowScaling.cpp        # LBM driver
├── CouetteFlowScaling.prm        # runtime parameters
├── CouetteFlowSweeps.py          # pystencils code generator
├── CouetteFlowMaterial.yaml      # MaterForge material spec
├── README.md, LICENSE, woody-sweepgen-requirements.txt
│
├── walberla/                     # git submodule
├── build/                        # CMake build output (gitignored)
│
├── scripts/                      # all SLURM + post-process helpers
│   ├── build_validation_binaries.sh
│   ├── run_validation.sh, run_validation_array.sh
│   ├── run_strong_scaling.sh
│   ├── run_perf_const.sh, run_perf_tempdep.sh
│   ├── run_perf_cache.sh, run_likwid_cache.sh, run_vtune_cache.sh
│   ├── extract_vtk_profiles.py
│   ├── generate_validation_plots.py, generate_performance_plots.py
│   ├── parse_scaling.py, plot_scaling.py
│   └── plot_material_demo.py
│
├── output/                       # all generated artifacts (gitignored)
│   ├── vtk/                      # waLBerla VTK output
│   ├── profiles/                 # cf_cpu_mf*_dat_*.csv  (z-axis profiles)
│   ├── data/                     # perf_data.csv, scaling_data.csv,
│   │                             # couette_validation_summary.csv
│   ├── plots/
│   │   ├── validation/           # couette_*.png, validation_const_nu*.png,
│   │   │                         # validation_tempdep_*.png
│   │   ├── performance/          # perf_{mlups,wall_time,timer*,bandwidth,
│   │   │                         #       overhead_decomposition,op_counts,summary}.png
│   │   ├── scaling/              # perf_scaling_{mlups,speedup,efficiency}.png
│   │   └── material/             # material_aisi304_properties.png
│   └── profiling/                # perf/, likwid/, vtune/ raw outputs
│
└── logs/                         # SLURM + cmake logs (gitignored)
    ├── build/                    # cmake configure + make logs
    ├── validation/               # run_validation_*.{log,err}
    ├── performance/              # run_perf_*.{log,err}, perf_cache.{log,err}
    ├── scaling/                  # scaling_*.{log,err}
    └── profiling/                # vtune_cache.*, likwid_cache.*
```

All generated artifacts live under `output/`, all logs live under `logs/`.
Both directories (plus `build/`) are gitignored as a whole - the contents
are reproduced by re-running the scripts.

## Prerequisites

```bash
# Initialise waLBerla submodule (one-time)
git submodule update --init --recursive

# Python dependencies for code generation
pip install lbmpy pystencils pystencilssfg sweepgen

# Activate the materforge virtualenv before any cmake/python invocation
source ~/.venvs/materforge/bin/activate     # adjust path to your venv
```

The CMake presets expect `$env{VIRTUAL_ENV}/bin/python3` - an active
virtualenv is mandatory.

---

## Configure

The build is driven by a small set of CMake cache variables. The
`CMakePresets.json` file provides sensible defaults for the woody cluster
and a typical workstation; pass extra `-D` flags after the preset name to
override individual options.

### Cache variables

| Variable           | Default     | Purpose                                                       |
|--------------------|-------------|---------------------------------------------------------------|
| `TARGET_PLATFORM`  | `GPU-HIP`   | `CPU`, `GPU-CUDA`, or `GPU-HIP`                               |
| `USE_MATERFORGE`   | `ON`        | `ON` -> MaterForge `nu(T)` baked in; `OFF` -> literal `CONST_NU` |
| `CONST_NU`         | `0.16667`   | Constant kinematic viscosity baked into the sweep when `USE_MATERFORGE=OFF` |
| `COLLISION_OP`     | `SRT`       | LBM collision operator: `SRT`, `TRT`, or `MRT`               |
| `WRITE_VISCOSITY`  | `OFF`       | `ON` -> write `nu` field each step for VTK (validation only); `OFF` -> omit it (performance). Sets the generated `StreamCollide` constructor arity, which `CouetteFlowScaling.cpp` matches automatically. |

### Available presets

| Preset                | Platform     | Build type | Notes                                |
|-----------------------|--------------|------------|--------------------------------------|
| `local-debug-cpu`     | CPU          | Debug      | local workstation                    |
| `local-release-cpu`   | CPU          | Release    | local workstation                    |
| `local-debug-gpu`     | GPU-CUDA     | Debug      | local workstation with CUDA          |
| `local-release-gpu`   | GPU-CUDA     | Release    | local workstation with CUDA          |
| `woody-release-cpu`   | CPU          | Release    | Xeon Gold 6326 (Ice Lake) on woody   |

For LUMI (HIP), no preset exists - pass `-DTARGET_PLATFORM=GPU-HIP`
directly to `cmake`.

### Examples

```bash
# Temperature-dependent (MaterForge ON), local CPU release
cmake --preset local-release-cpu

# Constant viscosity nu=0.1 baked at compile time
cmake --preset local-release-cpu -DUSE_MATERFORGE=OFF -DCONST_NU=0.1

# Debug build
cmake --preset local-debug-cpu

# woody release with default preset
cmake --preset woody-release-cpu
```

---

## Build

```bash
# Build the unified executable (the preset's *-build target)
cmake --build --preset local-release-cpu-build

# Or directly via make
make -C build/local-release-cpu CouetteFlowScaling -j8
```

The single binary works on the configured platform; the
`CouetteFlowSweeps.py` sweep is regenerated automatically when its inputs
change.

### Building all validation binaries (one-shot)

For the 10-case validation sweep you need nine constant-`nu` binaries plus
the tempdep binary. The helper script reconfigures and rebuilds the
project once per `nu` value, then copies each binary to a named file:

```bash
bash apps/scripts/build_validation_binaries.sh
```

Produces, under `apps/build/woody-release-cpu/`:
- `CouetteFlowScaling_const_0.04` ... `CouetteFlowScaling_const_1.0` (all `WRITE_VISCOSITY=OFF`)
- `CouetteFlowScaling_tempdep` — `WRITE_VISCOSITY=ON`, for validation/VTK runs
- `CouetteFlowScaling_tempdep_perf` — `WRITE_VISCOSITY=OFF`, the performance-benchmark binary used by `run_perf_srt.sh` (and `run_perf_tempdep.sh`)

cmake/make logs land in `apps/logs/build/{configure,build}_<case>.log`.

---

## Run

### Single binary

```bash
# CWD must be apps/ so the VTK output path "output/vtk" (relative)
# resolves into apps/output/vtk/.
cd apps/
./build/local-release-cpu/CouetteFlowScaling CouetteFlowScaling.prm
```

The binary prints **MLUPS per process** and **total MLUPS** at the end of
the timed run - the primary performance metric.

### Command-line overrides

`.prm` parameters can be overridden from the CLI using waLBerla's
`-<Block>.<Key>=<Value>` syntax (the `=` is required). Vector3 values must
use the angle-bracket form `<x,y,z>`:

```bash
./build/local-release-cpu/CouetteFlowScaling CouetteFlowScaling.prm \
    -DomainSetup.blocks="<4,1,1>" \
    -DomainSetup.cellsPerBlock="<64,32,32>" \
    -Parameters.timesteps=10000 \
    -Parameters.nu=0.08 \
    -Output.vtkWriteFrequency=200 \
    -Output.vtkOutputDir=output/vtk
```

`Output.vtkWriteFrequency=0` disables VTK output; any non-zero value
enables it with the given write frequency.

### .prm reference (`CouetteFlowScaling.prm`)

| Block         | Key                | Default     | Purpose                                        |
|---------------|--------------------|-------------|------------------------------------------------|
| `DomainSetup` | `blocks`           | `<1,2,2>`   | MPI process grid (Vector3)                     |
| `DomainSetup` | `cellsPerBlock`    | `<128,32,32>` | Cells per block (Vector3)                    |
| `DomainSetup` | `periodic`         | `<1,1,0>`   | Per-axis periodicity                           |
| `Parameters`  | `nu`               | `0.08`      | Initial fill / VTK label only (physics baked in) |
| `Parameters`  | `u_max`            | `0.025`     | Wall velocity (lattice units)                  |
| `Parameters`  | `timesteps`        | `60000`     | Total simulation steps                         |
| `Parameters`  | `errorThreshold`   | `1e-3`      | Convergence tolerance for steady-state check   |
| `Parameters`  | `T_bottom`         | `300.0`     | Bottom-wall temperature (K)                    |
| `Parameters`  | `T_top`            | `3000.0`    | Top-wall temperature (K)                       |
| `Output`      | `vtkWriteFrequency`| `0`         | Steps between VTK writes (0 = disabled)        |
| `Output`      | `vtkOutputDir`     | `output/vtk`| VTK output directory (relative to CWD)         |

> When `USE_MATERFORGE=OFF`, `Parameters.nu` must match the `CONST_NU` you
> built with - the `nu` field is only used for the VTK output filename and
> for the initial fill of the viscosity scalar field. The physics
> viscosity is the literal baked into the generated sweep.

### Validation sweep (SLURM)

10-case sweep (9 constant-`nu` + tempdep), one task per case:

```bash
sbatch apps/scripts/run_validation_array.sh
```

Per-task SBATCH parameters: 4 MPI ranks on a single node, ~35 min
walltime. Each task selects its `nu` / `timesteps` / `vtkWriteFrequency`
from a lookup table inside the script. Logs land in
`apps/logs/validation/run_validation_<task>.log`.

### Performance comparison (SLURM)

The const-vs-tempdep MLUPS overhead is < 10 %, which is smaller than the
node-to-node hardware variation across woody's Ice Lake pool. **Both cases must
therefore run on the same physical node**, otherwise the comparison is confounded
by node binning. Use the combined single-allocation job:

```bash
sbatch apps/scripts/run_perf_srt.sh        # const_0.08 + tempdep, same node, 5 trials each
```

It runs both binaries back-to-back inside one `--exclusive --constraint=icx`
allocation and redirects each case's output to `run_perf_const.log` /
`run_perf_tempdep.log` (the files `generate_performance_plots.py` parses), so the
post-processing is unchanged. The TRT counterpart is `run_perf_trt.sh`.

> The standalone `run_perf_const.sh` and `run_perf_tempdep.sh` are kept for
> single-case reruns, but submitting them as two independent jobs does **not**
> guarantee the same node — prefer `run_perf_srt.sh` for the headline comparison.

Logs land in `apps/logs/performance/`.

### Strong-scaling sweep (SLURM)

5-task array job, one task per rank count (1, 2, 4, 8, 16) inside a single
Xeon Gold 6326 node:

```bash
sbatch apps/scripts/run_strong_scaling.sh
```

Each task runs three trials for both binaries. Logs land in
`apps/logs/scaling/scaling_<task>.log`.

---

## Post-processing

### VTK extraction

Reads `.pvd` files from `apps/output/vtk/`, extracts a z-axis profile at
the domain midpoint for every timestep, and writes one CSV per simulation
into `apps/output/profiles/`:

```bash
# All simulations
python3 apps/scripts/extract_vtk_profiles.py

# One simulation by name fragment
python3 apps/scripts/extract_vtk_profiles.py "const_0.08"
python3 apps/scripts/extract_vtk_profiles.py "tempdep"

# Custom paths
python3 apps/scripts/extract_vtk_profiles.py \
    --vtk-dir path/to/vtk --out-dir path/to/csvs/
```

Output filenames:
- `cf_{cpu|gpu}_mfconst_<nu>_dat_<last_idx>_<last_step>.csv`
- `cf_{cpu|gpu}_mftempdep_dat_<last_idx>_<last_step>.csv`

### Validation plots

```bash
python3 apps/scripts/generate_validation_plots.py
```

Reads all `cf_*.csv` in `apps/output/profiles/`, produces per-case
4-panel + 2-panel + time-evolution plots in
`apps/output/plots/validation/`, and writes
`apps/output/data/couette_validation_summary.csv` (`L2`, `L_infinity`,
max/mean relative errors per case).

### Performance plots

```bash
python3 apps/scripts/generate_performance_plots.py
```

Parses the 5-trial logs from `apps/logs/performance/` (falls back to the
hardcoded `FALLBACK_DATA` block if logs are missing), produces:
- `perf_mlups.png`, `perf_wall_time.png` - throughput + wall time, 95 % CI
- `perf_timer_breakdown{,_pct}.png` - StreamCollide / Comms / BCs
- `perf_overhead_decomposition.png` - estimated vs measured overhead
- `perf_op_counts.png` - symbolic operation counts
- `perf_bandwidth.png` - algorithmic bandwidth vs hardware peak
- `perf_summary.png` - 6-panel summary
- `perf_data.csv` (under `apps/output/data/`) with Welch t-test summary

### Scaling plots

```bash
python3 apps/scripts/parse_scaling.py     # -> apps/output/data/scaling_data.csv
python3 apps/scripts/plot_scaling.py      # -> apps/output/plots/scaling/
```

Produces `perf_scaling_mlups.png`, `perf_scaling_speedup.png`,
`perf_scaling_efficiency.png` with 95 % CI error bars over the three
trials per (rank, configuration).

### Material demonstration figure

```bash
python3 apps/scripts/plot_material_demo.py
```

Loads `src/materforge/data/materials/1.4301.yaml` via the MaterForge API
and writes a 2x4 panel figure of AISI 304 properties to
`apps/output/plots/material/material_aisi304_properties.png`.

---

## Hardware-counter profiling (optional)

Three SLURM scripts collect raw hardware-counter data for the
cache/memory-bandwidth analysis used in the paper's
overhead-decomposition figure. All write raw output to
`apps/output/profiling/<tool>/` and SLURM logs to
`apps/logs/profiling/`.

```bash
# Linux perf-stat (cache + memory counters per MPI rank)
sbatch apps/scripts/run_perf_cache.sh
# -> apps/output/profiling/perf/{const_0.08,tempdep}/rank_{0..3}.txt

# LIKWID likwid-perfctr (CACHE + MEM groups)
sbatch apps/scripts/run_likwid_cache.sh
# -> apps/output/profiling/likwid/{const_0.08,tempdep}_{CACHE,MEM}.txt

# Intel VTune (opt-in; the launcher is gitignored - see below)
sbatch apps/scripts/run_vtune_cache.sh
# -> apps/output/profiling/vtune/...
```

`apps/scripts/run_vtune_cache.sh` is gitignored because the VTune launcher
reads setuid `/apps/` paths that are not part of the public reproduction
recipe. Add it back locally if you have a working VTune installation.

---

## Hardware

Performance numbers in the paper were measured on the **woody** cluster at
NHR@FAU:

- **Validation sweep:** mixed E3-1240 v6 / Xeon Gold 6326 nodes
  (4 ranks per task, single node).
- **Performance comparison & scaling:** Intel Xeon Gold 6326 (Ice Lake,
  32 physical cores, 2.9 GHz), DDR4-3200 8-channel
  (~204.8 GB/s/socket peak).

---

## End-to-end reproduction

```bash
# 1. Configure & build all validation binaries (~8-10 min on a login node)
source ~/.venvs/materforge/bin/activate
bash apps/scripts/build_validation_binaries.sh

# 2. Validation sweep (10 tasks, ~25 min walltime)
sbatch apps/scripts/run_validation_array.sh
# ... wait for completion ...
python3 apps/scripts/extract_vtk_profiles.py
python3 apps/scripts/generate_validation_plots.py

# 3. Performance comparison (const + tempdep, same node, 5 trials each, ~1.5 h walltime)
sbatch apps/scripts/run_perf_srt.sh
# ... wait for completion ...
python3 apps/scripts/generate_performance_plots.py

# 4. Strong-scaling sweep (~1 h walltime)
sbatch apps/scripts/run_strong_scaling.sh
# ... wait for completion ...
python3 apps/scripts/parse_scaling.py
python3 apps/scripts/plot_scaling.py

# 5. Material demonstration figure (~5 s, no SLURM needed)
python3 apps/scripts/plot_material_demo.py
```

All figures land under `apps/output/plots/` ready for inclusion in the
paper (`../paper/`).
