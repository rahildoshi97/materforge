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

## Reproducing the benchmark figures

Prerequisites: an active `materforge` Python virtualenv with `lbmpy`,
`pystencils`, `pystencilssfg`, `sweepgen` installed, and the waLBerla
submodule initialised (`git submodule update --init --recursive`).

```bash
# ── 1. Build all validation binaries (one-shot, on login node) ──────────────
source ~/.venvs/materforge/bin/activate
bash apps/scripts/build_validation_binaries.sh

# ── 2. Validation sweep (10 cases in parallel) ──────────────────────────────
sbatch apps/scripts/run_validation_array.sh

# When all tasks finish, post-process:
python3 apps/scripts/extract_vtk_profiles.py
python3 apps/scripts/generate_validation_plots.py
# Outputs land under apps/output/plots/validation/ and apps/output/data/.

# ── 3. Performance comparison (5 trials each, const vs tempdep) ────────────
sbatch apps/scripts/run_perf_const.sh
sbatch apps/scripts/run_perf_tempdep.sh
python3 apps/scripts/generate_performance_plots.py
# Outputs land under apps/output/plots/performance/ and apps/output/data/.

# ── 4. Strong-scaling sweep (1, 2, 4, 8, 16 MPI ranks on one icx node) ─────
sbatch apps/scripts/run_strong_scaling.sh
python3 apps/scripts/parse_scaling.py       # -> apps/output/data/scaling_data.csv
python3 apps/scripts/plot_scaling.py        # -> apps/output/plots/scaling/

# ── 5. Materials demonstration figure (no SLURM job needed) ─────────────────
python3 apps/scripts/plot_material_demo.py  # -> apps/output/plots/material/
```

For build/run details on a local workstation (CPU debug/release, GPU CUDA),
see `../CLAUDE.md`.

## Hardware

Performance numbers in the paper were measured on the **woody** cluster at
NHR@FAU:

- **Validation sweep:** mixed E3-1240 v6 / Xeon Gold 6326 nodes (4 ranks
  per task, single node).
- **Performance comparison & scaling:** Intel Xeon Gold 6326 (Ice Lake,
  32 physical cores, 2.9 GHz), DDR4-3200 8-channel
  (~204.8 GB/s/socket peak).
