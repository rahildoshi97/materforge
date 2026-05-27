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

## Layout

| File                            | Purpose                                          |
|---------------------------------|--------------------------------------------------|
| `CouetteFlowScaling.cpp`        | Main C++ application (LBM driver)                |
| `CouetteFlowSweeps.py`          | pystencils code-generation script (configure-time) |
| `CouetteFlowScaling.prm`        | Runtime parameters                               |
| `CouetteFlowMaterial.yaml`      | Material file consumed by MaterForge             |
| `build_validation_binaries.sh`  | Builds 10 binaries (9 const-nu + tempdep)        |
| `run_validation_array.sh`       | SLURM array job: 10-case validation sweep        |
| `run_validation.sh`             | Small 2-case validation job                      |
| `run_strong_scaling.sh`         | SLURM array job: strong-scaling sweep (1-16 ranks) |
| `run_const_0.08.sh`             | 5-trial performance run, constant nu             |
| `run_tempdep.sh`                | 5-trial performance run, temperature-dependent nu |
| `extract_vtk_profiles.py`       | Post-processing: VTK -> CSV                      |
| `generate_plots.py`             | Validation plots (per-case + summary)            |
| `generate_performance_plots.py` | Performance plots: MLUPS, timer breakdown, etc.  |
| `parse_scaling.py`              | Parses scaling-sweep logs -> CSV                 |
| `plot_scaling.py`               | Strong-scaling plots (MLUPS / speedup / efficiency) |
| `plot_material_demo.py`         | Demonstration figure for AISI 304 properties     |

## Reproducing the benchmark figures

Prerequisites: an active `materforge` Python virtualenv with `lbmpy`,
`pystencils`, `pystencilssfg`, `sweepgen` installed, and the waLBerla
submodule initialised (`git submodule update --init --recursive`).

```bash
# ── 1. Build all validation binaries (one-shot, on login node) ──────────────
source ~/.venvs/materforge/bin/activate
bash apps/build_validation_binaries.sh

# ── 2. Validation sweep (10 cases in parallel) ──────────────────────────────
sbatch apps/run_validation_array.sh

# When all tasks finish, post-process:
cd apps/
python3 extract_vtk_profiles.py
python3 generate_plots.py
# Outputs: couette_*.png, couette_validation_summary.csv

# ── 3. Performance comparison (5 trials each, const vs tempdep) ────────────
sbatch apps/run_const_0.08.sh
sbatch apps/run_tempdep.sh
python3 generate_performance_plots.py
# Outputs: perf_*.png, perf_data.csv

# ── 4. Strong-scaling sweep (1, 2, 4, 8, 16 MPI ranks on one icx node) ─────
sbatch apps/run_strong_scaling.sh
python3 parse_scaling.py     # -> scaling_data.csv
python3 plot_scaling.py      # -> perf_scaling_*.png

# ── 5. Materials demonstration figure (no SLURM job needed) ─────────────────
python3 plot_material_demo.py  # -> material_aisi304_properties.png
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
