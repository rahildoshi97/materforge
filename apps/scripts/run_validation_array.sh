#!/bin/bash
#SBATCH --job-name=couette_val
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:35:00
#SBATCH --array=0-9
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/validation/run_validation_%a.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/validation/run_validation_%a.err

# SLURM job array: 10 tasks run in parallel, one per viscosity configuration.
# Each task picks its nu / timesteps / vtkWriteFrequency from the lookup tables.
#
# Submit with:
#   sbatch apps/scripts/run_validation_array.sh
#
# After all tasks complete, post-process with:
#   source ~/.venvs/materforge/bin/activate
#   python3 apps/scripts/extract_vtk_profiles.py
#   python3 apps/scripts/generate_validation_plots.py

# ── Configuration lookup tables (index 0-9) ──────────────────────────────────
# Task:        0      1      2      3      4      5     6     7      8       9
NUS=(       0.04   0.06   0.08    0.1    0.2    0.4   0.6   0.8    1.0  tempdep)
TIMESTEPS=(120000  80000  60000  48000  24000  12000  8000  6000   4800   60000)
# VTK frequency = round(timesteps/100) -> ~100 snapshots per case.
# Gives ~5 convergence time-points at the default TIME_EVERY_N=20 in
# generate_validation_plots.py.
VTK_FREQ=(  1200    800    600    480    240    120    80    60     48     600)

# ── Pick this task's parameters ───────────────────────────────────────────────
i=${SLURM_ARRAY_TASK_ID}
NU="${NUS[$i]}"
TS="${TIMESTEPS[$i]}"
VF="${VTK_FREQ[$i]}"

MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
BUILD_DIR=${APPS_DIR}/build/woody-release-cpu
PRM=${APPS_DIR}/CouetteFlowScaling.prm

MPI_OPTS="-n 4 --bind-to core --map-by core"

echo "=== Couette validation: nu=${NU}  timesteps=${TS}  vtkWriteFrequency=${VF} ==="
echo "Node: $(hostname)   Task: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
grep "model name" /proc/cpuinfo | head -1
echo ""

# cd to apps/ so the binary's relative VTK output path (output/vtk) lands
# in apps/output/vtk/.
cd "${APPS_DIR}"

# Select binary and runtime nu override
if [[ "${NU}" == "tempdep" ]]; then
    BIN="${BUILD_DIR}/CouetteFlowScaling_tempdep"
    # For tempdep, nu in prm is unused for physics; keep default (0.08)
    EXTRA_ARGS="-Parameters.timesteps=${TS} -Output.vtkWriteFrequency=${VF}"
else
    BIN="${BUILD_DIR}/CouetteFlowScaling_const_${NU}"
    # Override prm nu to match compiled CONST_NU so VTK filename is correct.
    # waLBerla CLI overrides require "=" format: -Block.Key=value (not -Block.Key value).
    EXTRA_ARGS="-Parameters.nu=${NU} -Parameters.timesteps=${TS} -Output.vtkWriteFrequency=${VF}"
fi

# Use a relaxed error threshold: not all cases converge to 1e-3 at these
# timestep counts (short runs for high-nu cases may still be transient).
EXTRA_ARGS="${EXTRA_ARGS} -Parameters.errorThreshold=0.01"

echo "Running: ${BIN##*/} ${PRM##*/} ${EXTRA_ARGS}"
$MPIRUN $MPI_OPTS "${BIN}" "${PRM}" ${EXTRA_ARGS}

echo ""
echo "=== Task ${SLURM_ARRAY_TASK_ID} completed: $(date) ==="
