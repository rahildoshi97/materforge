#!/bin/bash
#SBATCH --job-name=cf_scaling
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --array=0-4
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/scaling/scaling_%a.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/scaling/scaling_%a.err

# Strong-scaling sweep on a single Xeon Gold 6326 (Ice Lake, 32 cores) node.
# Each array task pins the rank count from the RANKS table and runs 3 trials
# for both the constant-viscosity binary and the temperature-dependent
# (MaterForge) binary using the same fixed 128x64x64 domain.
#
# Submit:
#   sbatch apps/scripts/run_strong_scaling.sh
#
# After all tasks complete, post-process with:
#   python3 apps/scripts/parse_scaling.py
#   python3 apps/scripts/plot_scaling.py

# ── Configuration ────────────────────────────────────────────────────────────
# Vector3<uint_t> CLI overrides must use waLBerla's "<x,y,z>" syntax
# (parsed via operator>>( istream&, Vector3& ) which expects angle brackets).
RANKS=(    1            2            4            8            16     )
BLOCKS=( "<1,1,1>"    "<1,2,1>"    "<1,2,2>"    "<2,2,2>"    "<4,2,2>" )
CPB=(    "<128,64,64>" "<128,32,64>" "<128,32,32>" "<64,32,32>" "<32,32,32>" )

i=${SLURM_ARRAY_TASK_ID}
N="${RANKS[$i]}"
B="${BLOCKS[$i]}"
C="${CPB[$i]}"

TIMESTEPS=10000
TRIALS=3

MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
BUILD_DIR=${APPS_DIR}/build/woody-release-cpu
PRM=${APPS_DIR}/CouetteFlowScaling.prm
MPI_OPTS=(-n "${N}" --bind-to core --map-by core)

echo "=== Strong scaling: ranks=${N}  blocks=${B}  cellsPerBlock=${C} ==="
echo "Node: $(hostname)   Task: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
grep "model name" /proc/cpuinfo | head -1
echo ""

cd "${APPS_DIR}"

OVERRIDES=(
    "-DomainSetup.blocks=${B}"
    "-DomainSetup.cellsPerBlock=${C}"
    "-Parameters.timesteps=${TIMESTEPS}"
    "-Parameters.errorThreshold=0.01"
    "-Output.vtkWriteFrequency=0"
)

for CASE in const_0.08 tempdep; do
    BIN="${BUILD_DIR}/CouetteFlowScaling_${CASE}"
    EXTRA=()
    if [[ "${CASE}" != "tempdep" ]]; then
        EXTRA+=("-Parameters.nu=0.08")
    fi
    for ((t=1; t<=TRIALS; t++)); do
        echo ""
        echo "── ${CASE}  trial ${t}/${TRIALS}  ranks=${N} ──"
        "${MPIRUN}" "${MPI_OPTS[@]}" "${BIN}" "${PRM}" "${OVERRIDES[@]}" "${EXTRA[@]}"
    done
done

echo ""
echo "=== Task ${SLURM_ARRAY_TASK_ID} completed: $(date) ==="
