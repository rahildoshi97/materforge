#!/bin/bash
#SBATCH --job-name=couette_tempdep
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive          # sole occupant of the node - no shared-cache interference
#SBATCH --time=03:00:00      # 5 trials x ~11 min each + margin
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/performance/run_perf_tempdep.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/performance/run_perf_tempdep.err

APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
# Performance binary: WRITE_VISCOSITY=OFF (viscosity write omitted -> 344 B/cell).
# This isolates the temperature-field-read overhead vs the constant-viscosity case.
# The viscosity-write build (CouetteFlowScaling_tempdep) is for validation/VTK only.
BIN=${APPS_DIR}/build/woody-release-cpu/CouetteFlowScaling_tempdep_perf
PRM=${APPS_DIR}/CouetteFlowScaling.prm

echo "=== Couette Flow: Temperature-dependent viscosity (MaterForge) (4 MPI, icx, 128x64x64, 5 trials) ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
grep "model name" /proc/cpuinfo | head -1
echo ""

# cd to apps/ so any VTK output (if enabled) lands in apps/output/vtk/.
cd "${APPS_DIR}"

MPI_OPTS="-n 4 --bind-to core --map-by core --report-bindings"

for trial in 1 2 3 4 5; do
    echo "--- Trial $trial / 5  ($(date)) ---"
    $MPIRUN $MPI_OPTS $BIN $PRM
    echo ""
done

echo "=== Job completed: $(date) ==="
