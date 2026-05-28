#!/bin/bash
#SBATCH --job-name=couette_validation
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/validation/run_validation.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/validation/run_validation.err

# Validation run - VTK output enabled (vtkWriteFrequency > 0 in CouetteFlowScaling.prm).
# Runs const_0.08 then tempdep sequentially; both write to apps/output/vtk/ with
# distinct filenames (couette_flow_cpu_const_0.0800_128x64x64 and
# couette_flow_cpu_tempdep_128x64x64).  Not --exclusive: physics is
# deterministic and independent of co-resident jobs.

MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
BUILD_DIR=${APPS_DIR}/build/woody-release-cpu
PRM=${APPS_DIR}/CouetteFlowScaling.prm

MPI_OPTS="-n 4 --bind-to core --map-by core"

echo "=== Couette Flow Validation Run (VTK output) ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
grep "model name" /proc/cpuinfo | head -1
echo ""

# VTK output is written to apps/output/vtk/ (set in CouetteFlowScaling.prm).
# CWD must be apps/ so the relative path resolves correctly.
cd "${APPS_DIR}"

echo "--- Constant viscosity nu=0.08 ---"
$MPIRUN $MPI_OPTS ${BUILD_DIR}/CouetteFlowScaling_const_0.08 ${PRM}
echo ""

echo "--- Temperature-dependent viscosity (MaterForge) ---"
$MPIRUN $MPI_OPTS ${BUILD_DIR}/CouetteFlowScaling_tempdep ${PRM}
echo ""

echo "=== Validation run completed: $(date) ==="
