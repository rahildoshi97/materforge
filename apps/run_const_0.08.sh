#!/bin/bash
#SBATCH --job-name=couette_const
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive          # sole occupant of the node — no shared-cache interference
#SBATCH --time=02:30:00      # 5 trials × ~9 min each + margin
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/run_const_0.08.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/run_const_0.08.err

MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
BIN=/home/hpc/iwia/iwia133h/.local/materforge/apps/build/woody-release-cpu/CouetteFlowScaling_const_0.08
PRM=/home/hpc/iwia/iwia133h/.local/materforge/apps/CouetteFlowScaling.prm

echo "=== Couette Flow: Constant viscosity nu=0.08 (4 MPI, icx, 128x64x64, 5 trials) ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
grep "model name" /proc/cpuinfo | head -1
echo ""

# Bind each MPI rank to exactly one core in sequential order.
# With 4 ranks on a 32-core 2-socket Ice Lake node, --map-by core places
# ranks 0-3 on cores 0-3 of socket 0 — all within the same NUMA domain,
# consistent L3 reuse, deterministic memory placement.
MPI_OPTS="-n 4 --bind-to core --map-by core --report-bindings"

for trial in 1 2 3 4 5; do
    echo "--- Trial $trial / 5  ($(date)) ---"
    $MPIRUN $MPI_OPTS $BIN $PRM
    echo ""
done

echo "=== Job completed: $(date) ==="
