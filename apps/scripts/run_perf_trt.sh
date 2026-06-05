#!/bin/bash
#SBATCH --job-name=couette_trt
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive          # sole occupant of the node - no shared-cache interference
#SBATCH --time=03:30:00      # const + tempdep, 5 trials each x ~8-11 min + margin
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/performance/run_perf_trt.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/performance/run_perf_trt.err

# TRT-collision counterpart of run_perf_const.sh + run_perf_tempdep.sh.
# Both binaries are WRITE_VISCOSITY=OFF (344 B/cell tempdep, 336 B/cell const), so
# the const-vs-tempdep comparison is operator-matched to the SRT benchmark.
# Build them first with: bash apps/scripts/build_validation_binaries.sh  (SRT) and
# the TRT binaries via cmake -DCOLLISION_OP=TRT (see couette_performance_analysis.md).

APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
BUILD_DIR=${APPS_DIR}/build/woody-release-cpu
PRM=${APPS_DIR}/CouetteFlowScaling.prm

BIN_CONST=${BUILD_DIR}/CouetteFlowScaling_const_0.08_TRT
BIN_TEMPDEP=${BUILD_DIR}/CouetteFlowScaling_tempdep_perf_TRT

echo "=== Couette Flow TRT: const nu=0.08 vs temp-dependent (4 MPI, icx, 128x64x64, 5 trials each) ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
grep "model name" /proc/cpuinfo | head -1
echo ""

cd "${APPS_DIR}"
MPI_OPTS="-n 4 --bind-to core --map-by core --report-bindings"

for BIN in "${BIN_CONST}" "${BIN_TEMPDEP}"; do
    if [[ ! -x "${BIN}" ]]; then
        echo "ERROR: TRT binary not found: ${BIN}" >&2
        echo "       Build it with: cmake --preset woody-release-cpu -DCOLLISION_OP=TRT ... (see analysis doc)" >&2
        exit 1
    fi
    echo "######## Binary: ${BIN##*/} ########"
    for trial in 1 2 3 4 5; do
        echo "--- Trial $trial / 5  ($(date)) ---"
        $MPIRUN $MPI_OPTS "${BIN}" "${PRM}"
        echo ""
    done
done

echo "=== Job completed: $(date) ==="
