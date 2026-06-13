#!/bin/bash
#SBATCH --job-name=couette_srt
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive          # sole occupant of the node - no shared-cache interference
#SBATCH --time=02:30:00      # const + tempdep, 5 trials each x ~7-8 min + margin
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/performance/run_perf_srt.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/performance/run_perf_srt.err

# Same-node SRT performance comparison: const nu=0.08 vs MaterForge temp-dependent.
#
# WHY ONE JOB:  the const-vs-tempdep overhead is < 10 %, smaller than the node-to-node
# hardware variation on woody's Ice Lake pool (turbo/thermal/memory binning).  Submitting
# run_perf_const.sh and run_perf_tempdep.sh as two independent jobs does NOT pin them to the
# same node - SLURM is free to place them on different nodes under --constraint=icx, which
# confounds the comparison.  Running both binaries inside a single --exclusive allocation
# guarantees the same physical CPU for both cases, exactly like run_perf_trt.sh does for TRT.
#
# Each binary's full output is redirected to the per-case log that
# generate_performance_plots.py already parses (run_perf_const.log / run_perf_tempdep.log),
# so no downstream change is needed.  Both binaries are WRITE_VISCOSITY=OFF
# (344 B/cell tempdep, 336 B/cell const) so the comparison is traffic-matched.
#
# Build the binaries first with: bash apps/scripts/build_validation_binaries.sh

APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
BUILD_DIR=${APPS_DIR}/build/woody-release-cpu
LOG_DIR=${APPS_DIR}/logs/performance
PRM=${APPS_DIR}/CouetteFlowScaling.prm

# Performance binary: tempdep uses WRITE_VISCOSITY=OFF (the _perf build) so the
# const-vs-tempdep MLUPS comparison isolates the temperature-field-read overhead.
BIN_CONST=${BUILD_DIR}/CouetteFlowScaling_const_0.08
BIN_TEMPDEP=${BUILD_DIR}/CouetteFlowScaling_tempdep_perf

LOG_CONST=${LOG_DIR}/run_perf_const.log
LOG_TEMPDEP=${LOG_DIR}/run_perf_tempdep.log

mkdir -p "${LOG_DIR}"
cd "${APPS_DIR}"

# Guard: the perf binaries are produced by build_validation_binaries.sh. Fail with a
# clear message instead of an opaque mpirun error if one is missing.
for BIN in "${BIN_CONST}" "${BIN_TEMPDEP}"; do
    if [[ ! -x "${BIN}" ]]; then
        echo "ERROR: benchmark binary not found: ${BIN}" >&2
        echo "       Build it first with: bash apps/scripts/build_validation_binaries.sh" >&2
        exit 1
    fi
done

# Bind each MPI rank to one core on socket 0 (cores 0-3) - same NUMA domain, consistent L3.
MPI_OPTS="-n 4 --bind-to core --map-by core --report-bindings"

echo "=== Couette Flow SRT: const nu=0.08 vs temp-dependent (same node) ==="
echo "Node: $(hostname)   Job: ${SLURM_JOB_ID:-interactive}   Date: $(date)"
grep "model name" /proc/cpuinfo | head -1

# ── Constant viscosity nu=0.08 ───────────────────────────────────────────────
{
    echo "=== Couette Flow: Constant viscosity nu=0.08 (4 MPI, icx, 128x64x64, 5 trials) ==="
    echo "Node: $(hostname)   Job: ${SLURM_JOB_ID:-interactive}"
    echo "Date: $(date)"
    grep "model name" /proc/cpuinfo | head -1
    echo ""
    for trial in 1 2 3 4 5; do
        echo "--- Trial $trial / 5  ($(date)) ---"
        $MPIRUN $MPI_OPTS "${BIN_CONST}" "${PRM}"
        echo ""
    done
    echo "=== const job completed: $(date) ==="
} > "${LOG_CONST}" 2>&1

# ── Temperature-dependent viscosity (MaterForge) ─────────────────────────────
{
    echo "=== Couette Flow: Temperature-dependent viscosity (MaterForge) (4 MPI, icx, 128x64x64, 5 trials) ==="
    echo "Node: $(hostname)   Job: ${SLURM_JOB_ID:-interactive}"
    echo "Date: $(date)"
    grep "model name" /proc/cpuinfo | head -1
    echo ""
    for trial in 1 2 3 4 5; do
        echo "--- Trial $trial / 5  ($(date)) ---"
        $MPIRUN $MPI_OPTS "${BIN_TEMPDEP}" "${PRM}"
        echo ""
    done
    echo "=== tempdep job completed: $(date) ==="
} > "${LOG_TEMPDEP}" 2>&1

echo "=== Job completed: $(date) ==="
echo "Per-case logs: ${LOG_CONST}"
echo "               ${LOG_TEMPDEP}"
echo "Post-process:  python3 apps/scripts/generate_performance_plots.py"
