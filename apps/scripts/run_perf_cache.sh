#!/bin/bash
#SBATCH --job-name=couette_perf
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/profiling/perf_cache.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/profiling/perf_cache.err

# =============================================================================
#  perf stat cache counter analysis — Couette Flow
#
#  Third approach after VTune (paranoid=2) and LIKWID (MSR locked) both failed.
#  perf stat counting mode uses PERF_EVENT_OPEN via the kernel's own PMU driver,
#  which works at perf_event_paranoid=2 for user-space per-process events.
#
#  Per-rank measurement: each MPI rank runs its own perf stat instance and
#  writes to a separate file. All 4 ranks run the full simulation — we collect
#  rank 0's counters as representative (symmetric decomposition).
# =============================================================================

APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
BUILD_DIR=$APPS_DIR/build/woody-release-cpu
LOGS_DIR=$APPS_DIR/logs/profiling
RESULTS_DIR=$APPS_DIR/output/profiling/perf
PRM=$APPS_DIR/CouetteFlowScaling.prm

# Collision operator selector. COLLISION_LABEL="" profiles the SRT binaries
# (default); COLLISION_LABEL="_TRT" profiles the TRT binaries. Override at submit:
#   sbatch --export=ALL,COLLISION_LABEL=_TRT --output=.../perf_cache_trt.log run_perf_cache.sh
COLLISION_LABEL="${COLLISION_LABEL:-}"
# Profile the performance binary (WRITE_VISCOSITY=OFF, 344 B/cell) so the cache/IPC
# figures match the canonical const-vs-tempdep MLUPS benchmark (run_perf_tempdep.sh).
BIN_CONST=$BUILD_DIR/CouetteFlowScaling_const_0.08${COLLISION_LABEL}
BIN_TEMPDEP=$BUILD_DIR/CouetteFlowScaling_tempdep_perf${COLLISION_LABEL}

MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun

TIMESTEPS=5000

# ── Find perf ─────────────────────────────────────────────────────────────────
PERF=""
for candidate in \
    /usr/bin/perf \
    /usr/lib64/perf \
    /usr/local/bin/perf \
    /bin/perf; do
    if [[ -x "$candidate" ]]; then
        PERF="$candidate"
        break
    fi
done
if [[ -z "$PERF" ]]; then
    echo "ERROR: perf binary not found. Checked: /usr/bin/perf, /usr/lib64/perf"
    echo "  Install with: sudo yum install perf   (on AlmaLinux 8)"
    exit 1
fi

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

echo "=============================================================="
echo "  perf stat Cache Analysis — Couette Flow"
echo "  Node:     $(hostname)"
echo "  Date:     $(date)"
echo "  CPU:      $(grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)"
echo "  perf:     $PERF ($($PERF --version 2>/dev/null))"
echo "  paranoid: $(cat /proc/sys/kernel/perf_event_paranoid)"
echo "  Timesteps per run: $TIMESTEPS"
echo "=============================================================="

# ── Probe which hardware events are accessible ─────────────────────────────────
# At paranoid=2, generic hardware cache events typically work for user-space.
echo ""
echo "── Probing available hardware cache events ──────────────────"

probe_event() {
    local EVT=$1
    # Try to count the event for 0.5 seconds on an empty workload
    timeout 3 $PERF stat -e "$EVT" -- sleep 0.1 2>&1 | grep -qE "<not supported>|<not counted>|error|failed|unsupported" \
        && echo "SKIP $EVT" \
        || echo "OK   $EVT"
}

GENERIC_EVENTS=(
    "cache-references"
    "cache-misses"
    "L1-dcache-loads"
    "L1-dcache-load-misses"
    "L1-icache-loads"
    "L1-icache-load-misses"
    "LLC-loads"
    "LLC-load-misses"
    "LLC-stores"
    "dTLB-loads"
    "dTLB-load-misses"
    "cycles"
    "instructions"
)

AVAILABLE_EVENTS=()
for evt in "${GENERIC_EVENTS[@]}"; do
    result=$(probe_event "$evt")
    echo "  $result"
    if [[ "$result" == OK* ]]; then
        AVAILABLE_EVENTS+=("$evt")
    fi
done

if [[ ${#AVAILABLE_EVENTS[@]} -eq 0 ]]; then
    echo ""
    echo "ERROR: No hardware cache events available."
    echo "  All three profiling methods have been tried and failed on this system:"
    echo "  1. VTune  — perf_event_paranoid=2 blocks hardware-event sampling"
    echo "  2. LIKWID — MSR registers locked (perfmon_init error 2141)"
    echo "  3. perf   — No hardware events accessible"
    echo ""
    echo "  Resolution: Ask the NHR@FAU woody admins to either:"
    echo "    a) Set perf_event_paranoid=1 on icx nodes  (recommended)"
    echo "    b) Enable Intel VTune sampling drivers"
    echo "    c) Unlock MSR access for the LIKWID daemon"
    echo "  Contact: hpc-support@fau.de"
    exit 1
fi

EVENT_STR=$(IFS=,; echo "${AVAILABLE_EVENTS[*]}")
echo ""
echo "  Using events: $EVENT_STR"

# ── Helper: run perf stat for one case ────────────────────────────────────────
# Strategy: wrap each MPI rank's binary with perf stat.
# mpirun passes OMPI_COMM_WORLD_RANK per process; we use it to name output files.
# Only rank 0 runs under perf (others run naked) to avoid output interleaving
# and to get clean single-rank numbers (symmetric decomposition: rank 0 is representative).

run_perf_all_ranks() {
    local LABEL=$1
    local BIN=$2
    local OUTDIR=$RESULTS_DIR/$LABEL
    mkdir -p "$OUTDIR"

    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "  Case: $LABEL   (all 4 ranks under perf stat)"
    echo "  Output dir: $OUTDIR"
    echo "  Start: $(date)"
    echo "────────────────────────────────────────────────────────────"

    # Each rank runs under perf stat with its own output file.
    # The wrapper script uses OMPI_COMM_WORLD_RANK set by OpenMPI.
    local WRAPPER=$OUTDIR/perf_wrapper.sh
    cat > "$WRAPPER" << 'WRAPPER_EOF'
#!/bin/bash
RANK=${OMPI_COMM_WORLD_RANK:-0}
OUTDIR=$(dirname "$0")
PERF_BIN="${PERF_BIN:-perf}"
EVENTS="${PERF_EVENTS:-cache-references,cache-misses}"
exec "$PERF_BIN" stat \
    -e "$EVENTS" \
    --output "$OUTDIR/rank_${RANK}.txt" \
    -- "$@"
WRAPPER_EOF
    chmod +x "$WRAPPER"

    PERF_BIN=$PERF PERF_EVENTS=$EVENT_STR \
    $MPIRUN -n 4 --bind-to core --map-by core \
        "$WRAPPER" "$BIN" "$PRM" -Parameters.timesteps $TIMESTEPS

    local RC=$?
    echo "  mpirun exit code: $RC  ($(date))"

    # Show per-rank output
    for rank in 0 1 2 3; do
        local f=$OUTDIR/rank_${rank}.txt
        if [[ -f "$f" ]]; then
            echo ""
            echo "  ── Rank $rank ─────────────────────────────────────────────"
            cat "$f"
        fi
    done
}

# =============================================================================
#  Run both cases
# =============================================================================
run_perf_all_ranks "const_0.08${COLLISION_LABEL}" "$BIN_CONST"
run_perf_all_ranks "tempdep${COLLISION_LABEL}"    "$BIN_TEMPDEP"

# =============================================================================
#  Side-by-side comparison (rank 0 numbers)
# =============================================================================
echo ""
echo "=============================================================="
echo "  COMPARISON — rank 0 (representative for symmetric domain)"
echo "=============================================================="

FILE_A=$RESULTS_DIR/const_0.08${COLLISION_LABEL}/rank_0.txt
FILE_B=$RESULTS_DIR/tempdep${COLLISION_LABEL}/rank_0.txt

if [[ ! -f "$FILE_A" || ! -f "$FILE_B" ]]; then
    echo "  Result files missing — perf stat may have failed."
    exit 1
fi

echo ""
printf "  %-38s  %14s  %14s  %8s\n" "Event" "const ν=0.08" "temp-dep(MF)" "T/C ratio"
printf "  %-38s  %14s  %14s  %8s\n" \
    "$(printf '%.0s─' {1..38})" \
    "$(printf '%.0s─' {1..14})" \
    "$(printf '%.0s─' {1..14})" \
    "$(printf '%.0s─' {1..8})"

# Parse perf stat output: lines look like "  1,234,567,890   event-name"
extract_perf_count() {
    local FILE=$1
    local EVENT=$2
    grep -w "$EVENT" "$FILE" 2>/dev/null \
        | head -1 \
        | awk '{gsub(/,/,"",$1); print $1}'
}

EVENTS_TO_PRINT=(
    "cache-references"
    "cache-misses"
    "L1-dcache-loads"
    "L1-dcache-load-misses"
    "LLC-loads"
    "LLC-load-misses"
    "cycles"
    "instructions"
)

for evt in "${EVENTS_TO_PRINT[@]}"; do
    VA=$(extract_perf_count "$FILE_A" "$evt")
    VB=$(extract_perf_count "$FILE_B" "$evt")
    if [[ -n "$VA" || -n "$VB" ]]; then
        RATIO="-"
        if [[ "$VA" =~ ^[0-9]+$ && "$VB" =~ ^[0-9]+$ && "$VA" -gt 0 ]]; then
            RATIO=$(awk "BEGIN{printf \"%.3f\", $VB/$VA}")
        fi
        printf "  %-38s  %14s  %14s  %8s\n" "$evt" "${VA:--}" "${VB:--}" "$RATIO"
    fi
done

# Derived metrics
echo ""
echo "  ── Derived rates ──────────────────────────────────────────────────────────────"

derive_rate() {
    local FILE=$1
    local NUM_EVT=$2
    local DEN_EVT=$3
    local N=$(extract_perf_count "$FILE" "$NUM_EVT")
    local D=$(extract_perf_count "$FILE" "$DEN_EVT")
    if [[ "$D" =~ ^[0-9]+$ && "$D" -gt 0 && "$N" =~ ^[0-9]+$ ]]; then
        awk "BEGIN{printf \"%.4f\", $N/$D}"
    else
        echo "-"
    fi
}

L1_MISS_RATE_A=$(derive_rate "$FILE_A" "L1-dcache-load-misses" "L1-dcache-loads")
L1_MISS_RATE_B=$(derive_rate "$FILE_B" "L1-dcache-load-misses" "L1-dcache-loads")
LLC_MISS_RATE_A=$(derive_rate "$FILE_A" "LLC-load-misses" "LLC-loads")
LLC_MISS_RATE_B=$(derive_rate "$FILE_B" "LLC-load-misses" "LLC-loads")
LLC_MISS_RATE_A=$(derive_rate "$FILE_A" "LLC-load-misses" "LLC-loads")
CACHE_MISS_RATE_A=$(derive_rate "$FILE_A" "cache-misses" "cache-references")
CACHE_MISS_RATE_B=$(derive_rate "$FILE_B" "cache-misses" "cache-references")
IPC_A=$(derive_rate "$FILE_A" "instructions" "cycles")
IPC_B=$(derive_rate "$FILE_B" "instructions" "cycles")

printf "  %-38s  %14s  %14s\n" "L1D load miss rate" "$L1_MISS_RATE_A" "$L1_MISS_RATE_B"
printf "  %-38s  %14s  %14s\n" "LLC load miss rate" "$LLC_MISS_RATE_A" "$LLC_MISS_RATE_B"
printf "  %-38s  %14s  %14s\n" "Generic cache miss rate" "$CACHE_MISS_RATE_A" "$CACHE_MISS_RATE_B"
printf "  %-38s  %14s  %14s\n" "IPC (instr/cycle)" "$IPC_A" "$IPC_B"

echo ""
echo "=============================================================="
echo "  perf stat analysis completed: $(date)"
echo "=============================================================="
