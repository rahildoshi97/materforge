#!/bin/bash
#SBATCH --job-name=couette_likwid
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive            # sole occupant — no shared-cache/DRAM-bandwidth interference
#SBATCH --time=01:00:00        # 4 groups × 2 cases × ~45 s + finalization
#SBATCH --output=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/likwid_cache.log
#SBATCH --error=/home/hpc/iwia/iwia133h/.local/materforge/apps/logs/likwid_cache.err

# =============================================================================
#  LIKWID Cache Counter Analysis — Couette Flow (const ν=0.08 vs temp-dep)
#
#  Why LIKWID instead of VTune:
#    /apps/likwid/5.4.1/sbin/likwid-accessD is setuid root — it reads
#    CPU MSR registers directly, bypassing perf_event_paranoid=2 which
#    blocked VTune's hardware-event sampling mode.
#
#  Performance groups collected (2 runs per case):
#    CACHE  — L1/L2/L3 load miss rates, eviction counts
#    MEM    — DRAM read/write bandwidth per socket
#
#  All runs use 5000 timesteps (~42 s/run): enough for stable PMU statistics.
#  Each group requires its own run (hardware PMU has limited counter slots).
# =============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
APPS_DIR=/home/hpc/iwia/iwia133h/.local/materforge/apps
BUILD_DIR=$APPS_DIR/build/woody-release-cpu
LOGS_DIR=$APPS_DIR/logs
RESULTS_DIR=$APPS_DIR/likwid_results
PRM=$APPS_DIR/CouetteFlowScaling.prm

BIN_CONST=$BUILD_DIR/CouetteFlowScaling_const_0.08
BIN_TEMPDEP=$BUILD_DIR/CouetteFlowScaling_tempdep

MPIRUN=/apps/SPACK/0.18.1/opt/linux-almalinux8-skylake/gcc-12.1.0/openmpi-4.1.3-i76wyn4x5aobewottoyj2od33hhtyb2h/bin/mpirun
LIKWID=/apps/likwid/5.4.1/bin/likwid-perfctr

# Cores matching the MPI binding (4 ranks on cores 0-3 of socket 0)
CORES="0-3"
TIMESTEPS=5000

# ── Sanity checks ─────────────────────────────────────────────────────────────
echo "=============================================================="
echo "  LIKWID Cache Analysis — Couette Flow"
echo "  Node:     $(hostname)"
echo "  Date:     $(date)"
echo "  CPU:      $(grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)"
echo "  perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid)"
echo "  LIKWID:   $LIKWID"
echo "  Cores:    $CORES"
echo "  Timesteps per run: $TIMESTEPS"
echo "=============================================================="

for bin in "$BIN_CONST" "$BIN_TEMPDEP" "$MPIRUN" "$LIKWID"; do
    if [[ ! -x "$bin" ]]; then
        echo "ERROR: not found or not executable: $bin"
        exit 1
    fi
done

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

# ── Helper: run one LIKWID group ──────────────────────────────────────────────
# Usage: run_likwid LABEL GROUP BIN
#   LABEL  — short identifier (const_0.08 | tempdep)
#   GROUP  — LIKWID performance group (CACHE | MEM | ...)
#   BIN    — full path to the CouetteFlowScaling binary
run_likwid() {
    local LABEL=$1
    local GROUP=$2
    local BIN=$3
    local OUTFILE=$RESULTS_DIR/${LABEL}_${GROUP}.txt

    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "  Group: $GROUP   Case: $LABEL"
    echo "  Output: $OUTFILE"
    echo "  Start: $(date)"
    echo "────────────────────────────────────────────────────────────"

    # LIKWID wraps the whole mpirun.
    # -C 0-3   : measure hardware counters on exactly cores 0-3 (matches --map-by core)
    # -g GROUP : performance group
    # -f       : force (overwrite existing results)
    # -O       : machine-readable output (pipe-separated) — also written to OUTFILE below
    # MPI binding is left to mpirun --map-by core; LIKWID only measures, does not re-pin.
    $LIKWID \
        -C "$CORES" \
        -g "$GROUP" \
        -f \
        -- \
        $MPIRUN -n 4 --bind-to core --map-by core "$BIN" "$PRM" \
            -Parameters.timesteps $TIMESTEPS \
        2>&1 | tee "$OUTFILE"

    local RC=${PIPESTATUS[0]}
    echo ""
    echo "  LIKWID exit code: $RC  ($(date))"
    if [[ $RC -ne 0 ]]; then
        echo "  WARNING: non-zero exit — check $LOGS_DIR/likwid_cache.err"
    fi
}

# =============================================================================
#  PART 1 — CACHE group
#  Metrics: L1D load miss rate, L2 load miss rate, L3 load miss rate,
#           L2 eviction rate, L3 eviction rate
# =============================================================================
echo ""
echo "=============================================================="
echo "  PART 1: CACHE group (L1/L2/L3 miss rates)"
echo "=============================================================="

run_likwid "const_0.08" "CACHE" "$BIN_CONST"
run_likwid "tempdep"    "CACHE" "$BIN_TEMPDEP"

# =============================================================================
#  PART 2 — MEM group
#  Metrics: Memory read bandwidth [MB/s], Memory write bandwidth [MB/s],
#           Total memory bandwidth [MB/s]
#  This directly quantifies §3.3 of couette_performance_analysis.md:
#  the +5.3% memory traffic from the extra temperature/viscosity fields.
# =============================================================================
echo ""
echo "=============================================================="
echo "  PART 2: MEM group (DRAM bandwidth)"
echo "=============================================================="

run_likwid "const_0.08" "MEM" "$BIN_CONST"
run_likwid "tempdep"    "MEM" "$BIN_TEMPDEP"

# =============================================================================
#  PART 3 — Side-by-side comparison
#  Parse LIKWID output files and extract the "Sum" column (aggregate across
#  all 4 measured cores) for each metric, then print a comparison table.
# =============================================================================
echo ""
echo "=============================================================="
echo "  PART 3: Comparison summary"
echo "=============================================================="

# Extract the value from LIKWID's metric table.
# LIKWID metric lines look like:
#   |Metric name                | Core 0  | Core 1  | Core 2  | Core 3  | Sum     | Avg     |
# We want the "Sum" column (second-to-last before Avg).
extract_metric() {
    local FILE=$1
    local METRIC_PATTERN=$2
    # Find the line matching the metric, grab the "Sum" column
    # LIKWID table uses | as delimiter; Sum is the 6th field (after 4 core cols)
    grep -i "$METRIC_PATTERN" "$FILE" 2>/dev/null \
        | grep -v "^[+#-]" \
        | head -1 \
        | awk -F'|' '{
            # Strip whitespace from each field
            for(i=1;i<=NF;i++) gsub(/[[:space:]]/,"",$i);
            # Sum column is NF-2 (before Avg, before trailing empty)
            if(NF>=6) print $(NF-2);
            else print "-"
          }'
}

# Extract Avg column instead (per-core average, more intuitive for rates)
extract_avg() {
    local FILE=$1
    local METRIC_PATTERN=$2
    grep -i "$METRIC_PATTERN" "$FILE" 2>/dev/null \
        | grep -v "^[+#-]" \
        | head -1 \
        | awk -F'|' '{
            for(i=1;i<=NF;i++) gsub(/[[:space:]]/,"",$i);
            if(NF>=6) print $(NF-1);
            else print "-"
          }'
}

print_metric_row() {
    local LABEL=$1
    local METRIC_PATTERN=$2
    local FILE_A=$3
    local FILE_B=$4
    local USE_AVG=${5:-false}

    local VAL_A VAL_B
    if [[ "$USE_AVG" == "true" ]]; then
        VAL_A=$(extract_avg  "$FILE_A" "$METRIC_PATTERN")
        VAL_B=$(extract_avg  "$FILE_B" "$METRIC_PATTERN")
    else
        VAL_A=$(extract_metric "$FILE_A" "$METRIC_PATTERN")
        VAL_B=$(extract_metric "$FILE_B" "$METRIC_PATTERN")
    fi

    # Compute ratio if both are numbers
    local RATIO="-"
    if [[ "$VAL_A" =~ ^[0-9.eE+-]+$ && "$VAL_B" =~ ^[0-9.eE+-]+$ ]]; then
        RATIO=$(awk "BEGIN{if($VAL_A>0) printf \"%.3f\", $VAL_B/$VAL_A; else print \"-\"}")
    fi

    printf "  %-42s  %12s  %12s  %8s\n" "$LABEL" "${VAL_A:--}" "${VAL_B:--}" "$RATIO"
}

CACHE_A=$RESULTS_DIR/const_0.08_CACHE.txt
CACHE_B=$RESULTS_DIR/tempdep_CACHE.txt
MEM_A=$RESULTS_DIR/const_0.08_MEM.txt
MEM_B=$RESULTS_DIR/tempdep_MEM.txt

echo ""
printf "  %-42s  %12s  %12s  %8s\n" "Metric" "const ν=0.08" "temp-dep(MF)" "T/C ratio"
printf "  %-42s  %12s  %12s  %8s\n" \
    "$(printf '%.0s─' {1..42})" \
    "$(printf '%.0s─' {1..12})" \
    "$(printf '%.0s─' {1..12})" \
    "$(printf '%.0s─' {1..8})"

echo "  ── Cache miss rates (per-core avg) ─────────────────────────────────────────────"
print_metric_row "L1D load miss rate"            "L1D.*miss rate\|L1.*load.*miss rate"      "$CACHE_A" "$CACHE_B" true
print_metric_row "L2 load miss rate"             "L2.*miss rate\|L2.*load.*miss rate"       "$CACHE_A" "$CACHE_B" true
print_metric_row "L3 load miss rate"             "L3.*miss rate\|LLC.*miss rate"            "$CACHE_A" "$CACHE_B" true
print_metric_row "L2 eviction rate"              "L2.*evict"                                "$CACHE_A" "$CACHE_B" true
print_metric_row "L3 eviction rate"              "L3.*evict\|LLC.*evict"                   "$CACHE_A" "$CACHE_B" true

echo "  ── Cache event counts (sum across 4 cores) ──────────────────────────────────────"
print_metric_row "L1D load accesses (M)"         "L1D.*access\|L1.*load [^m]"              "$CACHE_A" "$CACHE_B" false
print_metric_row "L1D load misses (M)"           "L1D.*miss\|L1.*miss [^r]"                "$CACHE_A" "$CACHE_B" false
print_metric_row "L2 load accesses (M)"          "L2.*access\|L2.*request"                 "$CACHE_A" "$CACHE_B" false
print_metric_row "L2 load misses (M)"            "L2.*miss [^r]\|L2D.*miss"               "$CACHE_A" "$CACHE_B" false
print_metric_row "L3 load accesses (M)"          "L3.*access\|LLC.*access"                 "$CACHE_A" "$CACHE_B" false
print_metric_row "L3 load misses (M)"            "L3.*miss [^r]\|LLC.*miss [^r]"          "$CACHE_A" "$CACHE_B" false

echo "  ── DRAM bandwidth (MB/s, sum 4 cores = whole socket) ──────────────────────────"
print_metric_row "Memory read bandwidth (MB/s)"  "Memory read\|Read BW\|read bandwidth"    "$MEM_A"   "$MEM_B"   false
print_metric_row "Memory write bandwidth (MB/s)" "Memory write\|Write BW\|write bandwidth" "$MEM_A"   "$MEM_B"   false
print_metric_row "Total memory bandwidth (MB/s)" "Memory bandwidth\|Total.*BW\|total band" "$MEM_A"   "$MEM_B"   false

echo ""

# =============================================================================
#  Write summary to file
# =============================================================================
SUMMARY=$LOGS_DIR/likwid_cache_summary.txt
{
echo "LIKWID Cache Analysis Summary"
echo "Run date:  $(date)"
echo "Node:      $(hostname)"
echo "CPU:       $(grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)"
echo "Timesteps: $TIMESTEPS (full benchmark: 60000)"
echo "Domain:    128×64×64, 4 MPI ranks on cores 0-3 (socket 0)"
echo ""
echo "Raw result files:"
ls -lh "$RESULTS_DIR"/*.txt 2>/dev/null
echo ""
echo "See full LIKWID tables in each result file above."
echo ""
echo "Theoretical prediction (couette_performance_analysis.md §3.3):"
echo "  const  ν: 304 B/cell × 62.88 MLUPS = 19.1 GB/s"
echo "  temp-dep: 320 B/cell × 58.78 MLUPS = 18.8 GB/s  (+5.3% traffic, -1.6% throughput)"
echo "  → Expected MEM total BW ratio (T/C) ≈ 0.99 (slightly lower due to lower MLUPS)"
} > "$SUMMARY"

echo "  Summary written to: $SUMMARY"
echo ""
echo "  Raw LIKWID results:"
ls -lh "$RESULTS_DIR"/*.txt 2>/dev/null

echo ""
echo "=============================================================="
echo "  LIKWID cache analysis completed: $(date)"
echo "=============================================================="
