# Performance Analysis: Constant vs Temperature-Dependent Viscosity in a 3D LBM Couette Flow

**Benchmark setup** — 3D thermal Couette flow, D3Q19 lattice, generated with
lbmpy / pystencils / sweepgen and integrated with the waLBerla LBM framework. The
temperature-dependent viscosity is provided by **MaterForge** from a YAML material
specification and inlined into the generated stream-collide kernel as a SymPy
`Piecewise` expression.

This study reports **both collision operators** that the build supports without source
changes: **SRT** (single relaxation time, the default in `build_validation_binaries.sh`
and `run_perf_*.sh`) and **TRT** (two relaxation times). The two differ enormously in how
much extra *symbolic* arithmetic the temperature-dependent law introduces, which makes the
pair a clean test of whether `count_operations` predicts measured performance.

All performance binaries are built with **`WRITE_VISCOSITY=OFF`**: the viscosity field is a
visualisation-only output and has no analogue in the constant-ν kernel, so omitting it keeps
the const-vs-tempdep comparison operator- and traffic-matched (tempdep 344 B/cell, const
336 B/cell). The validation binary (`CouetteFlowScaling_tempdep`) uses `WRITE_VISCOSITY=ON`
and is *not* used for timing.

## 1. Experimental configuration

| Parameter                | Value                                                                |
|--------------------------|----------------------------------------------------------------------|
| Hardware                 | Intel Xeon Gold 6326 (Ice Lake, 16 cores/socket, 2 sockets @ 2.90 GHz) |
| Memory                   | DDR4-3200, 8 channels/socket (≈ 204.8 GB/s theoretical peak per socket) |
| Cluster                  | woody (NHR@FAU), Ice Lake nodes (`--constraint=icx`), AlmaLinux 8   |
| Compiler                 | GCC 13.3.0 (spack), OpenMPI 4.1.3                                  |
| Parallelism              | 4 MPI processes, each pinned to one core on socket 0 (cores 0–3)   |
| CPU pinning              | `--bind-to core --map-by core` (OpenMPI); verified via `--report-bindings` |
| Node allocation          | `--exclusive` — sole occupant of the node during all runs           |
| Warmup                   | 100 timesteps before timed region (cache, TLB, page-table warmup)  |
| Domain                   | 128 × 64 × 64 cells (524 288 total)                                 |
| Decomposition            | 1 × 2 × 2 blocks of 128 × 32 × 32 cells                            |
| Lattice / precision      | D3Q19, FP64                                                        |
| Collision operators      | **SRT** (default) and **TRT**, each const-vs-tempdep on one node   |
| Viscosity write          | **OFF** for all timed binaries (344 B/cell tempdep, 336 B/cell const) |
| Wall velocity `u_max`    | 0.025 (lattice units)                                               |
| Wall temperatures        | T_bottom = 300 K, T_top = 3000 K                                    |
| Timesteps                | 60 000                                                              |
| Trials per case          | **5** (serial, no overlap between cases)                            |
| Statistics               | mean ± sample std; 95 % CI from t-distribution (df = 4)             |

Two cases are compared on identical hardware, domain, and timestep count, for each operator:

1. **Constant viscosity** — ν = 0.08, baked into the generated kernel at compile time
   via `cmake -DUSE_MATERFORGE=OFF -DCONST_NU=0.08`. The viscosity field is not written.
2. **Temperature-dependent viscosity (MaterForge)** — `CouetteFlowMaterial.yaml` defines
   `dynamic_viscosity(T) = 0.16667·exp(−0.0005·(T − 300))` (a dimensionless lattice
   viscosity) over T ∈ [300, 3000] K. MaterForge fits it with a two-segment degree-2
   piecewise polynomial (break at T ≈ 1536 K), inlined as the `Piecewise` shown below. The
   temperature field is read each step; the viscosity field is **not** written (perf config).

```
ν(T) ≈
    1.502e-8·T² − 8.954e-5·T + 0.19194   if T < 1536.24 K
    8.116e-9·T² − 6.852e-5·T + 0.17593   otherwise
```

## 2. Code-generation metric — `count_operations`

`pystencils.sympyextensions.count_operations` reports the symbolic arithmetic operation
count of the stream-collide `AssignmentCollection` *before compilation*. Numbers below are
taken directly from the code-generation logs (`apps/logs/build/*.log`).

### 2.1 Operation counts per operator (`WRITE_VISCOSITY=OFF`)

| Operator | Case      | adds | muls | divs | **Total** | Algorithmic BW |
|----------|-----------|-----:|-----:|-----:|----------:|---------------:|
| **SRT**  | const     |  149 |  182 |    1 | **332**   | 336 B/cell     |
| **SRT**  | temp-dep  |  173 |  172 |    2 | **347**   | 344 B/cell     |
| **TRT**  | const     |  167 |  212 |    1 | **380**   | 336 B/cell     |
| **TRT**  | temp-dep  |  234 |  289 |    4 | **527**   | 344 B/cell     |

Bandwidth (both operators): base 336 B/cell = 19 PDFs read × 8 B + 19 PDFs write × 8 B +
density write 8 B + velocity write 24 B; temp-dep adds the 8 B/cell temperature read → 344.
(With `WRITE_VISCOSITY=ON` the tempdep kernel would add a further 8 B viscosity write → 352
B/cell and ≈ +10 ops; that build is for VTK only and is not benchmarked here.)

### 2.2 Symbolic deltas — the key contrast

| Operator | const → temp-dep ops | Δ ops | **Δ %**   | BW Δ        |
|----------|---------------------:|------:|----------:|------------:|
| **SRT**  | 332 → 347            | **+15**  | **+4.5 %**  | +8 B (+2.4 %) |
| **TRT**  | 380 → 527            | **+147** | **+38.7 %** | +8 B (+2.4 %) |

The temperature-dependent law inflates the symbolic op count **8.6× more for TRT than for
SRT**. For SRT, lbmpy common-subexpression-eliminates ω(T) so it is evaluated once per cell
in the subexpressions block and the collision main barely changes (+1 op on the rest PDF);
the +14 remaining ops are the ν(T) polynomial. For TRT, the symbolic factors `rr_0·(…)` and
`rr_1·(…)` cannot be folded into literal coefficients when ω is a per-cell symbol, so every
PDF update retains explicit relaxation factors → +147 ops. **§3 shows this 8.6× symbolic gap
all but vanishes in the measured wall-clock overhead.**

## 3. Runtime performance

Each case ran 60 000 timesteps on an exclusive node with CPU pinning, 5 independent trials,
const and temp-dep submitted on the **same node** for each operator (SRT on w2411, TRT on
w2502). Statistics use the two-sided t-distribution with df = 4.

### 3.1 Headline numbers (5-trial statistics)

| Metric                       | **SRT** const | **SRT** temp-dep | **TRT** const | **TRT** temp-dep |
|------------------------------|--------------:|-----------------:|--------------:|-----------------:|
| Total MLUPS (mean ± std)     | 73.43 ± 0.04  | 67.58 ± 0.15     | 67.38 ± 0.13  | 62.27 ± 0.13     |
| 95 % CI                      | [73.38, 73.48]| [67.40, 67.77]   | [67.22, 67.55]| [62.11, 62.43]   |
| Wall time (mean)             | 428.4 s       | 465.5 s          | 466.8 s       | 505.2 s          |
| **MLUPS overhead (T vs C)**  | —             | **7.97 %**       | —             | **7.59 %**       |
| **Wall-clock overhead**      | —             | **8.66 %**       | —             | **8.21 %**       |
| L∞ error (vs analytical)     | 1.51 × 10⁻⁷   | 1.76 × 10⁻⁶      | 1.51 × 10⁻⁷   | 1.76 × 10⁻⁶      |
| Convergence threshold passed | —             | yes (5/5)        | —             | yes (5/5)        |

**The central result:** SRT adds **+4.5 %** symbolic ops and TRT adds **+38.7 %** — an 8.6×
difference — yet the measured temperature-dependent overhead is **7.6–8.0 %** for *both*
operators (8.2–8.7 % wall-clock). The overhead is essentially independent of the symbolic
arithmetic increase, because both kernels are bandwidth-bound and the dominant marginal cost
is the +8 B/cell temperature read (+2.4 % algorithmic bandwidth), not the polynomial FLOPs.

SRT is the faster operator in absolute terms (≈ 73 vs ≈ 67 MLUPS for const), consistent with
its lower op count; TRT eliminates wall slip exactly for any ν/grid at ~8 % lower throughput.
Both reach the same L∞ accuracy here.

The 95 % confidence intervals for const and temp-dep do not overlap within either operator,
so the overhead is statistically significant (the per-operator std is < 0.25 % of the mean).

### 3.2 waLBerla timer breakdown (SRT, 60 000 steps, 4-rank sums, mean of trials)

| Timer              | const ν=0.08 | temp-dep   | Δ time    | const %  | temp-dep % |
|--------------------|-------------:|-----------:|----------:|---------:|-----------:|
| **StreamCollide**  | **1489.9 s** | **1638.4 s** | **+148.5 s** | **87.0 %** | **88.0 %** |
| LBM Communication  |    180.4 s   |    180.8 s   |   +0.4 s  |  10.5 %  |   9.7 %   |
| UBB boundary       |     31.9 s   |     31.8 s   |   ±0.0 s  |   1.9 %  |   1.7 %   |
| NoSlip boundary    |     10.3 s   |     10.1 s   |   ±0.0 s  |   0.6 %  |   0.5 %   |

The entire temp-dep overhead lands in **StreamCollide** (+10.0 % per-rank); communication and
boundary times are statistically identical between cases. TRT shows the same structure
(StreamCollide ≈ 89 % of wall time, all overhead in the collision kernel). `REDUCE_MAX /
REDUCE_AVG ≈ 1.0` for StreamCollide in all cases, confirming negligible rank imbalance with
4 well-pinned, equally-loaded processes.

### 3.3 Bandwidth-vs-compute regime

At the observed SRT const-ν throughput of 73.43 MLUPS:

  73.43 × 10⁶ cells/s × 336 B/cell = **24.7 GB/s of DRAM bandwidth**

against ≈ 204.8 GB/s peak per socket. The 4 pinned processes on socket 0 consume roughly
**12 % of one socket's peak bandwidth**. The temp-dep case (67.58 MLUPS × 344 B = 23.2 GB/s)
and both TRT cases (22.6 / 21.4 GB/s) sit in the same ≈ 11–12 % band. The kernel is far from
bandwidth-saturated, which is exactly why extra arithmetic is cheap (§4).

## 4. Why the overhead is independent of the symbolic FLOP ratio

| Quantity                          | SRT     | TRT     |
|-----------------------------------|--------:|--------:|
| Symbolic op ratio (count_operations) | 1.045 | 1.387   |
| Measured wall-time ratio          | 1.087   | 1.082   |
| Direction vs prediction           | measured **above** symbolic | measured **far below** symbolic |

The symbolic ratio is a wildly unreliable predictor in *both* directions: it under-predicts
the SRT overhead (1.045 symbolic vs 1.087 measured) and massively over-predicts the TRT
overhead (1.387 symbolic vs 1.082 measured). The reason the two operators land at the same
~8 % is structural, not coincidental:

1. **Bandwidth headroom hides arithmetic.** At ≈ 12 % of peak bandwidth, each memory access
   is serviced with short latency relative to compute throughput. The extra polynomial
   arithmetic executes in functional units that are otherwise idle waiting for the next
   cache line, so it is **absorbed inside memory-access stalls** by the Ice Lake
   out-of-order engine (512-entry reorder buffer). For TRT the extra `rr_0/rr_1` factors are
   far more numerous but still fit inside those stalls.

2. **AVX-512 execution width + FMA fusion.** The Xeon Gold 6326 executes 8 × FP64 lanes per
   instruction. Many symbolic adds/muls fuse into single FMA instructions, so the *retired
   instruction* increase is much smaller than the symbolic op increase (§8: SRT +14 %, TRT
   +20 % real instructions vs +4.5 % / +38.7 % symbolic).

3. **`Piecewise` becomes a branchless blend.** GCC emits the conditional as a vectorised
   select (`vblendvpd`/mask): both polynomial branches are evaluated and one is selected, so
   there is no branch-misprediction cost and the workload is divergence-free.

The marginal cost that *does* survive is the +8 B/cell temperature read (+2.4 % algorithmic
bandwidth) plus second-order memory-system effects, both operator-independent — hence the
common ~8 %.

A useful framing for a paper:

> *On a bandwidth-rich server CPU (Xeon Gold 6326, ≈ 12 % of peak bandwidth utilised), the
> MaterForge temperature-dependent viscosity kernel adds only 7.6–8.0 % MLUPS overhead
> (8.2–8.7 % wall-clock) for both SRT and TRT collision — despite the temperature-dependent
> law inflating the symbolic operation count by +4.5 % (SRT) and +38.7 % (TRT). The 8.6×
> difference in symbolic arithmetic does not appear in the measured overhead: the kernel is
> bandwidth-bound and the extra arithmetic is hidden inside memory-access stalls by the deep
> out-of-order engine and AVX-512 execution units. FLOP counts alone cannot predict LBM
> kernel performance.*

## 5. Why `count_operations` is not a sufficient performance metric

`count_operations` is a useful **code-generation diagnostic**: it cheaply exposes how complex
a generated kernel is and lets one compare material models without rebuilding. But the SRT/TRT
contrast above shows it fails as a *performance* metric on real hardware:

| Limitation                                   | Consequence on Xeon Gold 6326 (Ice Lake)                       |
|----------------------------------------------|----------------------------------------------------------------|
| Treats all symbolic ops identically          | over-predicts TRT overhead 5×; under-predicts SRT overhead     |
| Counts both `Piecewise` branches             | only one polynomial result is used; both are evaluated cheaply |
| Ignores compiler-pipeline effects            | misses AVX-512 FMA fusion (symbolic → retired-instruction gap) |
| Has no notion of memory traffic or bandwidth | misses the bandwidth headroom that absorbs extra compute       |
| Platform-agnostic by construction            | cannot predict how overhead scales with bandwidth utilisation  |

The right end-to-end metric is **MLUPS**, printed by the waLBerla driver after every run.
Intended usage: (1) use `count_operations` at code-gen time as a quick sanity check that the
AST is not pathologically large; (2) report MLUPS (with the StreamCollide-vs-Communication
breakdown) as the benchmark figure, ≥ 5 trials with 95 % CIs; (3) run on the *target* hardware
class — bandwidth utilisation, not raw FLOP count, governs how much extra arithmetic matters.

## 6. Summary for the publication

* **On woody (Xeon Gold 6326, Ice Lake, NHR@FAU, exclusive node, 5 trials, 60 000 steps,
  524 288-cell domain, 4 MPI ranks, `WRITE_VISCOSITY=OFF`):** MaterForge temperature-dependent
  viscosity adds **7.97 % MLUPS overhead for SRT** (8.66 % wall-clock) and **7.59 % for TRT**
  (8.21 % wall-clock) — even though the symbolic op count rises +4.5 % (SRT) vs +38.7 % (TRT).
* The kernel sustains ≈ 11–12 % of peak memory bandwidth; it is strongly **memory-bound**. The
  extra polynomial arithmetic is absorbed by the deep out-of-order pipeline and AVX-512 units.
  The +2.4 % bandwidth delta (temperature read) and memory-system effects — both
  operator-independent — account for essentially all of the measured overhead.
* The simulation **converges correctly and reproducibly** in all cases (L∞ identical across
  all 5 trials: 1.51 × 10⁻⁷ for constant ν; 1.76 × 10⁻⁶ for temp-dep, bounded by the degree-12
  analytical-reference polynomial fit, not by the LBM solution).
* **Recommendation for the paper:** report `count_operations` as a static complexity figure and
  report MLUPS with a StreamCollide breakdown (mean ± std, 95 % CI, ≥ 5 trials) as the
  wall-clock metric, for both collision operators. The SRT/TRT pair is itself the cleanest
  demonstration that symbolic FLOP counts do not predict bandwidth-bound LBM performance.

## 7. Reproducibility

### 7.1 Build environment (woody)

```bash
# Modules: cmake/3.30.5; toolchain GCC 13.3.0 + OpenMPI 4.1.3 (spack)
# Python virtualenv at $HOME/.venvs/materforge
git submodule update --init --recursive
source $HOME/.venvs/materforge/bin/activate
```

### 7.2 Build all benchmark binaries (SRT)

```bash
cd apps
bash scripts/build_validation_binaries.sh
# Produces (among others):
#   CouetteFlowScaling_const_0.08       (SRT, WRITE_VISCOSITY=OFF, 336 B/cell)
#   CouetteFlowScaling_tempdep_perf     (SRT, WRITE_VISCOSITY=OFF, 344 B/cell)  <- timed
#   CouetteFlowScaling_tempdep          (SRT, WRITE_VISCOSITY=ON, 352 B/cell)   <- VTK only
```

### 7.3 Build the TRT binaries

```bash
cd apps
export MODULEPATH=/apps/modules/data/tools
eval "$(MODULESHOME=/apps/modules /usr/bin/modulecmd bash load cmake/3.30.5)"
# const TRT
cmake --preset woody-release-cpu -DUSE_MATERFORGE=OFF -DCONST_NU=0.08 \
      -DCOLLISION_OP=TRT -DWRITE_VISCOSITY=OFF -S "$PWD" -B build/woody-release-cpu
make -C build/woody-release-cpu -j8 CouetteFlowScaling
cp build/woody-release-cpu/CouetteFlowScaling build/woody-release-cpu/CouetteFlowScaling_const_0.08_TRT
# temp-dep TRT (perf)
cmake --preset woody-release-cpu -DUSE_MATERFORGE=ON \
      -DCOLLISION_OP=TRT -DWRITE_VISCOSITY=OFF -S "$PWD" -B build/woody-release-cpu
make -C build/woody-release-cpu -j8 CouetteFlowScaling
cp build/woody-release-cpu/CouetteFlowScaling build/woody-release-cpu/CouetteFlowScaling_tempdep_perf_TRT
```

### 7.4 Run the benchmarks (SLURM)

```bash
# SRT — const + temp-dep in one --exclusive allocation, same node
sbatch scripts/run_perf_srt.sh
# TRT — const + temp-dep in one job, same node
sbatch scripts/run_perf_trt.sh
```

> **Same-node requirement.** The const→tempdep overhead is < 10 %, smaller than the
> node-to-node hardware variation across woody's Ice Lake pool. Both cases must run on the
> *same* physical node or the comparison is confounded. `run_perf_srt.sh`/`run_perf_trt.sh`
> enforce this by running both binaries inside a single allocation. Submitting the standalone
> `run_perf_const.sh` and `run_perf_tempdep.sh` as two independent jobs does **not** pin them
> to one node (SLURM may place them on different `icx` nodes) and must not be used for the
> headline comparison.

SLURM job IDs for this dataset (NHR@FAU woody, Ice Lake, exclusive):
* 11901976 — SRT const + temp-dep MLUPS, single same-node allocation, node w2411 (2026-06-13, `run_perf_srt.sh`)
* 11826268 — TRT const + temp-dep MLUPS (node w2502)
* 11844061 — validation array, 10 cases (per-case PASS, L∞ above)

Raw logs: `apps/logs/performance/run_perf_{const,tempdep,trt}.log`,
`apps/logs/build/*.log`, `apps/logs/validation/run_validation_*.log`.

## 8. Cache counter analysis (Linux `perf stat`)

### 8.1 Methodology

Direct hardware-counter profiling on woody required three attempts; only `perf stat` worked:

| Tool | Failure mode | Root cause |
|------|-------------|-----------|
| Intel VTune 2022.3 | `requires perf system-wide profiling` | `perf_event_paranoid = 2` blocks kernel-driver PMU sampling |
| LIKWID 5.4.1 (`CACHE`, `MEM`) | `Access to performance monitoring registers locked` | MSR-device access restricted cluster-wide |
| Linux `perf stat` (user-space counting) | **Success** | `PERF_EVENT_OPEN` counting mode allowed at `paranoid = 2` for per-process events |

`perf stat` ran with per-rank isolation (each MPI rank under its own `perf stat --output
rank_N.txt`), via `apps/scripts/run_perf_cache.sh` (set `COLLISION_LABEL=_TRT` to profile the
TRT binaries). **Collection scope:** 100 warmup + 5000 timed steps per rank — long enough for
the working set to reach steady state, so the cache-level *ratios* are representative; the
MLUPS/convergence lines printed by these short runs are timing artefacts and are ignored.
SLURM jobs: 11826261 (SRT) and 11826269 (TRT), node w2201. Raw outputs:
`apps/output/profiling/perf/{const_0.08,tempdep,const_0.08_TRT,tempdep_TRT}/rank_{0..3}.txt`.

### 8.2 Derived cache rates and IPC (4-rank totals)

| Metric                  | SRT const | SRT temp-dep | T/C (SRT) | TRT const | TRT temp-dep | T/C (TRT) |
|-------------------------|----------:|-------------:|----------:|----------:|-------------:|----------:|
| Retired instructions    | 25.28 G   | 28.89 G      | **1.142** | 26.95 G   | 32.24 G      | **1.196** |
| Cycles                  | 10.56 G   | 12.06 G      | 1.141     | 11.56 G   | 13.11 G      | 1.133     |
| **IPC**                 | 2.393     | 2.396        | **1.001** | 2.331     | 2.460        | **1.056** |
| L1-dcache-loads         | 7.98 G    | 9.28 G       | **1.163** | 8.97 G    | 10.60 G      | **1.182** |
| L1D load miss rate      | 4.06 %    | 3.58 %       | **0.882** | 3.61 %    | 3.12 %       | **0.866** |
| LLC load miss rate      | 62.4 %    | 62.5 %       | 1.002     | 61.0 %    | 63.4 %       | 1.039     |
| Generic cache miss rate | 28.3 %    | 27.8 %       | 0.982     | 27.8 %    | 27.6 %       | 0.994     |

### 8.3 Findings

**Finding 1 — Symbolic ops ≠ retired instructions.** The retired-instruction ratio is
**1.142 (SRT)** and **1.196 (TRT)** — almost identical at the hardware level — even though the
symbolic op ratios differ 8.6× (1.045 vs 1.387). For SRT the real instruction increase
*exceeds* the symbolic one (the ν(T) polynomial and its L1-resident intermediates are not
captured by the collision-main op count); for TRT it is *far below* the symbolic one (AVX-512
FMA fusion and the branchless `Piecewise` blend collapse the +147 symbolic ops). Either way,
`count_operations` mispredicts.

**Finding 2 — Extra L1 traffic, but L1 miss rate falls.** Both temp-dep kernels issue +16–18 %
more L1 data loads (polynomial coefficient broadcasts, intermediate products, piecewise
operands), yet the **L1D miss rate drops** (4.06 → 3.58 % SRT; 3.61 → 3.12 % TRT). The extra
arithmetic between consecutive PDF accesses gives the Ice Lake prefetcher more time to land
the next PDF cache line — better temporal locality despite a larger instruction footprint.

**Finding 3 — DRAM pressure essentially unchanged.** LLC load-miss rate moves by ≤ 2.4 points
in either direction; the generic cache-miss rate is flat. The +2.4 % algorithmic byte traffic
(temperature read) is well within L3 streaming throughput at 12 % bandwidth utilisation and
does not change the L3 hit/miss structure over a 5000-step window — consistent with §3.3. The
authoritative traffic signal is the full 60 000-step MLUPS measurement, not these short-run
cache counters.

**Finding 4 — IPC differs by operator.** SRT IPC is flat (2.393 → 2.396): the small extra
arithmetic neither helps nor hurts pipeline utilisation. TRT IPC *rises* (2.331 → 2.460,
+5.6 %): the many independent `rr_0/rr_1` multiply-add chains expose more instruction-level
parallelism than the const kernel's folded path, so the CPU retires the extra instructions
more efficiently. In both cases the net wall-clock overhead stays ~8 % because it is set by
memory traffic, not by the compute pipeline.

### 8.4 Summary table

| Metric (T/C ratio)              | SRT    | TRT    | Interpretation                                         |
|---------------------------------|-------:|-------:|--------------------------------------------------------|
| Symbolic FLOPs (count_operations) | 1.045 | 1.387  | code-generation complexity — **mispredicts both ways** |
| Retired instructions            | 1.142  | 1.196  | hardware instruction cost, similar for both operators  |
| L1-dcache-loads                 | 1.163  | 1.182  | L1-resident polynomial intermediates                   |
| L1D load miss rate              | 0.862  | 0.866  | miss rate *decreases* — prefetcher gains compute slack |
| IPC                             | 1.001  | 1.056  | flat (SRT) / improved (TRT) — extra ILP from rr-chains |
| **Wall-clock ratio (60 000 steps)** | **1.087** | **1.082** | **authoritative — ~8 % for both, bandwidth-bound** |

> **perf stat raw data:** `apps/output/profiling/perf/{const_0.08,tempdep,const_0.08_TRT,tempdep_TRT}/rank_{0..3}.txt`
> **Collection script:** `apps/scripts/run_perf_cache.sh` (`COLLISION_LABEL=_TRT` for TRT); SLURM 11826261 (SRT), 11826269 (TRT), node w2201.
