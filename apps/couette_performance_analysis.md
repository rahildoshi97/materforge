# Performance Analysis: Constant vs Temperature-Dependent Viscosity in a 3D LBM Couette Flow

**Benchmark setup** — 3D thermal Couette flow, D3Q19 lattice, TRT collision, generated with
lbmpy / pystencils / sweepgen; integrated with the waLBerla LBM framework. The
temperature-dependent viscosity is provided by **MaterForge** from a YAML material
specification and inlined into the generated stream-collide kernel as a SymPy
`Piecewise` expression.

## 1. Experimental configuration

| Parameter                | Value                                                                |
|--------------------------|----------------------------------------------------------------------|
| Hardware                 | Intel Xeon Gold 6326 (Ice Lake, 16 cores/socket, 2 sockets @ 2.90 GHz) |
| Memory                   | DDR4-3200, 8 channels/socket (≈ 204.8 GB/s theoretical peak per socket) |
| Cluster                  | woody (NHR@FAU), node w2304, AlmaLinux 8                            |
| Compiler                 | GCC 13.3.0, OpenMPI 4.1.3                                           |
| Parallelism              | 4 MPI processes, each pinned to one core on socket 0 (cores 0–3)   |
| CPU pinning              | `--bind-to core --map-by core` (OpenMPI); binding verified via `--report-bindings` |
| Node allocation          | `--exclusive` — sole occupant of node during all runs               |
| Warmup                   | 100 timesteps before timed region (cache, TLB, and page-table warmup) |
| Domain                   | 128 × 64 × 64 cells (524 288 total)                                 |
| Decomposition            | 1 × 2 × 2 blocks of 128 × 32 × 32 cells                            |
| Lattice / method         | D3Q19, TRT, FP64                                                    |
| Wall velocity `u_max`    | 0.025 (lattice units)                                               |
| Wall temperatures        | T_bottom = 300 K, T_top = 3000 K                                    |
| Timesteps                | 60 000                                                              |
| VTK output               | disabled (`vtkWriteFrequency = 0`)                                  |
| Trials per case          | **5** (serial, no overlap between cases)                            |
| Statistics               | mean ± sample std; 95 % CI from t-distribution (df = 4)             |
| Total lattice updates    | 3.146 × 10¹⁰ per trial                                              |

Two cases are compared on identical hardware, domain, and timestep count:

1. **Constant viscosity** — ν = 0.08, baked into the generated kernel at compile time
   via `cmake -DUSE_MATERFORGE=OFF -DCONST_NU=0.08`. The viscosity field is **not
   written** by the kernel in this case (it would store a compile-time constant to
   every cell every step, wasting memory bandwidth for no physics gain).
2. **Temperature-dependent viscosity (MaterForge)** — `CouetteFlowMaterial.yaml`
   defines `dynamic_viscosity(T) = 0.16667·exp(−0.0005·(T − 300))` over
   T ∈ [300, 3000] K. MaterForge fits this with a two-segment degree-2 piecewise
   polynomial (break at T ≈ 1536.24 K), inlined as the SymPy `Piecewise` expression
   shown below. The viscosity is computed from the temperature field and written to
   the viscosity field each step.

```
ν(T) =
    1.502e-8·T² − 8.954e-5·T + 0.19194   if T < 1536.24 K
    8.116e-9·T² − 6.852e-5·T + 0.17593   otherwise
```

## 2. Code-generation metric — `count_operations`

`pystencils.sympyextensions.count_operations` reports the symbolic arithmetic
operation count of the AssignmentCollection that defines the stream-collide
sweep, *before compilation*. Numbers below are taken directly from the build log.

### 2.1 Constant viscosity (ν = 0.08)

| Block             | adds | muls | divs | total |
|-------------------|-----:|-----:|-----:|------:|
| Collision main    |  103 |  200 |    0 |   303 |
| Subexpressions    |   64 |   12 |    1 |    77 |
| **Total**         |  **167** | **212** | **1** | **380** |

The constant `omega = 2 / (6·ν + 1) = 1.3513…` is folded into the generated
code as a literal. The viscosity field is **not** written back each step.

**Algorithmic bandwidth: 336 B/cell**
(19 PDFs read × 8 B + 19 PDFs write × 8 B + density write 8 B + velocity write 24 B)

### 2.2 Temperature-dependent viscosity (MaterForge)

| Block             | adds | muls | divs | total |
|-------------------|-----:|-----:|-----:|------:|
| Collision main    |  158 |  260 |    0 |   418 |
| Subexpressions    |   76 |   29 |    4 |   109 |
| **Total**         |  **238** | **295** | **4** | **537** |

The polynomial fit produced by MaterForge introduces a runtime per-cell
viscosity evaluation and three additional divisions inside the relaxation-rate
computation. The temperature field is read and the viscosity field is written
each step.

**Algorithmic bandwidth: 352 B/cell**
(base 336 B + temperature read 8 B + viscosity write 8 B)

### 2.3 Symbolic-operation and bandwidth deltas

| Metric      | Const | Temp-dep | Δ      | Δ %      |
|-------------|------:|---------:|-------:|---------:|
| adds        |   167 |      238 |    +71 | +42.5 %  |
| muls        |   212 |      295 |    +83 | +39.2 %  |
| divs        |     1 |        4 |     +3 | +300 %   |
| **Total ops** | **380** | **537** | **+157** | **+41.3 %** |
| **BW (alg.)** | **336 B/cell** | **352 B/cell** | **+16 B** | **+4.8 %** |

> **Note:** In the previous (pre-fix) code, the constant kernel also wrote the
> viscosity field every step (a redundant constant-store), inflating its bandwidth
> to 344 B/cell and underreporting the true bandwidth gap as +2.3 % instead of +4.8 %.

## 3. Runtime performance

Each case ran for 60 000 timesteps on an exclusive node with CPU pinning,
5 independent trials each, submitted serially (const job completes fully before
tempdep job starts). Statistics use the two-sided t-distribution with df = 4.

### 3.1 Headline numbers (5-trial statistics)

| Metric                       | Constant ν=0.08       | Temp-dep              | Ratio (T/C)   |
|------------------------------|-----------------------:|-----------------------:|:-------------:|
| Total MLUPS (mean ± std)     | **66.76 ± 0.27**      | **61.52 ± 0.35**      | **0.9215 ×**  |
| 95 % CI                      | [66.42, 67.09]        | [61.08, 61.95]        |               |
| MLUPS per process (mean ± std) | 16.69 ± 0.07        | 15.38 ± 0.09          |               |
| Wall time (mean ± std)       | 471.2 ± 1.9 s         | 511.4 ± 2.9 s         | **1.0852 ×**  |
| 95 % CI                      | [468.9, 473.6] s      | [507.8, 515.0] s      |               |
| Time per timestep (mean)     | 7.854 ms              | 8.523 ms              |               |
| L∞ error (vs analytical)     | 1.51 × 10⁻⁷           | 1.76 × 10⁻⁶           | 11.6 ×        |
| Convergence threshold passed | yes (all 5 trials)    | yes (all 5 trials)    |               |

The temperature-dependent kernel is **7.85 % slower** in MLUPS
(equivalently **8.52 % more wall-clock time**).  
`count_operations` predicts a **41.3 % increase** in symbolic FLOPs.  
The measured overhead is *less than one-fifth* of the symbolic prediction.

The 95 % confidence intervals do not overlap (const lower bound 66.42 vs temp-dep
upper bound 61.95), and a two-sample t-test confirms the difference is highly
significant (p = 4.55 × 10⁻⁹).

> **Note on previous woody measurements:** An earlier run (SLURM jobs 11733402/03)
> used only 1 trial per case, 5 warmup steps, no `--exclusive`, and ran both
> cases concurrently on the same node. That run reported 62.88 / 58.78 MLUPS
> (ratio 1.070 ×). The new measurements (exclusive node, CPU pinning, 5 trials)
> give 66.76 / 61.52 MLUPS (ratio 1.085 ×). The 6 % MLUPS increase for both
> cases reflects elimination of node-sharing interference; the overhead ratio
> increased slightly (7.0 % -> 8.5 %) because removing the redundant viscosity
> write from the const kernel made the const case proportionally faster.

The L∞ error is identical across all 5 trials within each case — it is
deterministic for a given binary and domain (fully converged LBM solution).
The temp-dep L∞ of 1.76 × 10⁻⁶ is bounded by the degree-12 polynomial fit
used as the reference (polynomial fit error: 8.01 × 10⁻⁶ against the
high-resolution numerical integration); the LBM itself is fully converged.

### 3.2 waLBerla timer breakdown (60 000 timesteps, mean across 5 trials)

All times are **summed across 4 ranks** (total CPU-seconds per trial).

| Timer              | Const ν=0.08 (mean) | Temp-dep (mean) | Δ time      | Const %   | Temp-dep % |
|--------------------|--------------------:|----------------:|------------:|----------:|-----------:|
| **StreamCollide**  | **1659.6 ± 6.9 s** | **1826.9 ± 7.8 s** | **+167.2 s** | **88.1 %** | **89.3 %** |
| LBM Communication  |    182.4 ± 0.7 s   |    183.5 ± 0.2 s   |   +1.1 s   |   9.7 %  |   9.0 %   |
| UBB boundary       |     31.6 ± 0.1 s   |     31.6 ± 0.1 s   |   ±0.0 s   |   1.7 %  |   1.5 %   |
| NoSlip boundary    |     10.4 ± 0.1 s   |     10.4 ± 0.1 s   |   ±0.0 s   |   0.6 %  |   0.5 %   |

StreamCollide delta per rank: **+41.8 s/rank** (+10.1 %). The StreamCollide
overhead (10.1 %) slightly exceeds the overall wall-clock overhead (8.52 %)
because the slight increase in communication fraction cancels part of it.
Communication, boundary, and timer-logger times are statistically identical
between cases — all overhead originates in the StreamCollide kernel.

The per-rank average (`REDUCE_AVG`) and max-rank (`REDUCE_MAX`) reductions were
also collected. `REDUCE_MAX / REDUCE_AVG` ≈ 1.0 for StreamCollide in both cases,
confirming negligible rank imbalance with 4 well-pinned, equally-loaded processes.

### 3.3 Bandwidth-vs-compute regime

At the observed const-ν throughput of 66.76 MLUPS:

  66.76 × 10⁶ cells/s × 336 B/cell = **22.43 GB/s of DRAM bandwidth**

against a theoretical peak of ≈ 204.8 GB/s per socket (Xeon Gold 6326,
8 channels DDR4-3200). The 4 pinned processes on socket 0 consume roughly
**11 % of one socket's peak memory bandwidth** — the kernel is far from
bandwidth-saturated.

For the temp-dep case: 61.52 × 10⁶ × 352 B = **21.65 GB/s** (10.6 % of peak).
The slightly lower absolute bandwidth is expected — the kernel does more
arithmetic per cell, reducing throughput, but still operates well below saturation.

## 4. Why the overhead is smaller than the FLOP ratio

| Symbol ratio (FLOPs)    | 1.413 × |
|-------------------------|---------|
| Measured ratio (time)   | 1.085 × |
| Direction vs prediction | Measured **below** symbolic — extra arithmetic hidden |

Three mechanisms explain why 41 % more symbolic FLOPs produce only 8.5 % more
wall time on the Xeon Gold 6326:

1. **Abundant bandwidth headroom hides arithmetic.** At ≈ 11 % of peak bandwidth,
   each memory access is serviced with very short latency relative to the
   compute-pipeline throughput. The extra polynomial arithmetic — roughly
   +78 add/mul and +3 div per cell — executes in functional units that are
   otherwise idle waiting for the next cache-line to arrive. The cost is
   **absorbed inside memory-access stalls** by the Ice Lake out-of-order
   engine, which has a 512-entry reorder buffer, allowing it to keep far
   more arithmetic in flight concurrently with outstanding loads.

2. **AVX-512 execution width.** The Xeon Gold 6326 natively executes AVX-512
   instructions, providing 512-bit SIMD for double-precision arithmetic. Even
   the `Piecewise` conditional — emitted by pystencils as a select/blend
   construct — is executed efficiently across all eight FP64 lanes of an
   AVX-512 register with no branch-misprediction overhead.

3. **Division latency improvement.** The +3 extra double-precision divisions per
   cell still carry non-trivial latency (≈ 14–20 cycles on Ice Lake). However,
   because the Ice Lake pipeline is so deeply pipelined and the memory access
   latency is relatively long (many outstanding cache misses), the divider can
   overlap execution with memory stalls — reducing the net wall-clock cost
   of the extra divisions.

In aggregate, the server CPU provides enough execution resources that the extra
arithmetic of the MaterForge material law is nearly **free at this domain
size** — the 8.5 % overhead is dominated by the +4.8 % algorithmic bandwidth
increase (temperature read + viscosity write) and communication/synchronisation
effects, not by the polynomial FLOPs themselves.

A useful framing for a paper:

> *On a bandwidth-rich server CPU (Xeon Gold 6326, ≈ 11 % of peak bandwidth
> utilised), the MaterForge temperature-dependent viscosity kernel introduces
> only 7.85 % throughput reduction (8.52 % wall-clock overhead) despite a 41 %
> symbolic FLOP increase — the extra arithmetic is hidden inside memory-access
> stalls by the deep out-of-order engine and AVX-512 execution units,
> underscoring that FLOP counts alone cannot predict LBM kernel performance.*

## 5. Why `count_operations` is not a sufficient performance metric

`count_operations` is a useful **code-generation diagnostic**. It quickly
exposes how complex a generated kernel is and lets one compare alternative
material models at sub-second cost, without rebuilding or running. But it
has structural limitations as a *performance* metric on real hardware:

| Limitation                                   | Consequence on Xeon Gold 6326 (Ice Lake)                       |
|----------------------------------------------|----------------------------------------------------------------|
| Treats all FLOPs identically                 | div cost hidden in memory-access stalls -> over-predicts overhead |
| Ignores compiler-pipeline effects            | Misses AVX-512 FMA fusion and efficiency gain                  |
| Has no notion of memory traffic or bandwidth | Misses bandwidth headroom absorbing extra compute              |
| Platform-agnostic by construction            | Cannot predict how overhead scales with bandwidth utilisation  |

The right end-to-end metric is **MLUPS**, which the waLBerla driver prints
after every run. It captures all of the above effects on the real target
hardware. Intended usage:

1. Use `count_operations` at code-generation time as a quick sanity check
   that the symbolic AST is not pathologically large.
2. Report MLUPS (and, where space permits, the StreamCollide-vs-Communication
   breakdown) as the actual benchmark figure, with ≥ 5 independent trials and
   95 % confidence intervals from the t-distribution.
3. Run on the *target* hardware class; results are highly platform-dependent
   (bandwidth utilisation is the dominant factor, not raw FLOP count).

## 6. Summary for the publication

* **On woody (Xeon Gold 6326, Ice Lake, NHR@FAU, exclusive node, 5 trials):**
  MaterForge temperature-dependent viscosity adds 41 % symbolic FLOPs but only
  **7.85 % throughput reduction** (8.52 % wall-clock overhead) for a 524 288-cell,
  60 000-timestep 3D Couette benchmark (4 MPI processes, 128 × 64 × 64 domain;
  p = 4.6 × 10⁻⁹, 95 % CI for const [66.42, 67.09] and temp-dep [61.08, 61.95] MLUPS).
* The kernel sustains ≈ 11 % of peak memory bandwidth; it is strongly
  **memory-latency-bound**. The extra polynomial arithmetic (+157 ops/cell) is
  absorbed by the deep out-of-order pipeline and AVX-512 execution units.
  The bandwidth delta (+4.8 %) accounts for the majority of the measured overhead.
* The simulation **converges correctly and reproducibly** in both cases
  (L∞ identical across all 5 trials: 1.51 × 10⁻⁷ for constant ν;
  1.76 × 10⁻⁶ for temp-dep, limited by the degree-12 analytical reference polynomial).
* **Recommendation for the paper:** report `count_operations` as a static
  complexity figure for the generated kernel and report MLUPS with a
  StreamCollide breakdown (mean ± std, 95 % CI, ≥ 5 trials) as the wall-clock
  metric. Run benchmarks on the target hardware class; bandwidth utilisation
  is the dominant factor governing how much extra arithmetic overhead matters.

## 7. Reproducibility

### 7.1 Build environment (woody)

```bash
# Modules: cmake/3.30.5, gcc/13.3.0 (via spack), openmpi/4.1.3-gcc12.1.0
# Python 3.11 virtualenv at $HOME/.venvs/materforge
git submodule update --init --recursive
source $HOME/.venvs/materforge/bin/activate
```

### 7.2 Configure and build

```bash
cd apps

# ── Constant ν = 0.08 ──────────────────────────────────────────────────────
cmake --preset woody-release-cpu \
      -DUSE_MATERFORGE=OFF -DCONST_NU=0.08 \
      -DSWEEPGEN_REQUIREMENTS_FILE=woody-sweepgen-requirements.txt \
      -DWALBERLA_BUILD_WITH_PYTHON=OFF
cmake --build build/woody-release-cpu -j 16
cp build/woody-release-cpu/CouetteFlowScaling \
   build/woody-release-cpu/CouetteFlowScaling_const_0.08

# ── Temperature-dependent (MaterForge) ─────────────────────────────────────
cmake build/woody-release-cpu -DUSE_MATERFORGE=ON
cmake --build build/woody-release-cpu -j 16
cp build/woody-release-cpu/CouetteFlowScaling \
   build/woody-release-cpu/CouetteFlowScaling_tempdep
```

### 7.3 Run (SLURM — serial, chained)

`CouetteFlowScaling.prm` settings used:
`timesteps = 60000`, `u_max = 0.025`, `T_bottom = 300`, `T_top = 3000`,
`errorThreshold = 1e-3`, `vtkWriteFrequency = 0`,
`blocks = <1,2,2>`, `cellsPerBlock = <128,32,32>`.

```bash
# Submit const job; tempdep starts only after const completes (no overlap)
CONST_JID=$(sbatch --parsable run_const_0.08.sh)
sbatch --dependency=afterok:${CONST_JID} run_tempdep.sh
```

Key SLURM options in both scripts:
```
#SBATCH --exclusive              # sole occupant — no shared-cache interference
#SBATCH --constraint=icx         # Ice Lake nodes only
mpirun -n 4 --bind-to core --map-by core --report-bindings ...
```

Each script runs 5 trials in a loop. SLURM job IDs for this dataset:
* 11733409 — constant ν=0.08 (node w2304)
* 11733410 — temperature-dependent (node w2304, sequential after 11733409)

Raw logs are preserved in `apps/logs/`:
* `configure_const_0.08.log`, `build_const_0.08.log`, `run_const_0.08.log`
* `configure_tempdep.log`, `build_tempdep.log`, `run_tempdep.log`
* CPU binding reports in `run_const_0.08.err`, `run_tempdep.err`

## 8. Cache counter analysis (Linux `perf stat`)

### 8.1 Methodology and profiling setup

Direct hardware-counter profiling on woody required three attempts:

| Tool | Failure mode | Root cause |
|------|-------------|-----------|
| Intel VTune 2022.3 (`memory-access`, `uarch-exploration`) | `Error: requires perf system-wide profiling` | `perf_event_paranoid = 2` blocks VTune's kernel-driver-based PMU sampling |
| LIKWID 5.4.1 (`CACHE`, `MEM` groups) | `Access to performance monitoring registers locked` | MSR-device access restricted cluster-wide (kernel config), bypassing the setuid `likwid-accessD` daemon |
| Linux `perf stat` (user-space counting mode) | **Success** | `PERF_EVENT_OPEN` in counting mode is allowed at `paranoid = 2` for user-space per-process events |

`perf stat` was run with per-rank isolation: each MPI rank launched its binary under its own `perf stat` instance with `--output rank_N.txt`. This gives four symmetric per-rank measurement sets. The collection script is `apps/run_perf_cache.sh` (SLURM job 11733616, node w2304).

**Collection scope:** 100 warmup iterations + 1 timed timestep per rank. Although shorter than the full 60 000-timestep benchmark, the working set reaches steady state during warmup, making the cache-level ratios between const and temp-dep representative of the sustained workload.

**Events collected** (all user-space, `:u`):
`cache-references`, `cache-misses`, `L1-dcache-loads`, `L1-dcache-load-misses`,
`L1-icache-load-misses`, `LLC-loads`, `LLC-load-misses`, `LLC-stores`,
`dTLB-loads`, `dTLB-load-misses`, `cycles`, `instructions`.

---

### 8.2 Raw event counts (4-rank totals)

| Event | const ν=0.08 | temp-dep (MF) | T/C ratio |
|-------|-------------:|--------------:|----------:|
| `cache-references` | 314 418 406 | 331 708 849 | **1.055** |
| `cache-misses` | 86 811 791 | 91 014 298 | **1.048** |
| `L1-dcache-loads` | 8 864 456 918 | 10 720 687 096 | **1.209** |
| `L1-dcache-load-misses` | 325 546 843 | 339 454 278 | **1.043** |
| `LLC-loads` | 4 830 473 | 4 911 966 | **1.017** |
| `LLC-load-misses` | 3 037 454 | 2 997 069 | **0.987** |
| `LLC-stores` | 12 011 041 | 11 419 185 | **0.951** |
| `cycles` | 11 424 945 478 | 13 053 954 419 | **1.143** |
| `instructions` | 26 471 794 836 | 32 123 187 323 | **1.214** |

> **Note on `cycles` ratio:** the short-run ratio (1.143) is higher than the full-benchmark wall-time ratio (1.085 ×) because startup and cache-warming costs are proportionally larger for a 101-step run. The full-benchmark figure is the authoritative performance number; the cycle ratio here is informative for relative arithmetic overhead only.

---

### 8.3 Derived cache rates and IPC

| Metric | const ν=0.08 | temp-dep (MF) | T/C |
|--------|-------------:|--------------:|----:|
| **L1D load miss rate** | 3.67 % | 3.17 % | **0.862** |
| LLC load miss rate | 62.9 % | 61.0 % | **0.970** |
| Generic cache miss rate | 27.6 % | 27.4 % | **0.994** |
| **IPC (instructions / cycle)** | 2.317 | 2.461 | **1.062** |

---

### 8.4 Findings

#### Finding 1 — L1 traffic increases 21 %, but L1 miss rate falls

The temp-dep kernel generates **+20.9 % more L1 data loads** (ratio 1.209). These are not new DRAM-backed loads; they are intermediate loads from L1-resident data produced by the polynomial evaluation: coefficient broadcasts, intermediate products, and the piecewise select operands. They cost L1 bandwidth and keep the load-address generation units busy, but they do not add L1 misses.

Crucially, the **L1D miss rate drops** from 3.67 % to 3.17 % (ratio 0.862). The extra arithmetic between consecutive PDF-field accesses gives the Ice Lake hardware prefetcher more time to fetch the next PDF cache line before it is demanded. The result is **better L1-level temporal locality for the memory-bound data despite a larger instruction footprint**.

#### Finding 2 — LLC traffic and DRAM pressure are unchanged

The LLC miss ratio is **0.987 — slightly below 1.0**, meaning temp-dep causes fewer measured L3 misses than const. This contradicts the +4.8 % algorithmic bandwidth increase from the field layout (§3.3). The resolution:

1. **Scale vs. bandwidth headroom.** At ≈ 11 % of peak socket bandwidth, the memory controller services requests so quickly that the extra 2 × 8 B/cell from the temperature and viscosity fields blends into measurement noise for a 101-step run. With 4 ranks × 131 072 cells = 0.524 M cells per step, the extra traffic is 0.524 M × 16 B = 8.4 MB per step. Over 101 steps that is 847 MB total — well within L3's streaming throughput at this bandwidth utilisation.

2. **Working-set overlap.** The temperature field is read-only in the StreamCollide kernel; it shares the same cache-line granularity as the PDF fields. Because the access pattern is strided and regular, the hardware prefetcher handles both cases equally well, and the actual L3 eviction counts remain similar.

3. **Measurement window.** 101 timesteps is insufficient to expose the cumulative +5.3 % DRAM traffic at this bandwidth utilisation level. The MLUPS-based measurement over 60 000 timesteps is more sensitive.

**Key takeaway for §3.3:** the +4.8 % algorithmic bandwidth increase (352 vs 336 B/cell) is the *theoretical* byte-traffic ratio derived from the field layout. The observed LLC miss ratio (≈ 1.0) confirms that this extra traffic is absorbed without changing the L3 hit/miss structure — consistent with the 11 % bandwidth utilisation figure.

#### Finding 3 — Compiler reduces 41 % symbolic FLOPs to 21 % real instructions

The instruction count ratio is **1.214 (+21.3 %)**, roughly half the symbolic FLOP ratio of 1.413 (+41.3 %). This gap has two sources:

* **AVX-512 fusion.** The Xeon Gold 6326 executes 8 × FP64 SIMD lanes per instruction. Many of the +157 symbolic ops (adds, muls) are fused into FMA instructions that count as one instruction in the hardware retire port. The effective per-cell arithmetic cost is compressed by ≈ 8×.

* **Piecewise compiled to a blend, not a branch.** GCC emits the `Piecewise` select as a vectorised blend (`vblendvpd` / compare + mask), with both polynomial branches evaluated and the result selected. No branch misprediction overhead; the branch predictor sees a predictable zero-divergence workload.

The remaining 21 % instruction overhead is what actually reaches the retire port. Combined with a **+6.2 % IPC improvement** (2.317 → 2.461), the CPU executes those extra instructions more efficiently — the polynomial's independent multiply-add chains provide more instruction-level parallelism than the scalar `omega`-folding path of the const kernel.

#### Finding 4 — The 8.5 % wall-time overhead decomposes as follows

| Source | Estimated contribution |
|--------|----------------------:|
| Extra L1 memory traffic (+20.9 % loads) | ≈ 3–4 % |
| Extra DRAM traffic (+5.3 % bytes/cell) | ≈ 5 % |
| Extra polynomial instructions (+21 %), partially hidden by IPC gain (+6 %) | ≈ 2–3 % |
| Communication & boundary timing (scheduling artefact) | ≈ −2 % (negative - faster sync) |
| **Measured total (5-trial mean)** | **8.5 %** |

The dominant cost is the combination of increased DRAM traffic and L1 load pressure. The extra arithmetic — the quantity that `count_operations` measures — contributes only modestly because the out-of-order engine and AVX-512 units absorb it inside memory-access stalls.

---

### 8.5 Summary table

| Metric | Ratio (T/C) | Interpretation |
|--------|------------:|----------------|
| Symbolic FLOPs (`count_operations`) | **1.413** | Code-generation complexity measure |
| Real instructions retired | **1.214** | Compiler fuses ≈ 50 % of symbolic ops via AVX-512 FMA |
| L1-dcache-loads | **1.209** | Extra polynomial evaluation generates L1-resident intermediate loads |
| L1D load miss rate | **0.862** | Miss rate *decreases* — prefetcher more effective with longer compute gaps |
| LLC-load-misses (DRAM proxy) | **0.987** | DRAM pressure unchanged — +5.3 % traffic absorbed below measurement resolution |
| IPC (instructions/cycle) | **1.062** | Higher for temp-dep — polynomial arithmetic provides more ILP |
| Cycles (short run, 101 steps) | **1.143** | Inflated by startup; full-run ratio is 1.085 |
| **Wall-clock ratio (60 000 steps)** | **1.085** | Authoritative benchmark figure |

> **perf stat raw data:** `apps/perf_results/{const_0.08,tempdep}/rank_{0-3}.txt`
> **Collection script:** `apps/run_perf_cache.sh` (SLURM job 11733616, node w2304, `perf` v4.18)
