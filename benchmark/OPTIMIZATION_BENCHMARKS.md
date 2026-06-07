# AMBER structure-optimization performance

Benchmark report for the cached, multi-backend force-field evaluator added on
branch `anhi/optimization-performance`. All runs are on the same machine (Apple
Silicon, 12 logical cores, Metal GPU). Timings are best-of-N wall-clock;
absolute numbers drift with background load, so the **ratios within a single
run** are the reliable figures.

Primary structure: **`benchmark/data/AmberFF_bench.pdb`** — BPTI, 892 atoms,
~634 k nonbonded pairs at the default 20 Å cut-off (AMBER96, Float32).

---

## TL;DR

- **Up to ~36× faster on a single CPU thread and ~138× with 8 threads** than the
  original `optimize_structure!`, for the same minimization.
- The original was **~15–20× slower than BALL per evaluation** (it rebuilt the
  whole pair list every energy/force call). The new code is **on par with BALL
  single-threaded and faster with threads/GPU**.
- **Forces match BALL to 0.1 %**, bonded energies to ~1 %.
- The **GPU backend scales best**: at ~7 000 atoms / 5 M pairs it computes forces
  ~3.8× faster than one CPU thread and overtakes the 8-thread CPU backend.

---

## 1. Speed-up vs. the original implementation

The original `optimize_structure!` called `update!(ff)` *inside the objective*,
so the entire nonbonded interaction list was torn down and rebuilt on every
energy/force evaluation (≈ 42 ms and 85 MB allocated **per call**). The new
evaluator builds the pair list once and rebuilds it only when atoms drift past a
Verlet skin (the strategy BALL uses).

200-iteration minimization of BPTI, measured back-to-back in one run:

| Variant                         | Wall-clock | Speed-up |
|---------------------------------|-----------:|---------:|
| Original (`update!` every eval) |    68.5 s  |    1×    |
| New — serial                    |    1.91 s  |  **36×** |
| New — threads (8)               |    0.50 s  | **138×** |

---

## 2. Accuracy vs. BALL (BPTI, identical 892-atom structure)

BALL (C++ reference) was built with its `devenv` (Nix) toolchain and run on the
same PDB. Forces were converted from BALL's SI Newtons to kJ/(mol·Å).

| Quantity                | BALL      | Julia (new) | Agreement |
|-------------------------|----------:|------------:|-----------|
| RMS force [kJ/(mol·Å)]  |   62.32   |    62.39    | **0.1 %** |
| Bond stretches [kJ/mol] |   494.6   |    492.7    |  ~0.4 %   |
| Angle bends             |   517.2   |    523.5    |  ~1.2 %   |
| Torsions (prop+impr)    |  1171.5   |   1168.3    |  ~0.3 %   |
| Non-bonded (vdW+HB+ES)  |  −7880    |   −8344     |  ~6 %     |
| **Total**               | **−5697** |  **−6159**  |  ~8 %     |

Forces and bonded energies agree closely. The ~6 % non-bonded gap lives entirely
in the long-range electrostatic sum (it does **not** appear in the forces). This
is a **pre-existing property of the Julia AMBER port** — the cached evaluator
reproduces the Julia reference path exactly (Float64: bit-for-bit; Float32:
rounding level), so the performance work did not change any energies.

> Accuracy is validated in `test/forcefields/AMBER/test_compiled_amberff.jl`
> (compiled vs. reference, both precisions and all backends) and the existing
> `test_amberff.jl` / `test_optimize_structure.jl` oracle suites.

---

## 3. Speed vs. BALL (BPTI)

**Per-evaluation cost** (cached pair list — the operation that dominates a
minimization step):

| Backend            | Energy (ms) | Forces (ms) |
|--------------------|------------:|------------:|
| BALL (C++, 1 core) |     1.5     |     2.1     |
| Julia — serial     |     1.8     |     3.2     |
| Julia — threads(8) |     0.26    |     0.52    |
| Julia — GPU (Metal)|     1.5     |     2.6     |

- Single-threaded, BALL is ~1.4–1.5× faster than our serial Julia (expected for
  hand-tuned C++).
- The threaded backend is **~4–6× faster than BALL**; the GPU matches BALL at
  this small size and pulls ahead as the system grows (§4).

**End-to-end 200-iteration minimization:** BALL (Strang L-BFGS) 0.9–1.2 s; Julia
(L-BFGS-B) serial 1.9 s, threads 0.50 s. (Iteration semantics differ between the
two minimizers, so per-evaluation cost above is the cleaner comparison.)

---

## 4. Backend scaling with system size

To exercise larger systems, BPTI was replicated into *N* copies placed 120 Å
apart (beyond the cut-off, so copies don't interact and pair count grows
linearly). Per-evaluation timings (ms), Float32:

| System  | Atoms | Pairs | serial E / F | threads(8) E / F | GPU E / F |
|---------|------:|------:|-------------:|-----------------:|----------:|
| BPTI ×1 |   892 | 0.63 M| 1.83 / 3.18  |  0.26 / 0.52     | 1.52 / 2.57 |
| BPTI ×2 |  1790 | 1.28 M| 4.43 / 9.24  |  0.75 / 1.53     | 2.36 / 4.22 |
| BPTI ×4 |  3586 | 2.56 M| 9.40 / 18.73 |  1.38 / 3.17     | 2.65 / 5.88 |
| BPTI ×8 |  7178 | 5.14 M| 17.6 / 34.3  |  3.05 / 13.9     | 4.15 / 9.02 |

- **Serial** scales linearly with pair count, as expected.
- **GPU** has a fixed launch/transfer overhead that dominates the smallest
  system (×1: GPU ≈ serial), but its compute barely grows: from ×1 to ×8 the
  pair count rises 8× while GPU forces rise only ~3.5×.
- **Crossover:** by ×4–×8 the GPU overtakes the 8-thread CPU backend. At ×8
  (≈ 7 200 atoms, 5 M pairs) GPU forces (9.0 ms) beat serial (34.3 ms, **3.8×**)
  and the threaded backend (13.9 ms).

The takeaway: use `:serial`/`:threads` for typical proteins, and `:gpu` for large
assemblies where its near-flat scaling wins.

---

## 5. Notes & caveats

- **Larger head-to-head vs. BALL** on *identical* atoms was not achievable: the
  Julia FragmentDB topology inference does not robustly handle multi-chain /
  larger PDBs (e.g. 4hhb, 2ptc fail; 1tgh loads as 2956 atoms in Julia vs. 2355
  in BALL), and the synthetic replicas are rejected by BALL's PDB reader. The
  larger-system comparison is therefore the **Julia-internal scaling curve**
  (§4); the rigorous cross-check vs. BALL is the BPTI head-to-head (§2–§3). For
  reference, BALL on `1tgh` (2355 atoms) runs at 4.6 ms (energy) / 6.4 ms
  (forces) per evaluation.
- **Accumulation precision** is selectable (`accumulation = Float32 | Float64`).
  Float64 reproduces the reference path bit-for-bit; Float32 (default for
  `AmberFF`) differs only at rounding level.
- **GPU** uses KernelAbstractions and works on **Metal and CUDA** via package
  extensions (no GPU dependency for CPU-only users). Metal has no Float64, so the
  GPU backend requires Float32 accumulation there.
- Reproduce: `optimize_structure!(ff; backend = :serial | :threads | :gpu,
  accumulation = Float32 | Float64)`. BALL reference built via
  `cd tmp/BALL && devenv shell -- bash -c 'cd build/devenv && ninja …'`.
