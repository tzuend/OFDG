# Comparative DG, OFDG, and OEDG study in MFEM

This repository contains the implementation and reproducible numerical study
for the thesis *Oscillation-Free Discontinuous Galerkin Stabilization in MFEM*.
It compares four configurations:

- `dg`: conventional RKDG without an oscillation-elimination step;
- `ofdg-kxrcf`: the adapted OFDG filter restricted to KXRCF troubled cells
  (the primary project method);
- `oedg`: the 2024 OEDG method of Peng, Sun, and Wu;
- `ofdg`: the adapted OFDG filter on every cell (an ablation).

The implementation favors readable, shared building blocks over duplicated
drivers. MFEM supplies meshes, DG spaces, flux assembly, and the normal explicit
solvers. The scalar and Euler examples share face-local physics, filter
dispatch, MPI ownership rules, diagnostics, and study output.

## Requirements

- C++17 and an MPI compiler wrapper;
- MFEM built with MPI, Hypre, and METIS;
- `make` and Python 3 with NumPy and Matplotlib;
- `latexmk` for the thesis report.

The default layout expects MFEM at `../mfem`. Override it when needed:

```sh
make MFEM_DIR=/path/to/mfem test
```

## Build and validation

```sh
# Address/undefined-behavior sanitizers plus the two-rank MPI test
make test

# Optimized unit tests and independent-study validation
make test-release
make test-study

# Optimized example programs and the WENO5 reference solver
make release reference-solver
```

The suite covers frozen OFDG fingerprints, affine and mixed mesh geometries,
mean preservation, OEDG attenuation and scale/evolution invariance, KXRCF
component pooling, Euler normal wave speeds, positivity limiting, rejected RK
steps, synthetic rate/plot calculations, a smooth WENO convergence check, an
exact Sod Riemann check, and shared-face consistency on two MPI ranks.
The MPI target also runs bounded one- versus two-rank comparisons through all
maintained example drivers, so a mismatched collective fails in CI instead of
hanging indefinitely.

## Common command-line interface

Every example accepts:

```text
-method dg|ofdg|ofdg-kxrcf|oedg
-cadence auto|step|stage
```

`auto` uses one post-step filter application with MFEM's standard solver for
OFDG and stage-wise filtering for the paper-faithful OEDG comparison. The small
experiment RK driver is also used when positivity checks must occur at every
stage. The older `-stab` and `-kxrcf` flags remain compatibility aliases.

Representative runs are:

```sh
# 1D smooth advection and the Gauss--Radau superconvergence experiment
build/release/examples/advection -d 1 -n 80 -o 2 -tf 1.1 -method oedg
build/release/examples/advection -d 1 -n 80 -o 2 -tf 1.1 -radau-init -method dg

# 2D periodic advection on triangles
mpirun -np 2 build/release/examples/advection -d 2 -n 40 -o 2 \
  -triangular -method ofdg-kxrcf

# Smooth and post-shock Burgers problems
build/release/examples/burgers -p 1 -n 160 -o 3 -tf 0.6 -method dg
build/release/examples/burgers -p 2 -n 160 -o 2 -tf 1.5 -method oedg

# Transported Euler vortex and the first 2D Euler Riemann configuration
build/release/examples/euler -p 1 -r 1 -o 2 -tf 1 -method ofdg-kxrcf -no-vis
mpirun -np 4 build/release/examples/euler -p 6 -o 1 -tf 0.25 \
  -method oedg -positivity -no-vis
```

Euler problems are: transported vortices (1--2), a periodic 2D density wave
(3), the smooth 1D wave from the 2024 paper (4), Woodward--Colella blast data
(5), the first 2D Riemann problem (6), scaled Lax data (7), and Shu--Osher (8).
Euler face damping uses the local normal spectral radius
`abs(velocity · normal) + sound_speed`; it never assumes unit wave speed.

For hard Euler tests, `-positivity` enables the common cell-average-preserving
scaling limiter. It checks volume and face evaluation states, reports limited
cells and rejected steps, and rejects a step if its cell mean is inadmissible.
Unassisted failure runs use `-no-positivity` and are recorded separately.

## Study-design phase and reproducibility guard

The JSON manifest in `experiments/study_manifest.json` is the source of truth.
It expands into commands and records raw logs, rank-local profiles, long-form
CSV, summarized CSV, generated LaTeX, and PDF/PNG figures.

The thesis-final matrix is intentionally still marked `draft`.  At this stage,
inspect the expanded experiment design without launching a solver:

```sh
# Expand all proposed rows and estimate their runtime; run no measurements
make study-plan

# Refuses to build or execute while the full profile is not approved
make study-full

# Regenerate publication-style draft assets from the existing quick checkpoint
make figures-draft
```

`measurements/study/plans/full.csv` is the approval-ready matrix.  Every row
records the scientific question, equations and domain, methods, numerical
settings, metrics, independent-reference check, priority, intended artifact,
estimated runtime, and proposed command.  Changing the full-profile status to
`approved` is a deliberate later step after this design has been reviewed.

`make study-quick` remains available as a development and CI checkpoint.  It
does execute a reduced matrix, but its values are never thesis-final evidence.
Its analysis can write only to `report/generated/draft`; the generator rejects
any attempt to write quick data into `report/generated/final`.

When it is eventually approved, the full profile will perform one unrecorded
warm-up and five recorded timing repetitions. Every recorded row includes the
actual DOFs, degree, CFL/time step, RK method, cadence, ranks, three error norms,
conservation drift, extrema, reference-profile error, total variation,
filter/limiter activity, failure information, runtime, compiler/MFEM revisions,
commit, and host. The runner checkpoints the CSV after each completed row.

The independent component-wise WENO5/local-Lax--Friedrichs/SSPRK3 solver in
`benchmarks/euler_reference_weno.cpp` supplies 1D discontinuous references; no
compared DG method is used as its own reference. The optional `hpc` manifest
contains workloads that are intentionally excluded from routine local runs.

## Visualization

Study figures are generated headlessly with common plot scales. All maintained
example drivers (`advection`, `burgers`, and `euler`) also support optional
live GLVis output. Start the GLVis server in one terminal and pass `-vis` to an
example in another; `-vs` controls the timestep interval between updates:

```sh
../glvis-4.5/glvis
build/release/examples/advection -d 2 -n 40 -o 2 -vis -vs 10
build/release/examples/burgers -p 2 -n 160 -o 2 -vis -vs 10
mpirun -np 2 build/release/examples/euler -p 1 -r 1 -o 2 -vis -vs 10
```

Visualization is disabled by default, so studies and CI remain headless. If a
requested GLVis connection cannot be established on every MPI rank, the
example reports that fact and continues without visualization.

Euler can additionally write rank-local MFEM field files with `-save` for
later inspection. A compatible checkout is available at either
`../glvis-4.5` or `/Users/timozund/Documents/master_thesis/glvis-4.5`; for
example:

```sh
../glvis-4.5/glvis -m euler-mesh.000000 -g euler-0-final.000000
```

GLVis is useful for interactive spot checks, while the thesis figures remain
script-generated so regeneration does not depend on a display server.

## Report modes and principal files

```sh
# Draft report with unmistakably preliminary quick-checkpoint assets
make report-draft

# Refuses to build until the full manifest is approved and final generated
# Results, Discussion, and Conclusion files all exist
make report-final
```

`make report` is an alias for `make report-draft`; its output is
`report/thesis.pdf`.  The draft Abstract makes no numerical claim, Results
contains a prominent preliminary-data warning, and Discussion and Conclusion
contain delimited finding placeholders.

The convergence generator follows the reference-paper layout: separate tables
for RKDG, OFDG--KXRCF, and OEDG; multirow degree blocks; all mesh levels; paired
errors and rates; scientific notation; and no compressed `resizebox`.  Draft
figures separate polynomial degrees and problems, limit each axis to the three
headline methods plus one reference, and split 2D contours from line cuts.

Important implementation files are:

- `src/ofdg_serial_optimized.hpp`: adapted OFDG, immutable caches, MPI halos,
  and exact modal decay;
- `src/oedg_2024.hpp`: fixed 2024 OEDG choices (interface L1 jumps,
  trapezoidal 2D faces, face heights, component pooling, and scale invariance);
- `src/kxrcf.hpp`: troubled-cell indicator and profiling;
- `src/euler_positivity.hpp`: conservative Euler positivity limiter;
- `src/face_physics.hpp`: advection, Burgers, and Euler propagation policies;
- `scripts/run_study.py` and `scripts/analyze_study.py`: reproducible runs,
  metrics, tables, and figures;
- `tests/`: deterministic serial, MPI, and study-level validation.

The full local study is intentionally smaller than the published $1280^2$
suite. Double-Mach reflection, the Mach-2000 jet, and the largest grids are HPC
workloads, not prerequisites for rebuilding this thesis checkpoint.
