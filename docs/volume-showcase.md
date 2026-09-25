# Small 3D volume OFDG showcase

This isolated benchmark links the existing production filter library. No production source was changed. It demonstrates affine hexahedral Q2 scalar advection, not surface transport, curved 3D, Euler, or all supported element families.

## Reproduce

From codex_dev, with the existing MFEM installation:

```sh
make build/release/benchmarks/volume_showcase -j1
python3 scripts/run_volume_showcase.py
python3 scripts/analyze_volume_showcase.py
python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner resumes matching completed data and rejects changed source fingerprints. Preserve the output directory before an intentional new run. Plotting requires NumPy/Matplotlib; PDF generation uses ReportLab. Native initial/final MFEM mesh and grid functions can be viewed with GLVis, e.g.:

```sh
glvis -m measurements/volume-showcase/visualization/n10-ofdg-p1-dt0.003-final.mesh -g measurements/volume-showcase/visualization/n10-ofdg-p1-dt0.003-final.gf
```

## Problem and numerics

Periodic unit cube; velocity (0.5,0.3,0.2); T=0.3. Pulse: indicator of the radius-0.2 ball centered at (0.30,0.35,0.40), with periodic distance. Exact solution is translated initial data. Smooth control: 1+0.2 sin(2 pi x) sin(2 pi y) sin(2 pi z). Constant control: 1.

Q2 tensor-product DG on 6^3 and 10^3 affine hexes, 27 unknowns/cell. MFEM Rusanov flux and classical RK4; fixed dt=0.003 and OFDG after each full step. Identical nodal initialization; no positivity limiter or KXRCF. Finest pulse OFDG repeats at dt=0.0015. The existing StudyFilter uses the production derivatives and nested damping unchanged.

Order-16 positive tensor quadrature measures errors and sampled extrema. Pulse L2 is approximate because the exact solution is discontinuous; no quadrature-sensitivity certification is claimed. Conservation compares each discrete final mass to its discrete initial mass, not the analytic ball volume. Five interior samples per axis per cell are exported in physical coordinates. Plots use actual polynomial values, with line segments broken at cell boundaries. The 3D threshold point cloud is qualitative; the common-scale slice and enlarged cut reveal excursions it hides.

## Evidence and limits

See results.csv, run.json (exact commands, source/executable hashes, memory and times), summary.json, raw/, visualization/, figures/, source/. The runner checks available logical CPUs and uses just one rank and one numerical-library thread. Sequential checks are killed above five minutes or 4 GiB aggregate sampled RSS, and launching stops at 45 minutes accumulated compute. All ten runs passed in 52.12 seconds, peak sampled RSS 80.92 MiB.

OFDG suppresses excursions but increases pulse L2 error through smearing. Smooth errors fall on refinement; two grids do not establish asymptotic order. No performance ranking follows from these single runs. This is a compact capability showcase, not a new OFDG formulation or a general 3D validation suite.
