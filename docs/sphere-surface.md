# Screened Poisson on an embedded spherical surface

This is an isolated surface-FEM example for the supervisor brief, based on the installed `../mfem/examples/ex7.cpp`. It solves on a two-dimensional surface embedded in R^3, not on a three-dimensional volume. The aim is to demonstrate a verified solve on a curved surface. Production DG/OFDG, all previous experiments and the thesis are unchanged.

An earlier supervisor PDF appended two pages (15-16) containing the sphere solution/error rendering, the surface mathematical construction, convergence and geometry checks. The first fourteen pages are unchanged. The earlier FEM-containing PDF is archived under `tmp/before-surface-ofdg/`. Those pages have now been replaced by the surface-DG/OFDG prototype; this FEM dataset is preserved.

## Problem and discretization

On the unit sphere Gamma, solve

    -Delta_Gamma u + u = f,    u(x,y,z) = xy,    f(x,y,z) = 7xy.

This is the same manufactured solution as MFEM Example 7. Since xy is a degree-two spherical harmonic, Delta_Gamma(xy)=-6xy. The sphere is closed and has no boundary. The reaction term makes the problem coercive without a mean-zero constraint.

Mesh geometry is polynomial. For triangles, start from an octahedron with eight faces; for quadrilaterals, start from the six cube faces. Set geometry degree g, uniformly refine, then snap all final nodal coordinates to the unit sphere. This matches Example 7's default snap-at-the-end construction. The curved surface Gamma_h only approximates Gamma between the nodes. Solution and geometry degrees may be set independently.

For a reference-to-surface map F, J=dF/dxi is a 3-by-2 matrix and G=J^T J is a 2-by-2 metric. The surface gradient and measure are

    grad_Gamma_h v = J G^{-1} grad_reference(v_hat),
    dS = sqrt(det G) dxi.

The discrete variational equation is

    integral_Gamma_h [grad u_h . grad v_h + u_h v_h] dS
      = integral_Gamma_h f_ext v_h dS.

The solution space is conforming H1, not DG. MFEM's DiffusionIntegrator, MassIntegrator and DomainLFIntegrator assemble this problem. The solve uses preconditioned conjugate gradients and a symmetric Gauss-Seidel preconditioner, relative tolerance 1e-12, absolute tolerance 1e-14, maximum 1000 iterations. Relative algebraic residuals are recomputed independently from A*x-b.

The radial extensions to Gamma_h are

    u_ext(x,y,z) = xy / (x^2+y^2+z^2),    f_ext = 7u_ext.

The continuum forcing identity applies on the exact unit sphere. On the approximate surface, geometry introduces a consistency error. Consequently the computed error contains both approximation and geometry effects.

## Configurations

28 sequential solves:

- Both triangles and quadrilaterals; p=g=1,2,3; refinement levels 1,2,3.
- One extra refinement (level 4) for p=g=2 on both families.
- Hold p=2, level 3 fixed and compare geometry g=1,2,3.
- At p=g=2, level 3, increase both assembly and evaluation quadrature orders by six.
- Constant control u=f=1 at p=g=2, level 2 on both families.

Element counts are 8*4^level for triangles and 6*4^level for quadrilaterals. Assembly quadrature order is 2p+2g+6 (+extra); error quadrature is eight orders higher. Full configurations and exact commands are in `measurements/sphere-surface/run.json`.

## Errors and results

Volume here means surface area; all error integrals are over Gamma_h. We record absolute and relative surface L2 error against u_ext, tangential H1-seminorm error against grad_Gamma_h u_ext, energy error, sampled maximum error, surface area error against 4*pi, radial RMS and sampled maximum distance to the unit sphere, algebraic residuals, and the reaction balance defect abs(integral_Gamma_h(u_h-f_ext)). The latter checks the constant-test-function balance, not a time-dependent mass-conservation claim.

The reference gradient is calculated from the analytic radial extension and projected onto the actual discrete element tangent plane. The reported norms are not integrals of a lifted error over the exact sphere. Numerical samples do not certify a continuous maximum norm.

| Family | Degree p=g | Finest elements | Unknowns | Relative L2 error | L2 rate |
|---|---:|---:|---:|---:|---:|
| Triangles | 1 | 512 | 258 | 3.136% | 1.90 |
| Triangles | 2 | 2048 | 4098 | 0.01084% | 2.96 |
| Triangles | 3 | 512 | 2306 | 0.005716% | 4.12 |
| Quads | 1 | 384 | 386 | 1.596% | 1.93 |
| Quads | 2 | 1536 | 6146 | 0.005587% | 2.98 |
| Quads | 3 | 384 | 3458 | 0.002732% | 4.18 |

Rates use log2(E_coarse/E_fine) for absolute L2 error on the final two refinements. Observed H1-seminorm rates are 0.99/1.00 for p=1, 1.98/1.99 for p=2, and 3.07/3.14 for p=3 (triangles/quads). Different families and degrees have different unknown counts; the table is not a matched-cost ranking.

At level 3, p=2:

| Family | Relative L2, g=1 | Relative L2, g=2 | Relative L2, g=3 |
|---|---:|---:|---:|
| Triangles | 1.253% | 0.08434% | 0.1143% |
| Quads | 0.6944% | 0.04410% | 0.07463% |

The planar geometry is substantially less accurate. Geometry g=3 has lower area error than g=2, but its solution error is higher in these cases. Geometry and solution approximation errors can interact or cancel; the data do not establish which contribution dominates each case. Increasing geometry degree is not a guarantee of monotonically decreasing solution error for every fixed discretization.

## Controls and resources

All 28 cases passed. Constant relative L2 error is at most 4.99e-14. Increasing quadrature changes absolute L2 error by at most 1.02e-13 (the max(1,error) scaling is one for these cases). Maximum algebraic relative residual is 1.38e-12. Maximum reaction balance defect is 2.00e-12. Geometry measure is positive at all evaluated points.

Eight logical processors were confirmed. Runs used one MPI rank and one numerical-library thread, sequentially, with no oversubscription. The runner enforces five minutes/check, 4 GB sampled process-tree RSS, and stops launching after 45 minutes accumulated numerical execution. All 28 cases completed in 5.61 seconds monitored wall time, with maximum sampled RSS 23.1 MiB. This includes visualization export for the two representative runs. Compilation, plotting and PDF generation are excluded. These short runs are not a performance benchmark.

## Visualizations and reproduction

`measurements/sphere-surface/visualization/` contains `tri.mesh`, `tri.gf`, `quad.mesh`, `quad.gf`, and their sampled-surface CSV files, all at p=g=2 and refinement level 3. These are native MFEM files usable with GLVis, for example:

```sh
glvis -m measurements/sphere-surface/visualization/tri.mesh -g measurements/sphere-surface/visualization/tri.gf
```

The PDF renderings were produced with Matplotlib from samples of the actual curved MFEM surface, not from ideal-sphere coordinates. Each element is subdivided for display only. Patch colors are means of the three sampled vertex values. The displayed mesh edges belong to original elements, not to the rendering subdivision. The right panel uses its own symmetric signed-error scale. This stationary problem has one computed state rather than an initial/final pair.

Reproduce from codex_dev using the existing MFEM installation:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make build/release/benchmarks/sphere_surface
/opt/homebrew/bin/python3 scripts/run_sphere_surface.py
MPLCONFIGDIR=/tmp/ofdg-matplotlib /opt/homebrew/bin/python3 scripts/analyze_sphere_surface.py
/Users/timozund/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner resumes only when source/executable fingerprints match. Preserve/rename the measurements directory before a fresh rerun with changed code. Raw per-case output, logs, commands, resources, aggregate CSV, JSON summaries and source snapshots are retained.

## What this establishes

The pinned MFEM can assemble and solve this curved embedded-surface elliptic problem accurately at small cost. This does not extend our production OFDG implementation to surfaces. That would require choosing appropriate tangential derivative/sensor definitions and surface flux geometry; the square inverse-map construction for volume cells cannot be used unchanged when J is rectangular. No DG, OFDG, KXRCF, shock, transport or time-integration comparison is claimed here.
