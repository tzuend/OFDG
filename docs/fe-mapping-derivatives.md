# Derivatives from finite element geometry

This is an isolated experiment. Production sources, existing benchmark sources, and thesis files are unchanged. The supervisor preview retains its original 19 pages and appends the experiment on pages 20–22. The previous PDF and preservation fingerprints are in `tmp/before-fe-mapping-derivatives/`.

## What is being compared

All main-suite methods receive identical DG coefficients, the same represented mesh, and the same physical evaluation points and quadrature weights:

1. `recursive`: the unchanged production `ProjectedPhysicalDerivatives` matrices, applied in canonical chronological x/y/z order.
2. `direct`: direct differentiation of the mapped solution with analytic fixture geometry.
3. `fe_geometry`: direct differentiation using only geometry FE coefficients and its basis, with no analytic mapping formula available to this evaluator.

`IsoparametricTransformation::GetFE()` supplies the geometry basis, and `GetPointMat()` supplies its vector coordinate field. Those coefficients are reused, not projected again. Geometry and solution degree need not match. A cached change of polynomial basis makes the same reference evaluator available to each coordinate component and to the scalar solution. No physical derivative is reprojected along either direct path.

The new 2D benchmark includes the preserved earlier direct benchmark with its entry point renamed, reusing its definitions without editing it. The generic three-dimensional helper and the high-order helper live only in `benchmarks/`. Existing experiments and production interfaces are unchanged.

## Construction and arbitrary finite order

For reference-to-physical map F, let xi0 be the evaluation point and h a nominal reference-grid scale. Form the inverse Taylor series R(z) satisfying

    F(R(z)) = F(xi0) + h z.

Starting at R=xi0, update

    R <- R - J(xi0)^(-1) [F(R) - F(xi0) - h z].

With all arithmetic truncated after total degree m, each correction resolves another degree. We use m+1 corrections and record the remaining composition residual. Compose the original solution polynomial with R. Multiplying coefficient alpha by alpha! / h^|alpha| gives its physical derivative. Geometry derivatives beyond the polynomial geometry degree vanish, but derivatives of its inverse need not vanish.

The main helper retains degrees through three. `high_order::Jet<D,O>` provides a compile-time finite-order evaluator, instantiated at O=5 and O=10 for this experiment; the CLI accepts those two capacities. It has binomial(O+D,D) coefficients and caches valid multiplication index triples. The order-five specialization avoids carrying tenth-order terms for ordinary use. This is an algebraically extensible prototype, not a production arbitrary-order API or a guarantee of stable high-order accuracy.

## Main suite: 59 configurations

The 50 2D configurations exactly repeat the earlier direct-derivative study: triangles and quadrilaterals, P/Q2–P/Q4 solution spaces, degree-four geometry, derivatives through min(p,3), n=4,8 and selected n=16, both strong and mild coupled/separable maps, affine controls, Gauss-Lobatto/Gauss-Legendre/Bernstein solution bases, and extra-quadrature controls. Input projection, production matrix assembly, and physical evaluation quadrature are unchanged from that study.

- Coupled: x=s+b, y=t+0.7b, b=a s(1-s)t(1-t), a=0,0.4,1.2. Fields are the earlier sine, shortwave sine, exponential and constant.
- Separable: x=s+a s(1-s), y=t+0.7a t(1-t), a=0.3,0.7. The pullback s²+st+t² is exactly representable. The affine suite also includes a physical quadratic.

The nine 3D configurations use Q2 geometry with

    x=s+0.2 t(1-t), y=t+0.15 r(1-r), z=r.

The map has determinant one; local reference-cell determinants also contain h³. It has an explicit inverse and curved faces. Compare Q2/Q3 solutions on 2³ and 4³ cells, the corresponding four affine controls, and a finer Q3 curved repeat with eight extra evaluation quadrature orders. Fields are 1, the exactly represented pullback s²+t²+r²+str, and projected sin(2x+3y+z). All 19 nonconstant multi-indices through order three are evaluated, including third derivatives of Q2 solutions.

## High-order extension: 14 local-element configurations

These are representative single-element probes, not a global refinement convergence study. Use h=1/2,1/8,1/32, original-coordinate origins (0.2,0.3) or (0.2,0.3,0.15), and three local points: (0.23,0.31,0.37), (0.61,0.47,0.53), and (0,0.4,0.7), ignoring the third coordinate in 2D. The third point tests a face trace.

Use the 2D strong separable map (Q4 geometry) or the 3D shear map (Q2 geometry), plus affine controls. Q5 Gauss-Lobatto input coefficients are nodal values of a constant or the exactly representable original-coordinate polynomial

    2D: s^5+t^5+s²t²
    3D: s^5+t^5+r^5+s²t²+str.

The order-ten runs cover every pure and mixed multi-index through ten for all three methods. Two curved h=1/8 cases repeat at capacity five to check the specialization against the same low-order coefficients in the order-ten implementation. Initialization differs deliberately from the main suite: this is nodal interpolation of a representable field, not an L2 solve. Floating-point nodal values and basis conversion can still introduce errors even though the function belongs to the space mathematically.

`fe_mapping_reference.py` supplies independent 90-digit references: explicit quadratic inverse derivatives with the Leibniz rule in 2D, and expanded explicit-inverse physical polynomials in 3D. Selected fifth/tenth-order references are repeated at 120 digits. No C++ inverse-series code is used by these references. All 39,150 high-order rows, including failures of accuracy, are retained.

## Metrics, controls, and interpretation

Main CSV rows contain absolute/relative physical volume L2 error, sampled maximum error, one-sided face RMS/relative error, interior jump RMS/mean absolute value, input error, and discrepancy from the analytic direct derivative of the same computed solution. Face norms count both element traces; interior jump norms use physical interior-face measure. These are derivative diagnostics, not complete OFDG damping rates. Sampled maxima are not certified supremum norms. A vanishing analytic derivative gets an undefined relative error rather than an invented percentage.

The FE/direct agreement control uses Taylor derivatives scaled by h^order and divides by max(1, magnitude of the corresponding analytic direct scaled derivative). Geometry-series checks divide by h; first derivatives also compare independently against MFEM physical basis gradients. Geometry positions, Jacobians, inverse residuals and shared-face positions use a 1e-10 threshold; direct-method/first-derivative agreement uses 1e-8. The 2D geometry-series control compares ordinary Taylor coefficients; the 3D control includes factorials. All these tolerances were retained as planned.

The independent validation includes 18 oblique-map 2D probes on coupled and separable maps, using the earlier 65-digit Newton inversion and physical finite-difference reference; another 28 pure/mixed probes cover the high-order evaluator through third order in 2D/3D, with 90-digit explicit-inverse physical finite differences. Steps 1e-5 and 1e-6 check reference sensitivity. Both reference agreement and step sensitivity must remain below 1e-7 on a max(1, reference magnitude) scale. Affine polynomial and constant controls use a 1e-6 absolute L2 ceiling. These are validation gates, not tolerances selected after seeing failures.

High-order tables retain absolute errors and nonzero relative errors. Their plotted diagnostic is max |error|/max(1, |reference|) over the sampled points and multi-indices. It becomes absolute error for small reference derivatives and must not be described uniformly as a relative error. Near-zero references can produce enormous relative errors; the absolute/scaled measures are essential. High-order tests are diagnostic, with no silent accuracy-based exclusion or relaxation of the main acceptance gates.

The direct methods share inverse-series arithmetic, so their agreement alone is not independent validation. Mixed partials are represented by one Taylor coefficient; the experiment does not independently evaluate every differentiation ordering. Better differentiation of u_h does not imply a uniformly better approximation of derivatives of the underlying analytic u.

## Reproduce

From the project root with its existing MFEM installation:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make -j1 build/release/benchmarks/fe_mapping_derivatives \
  build/release/benchmarks/fe_mapping_derivatives_3d \
  build/release/benchmarks/fe_mapping_high_order \
  build/release/benchmarks/fe_mapping_controls
python3 scripts/run_fe_mapping_derivatives.py
python3 scripts/run_fe_mapping_high_order.py
build/release/benchmarks/fe_mapping_controls > measurements/fe-mapping-derivatives/control-probes.csv
MPLCONFIGDIR=/tmp/ofdg-fe-matplotlib python3 scripts/analyze_fe_mapping_derivatives.py
/Users/timozund/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  scripts/build_preview_pdf.py --input measurements/study/preview
```

The analysis Python needs NumPy/Matplotlib; references otherwise use only the standard library. The bundled Python provides ReportLab for the PDF. The main benchmark CLI retains `family n p amplitude basis map extra`; the 3D CLI is `n p curved extra`; the high-order CLI is `dimension n curved capacity` with capacity 5 or 10.

Measured cases run sequentially, one MPI rank and one numerical-library thread. Runners monitor process-tree sampled RSS, stop a case after 300 seconds or 4 GiB, and refuse new work after 2700 seconds combined compute. Resumption requires identical source and executable fingerprints; preserve existing measurement directories before intentionally rerunning changed code. Failed or interrupted outputs are not silently overwritten. Raw CSVs/logs, commands, source hashes, source snapshots, analysis tables, figures, and verification provenance remain under `measurements/fe-mapping-derivatives/`.

The CSV initialization/evaluation timers cover instrumented stages and are diagnostic, not a production performance benchmark or a complete partition of elapsed time. Process sampling may miss short-lived memory peaks. Independent references, smoke checks, compilation and PDF generation are outside the reported monitored case time.

## Measured findings

All 59 main configurations passed; 5,505 rows. All 14 high-order jobs completed, preserving 39,150 diagnostic rows (completion does not mean their high-order errors are small). The 46 independent probes have maximum scaled error 1.06e-10. All 2,644 preserved baseline rows reproduce the tested norms exactly. The largest element-scaled FE/direct discrepancy is 1.25e-09.

At n=8, the represented-quadratic Q3 worst third-derivative relative L2 errors are 0.0846 (recursive), 3.05e-10 (analytic direct), and 3.76e-10 (FE geometry). For P3 they are 0.0722, 7.89e-10, and 4.62e-9. Projected sine errors remain approximately the same between the two direct methods and are not uniformly improved over recursion.

In the curved 3D Q3 n=4 case, represented-input worst nonzero third-derivative relative errors are 1.02e-6 (recursive) and 7.21e-10 (both direct methods). Projected sine errors are about 0.471 and 0.472, respectively.

The higher-order results are intentionally not summarized as an accuracy success:

| Curved geometry | h | Order-5 scaled error | Order-10 scaled error |
|---|---:|---:|---:|
| 2D | 1/2 | 4.945e-11 | 2.268e-09 |
| 2D | 1/8 | 5.703e-09 | 7.693e-05 |
| 2D | 1/32 | 2.45e-06 | 28.61 |
| 3D | 1/2 | 3.431e-08 | 0.004261 |
| 3D | 1/8 | 3.633e-06 | 363.8 |
| 3D | 1/32 | 0.001514 | 3.245e+08 |

At the 3D h=1/32 probe, analytic direct evaluation has essentially the same order-5 and order-10 errors as FE geometry (0.00151 and 3.24e8); recursive errors are 0.000649 and 1.15e9. This is not solely an error introduced by reading the geometry as an FE field. The observations are consistent with amplification of initialization, polynomial-conversion and arithmetic errors in the existing representation; the experiment does not separate those contributions quantitatively. A stable basis-evaluation study is warranted before claiming uniformly accurate fifth derivatives. No high-order superiority or precision guarantee is claimed.

Order-five and order-ten implementations agree through order five within 7.88e-16 on a max(1,magnitude) scale. Combined monitored compute was 319.43 seconds, with 27.28 MiB peak sampled process-tree memory. Source/binary fingerprints and complete commands accompany the results.
