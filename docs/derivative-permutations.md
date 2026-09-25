# Isolated curved-derivative permutation study

This experiment changes no solver or thesis source. It uses the production
`ofdg::detail::ProjectedPhysicalDerivatives` assembler directly, then applies
its matrices in every distinct permutation of `xy`, `xxy`, and `xyy`.
Words are chronological: `xxy` means x, then x, then y, i.e. `Dy Dx Dx U`.
The sorted word is the sequence used by the production derivative graph.

## Reproduce

```sh
make build/release/benchmarks/derivative_permutations
python3 scripts/run_derivative_permutations.py
python3 scripts/analyze_derivative_permutations.py
python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner refuses to overwrite an existing `run.json`; archive the existing
`measurements/derivative-permutations` directory before a fresh run. Analysis
requires NumPy and Matplotlib; PDF generation requires ReportLab. The two steps
may use different Python environments. A single case is also runnable as:

```sh
build/release/benchmarks/derivative_permutations quad 8 3 1.2 0
```

Arguments: family (`quad` or `tri`), n, solution order, curvature amplitude,
additional derivative-assembly quadrature order. Output is CSV on stdout.

## Geometry and fields

The fixed domain is the unit square. On the undeformed coordinates (s,t), set
b = a s(1-s)t(1-t), x = s+b, y = t+0.7b. MFEM stores this polynomial map at
geometry order 4. Amplitudes are 0, 0.4, and 1.2. The boundary stays fixed and
its global Jacobian determinant is
1 + a[(1-2s)t(1-t) + 0.7s(1-s)(1-2t)], bounded below by 1-0.425a.
Local determinants in the CSV also contain the reference-to-grid factor 1/n^2.

Meshes use quadrilaterals (n^2 cells) or triangles (2n^2 cells), with n=2,4,8;
strong curvature additionally uses n=16. Solution spaces are Q3-Q5 and P3-P5,
respectively, with Gauss-Lobatto coefficient bases. These are not equal-DOF
comparisons. Only 2D static polynomial geometry is covered.

Fields are sin(2x+3y), exp(x+y/2), and x^2 y + x y^2. Their analytic mixed
partial derivatives are evaluated directly. The cubic has xxy = xyy = 2 and
is exactly representable on every affine space tested. Input coefficients are
an independently assembled physical L2 projection, not nodal interpolation.

## Accuracy and sensitivity metrics

For each multi-index, compute every distinct sequence and their arithmetic
mean as coefficient vectors (equivalently, as functions in the same space).

- `error_l2`: physical L2 norm of sequence minus the exact derivative of the
  original analytic field. This is the accuracy metric; the permutation mean
  is not the reference solution.
- `path_minus_mean_l2`: ordering sensitivity relative to the mean, separately
  from accuracy. A sizeable change need not yield a sizeable improvement.
- `exact_l2`: normalize errors by this value for a relative error. All chosen
  mixed derivatives have nonzero norms.
- `best_projection_error_l2`: independently project the exact derivative into
  the same space and measure its error. This is the approximation floor under
  the integration rule. Recursive differentiation need not attain it.
- `input_error_l2`: error of the initial projection. Reported higher derivative
  errors include its differentiated effect, as well as intermediate projection
  errors. We do not claim to isolate the exact higher derivative of the mapped
  DG function itself.

`comparisons.csv` gives canonical and mean errors, their ratio, and the change
relative to canonical error. `summary.json` and figures are derived entirely
from raw CSVs. The main comparison summary restricts to n=8, nonzero curvature,
and third derivatives, giving 72 field/space/geometry/multi-index combinations.

## Independent controls and limits

The test assembles a separate physical mass matrix and projects the analytic
input and exact derivative using quadrature 16 orders above the production
rule. It checks the weak first-derivative residual using MFEM physical basis
gradients integrated directly, rather than multiplying the production
assembly matrix. Constants must differentiate to zero and sampled Jacobian
determinants must be positive. Four strong-curvature n=8 cases repeat the
production derivative assembly at quadrature order +8 (both families, p=3,5).
The affine runs control permutation agreement and cubic exactness.

Checks run sequentially on one rank with numerical-library threads set to one.
Each process tree is sampled every 0.1 seconds with ceilings of 300 seconds and
4 GiB. No new case starts after 45 minutes of accumulated numerical time.
Statuses, failures, commands, measured memory and wall time, source/executable
SHA-256 hashes, and the MFEM configuration are retained in `run.json` and `raw/`.
This is an operator-accuracy study, not a performance benchmark.

The supervisor PDF appends this study when its `summary.json` exists. Its
original six pages still describe the earlier PDE dataset; the operator study
has its own resource accounting. No evidence here establishes an improvement
in OFDG face-jump sensors or PDE solutions. Those are separate follow-up tests.

## Measured findings (64 completed checks)

All checks passed in 277.78 seconds (4.63 minutes); sampled peak process-tree
RSS was 19.9 MiB. The affine controls use the same physical projection assembler
on affine maps, rather than the production affine fast path. This isolates the
geometric effect; it is not a replacement for the affine regression suite.

Across the 72 curved n=8 third-derivative comparisons, the mean improves 41
and worsens 31. Mean/canonical error ratios are 0.712 to 2.220, median 0.999456.
The largest worsening is Q3, a=0.4, cubic field, xyy: 0.001092 to 0.002424.
Since the exact derivative is 2, these are relative errors 0.0546% and 0.1212%.
The best reduction is Q3, a=0.4, exponential field, xxy: 0.006115 to 0.004354.
These outcomes rule out treating the permutation mean as an accuracy reference.

For sine, a=1.2, n=8, xxy:

| Space | Canonical L2 error | Mean L2 error | Error change |
|---|---:|---:|---:|
| Q3 | 0.151500 | 0.144589 | -4.56% |
| Q4 | 0.008933 | 0.008491 | -4.95% |
| Q5 | 0.00043967 | 0.00042827 | -2.59% |
| P3 | 2.649333 | 2.648790 | -0.02% |
| P4 | 0.428397 | 0.427788 | -0.14% |
| P5 | 0.049324 | 0.049219 | -0.21% |

For n=8 to 16 on this field and geometry, measured log2 error-reduction rates
for canonical/mean are Q3 1.263/1.220, Q4 2.224/2.174, Q5 3.120/3.093;
P3 0.997/0.997, P4 2.003/2.001, P5 3.006/3.004. There is no observed order
improvement from averaging on this sequence. These are empirical rates on a
small refinement range, not asymptotic convergence proofs.

Affine canonical-minus-mean L2 differences are at most 5.57e-9; affine cubic
errors at most 3.25e-9 after repeated differentiation. Independent first weak
residuals are at most 3.37e-15 and constant derivative coefficient defects
8.54e-12. The quadrature repeats change error norms by at most 2.41e-5 relative.
All 2,112 reported numerical rows are finite and respect the independently
computed projection-error floor (allowing 1e-10 for roundoff).

Recommendation: retain the current canonical order pending a targeted
comparison of face jumps, damping rates, and a complete smooth-advection run.
Permutation averaging is a different operator with no guaranteed accuracy gain;
this volume-norm experiment does not establish whether changing it helps OFDG.
