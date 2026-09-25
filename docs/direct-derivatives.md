# Direct differentiation of the original mapped DG function

This isolated study preserves the production OFDG derivative recursion. It compares that recursion with direct physical derivatives of the same original mapped DG function. It does not reconstruct a new physical polynomial, project any derivative, evolve a PDE, or alter the damping implementation.

The findings replace the physical-polynomial reconstruction pages in `output/pdf/supervisor-preview.pdf` (pages 12-14). That earlier study remains in `measurements/physical-reconstruction/` and `docs/physical-reconstruction.md`; its PDF is archived as `tmp/before-direct-derivatives/supervisor-preview.pdf`. It is not included in the new brief's appendix because it answers a different question. The earlier eleven pages are unchanged. The thesis is unchanged by this study.

## Mathematical construction

Let F map an element's reference coordinates to physical coordinates and let u_h = u_hat composed with F^{-1}. At a physical evaluation point x_q = F(xi_q), expand the local inverse map as a bivariate Taylor polynomial in z:

    xi(z) = F^{-1}(x_q + h z),    h = 1/n.

Coefficients are ordinary Taylor coefficients (derivatives divided by multi-index factorials). Retain all terms of total degree at most three. Composition with the original reference polynomial gives

    u_hat(xi(z)) = sum_{a+b<=3} t_ab z_x^a z_y^b + higher terms,
    partial_x^a partial_y^b u_h(x_q) = a! b! t_ab / h^(a+b).

The truncated arithmetic is sufficient for exact derivatives through order three in exact arithmetic; it is not a finite-difference approximation. The implementation constructs the inverse series with

    xi_{m+1}(z) = xi_m(z) - J(xi_q)^{-1} [F(xi_m(z)) - x_q - h z],
    xi_0(z) = xi_q.

Because the residual iteration has zero linearization at the base point, each correction resolves another Taylor degree. Three corrections suffice through degree three in exact arithmetic; the implementation uses four and checks the composition residual. This is formal local inversion, not a nonlinear physical-point search with a finite stopping tolerance. The input reference polynomial is expressed in centred local monomials by interpolation. For the bounded degrees in this experiment, that is a change of basis, not an L2 projection or an approximation-space change. Basis conversion and input initialization incur floating-point error, which high derivatives amplify.

The direct evaluator uses the known analytic polynomial geometry of these fixtures, composed with each original affine Cartesian element map. MFEM stores degree-four nodal geometry representing that same map. Position agreement is checked after scaling by h, and first physical solution derivatives are independently checked against MFEM CalcPhysDShape. This is not a general implementation of higher geometry derivatives for arbitrary MFEM elements.

Each mixed derivative is one Taylor coefficient. Mathematical commutation is built into this construction; the study does not independently compute all differentiation orderings. First through third derivatives are independently verified as described below.

## Configurations and references

50 sequential configurations, 2644 derivative/method rows:

- Triangles and quadrilaterals, P/Q2-P/Q4, canonical physical derivatives through min(k,3).
- Strong-map refinement n=4,8 for all orders; n=16 for k=3. Quad and triangle element counts are n^2 and 2n^2, respectively.
- Coupled map: x=s+b, y=t+0.7b, b=a*s*(1-s)*t*(1-t). Strong a=1.2; mild a=0.4. Fields sin(2x+3y), sin(8x+6y), exp(x+y/2), and a constant.
- Separable map: x=s+a*s*(1-s), y=t+0.7a*t*(1-t). Strong a=0.7; mild a=0.3. Field s(x)^2+s(x)*t(y)+t(y)^2 and a constant. This field is exactly representable for every tested order. The quadrilateral edges remain straight, despite the non-affine map.
- Affine controls at n=4 for every order include a physical quadratic.
- Mild maps use n=8,k=3. Basis controls use n=4,k=3 and Gauss-Lobatto, Gauss-Legendre, Bernstein on both strong maps/families. Evaluation-quadrature controls add eight orders in those configurations.

Initialization and production projected matrices exactly reuse the prior study's quadrature and implementation. Input projection uses production order +16; derivative matrices use the unchanged production order. Evaluation uses production order +16 (+8 for sensitivity controls), positive physical quadrature, and physical face measures. Both methods receive identical DG coefficients.

For the separable map, inverse coordinates and analytic derivatives follow from s(x)=2x/(1+a+sqrt((1+a)^2-4*a*x)). Reference solution derivatives are never used to compute direct numerical derivatives: the direct method differentiates the original computed coefficient vector.

## Measures and interpretation

`results.csv` contains absolute volume L2 error, relative volume L2 error, sampled maximum error, one-sided face RMS/relative errors, interior jump RMS/mean absolute value, input approximation error, and the current recursion's L2 discrepancy from direct differentiation of u_h (`operator_error_l2`, plus normalization by the direct derivative norm). Errors relative to an analytically zero derivative are NaN rather than a misleading percentage. The direct method's operator discrepancy is zero by definition; it is not an independent correctness check.

Face RMS integrates both traces and divides by the sum of physical element-boundary measures. Jump RMS divides the interior-face jump integral by total interior measure. All exact fields are smooth, so their analytic derivative jumps vanish. These are not complete OFDG damping rates: no wave-speed, solution-amplitude, factorial, or element-length sensor factors are applied.

For n=8,k=3, strong maps, matching the xxx derivative:

| Field | Space | Recursive relative error | Direct relative error | Recursive absolute L2 | Direct absolute L2 |
|---|---|---:|---:|---:|---:|
| Exactly represented quadratic | Q3 | 8.4643e-2 | 1.8334e-12 | 52.516 | 1.1375e-09 |
| Exactly represented quadratic | P3 | 7.2244e-2 | 1.3763e-11 | 44.823 | 8.5393e-09 |
| Projected sine | Q3 | 1.9819e-1 | 1.9841e-1 | 1.1245 | 1.1258 |
| Projected sine | P3 | 4.2834e-1 | 4.3468e-1 | 2.4304 | 2.4663 |

These relative values are fractions, not percentages. Plot envelopes use the worst relative error across xxx,xxy,xyy,yyy, which need not select xxx for direct differentiation of the represented input.

The exact-input comparison removes input approximation error in exact arithmetic. Direct differentiation therefore improves dramatically, and its remaining errors are consistent with amplified floating-point initialization/basis-conversion errors. Across all exact-input cases its largest relative error is 6.16e-9. Refining at fixed order can increase that roundoff contribution; it is not evidence of divergent truncation error.

For smooth projected fields on strong n=8 maps, 138 matched comparisons have median direct/recursive volume error ratio 1.00089 (range 0.98277-1.37311), face ratio 1.00170, and jump ratio 1.00012 (range 0.565-1.99626). Only 13/138 volume errors improve. For the represented quadratic, all 46 comparisons improve in volume, face and jump errors. Ratios exclude recursive norms below 1e-9. The sample is a designed diagnostic set, not a statistical distribution of arbitrary meshes.

The difference is explained by

    direct derivative error = partial^alpha(u_h-u),
    recursive derivative error = [recursive(u_h)-partial^alpha u_h] + partial^alpha(u_h-u).

Direct differentiation eliminates the first bracket. The two terms can reinforce or cancel, so removing it does not guarantee a smaller error against the analytic solution. Intermediate projection can also suppress components of input error. The study supports further sensor/PDE testing, not replacement of the production method or claims of better stability.

## Independent validation

A separate probe harness tests a polynomial r^3*s + 0.4*r*s^2 + 0.7*r^2 - 0.2*s on an oblique affine reference-to-original-coordinate map followed by each curved map. The reference implementation uses Python Decimal at 65 digits, independent physical Newton inversion, and tensor-product central finite differences in physical x,y. All 18 first-to-third derivative probes agree with the jet implementation within 1.055e-10 after division by max(1, reference magnitude). Changing the reference step from 1e-5 to 1e-6 changes it by at most 1.045e-8 on that scale. Thus the verification does not differentiate using the same Taylor arithmetic.

Other controls:

- 934 matched recursive volume error norms reproduce the preserved accuracy study exactly in the saved CSV values.
- Inverse-series residual / h: <=1.07e-14; geometry position discrepancy / h: <=1.43e-13.
- Recovered reference polynomial vs MFEM evaluation: <=7.33e-14 absolute.
- First physical derivatives vs MFEM CalcPhysDShape: <=1.76e-11 scaled by max(1, magnitude).
- Basis changes: <=1.15e-9 in error norms scaled by max(1, analytic derivative L2 norm).
- Eight extra evaluation quadrature orders: <=2.01e-14 on that norm scale.
- Affine quadratic direct derivative errors: <=3.12e-10 absolute L2.
- Constant-input derivative errors: <=2.06e-8 absolute L2. The prototype does not enforce exact zero for constants; this would need attention before production integration.

## Resources and reproducibility

Eight logical processors were confirmed with sysctl. The study used one MPI rank and one numerical-library thread, sequentially, with no oversubscription. The runner monitors the process tree, kills checks exceeding five minutes or 4 GB sampled aggregate RSS, and stops launching after 45 minutes accumulated numerical time. The 50 monitored configurations used 52.54 seconds total, maximum single check 11.13 seconds, and maximum sampled RSS 21.0 MiB. Independent probe/Decimal checks add less than one second of execution; compilation and report generation are not numerical solver timing. These figures are not a production performance benchmark.

From codex_dev with the existing MFEM installation:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make build/release/benchmarks/direct_derivatives build/release/benchmarks/direct_derivatives_control
python3 scripts/run_direct_derivatives.py
build/release/benchmarks/direct_derivatives_control > measurements/direct-derivatives/control-probes.csv
MPLCONFIGDIR=/tmp/ofdg-matplotlib /opt/homebrew/bin/python3 scripts/analyze_direct_derivatives.py
/Users/timozund/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner resumes completed cases. To repeat from scratch, preserve/rename the existing measurements/direct-derivatives directory first. `run.json` stores exact commands, resource measurements and source/binary fingerprints. `source/` retains the experiment sources. `independent-controls.csv`, `comparisons.csv`, `representatives.csv` and `summary.json` contain the analysis and verification outputs.

Remaining support boundary: two explicit polynomial geometry fixtures in 2D, derivative order at most three. Higher geometry derivatives for arbitrary MFEM maps, NURBS, 3D, cached evaluation matrices, actual damping rates and PDE evolution are not implemented here. No production solver file was changed.
