# Curvilinear derivative accuracy (2026-09-10)

135 completed checks, 3072 rows, one rank/thread, 14.69 minutes charged compute and 21.33 MiB sampled peak process-tree RSS. Eight logical processors were confirmed; each check stayed below five minutes and 4 GB. One interrupted check was retained separately and charged a conservative 60 seconds when reducing the original schedule. No production solver code was changed.

## Reproduction

From codex_dev:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make build/release/benchmarks/derivative_accuracy
python3 scripts/run_derivative_accuracy.py
python3 scripts/check_derivative_accuracy_reference.py
MPLCONFIGDIR=/tmp/ofdg-matplotlib /opt/homebrew/bin/python3 scripts/analyze_derivative_accuracy.py
/Users/timozund/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner resumes existing checkpoints. Preserve/rename measurements/derivative-accuracy before a fresh run. run.json contains exact commands, source and executable hashes, quadrature increments, timings, failures and interruption accounting. The analysis requires NumPy/Matplotlib; PDF building requires ReportLab. Benchmark arguments are: family n p amplitude basis map assembly_increment evaluation_increment.

## Fields and maps

A: x=s+b, y=t+0.7b, b=a*s*(1-s)*t*(1-t), a=0,.4,1.2. Project sin(2x+3y), sin(8x+6y), exp(x+y/2), and a constant using physical L2 initialization. This measures input approximation plus derivative recursion error. The global determinant is at least 1-.425a, hence .49 at strongest curvature.

B: x=s+a*s*(1-s), y=t+.7a*t*(1-t), a=.3,.7. Use u=s(x)^2+s(x)*t(y)+t(y)^2. Its reference pullback is quadratic and is exactly represented by P/Q2 and higher. The inverse s=2x/(1+a+sqrt((1+a)^2-4ax)) has derivatives 1/j, 2a/j^3, 12a^2/j^5, where j=1+a-2as. Combining these through the chain rule gives exact physical derivatives. The global determinant is at least (1-a)(1-.7a), hence .153. Quad edges remain straight with non-affine parameterization; triangle diagonals may curve. Study A supplies curved quad faces. Degree-4 geometry represents both maps exactly.

All canonical derivatives through min(p,3) are tested: x,y; xx,xy,yy; xxx,xxy,xyy,yyy. Words are chronological, e.g. xxy=Dy Dx Dx. P/Q1-P/Q4 in A; P/Q2-P/Q4 in B. Gauss-Lobatto and Gauss-Legendre nodal bases and positive Bernstein basis are compared. The study uses production ProjectedPhysicalDerivatives, with no filtering/PDE evolution. Reference inverse derivative formulas passed 135 independent 70-digit centered finite-difference checks, maximum scaled discrepancy 1.26e-16. See analytic_controls.json.

## Error measures and controls

Absolute L2 error E=||computed derivative - analytic derivative|| over physical volume; relative R=E/||analytic derivative||. Relative errors are undefined (NaN) for zero derivatives. L1, sampled maximum absolute/relative error, input projection error, exact derivative norm, direct best-projection error, and excess over that projection are also exported. Squared total error equals squared best-projection error plus squared excess up to quadrature and roundoff. For study A the excess includes differentiated input error; study B removes it (input error below 1.27e-14).

Face RMS sums one-sided errors over all element boundaries, counting interior faces twice, and divides by total physical boundary measure before taking the square root. Face relative L2 uses the matching analytic trace norm. These are trace diagnostics, not actual OFDG jumps or damping rates. Sampled maxima are not certified suprema. scaled_l2 is E*(1/n)^order/||u||; 1/n is a nominal reference scale, not OFDG's physical cellwise length.

Independent initialization/evaluation uses production quadrature order +16. Selected n=8 strong P/Q3 cases separately raise assembly or evaluation order by eight. Weak first-derivative checks apply the chain rule independently from reference basis gradients. Basis comparisons concern physical error norms, not coefficient or pointwise equality.

## Assessment

Strong coupled curvature, n=8, sine: Q3 worst relative errors by derivative order are .036%, 1.06%, 19.82%; P3 .140%, 3.64%, 42.83%. Q4 gives .00149%, .0641%, 1.81%; P4 .0108%, .450%, 8.46%. At n=16, Q3/P3 third-derivative errors decrease to 9.92%/21.34%. Highest derivatives converge slowly at low order.

Exactly represented quadratic, strong separable map, n=8: Q3 gives .0173%, .831%, 8.46%; P3 .0136%, .655%, 7.22%. Q4 gives .00196%, .130%, 1.98%; P4 .00153%, .0997%, 1.61%. The raw tables contain absolute errors, all derivatives, milder curvature and shorter-wave cases.

Changing basis changes nonconstant-field error norms by at most 8.39e-10 scaled by max(analytic derivative norm,input norm). Constant-field high derivatives reveal amplified roundoff, reaching 1.84e-6 at sampled points after repeated differentiation of initial errors around 1e-15 to 1e-14. This benchmark does not exercise OFDG's amplitude cutoff or cell-mean preservation.

This is a sensible construction to retain as an OFDG roughness sensor, but its highest derivatives are not uniformly accurate. Refinement/order help; changing basis does not cure approximation error. Face errors can be larger than volume errors. These norms do not prove acceptable damping, smooth-solution convergence or shock robustness. The next useful validation is an independent face-jump/damping comparison followed by smooth-advection convergence. No universal adequacy claim follows from this study alone.

The supervisor brief includes the study; the thesis separately has an expanded step-by-step derivation. Backups are in tmp/before-derivative-accuracy-study/. Raw data, plots, controls, convergence tables and summaries are in measurements/derivative-accuracy/.
