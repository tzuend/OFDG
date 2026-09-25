# Physical-polynomial reconstruction versus recursive projected derivatives

This isolated experiment preserves the production OFDG implementation, the derivative-accuracy and derivative-permutation studies, and all previous supervisor-report pages. It changes only how an auxiliary derivative estimate is formed in a benchmark. No PDE or damping step is advanced. The pre-comparison supervisor PDF is retained in tmp/before-physical-reconstruction/.

## Definitions

Current method: production recursive D_j = M^{-1} G_j in the mapped DG space; the canonical chronological word xxy means Dy Dx Dx. Alternative: physically L2-project the exact same initialized DG function u_h into total-degree physical P_q(K), q=k,k+1,k+2, once, and evaluate analytic derivatives of that physical polynomial. Do not reproject derivatives into the mapped space. The polynomial basis consists of products of Legendre polynomials in centred/scaled physical Cartesian coordinates with total-degree indices. Local scales use a sampled element bounding box; they change the representation, not the spanned physical P_q space.

The oracle_k2 rows directly project the analytic field into physical P_(k+2). They are explicitly an unavailable diagnostic reference, excluded from feasible-method improvement counts. This reveals the limitation of reconstructing an already-approximated input. It is not an alternative that can be used by the solver without knowing the exact solution.

Mapped Q_k has (k+1)^2 coefficients in 2D; physical P_q has (q+1)(q+2)/2. Thus q=k is not dimension matched for quadrilaterals, and may discard content even on affine meshes. Taking l derivatives reduces physical degree to q-l, whereas current projected derivative states retain the full mapped order-k space. In particular q=k+2 leaves only degree k-1 after three differentiations. Commutation alone does not guarantee accuracy.

## Inputs and coverage

Geometry and analytic fields match the preceding derivative-accuracy study. Coupled map: x=s+b, y=t+.7b, b=a*s*(1-s)*t*(1-t), a=0,.4,1.2. Use sin(2x+3y), sin(8x+6y), exp(x+y/2), and a constant. Separable map: x=s+a*s*(1-s), y=t+.7a*t*(1-t), a=.3,.7. Use the exactly representable mapped quadratic s^2+st+t^2 and a constant. Both are represented by degree-4 geometry. The non-affine separable quad map retains straight edges; the coupled map has curved quad faces. Analytic inverse derivatives and geometry nonsingularity are documented in docs/derivative-accuracy.md.

P/Q2-P/Q4; all coordinate derivatives through min(k,3). Strong-map refinement uses n=4,8,16 for quads and n=2,4,8 for triangles. Mild maps use n=8 (all orders for the coupled map, k=3 for the separable map). Affine n=4 controls include an exact physical quadratic. Selected n=4, k=3 cases compare Gauss-Lobatto, Gauss-Legendre and Bernstein input bases and increase reconstruction/evaluation quadrature by eight. Exact configurations and commands are in run.json. This is a configured comparison, not a representative sample of all meshes.

Initialization is unchanged from the earlier study: physical projection with its production quadrature order +16. Production derivative matrices use their existing quadrature. Reconstruction uses positive quadrature order 8*(k+2)+8, reflecting degree-4 geometry; independent error evaluation uses eight more. Thus the original input is identical for all five methods. The QR construction solves the weighted sampling least-squares problem representing the physical L2 projection.

## Numerical stability and retained first attempt

Initially the physical reconstruction was obtained through the mass matrix (normal equations). An affine polynomial-reproduction control failed: degree-6 reconstruction on triangles produced a spurious third derivative of about 2.3e-5 in L2. This attempt and the precise source are retained in measurements/physical-reconstruction-normal-equations/.

The reported results use weighted QR with two-pass modified Gram-Schmidt and a triangular solve. This avoids squaring the evaluation matrix condition number. The analytic polynomial control and the nested-space, cell-integral, basis and quadrature checks are evaluated in scripts/analyze_physical_reconstruction.py. The first attempt's 205.75 seconds are charged to the final compute total; neither attempt is hidden or overwritten.

## Error measures

Absolute physical volume L2 error and relative error against analytic physical derivatives; sampled volume maximum error; one-sided face RMS and relative L2 errors; and interior derivative-jump RMS/mean absolute jump. Relative derivative errors are undefined for zero analytic derivatives. Interior exact jumps vanish for all these smooth fields. Numerical jumps therefore diagnose sensor contamination, but they are not the full OFDG rate: wave speeds, amplitudes, factorials and element lengths are not included. A small jump may also reflect loss of derivative information.

Input_error_l2 measures the original DG approximation. Reconstruction_l2 measures the physical polynomial minus the actual DG input, including for the oracle; it is zero by definition for recursive rows. Cell-integral defects check projection mean preservation, and normal-equation residuals independently check the QR solve. The benchmark is not a test of final filtered solution conservation.

comparisons.csv contains matched-direction ratios. representatives.csv takes the largest relative error over all directions of a specified derivative order (the absolute error column belongs to that selected direction). The heatmap ratios compare these order-wise envelopes; worst directions can change. aggregate.csv summarizes strong n=8, P/Q2-P/Q4 cases, excluding constants and baseline errors below 1e-9 or target norms below 1e-10. Oracle results never enter win counts. Basis and quadrature controls compare error norms, not pointwise equality.

## Reproduction

From codex_dev, using the existing MFEM build:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make build/release/benchmarks/physical_reconstruction
python3 scripts/run_physical_reconstruction.py
MPLCONFIGDIR=/tmp/ofdg-matplotlib /opt/homebrew/bin/python3 scripts/analyze_physical_reconstruction.py
/Users/timozund/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner resumes an existing checkpoint. For a clean rerun, preserve/rename measurements/physical-reconstruction first. All cases are sequential, one rank/thread, on a laptop with eight logical processors. Each check has a five-minute/4-GB process-tree ceiling; no new case is launched after 45 charged compute minutes. The runner currently includes the retained first-attempt charge. run.json supplies timings, memory samples, failure status and hashes. Summary, raw CSV/logs, comparison tables and figures are in measurements/physical-reconstruction/.

## Interpretation

Retain the current recursive projection as the production approach. Physical reconstruction is a coherent auxiliary sensor and removes permutation dependence, but these degrees do not show a universal accuracy gain. In some smooth cases q=k improves the errors, in others it discards useful content; enrichment often approaches the current answer. On the exactly represented strongly distorted k=3 input the current approximation is more accurate than the displayed physical reconstructions.

This does not rule out physical-space DG, which would evolve a different solution. A useful next isolated comparison would account for derivative degree loss, e.g. q=k+3 for third derivatives, or investigate a derivative-aware reconstruction objective. No final stability, convergence or shock-robustness claim follows without a PDE comparison.

## Completed run

62 QR checks passed, 8010 rows, 7.09 charged compute minutes including the retained first attempt; QR peak sampled RSS 42.6 MiB. The affine quadratic L2 derivative error is at most 3.82e-10; baseline norm agreement with the previous study is within 2.82e-14 on the comparison scale. The supervisor brief has 14 pages, with the previous 11 pages text-identical. The new pages were visually checked. See validation.json and summary.json for all control values.
