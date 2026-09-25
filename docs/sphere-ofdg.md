# Isolated first-order OFDG on a spherical surface

This experiment addresses surface DG with OFDG, rather than the earlier conforming surface-FEM screened-Poisson example. The latter remains saved in `measurements/sphere-surface/` and `docs/sphere-surface.md`; its PDF is archived under `tmp/before-surface-ofdg/`. The supervisor brief now has 17 pages, with the surface-DG/OFDG experiment on pages 15-17. The first fourteen pages are unchanged.

No production OFDG, KXRCF or volume-DG source is changed. The surface method is an explicit, isolated P1/Q1 adaptation in `benchmarks/sphere_ofdg.cpp`, not general surface support through the production API. Further changes to the OFDG formula or production integration must be explained to the user before implementation.

## PDE and exact solution

On the unit sphere Gamma:

    u_t + div_Gamma(v u) = 0,    v=(-y,x,0),    T=pi/2.

The velocity is tangential and surface-divergence-free. Its solution is rigid rotation:

    u(x,y,z,t) = u_0(x cos t+y sin t, -x sin t+y cos t, z).

Inputs are smooth u_0=1+0.5x, a discontinuous cap u_0=1 for x>0.5 and zero otherwise, and a constant control u_0=1. The cap is initialized by the same unfiltered physical L2 projection for DG and OFDG; this projection already overshoots. No initial clipping or positivity limiter is used.

## Surface mesh and transport operator

Start from MFEM's octahedral triangle or cubed quadrilateral construction, refined uniformly and with vertices on the unit sphere. At each quadrature point radially lift the polynomial base map F_h:

    X(xi)=F_h(xi)/|F_h(xi)|,
    J = (I-X X^T) J_base / |F_h|.

This is a custom exact-sphere surface map in the prototype. The mesh stored by MFEM supplies topology and reference coordinates; the pointwise radial map defines the actual geometry used for assembly, projection, measurement and rendering. It should not be mistaken for simply assembling on MFEM's planar geometry.

J is 3-by-2. Let G=J^T J. Then

    grad_Gamma u_h = J G^{-1} grad_reference(u_hat),
    dS = sqrt(det G) dxi.

Assemble the weak conservative advection form with an upwind numerical flux. An edge flux uses velocity dotted with the tangential outward co-normal, multiplied by physical edge length. Opposite-side factors are checked, then a single common flux is added with opposite signs to the two cells. This is a closed surface; there are no boundary conditions or boundary fluxes.

The spaces are discontinuous P1 on triangles and Q1 on quadrilaterals. Levels 1,2,3 give 32/128/512 triangles and 24/96/384 quads, respectively, with equal scalar unknown counts 96/384/1536 across the families. This equal count does not imply identical geometry, operator cost or accuracy.

## Exactly which OFDG adaptation is being tested

For each cell K and Cartesian component a=1,2,3, form the physical surface L2 projection

    g_{K,a}=Pi_K[(grad_Gamma u_h)_a].

These matrices have the same mass-inverse-times-weak-gradient construction as the current projected first derivative. They are assembled separately in the prototype because a surface has two reference coordinates and three tangential-gradient components in ambient Cartesian coordinates.

For each interior edge F define

    J_F,0 = (1/|F|) integral_F |[u_h]| ds,
    J_F,1 = (1/|F|) integral_F sum_{a=1}^3 |[g_a]| ds.

The componentwise absolute sum retains the Cartesian convention of the existing sensor. The individual projected vector components need not reconstruct a pointwise tangent vector; the definition is the jump of those projected components. No post-projection tangential correction, parallel transport or rotational invariance is claimed. For order one, no repeated or higher surface derivative is needed.

Let h_K=sqrt(|K|/|reference K|), beta_F be the sampled maximum absolute co-normal transport speed, and A the sampled global maximum of |u_h-global surface mean| over volume and edge points. Then

    lambda_K = sum_F beta_F [J_F,0/(2 h_K) + 3 J_F,1/2] / A,
    u_h^+ = mean_K + exp(-dt lambda_K) (u_h^- - mean_K).

The 1/2 and 3/2 prefactors and damping toward the physical cell mean are the existing order-one OFDG construction. The reference square includes all Q1 nonconstant modes in the damped tail. A relative near-constant guard disables damping when A <= 1e-12*max(1,sampled volume magnitude). There is no KXRCF gate and no positivity limiter.

Necessary surface changes are the tangential first derivative, physical surface measures, edge co-normal speeds, and the three-component representation. This is an adaptation, not a proof that a Euclidean-domain OFDG theorem extends unchanged to the sphere. Higher order would require a deliberate definition of higher surface derivatives and must not be enabled merely by changing the polynomial-order flag.

## Runs and numerical controls

Both methods use identical initial physical L2 projections and transport operators, classical RK4, and dt <= 0.12*min(h_K)/3, adjusted to end at pi/2. OFDG acts once after each full RK4 step.

32 final configurations:

- Two element families, three refinement levels, DG/OFDG, smooth/cap inputs: 24 runs.
- Constant controls with OFDG on both families at level 2: 2 runs.
- Halved CFL=0.06 on finest triangular smooth DG/OFDG: 2 runs.
- Smooth triangular level-2 quadrature 20 versus 26: 2 additional runs.
- Finest triangular cap quadrature 20 versus 30: 2 additional runs.

Volume quadrature is positive tensor/Duffy quadrature from the existing `CurvedRule` helper; edge quadrature is positive Gauss quadrature. Base assembly/input/sensor quadrature order is 20. Smooth error evaluation uses order 28, or assembly order +8 for sensitivity checks. Cap error evaluation uses order 60. The cap discontinuity is integrated numerically; no certified continuous norm is claimed. The maximum is sampled at error quadrature points and both one-sided edge samples, not certified between samples.

### Retained unsuccessful first measurement attempt

The first attempt used MFEM's generic high-order triangle rules. Negative quadrature weights produced a negative squared-error integral in the constant control and unreliable cap error norms. The PDE trajectories remained finite, but that did not make the error measurement valid. The first source and all its measurements are archived in `measurements/sphere-ofdg-initial-quadrature/`. These results are excluded from report tables.

The final rerun uses positive quadrature and explicitly rejects negative/nonfinite squared-error integrals. The analysis also asserts finite nonnegative errors. No change to the OFDG formula was needed. The first attempt's 9.44 seconds are charged to the final total.

## Findings

Finest level, 1536 scalar unknowns:

| Mesh/method | Smooth relative L2 | Cap relative L2 | Cap sampled min | Cap sampled max |
|---|---:|---:|---:|---:|
| Triangles / DG | 0.08919% | 22.30% | -0.1778 | 1.1728 |
| Triangles / OFDG | 0.1554% | 27.40% | -0.0440 | 1.0269 |
| Quads / DG | 0.08353% | 21.86% | -0.1150 | 1.1222 |
| Quads / OFDG | 0.2026% | 27.80% | -0.0222 | 1.0194 |

Relative errors divide by the full exact solution norm (including the offset in the smooth case). Smooth finest-pair L2 rates are 1.90/2.49 for triangular DG/OFDG and 1.93/2.63 for quadrilateral DG/OFDG. The short sequence is not evidence for superconvergence. OFDG reduces cap oscillations while increasing smearing and L2 error. It does not make the solution positive or strictly bound-preserving. The equatorial cut shows the tradeoff but does not contain all global extrema.

Validation maxima over the final suite:

- Mass drift from each run's discrete initial surface integral: 7.82e-14.
- Physical cell-mean change in one filter application: 7.78e-16.
- Final constant-field absolute L2 error: 6.67e-15.
- Semidiscrete constant-state RHS infinity norm: 3.95e-11.
- Opposite-side edge co-normal flux-density factors: difference <=3.34e-16.
- Integrated sphere area error against 4*pi: <=6.04e-14.
- Halving the finest triangular smooth timestep changes L2 error by 2.72e-9% for DG and 0.160% for OFDG.
- Increasing smooth quadrature changes L2 errors by 0.0111% DG / 0.0267% OFDG.
- Increasing cap assembly/input quadrature changes L2 errors by 0.0528% DG / 0.0423% OFDG, with the same order-60 error rule.

Mass conservation is measured relative to numerical initialization, not the exact cap area; the initial discontinuous projection has quadrature error. Changing the timestep also changes filter application cadence, so the comparison measures the combined time-discrete scheme.

## Resources, files and reproduction

One MPI rank and one numerical-library thread; sequential execution on the laptop with eight confirmed logical processors. Every check is monitored against five minutes and 4 GB process-tree RSS; no new check launches after 45 minutes charged execution. The 32 final runs plus the retained first attempt total 49.26 seconds monitored numerical time. Peak sampled memory is 389.5 MiB. Compilation, analysis and PDF generation are excluded. These prototype timings are not a production overhead benchmark.

`measurements/sphere-ofdg/` contains the commands, per-case logs, aggregate results, source fingerprints, summary and initial/final surface samples. The visualization files retain separate element traces. The sphere plot uses actual radial-map points and a common color scale. The initial sphere panel shows the exact cap; both numerical final panels start from the same projected coefficients. The equatorial cuts take northern one-sided traces and draw each element separately, with no smoothing across DG interfaces.

From codex_dev:

```sh
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
make build/release/benchmarks/sphere_ofdg
/opt/homebrew/bin/python3 scripts/run_sphere_ofdg.py
MPLCONFIGDIR=/tmp/ofdg-matplotlib /opt/homebrew/bin/python3 scripts/analyze_sphere_ofdg.py
/Users/timozund/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner resumes only matching source/executable fingerprints. Preserve the measurements directory before rerunning changed code. The historical first-attempt charge is explicitly retained in the run manifest and runner.

## Remaining support boundary

This is a working isolated P1/Q1 OFDG adaptation for rotational transport on a radially mapped sphere. It does not modify or enable surface meshes in the production OFDG API. Arbitrary embedded surfaces, higher-order covariant derivatives, surface KXRCF, production caching and MPI surface-neighbor behavior remain outside the prototype. Explain and discuss any additional OFDG change before implementing it.
