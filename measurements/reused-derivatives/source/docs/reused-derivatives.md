# Curved derivatives from the existing affine reference machinery

This is an isolated experiment. Production code and all older benchmarks remain unchanged. Results live in `measurements/reused-derivatives/`. The supervisor preview appends the measured conclusions to the earlier FE-geometry study.

## What is actually reused

`benchmarks/reused_derivatives.hpp` links the unchanged production `ofdg::detail::AssembleReferenceDerivatives` from `libofdg.a`. The assembler uses only `CalcShape`, `CalcDShape` and an unisolvent reference sampling matrix: for each reference direction, `D_a = V^{-1} G_a`. There is no new monomial conversion or analytic mapping formula in the new evaluator. Geometry coefficients come from MFEM's `GetFE()` and `GetPointMat()`.

The small experiment-local coefficient graph follows the production graph's canonical parent rule and multiplication operation. It extends the graph to the requested order and inserts structural polynomial zeros. It does **not** instantiate the production OFDG operator, which would also assemble unrelated damping/projection data. Geometry and solution have separate operators and coefficient arrays. The polynomial space must be closed under reference differentiation. This benchmark supports scalar P spaces on simplices and Q spaces on tensor-product cells, with vector-valued coordinate coefficients; it does not claim arbitrary non-polynomial, surface or NURBS support.

## Pointwise recurrence

Let `h` be a positive, constant element scale, `J = dF/dxi`, `K = J/h`, and `B = K^{-1} = h J^{-1}`. Coordinate FE derivatives provide every required `K_beta`. Differentiate `K B = I` to obtain

```
B_0 = inverse(K_0)
B_beta = -B_0 sum_{0 < gamma <= beta} binomial(beta,gamma) K_gamma B_(beta-gamma).
```

The matrix product order matters. Let `G_(alpha,beta)` denote the reference derivative `beta` of the scaled physical derivative `h^|alpha| D_x^alpha u_h`, pulled back to the element. Initialize `G_(0,beta)` from the same reference derivative graph. Then

```
G_(alpha+e_i,beta) = sum_a sum_{gamma <= beta}
    binomial(beta,gamma) (B_gamma)_(a,i) G_(alpha,beta-gamma+e_a).
D_x^alpha u_h = h^(-|alpha|) G_(alpha,0).
```

Only states with `|alpha|+|beta| <= m` are used. No intermediate physical derivative is projected into an FE space. No inverse-map Taylor series is used by this new method. An order-five request needs geometry derivatives through order five and inverse-Jacobian derivatives through order four. Polynomial reference derivatives can vanish while arbitrarily high physical derivatives remain nonzero on curved cells. The `h^{-m}` amplification remains unavoidable.

The graph is cached per basis/order and derivative coefficients per element; pointwise inverse-Jacobian and chain-rule states are computed at evaluation points. The current dense state storage favors a transparent experiment over production optimization.

## Centering is a separate comparison

`reused` applies the matrices to the original coefficients. `reused_centered` subtracts a constant before differentiating: `U - c q`, where `V q = 1` and `c` is the field value at the reference center. Geometry coordinates are centered independently. Values at order zero still use the original coefficients. This avoids assuming that constant coefficients are all ones. In exact arithmetic centering cannot change a derivative; in floating point it changes the subtraction and multiplication errors and is not guaranteed to improve every result.

## Cases, references and measures

The main suite repeats the previous 50 two-dimensional and nine three-dimensional configurations with all five methods and identical original solution coefficients, mesh geometry, quadrature points and weights. It retains volume L2, relative volume L2, sampled maximum, one-sided face RMS, relative face L2, interior jump RMS/mean absolute, and input approximation error. Operator-error columns use analytic direct differentiation of the represented polynomial as a **comparison**, not an independent oracle. Analytic direct and the earlier FE method share inverse-series arithmetic.

Fourteen local Q5 stress configurations test all pure and mixed derivatives through order ten, affine/curved maps, 2D/3D, and `h=1/2,1/8,1/32`, including order-five capacity controls. These are three-point samples, not integrated volume norms. Independent 90-digit explicit-map formulas supply the ideal-field references.

New error measures are:

- Absolute error, including all zero-reference derivatives.
- Relative error where the reference is nonzero; near-zero references can make it misleading.
- Order-group discrete relative L2: `sqrt(sum(error^2)/sum(reference^2))` across the sampled points and multi-indices. This is not the physical-volume L2 and is coordinate dependent.
- `max |error|/max(1,|reference|)` for continuity with the earlier stress plots. This mixed absolute/relative scale assumes the dimensionless fixture units; it is not a universal relative error.
- Element-scaled absolute error `h^m |error|` (and `h^m` times the volume L2/operator L2 error). This exposes the reference-scale error that physical differentiation amplifies; it must be read beside physical error, not used to conceal it.
- Independent **evaluation error** against a high-precision derivative of the same stored FE geometry and solution. This reference reconstructs tensor Lagrange polynomials from exact binary64 nodes and coefficients, performs 90-digit Newton inversion, and applies physical finite-difference stencils at steps `h*1e-4` and `h*5e-5`. It shares neither new recurrence nor old inverse-series code. Twenty selected pure/mixed probes cover orders one, three and five in 2D/3D at `h=1/8,1/32`.
- **Representation error**: stored-function reference minus the ideal analytic field derivative at that same physical point. For these exactly representable fixtures this measures initialization/geometry coefficient perturbations, not truncation from insufficient degree. Signed evaluation and representation errors add to signed total error; their absolute magnitudes need not add because they can cancel.

Additional basis controls reproduce fifth-order polynomial derivatives and constants in Gauss-Lobatto, Gauss-Legendre and Bernstein bases, and in independently scaled versions of those bases whose constant vector is not all ones.

## Reproduction and limits

```
make -j1 build/release/benchmarks/reused_mapping_2d build/release/benchmarks/reused_mapping_3d build/release/benchmarks/reused_high_order build/release/benchmarks/reused_basis_controls
python3 scripts/run_reused_derivatives.py
OMP_NUM_THREADS=1 build/release/benchmarks/reused_basis_controls > measurements/reused-derivatives/basis-controls.csv
python3 scripts/reused_stored_reference.py
python3 scripts/analyze_reused_derivatives.py
python3 scripts/build_preview_pdf.py --input measurements/study/preview
```

The runner records exact commands, source/binary fingerprints, stdout/stderr, return codes, sampled process-tree memory and elapsed time. It runs sequentially, one rank/thread, with five minutes and 4 GiB per case and a 45-minute total compute ceiling. Resume requires matching fingerprints; interrupted outputs must be preserved before retry. Initialization/evaluation columns time the combined five-method diagnostic program, not individual implementations. No production performance claim follows from them. Reference analysis and compilation are separate from monitored benchmark compute.

Acceptance controls use `1e-10` for geometry/inverse identities, `1e-8` for MFEM first derivatives, and `1e-7` for independent finite-difference step sensitivity and fifth-order reference-basis reproduction. High-order physical errors are measured outcomes, not presumed to satisfy these validation thresholds. Old result rows are checked for exact reproduction. No solver evolution or OFDG damping changes are included.

## Measured results

All 73 monitored configurations completed: 9,175 main rows and 65,250 high-order rows, 409.2 seconds total benchmark compute, 34.0 MiB peak sampled memory. All 5,505 older main rows and 39,150 older high-order values reproduced exactly in the checked error columns/derivative values. The 24 basis controls passed. Worst new MFEM Jacobian discrepancy was 3.65e-12, first-derivative discrepancy 1.99e-11, and differentiated inverse identity residual 4.87e-15.

For curved Q5 probes at h=1/32, the maximum `|error|/max(1,|ideal derivative|)` at order five was:

| Dimension | Recursive | Analytic direct | Earlier FE series | Reused | Reused + centered |
|---|---:|---:|---:|---:|---:|
| 2D | 1.63e-3 | 2.11e-6 | 2.45e-6 | 1.84e-7 | 3.64e-7 |
| 3D | 6.49e-4 | 1.51e-3 | 1.51e-3 | 1.89e-5 | 1.66e-5 |

This improves on the earlier FE-series method by about 13 times in 2D and 80 times in 3D without centering. Centering is not uniformly better. At order ten the centered errors reach 27.0 in 2D and 1.40e7 in 3D on this same scale; those fine-element results are not accurate. Capacity-five and capacity-ten outputs agree through fifth order to 7.88e-16 on a max(1,magnitude) scale.

The independent stored-function reference has at most 1.80e-17 step sensitivity across 20 probes. At h=1/32, fifth-order maximum absolute evaluation errors (selected directions) are 1.16e-6 versus 8.88e-8 for earlier FE series versus centered reuse in 2D, and 8.45e-5 versus 7.64e-7 in 3D. The corresponding maximum absolute representation errors are 1.86e-7 and 1.37e-6. These maxima can occur in different directions. Raw per-direction signed values are saved so cancellation and the exact decomposition can be inspected.

The initial 1e-7 accuracy aspiration is **not uniformly achieved** by fifth-order physical evaluation: selected centered evaluation errors reach 7.64e-7 after max(1,reference) normalization. Stable independent steps rule out stencil uncertainty at that magnitude. Reference-basis errors, coefficient perturbations, inverse-Jacobian products and h^-m amplification remain. No acceptance threshold was relaxed to conceal these outcomes. Broad first/third-order controls pass; fifth/tenth-order stress results characterize the limits instead of constituting a blanket success claim.

The projected sine remains dominated by input approximation across the direct approaches. Improved evaluation of a represented FE polynomial does not reconstruct missing derivatives of the underlying analytic field. Full absolute/relative volume and face errors, jumps, size-scaled errors, discrete high-order relative norms, and independent error decomposition are retained in CSV files.
