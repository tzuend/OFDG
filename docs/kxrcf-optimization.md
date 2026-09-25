# KXRCF optimization check

KXRCF computes no solution derivatives. The earlier combined filtering share included both the indicator and the remaining OFDG work. The indicator itself still evaluates cell norms and face traces at quadrature points and determines inflow. Curved overintegration increases those costs.

The implementation previously requested MFEM face transformations on every indicator call. That reconstructs element and face geometry from the mesh. The change retains owned transformations alongside the existing static face cache, using MFEM’s supported overload. Integration points, solution traces, and velocity/physics are still evaluated on each call. Public interfaces, quadrature, thresholds, and arithmetic of the indicator are unchanged. Reconstruct the indicator after mesh changes, as before. Shared MPI faces retain their existing path.

## Before/after measurements

One rank, one library thread, five measured repetitions plus warm-up; unchanged RK4 benchmark at T=0.005. The order of old/new binaries alternates between configurations. Reported times are medians per step. Total speed changes include run-to-run processor/cache variation; direct indicator timing is the more targeted measurement. Original results remain separate.

| Geometry / degree / grid / data | KX before → after (ms) | KX reduction | Total gated before → after (ms) | Gated / DG before → after |
|---|---:|---:|---:|---:|
| Affine / P2 / 32² / smooth | 2.66 → 1.70 | 36.0% | 15.01 → 15.02 | 1.393× → 1.353× |
| Curved / P2 / 32² / smooth | 7.53 → 5.87 | 22.0% | 23.53 → 21.99 | 1.782× → 1.639× |
| Affine / P2 / 64² / smooth | 10.66 → 7.49 | 29.7% | 60.00 → 57.44 | 1.390× → 1.308× |
| Curved / P2 / 64² / smooth | 29.50 → 23.54 | 20.2% | 93.83 → 89.17 | 1.748× → 1.624× |
| Curved / P1 / 64² / smooth | 28.42 → 21.07 | 25.9% | 83.30 → 73.52 | 2.046× → 1.791× |
| Curved / P2 / 64² / jump | 28.29 → 22.90 | 19.0% | 96.47 → 92.86 | 1.843× → 1.689× |

## Checks and limits

All 180 measured samples agree with their matching old/new counterpart in final L2 error and mass drift within 1e-12, and average activation counts within 1e-12. Existing KXRCF fingerprints, mixed-element tests, and curved geometry tests in serial and two ranks pass. A dedicated regression checks spatially varying, changing velocity and reuse after unrelated mesh transformation lookups; it passes in both release and AddressSanitizer/UndefinedBehaviorSanitizer builds.

Numerical benchmark elapsed: 3.47 minutes. Maximum sampled aggregate RSS: 261.3 MiB. Every process completed within 300 seconds and 4 GB. The CSV includes peak RSS for each old/new run: these are whole-process peaks, not an isolated cache allocation measurement.

The tradeoff is additional per-face geometry storage. No change to OFDG derivative preparation was made. Avoiding derivatives on inactive cells and their unneeded neighbors is a possible next optimization, but requires careful local/MPI dependency handling. Reducing curved quadrature or caching state-dependent inflow would affect numerical behaviour and is not part of this change.

Raw logs, exact commands, source snapshots, executable hashes, and timing CSV are in measurements/performance/kxrcf-optimization. The original benchmark brief describes the pre-optimization implementation; this document records the subsequent change.
