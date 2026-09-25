# Core readability refactor

This pass simplifies the existing files without changing the public headers,
method definitions, derivative ordering, or supported geometry families.

## Changes

- OFDG core: 1,796 -> 1,660 lines. KXRCF: 811 -> 784 lines.
- The private implementation classes now contain their actual implementation,
  instead of forwarding to an additional implementation object. Public OFDG and
  OEDG façades continue sharing the same numerical kernel.
- A single `AddFaceRates` routine implements the mathematical contribution of
  a face to one element, used by both local and MPI traversal.
- A single `ApplyDerivativeGraph` routine evaluates local and imported states.
  The parent choice and direction ordering are unchanged.
- Reference derivative construction is an ordinary function rather than a
  one-pointer builder class. It fills the Vandermonde and derivative-value
  matrices together and retains the coefficient-basis-independent solve.
- Projection construction reuses its already assembled high-order mass matrix.
  Only the remainder matrices actually consumed by damping are built. The
  pre-existing constant-space path retains its zero-order remainder.
- Removed unused internal constructors, the unused element vector-DOF array,
  and the derivative-pair wrapper. Corrected the imported-face comment and
  renamed `ComputeJumps` to `ComputeDampingRates` to describe its actual role.

The reference/curved operator-cache representation, ownership of geometry and
scratch, affine physical matrix construction, and sensor option API remain as
before. This is a focused simplification, not a wholesale rewrite. No core
implementation was moved into additional files. No performance gain is claimed.

## Before/after verification

`tests/test_refactor_outputs.cpp` records every coefficient in hexadecimal.
It covers intervals, triangles, quadrilaterals, tetrahedra, hexahedra, and prisms;
affine and curved maps (curved dimensions 2/3); Gauss-Lobatto and Gauss-Legendre
bases; two solution components; orders 1-3 in 1D/2D and order 2 in 3D. It emits
KXRCF pooled/component values, OFDG stabilization and decay, ungated and gated
outputs, a deterministic partial mask, an all-inactive mask, and affine OEDG.

The saved before and after runs use one rank and two ranks, with all numerical
library threads restricted to one. Across 69,824 recorded values:

- Serial outputs are byte-for-byte identical.
- Two-rank outputs differ by at most 1.1102230246251565e-16, both absolutely and
  scaled by max(1, |before|, |after|).
- 2,211 hexadecimal values differ across the MPI files. The shared rate helper
  uses multiplication by inverse amplitude, as the original local path did;
  the old MPI path used division. This is an arithmetic rounding difference,
  not bitwise identity of MPI results. The comparison tolerance is 1e-12.

The existing affine, mixed-element, KXRCF, OFDG fingerprint, OEDG, positivity,
face-physics, and RK checks pass. Curved geometry and Euler checks pass on both
one and two ranks. Public-header checks pass. An added regression checks that
the constant-space stabilization and decay leave piecewise constants unchanged.
The OFDG regression is also checked with address and undefined-behavior
sanitizers. Full logs and the final status manifest are retained with the data.

## Evidence and reproduction

`measurements/refactor-readability/` contains the pre-edit source snapshot and
baseline executable, full hexadecimal outputs, comparison summary, and logs.
`after/regressions.json` lists the existing numerical checks and their commands.
The checks run sequentially with a five-minute and 2-GiB process-tree ceiling,
below the allowed 4 GiB. The machine reports eight logical CPUs; two MPI ranks
are within its half-CPU cap.

```sh
make build/release/tests/test_refactor_outputs
python3 tests/run_bounded.py build/release/tests/test_refactor_outputs measurements/refactor-readability/after/serial
python3 tests/run_bounded.py mpirun -np 2 build/release/tests/test_refactor_outputs measurements/refactor-readability/after/mpi
python3 tests/compare_refactor_outputs.py measurements/refactor-readability/before measurements/refactor-readability/after
```

Use the recorded one-thread environment when reproducing these commands. This
verification samples the listed configurations; it is not a claim of bitwise
identity for arbitrary simulations or a new proof of numerical convergence.
