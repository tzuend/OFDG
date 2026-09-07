# Curvilinear OFDG: implementation and validation

## Support boundary

OFDG and KXRCF support fixed polynomial curvilinear meshes with uniform scalar
DG solution order, repeated over the component dimension. Validation covers
triangles, quadrilaterals, tetrahedra, hexahedra, and prisms at P1--P3, plus cubic
polynomial projections of MFEM's disc and ball NURBS fixtures at P2. Existing
affine mixed-element and OEDG regressions are retained. Curved OEDG is not a
validated paper-faithful extension.

**Exact native NURBS support remains blocked in the pinned MFEM face path.**
No mesh is silently converted. OFDG, KXRCF, and the positivity limiter reject
native NURBS geometry before entering that path. An application can explicitly
call `mesh.SetCurvature(3)` before constructing its solution space, accepting a
polynomial approximation of the rational geometry. This is what the projected
NURBS tests exercise.

Moving meshes, variable solution order, pyramids, and embedded surfaces are
outside this verification. Reconstruct the filters after any topology, nodal
geometry, or solution-space change.

## Exact native NURBS failure

Dependency: MFEM commit `a1ce49fb5742ad3c46e3e3732e5f9b04e4ca6e18`.

The first native disc test failed during OFDG face-cache construction:

1. Load MFEM `data/disc-nurbs.mesh` without converting its geometry.
2. Request `Mesh::GetFaceElementTransformations` for an interior face, including
   both adjacent element transformations.
3. Evaluate the physical face measure with `Face->Weight()`.
4. The release build terminates with `SIGSEGV` (signal 11). The observed stack
   enters `ElementTransformation::EvalWeight` and
   `IsoparametricTransformation::EvalJacobian`.

`tests/probe_native_nurbs.cpp` reproduces the failure using MFEM alone. It does
not include or link the OFDG library. To rebuild and inspect its process status:

```sh
make build/release/tests/probe_native_nurbs
python3 -c 'import subprocess; p = subprocess.run(["build/release/tests/probe_native_nurbs"]); print("returncode:", p.returncode)'
```

The observed Python return code is `-11`. The probe is deliberately excluded
from passing test targets.

The corresponding MFEM source explicitly guards this operation in
`mesh/mesh.cpp`, lines 1141--1144: under `MFEM_DEBUG`, a NURBS mesh with both
neighbor transformations requested aborts with `NURBS mesh not supported!`.
That check is absent in the release build. This establishes an unsupported
**two-sided native NURBS face operation**, not an absence of NURBS support
throughout MFEM and not a mathematical impossibility for OFDG.

Completing exact rational support requires a separately validated replacement
or repair of the native two-sided face geometry path, including independent
neighbor transformations and MPI geometry exchange. The DG evolution face
integrators also use that path. Changing OFDG's volume quadrature alone cannot
resolve this failure. No MFEM source files were modified by this implementation.

## Numerical construction

- Reference derivatives, derivative graphs, and sensor constants are shared per
  finite-element signature. Curved elements and MPI neighbors cache their own
  physical projectors and derivative matrices.
- Curved mass and cross-mass matrices use physical volume weights. Lower-order
  mapped reference spaces are nested, giving conservative shell attenuation.
- First derivatives satisfy `M D_j = G_j`, where
  `G_j[a,b] = integral(phi_a * physical_derivative_j(phi_b))`. Repeated matrix
  applications are successively projected derivatives, not exact higher mapped
  derivatives. Mixed projected derivatives need not commute; the existing
  derivative graph fixes their application order.
- Curved quadrature order is `2*p + 2*max(0, transformation.OrderW()) + 8`.
  Positive tensor/Duffy rules avoid cancellation in high-order simplex rules.
  All six barycentric permutations are included on triangular faces so
  reorientation across MPI ranks cannot change sampling.
- Face weights are normalized by integrated physical face measure; normal
  speeds use normals evaluated at the same quadrature points.
- The curved element scale is `(physical_volume/reference_volume)^(1/d)`.
  Face heights use the corresponding volume/face-measure ratio. KXRCF radius
  sampling includes vertices and geometry-dependent boundary quadrature.
- KXRCF treats transport within `128*machine_epsilon*speed_scale` of zero as
  tangential. This prevents normal roundoff from changing the inflow area while
  preserving velocity scaling; spectral radii are unchanged.

Two implementation failures were resolved during validation: high-order
simplex quadrature initially produced a constant-state defect of about
`2.0e-11`, and non-symmetric triangular sampling produced MPI indicator
mismatches of order `1e-5`--`1e-4`. Positive and orientation-symmetric quadrature
removed these failures. The projected ball then exposed tangential-flow
classification sensitivity, resolved by the relative transport tolerance.

## Reproduction and measured accuracy

```sh
make test-release
make test
make test-curved
```

`test-curved` checks conservation of every component, constants, unchanged
inactive cells, coefficient-basis independence, physical face moments against
higher-order quadrature, continuous physical linear states, and serial/two-rank
OFDG and KXRCF agreement. It includes positivity tests on curved and projected
NURBS geometry, a small periodic Euler evolution, and smooth advection.

For a P4 representation of `sin(x)` on a fixed smooth deformation of the unit
square, the physical L2 derivative errors are:

| Cells per direction | First derivative | Second derivative | Third derivative |
|---:|---:|---:|---:|
| 2 | 1.37463e-5 | 5.49751e-4 | 1.38358e-2 |
| 4 | 8.70317e-7 | 6.94008e-5 | 3.46324e-3 |
| 8 | 5.46137e-8 | 8.69929e-6 | 8.66081e-4 |

These approach orders four, three, and two. The first derivative additionally
satisfies an independently assembled weak projection residual, and increasing
quadrature by eight orders leaves the derivative matrix within the test tolerance.

Periodic curved P2 advection, final time 0.05 and time step 0.0005:

| Cells per direction | OFDG L2 error | OFDG--KXRCF L2 error |
|---:|---:|---:|
| 8 | 3.826409855e-3 | 2.044647431e-3 |
| 16 | 4.948469787e-4 | 3.277395723e-4 |
| 32 | 6.316891208e-5 | 4.811246351e-5 |

OFDG rates are 2.95 and 2.97; gated rates are 2.64 and 2.77. Both agree between
one and two ranks. The coarser 4-to-8 gated rate was 2.31, so the retained small
convergence test starts at 8 cells rather than asserting asymptotic behavior on
the coarsest grid. These are implementation checks, not thesis-final study data.

The Euler smoke test uses a smooth density wave on a periodic 4-by-4 curved
mesh for 20 steps with OFDG--KXRCF and positivity limiting. Conservation is
checked independently for all four conservative variables.

Equivalent-parameterization coverage is deliberately limited: an interior
bubble changes the square-cell map while preserving its physical boundary;
both representations preserve an exactly representable physical linear field.
General mapped polynomial spaces change under nonlinear reparameterization, so
arbitrary filtered states are not claimed to be parameterization-invariant.

Each `test-curved` check is monitored by `tests/run_bounded.py`, with a
five-minute deadline and a 2-GiB sampled resident-memory ceiling for its process
tree. Checks run sequentially and use no more than two MPI ranks. The final
validation batch must remain within the agreed 30-minute budget.

## Final local verification (2026-09-07)

- `make test-release`: passed, including affine OFDG fingerprints and OEDG regressions.
- `make test` under the bounded runner: passed; 126.67 seconds, sampled peak
  process-tree resident memory 1657.0 MiB (below 2 GiB). This includes sanitizer
  checks, existing two-rank mixed geometry/positivity tests, and five driver
  comparisons between one and two ranks.
- `make test-curved` under the bounded runner: passed; 56.17 seconds, sampled peak
  402.9 MiB including the build. Its serial geometry check, two-rank geometry
  check, serial/two-rank Euler smoke tests, and advection convergence checks all
  passed.
- The standalone MFEM-only native NURBS probe returned `-11` (`SIGSEGV`), as
  expected for the documented upstream limitation.

Local logs are retained under `build/curved-validation/` (ignored build output).
No thesis study measurements were replaced or marked final.
