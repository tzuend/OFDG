# Supported-space contract

This page distinguishes the portability of the mathematical construction from
the behavior verified in the present MFEM implementation. “Not yet supported”
means that construction fails immediately with a diagnostic; it does not mean
that the method is mathematically restricted to that case.

| Capability | Current implementation | Verification |
|---|---|---|
| Solution coefficient basis | No modal, orthogonal, or Legendre storage basis is required | Gauss–Legendre and Gauss–Lobatto fingerprints agree |
| Equations | Scalar fields and systems represented as repeated scalar DG spaces | Advection, Burgers, and compressible Euler tests |
| Dimensions | Full-dimensional meshes in 1D, 2D, or 3D are accepted | Routine regression coverage is strongest in 1D and 2D |
| Element families | One geometry family per finite-element space | Segments, triangles, rectangles, and skewed quadrilaterals are tested independently |
| Mixed element geometries | Not yet supported | A mixed triangle–quadrilateral mesh is rejected explicitly |
| Polynomial order | One order and local DOF layout throughout the space | Variable-order spaces are rejected explicitly |
| Geometry mapping | Affine mappings only | Nonuniform, translated, rotated, and skewed affine elements pass; a NURBS mesh is rejected |
| Mesh evolution | Geometry and topology remain fixed for a filter's lifetime | Sequence changes are detected before cached data are accessed |
| Parallel execution | Distributed affine meshes with shared-face exchange | One- and two-rank regression tests |
| Live visualization | Optional GLVis output in all maintained drivers | Serial and MPI connection behavior is covered by driver checks |

## Claim boundaries

The projection and exact-decay operators are represented in the coefficient
basis supplied by MFEM. This is a basis-agnostic realization, not a claim that
every finite-element range or mapping is supported. The current code requires
scalar, `VALUE`-mapped DG elements repeated over the vector dimension.

KXRCF gating prevents OFDG from changing inactive elements. Their numerical
update is therefore the unfiltered DG update, apart from roundoff in shared
driver operations. The KXRCF detector still consumes time, so “restores DG
runtime” is not a valid claim without measurements. Its supported claim is
reduced dissipation and reduced filtering work relative to all-cell OFDG.

The 2024 OEDG paper can be formulated with orthogonal modal coefficients, but
the damping operator can also be represented through basis-independent
projectors. The contribution here is a portable MFEM realization; it should not
be described as proof that other OEDG implementations intrinsically require a
Legendre solution basis.

See [geometry-design.md](geometry-design.md) for the ordered extension path.
