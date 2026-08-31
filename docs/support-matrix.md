# Verified support and preconditions

This table separates the mathematical portability of the construction from the
cases exercised by the present MFEM implementation. Unsupported inputs and
post-construction mesh changes are precondition violations; the core filters do
not promise runtime diagnostics for them.

| Capability | OFDG and OEDG | KXRCF | Verification |
|---|---|---|---|
| Solution coefficient basis | No modal, orthogonal, or Legendre storage basis is required | Same | Gauss--Legendre and Gauss--Lobatto OFDG fingerprints agree |
| Equations | Scalar fields and systems represented as repeated scalar DG spaces | Same | Advection, Burgers, and compressible Euler tests |
| Dimensions | Full-dimensional affine meshes in 1D--3D | Same | Serial and MPI tests in 1D--3D |
| Homogeneous element families | Operators are constructed from each MFEM finite element | Same | Segments, triangles, quadrilaterals, tetrahedra, hexahedra, and prisms |
| Mixed element geometries | Multiple affine signatures at one uniform order | Same | Triangle--quadrilateral and tetrahedron--hexahedron--prism meshes, P1--P3 |
| Polynomial order | One order throughout the space | One order throughout the space | Variable-order spaces are not verified |
| Geometry mapping | Affine mappings only | Affine mappings only | Nonuniform and skewed affine meshes |
| Mesh evolution | Geometry and topology remain fixed after construction | Same | Reconstruct the method object after any mesh change |
| Parallel execution | Independently sized local and remote face sides | Same | Mixed 2D and 3D shared faces across two ranks |
| Live visualization | Provided by the maintained drivers | Provided by the maintained drivers | Serial and MPI driver checks |

## Claim boundaries

The projection and decay operators are represented in the coefficient basis
supplied by MFEM. This is a basis-agnostic realization, not a claim that every
finite-element range or mapping is supported. The current code requires scalar,
`VALUE`-mapped DG elements repeated over the vector dimension.

The implementation has no geometry allow-list. Other compatible affine MFEM
element combinations follow the same per-signature path, but are unverified
until a representative serial and MPI case is added. Pyramid and curvilinear
elements are not included in the present affine claim.

KXRCF gating prevents OFDG from changing inactive elements. Their numerical
update is therefore the unfiltered DG update, apart from roundoff in shared
driver operations. Detection still has a cost, so runtime equivalence with pure
DG is not claimed. Reduced filtering work and reduced dissipation relative to
all-cell OFDG are the mechanisms to be measured in the final study.

The 2024 OEDG paper can be formulated with orthogonal modal coefficients, but
the damping operator can also be represented through basis-independent
projectors. The contribution here is a portable MFEM realization; this does not
show that competing implementations intrinsically require Legendre solution
coefficients.

See [geometry-design.md](geometry-design.md) for the remaining extension path.
