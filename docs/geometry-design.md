# Mixed-element and curvilinear design decision

## Decision

KXRCF now supports the first two geometry milestones: affine 2D meshes mixing
triangles and quadrilaterals, including MPI faces with unequal side-specific
DOF counts. The remaining production order is:

1. add the same mixed affine 2D support to OFDG and OEDG;
2. add further affine three-dimensional element families;
3. add static curvilinear mappings.

The OFDG mixed-element refactor must remain separate from curved-map changes. Every stage retains the current affine fingerprints and public
`ofdg::` interface.

For curvilinear elements, repeated physical derivatives will be defined through
recursive element-local (L^2) projection of first physical derivatives. This
is the selected portable extension:

\[
  M_K D_{K,d}=G_{K,d},\qquad
  (M_K)_{ij}=\int_K\phi_i\phi_j\,dx,\qquad
  (G_{K,d})_{ij}=\int_K\phi_i\,\partial_{x_d}\phi_j\,dx.
\]

Applying the matrices \(D_{K,d}\) along the derivative DAG produces the
coefficient representation used for face traces. It requires basis values,
first physical basis derivatives, quadrature, and local linear algebra. On an
affine polynomial element the derivative remains in the discrete space, so the
construction reduces to the current exact coefficient derivative up to
roundoff. On a curved mapped element it represents successively projected
physical derivatives, not exact higher derivatives of the mapped function.
That distinction must be stated in the report and API documentation.

Exact higher physical derivatives were rejected as the default because they
require mapping and basis derivatives beyond first order, are not uniformly
available for all requested polynomial degrees, and would weaken the intended
framework portability. They may later be implemented only as a validation
oracle for element families where MFEM exposes the required Hessians.

## Remaining OFDG/OEDG mixed-element architecture

Replace the single operator built from `GetFE(0)` by an immutable repository
keyed by an element signature containing geometry, polynomial order, scalar DOF
count, map type, and coefficient-basis identifier. Each local element stores a
handle to its operator set rather than assuming a common matrix size.

Face cache entries own independently sized evaluation matrices for element 1
and element 2. Runtime scratch storage grows to the largest local signature but
views retain the actual side-specific dimensions. No padding participates in a
matrix operation.

For shared MPI faces, initialize face-neighbour metadata before building the
cache. Obtain the remote finite element with
`ParFiniteElementSpace::GetFaceNbrFE`, obtain its vector DOFs with
`GetFaceNbrElementVDofs`, and construct the remote derivative state from that
finite element's signature and transformation. Ownership remains unchanged:
each rank accumulates damping only for its local element.

KXRCF is already tested on an in-memory P2 mesh with two triangles and two
quadrilaterals. Its explicit two-rank partition places triangle--quadrilateral
faces across ranks and compares active count, indicator sum, and maximum with
the serial result. The OFDG/OEDG refactor will reuse this test geometry.

## Curvilinear geometry data

Curved support makes the following quantities element- or quadrature-specific:

- mass, cross-mass, shell-projector, and projected-derivative matrices;
- inverse Jacobians and physical volume weights at every volume point;
- physical face weights and unit normals at every face point;
- normalized face weights formed with the integrated physical face measure;
- physical element volumes and reference volumes;
- element and face characteristic lengths.

Use

\[
  h_K=\left(\frac{|K|}{|\widehat K|}\right)^{1/d}
\]

for the OFDG element scale. It equals the current determinant-based value on an
affine map. For a face \(F\), use the volume-to-face ratio corrected by the
reference measures,

\[
  h_{K,F}=\frac{|K|}{|F|}\frac{|\widehat F|}{|\widehat K|}.
\]

This gives the normal height for affine simplices and tensor-product elements
and remains meaningful for curved faces. KXRCF retains its radius meaning by
approximating

\[
  \max_{\xi\in\partial\widehat K}
  \|T_K(\xi)-T_K(\xi_c)\|
\]

with a boundary rule chosen from the geometry order; affine convex elements
therefore reproduce the current vertex result.

Caches remain immutable. Mesh topology, nodal geometry, polynomial order, or
finite-element-space sequence changes invalidate every cache and require
reconstruction. Arbitrary Lagrangian–Eulerian mesh motion is explicitly out of
scope for the first curved implementation.

## Required validation before claiming support

- Preserve constants and every component's cell average to roundoff.
- Reproduce all existing affine fingerprints on homogeneous meshes.
- Match direct physical first derivatives for manufactured functions.
- Demonstrate the expected refinement convergence of recursively projected
  second- and third-order derivatives on curved meshes.
- Show invariance under a change of coefficient basis and test sensitivity to
  equivalent parameterizations of the same physical geometry.
- Normalize curved face jumps with physical face measure and verify equal and
  opposite two-sided traces for a continuous manufactured state.
- Reuse the verified KXRCF mixed local and MPI cases for OFDG/OEDG.
- Run OFDG, OFDG--KXRCF, and OEDG checks on mixed meshes, then repeat the
  relevant filter and positivity checks on MFEM NURBS meshes.
- Document that filters must be reconstructed after mesh or nodal-coordinate
  changes.

Curvilinear support is complete only when the serial, MPI, conservation, and
manufactured-convergence tests all pass. Merely accepting a NURBS mesh is not
sufficient.
