# Mixed-element and curvilinear design decision

## Decision

OFDG, OEDG, and KXRCF now use the finite element and actual DOF count from each
side of every face. Their affine path has no geometry allow-list. Verification
covers mixed triangle--quadrilateral meshes and MFEM's mixed
tetrahedron--hexahedron--prism mesh, including unequal-signature MPI faces.
The remaining production geometry task is static curvilinear mappings.

The completed mixed-element refactor remains separate from curved-map changes
and preserves the public `ofdg::` interfaces and homogeneous fingerprints.

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

## Implemented affine operator architecture

An immutable repository replaces the former operator built from `GetFE(0)` and is
keyed by an element signature containing geometry, polynomial order, scalar DOF
count, map type, and coefficient-basis identifier. Each local element stores a
handle to its operator set rather than assuming a common matrix size.

Face cache entries own independently sized evaluation matrices for element 1
and element 2. Runtime scratch storage grows to the largest local signature but
views retain the actual side-specific dimensions. No padding participates in a
matrix operation.

For shared MPI faces, face-neighbour metadata is initialized before the
repository is frozen. The remote finite element is obtained with
`ParFiniteElementSpace::GetFaceNbrFE`, obtain its vector DOFs with
`GetFaceNbrElementVDofs`, and construct the remote derivative state from that
finite element's signature and transformation. Ownership remains unchanged:
each rank accumulates damping only for its local element.

The in-memory affine 2D fixture contains a parallelogram and two triangles; its
explicit two-rank partition places a triangle--quadrilateral face across ranks.
The pinned MFEM `fichera-mixed.mesh` fixture similarly exercises tetrahedra,
hexahedra, prisms, triangular and quadrilateral faces in serial and across two
ranks. All three methods use these fixtures.

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
- Retain the verified mixed local and MPI cases for all three filters.
- Repeat the relevant filter and positivity checks on MFEM NURBS meshes.
- Document that filters must be reconstructed after mesh or nodal-coordinate
  changes.

Curvilinear support is complete only when the serial, MPI, conservation, and
manufactured-convergence tests all pass. Merely accepting a NURBS mesh is not
sufficient.
