#include "mfem.hpp"

#include "../src/kxrcf.hpp"
#include "../src/oedg_2024.hpp"
#include "../src/ofdg.hpp"

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using namespace mfem;
using namespace ofdg;

namespace
{

void Require(bool condition, const std::string &message)
{
   if (!condition) { throw std::runtime_error(message); }
}

Mesh MakeMixedTwoDimensionalMesh()
{
   Mesh mesh(2, 6, 3);
   mesh.AddVertex(0.0, 0.0);
   mesh.AddVertex(1.0, 0.0);
   mesh.AddVertex(2.0, 0.0);
   mesh.AddVertex(0.0, 1.0);
   mesh.AddVertex(1.0, 1.0);
   mesh.AddVertex(2.0, 1.0);
   mesh.AddQuad(0, 1, 4, 3);
   mesh.AddTriangle(1, 2, 5);
   mesh.AddTriangle(1, 5, 4);
   mesh.FinalizeTopology();
   return mesh;
}

Mesh MakePrismMesh()
{
   Mesh mesh(3, 8, 2);
   mesh.AddVertex(0.0, 0.0, 0.0);
   mesh.AddVertex(1.0, 0.0, 0.0);
   mesh.AddVertex(1.0, 1.0, 0.0);
   mesh.AddVertex(0.0, 1.0, 0.0);
   mesh.AddVertex(0.0, 0.0, 1.0);
   mesh.AddVertex(1.0, 0.0, 1.0);
   mesh.AddVertex(1.0, 1.0, 1.0);
   mesh.AddVertex(0.0, 1.0, 1.0);
   int vertices[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   mesh.AddHexAsWedges(vertices);
   mesh.FinalizeTopology();
   return mesh;
}

void ContinuousPolynomial(const Vector &point, Vector &value)
{
   value.SetSize(2);
   value(0) = 1.0 + 0.2 * point(0) - 0.3 * point(1) +
              (point.Size() == 3 ? 0.4 * point(2) : 0.0);
   value(1) = -0.2 + 0.5 * point(0) + 0.1 * point(1) -
              (point.Size() == 3 ? 0.2 * point(2) : 0.0);
}

real_t ElementIntegral(const GridFunction &state, int element, int component)
{
   const FiniteElementSpace *space = state.FESpace();
   const FiniteElement *finite_element = space->GetFE(element);
   ElementTransformation *transformation =
      space->GetMesh()->GetElementTransformation(element);
   const IntegrationRule &rule = IntRules.Get(
      finite_element->GetGeomType(), 2 * finite_element->GetOrder() + 2);
   Vector value;
   real_t integral = 0.0;
   for (int q = 0; q < rule.GetNPoints(); ++q)
   {
      const IntegrationPoint &point = rule.IntPoint(q);
      transformation->SetIntPoint(&point);
      state.GetVectorValue(element, point, value);
      integral += point.weight * transformation->Weight() * value(component);
   }
   return integral;
}

void ConstantVelocity(const Vector &point, Vector &velocity)
{
   velocity.SetSize(point.Size());
   velocity = 1.0;
}

void ValidateMixedMesh(Mesh &mesh, int order, int ordering,
                       const std::string &name,
                       bool require_multiple_signatures = true)
{
   const int dimension = mesh.Dimension();
   DG_FECollection collection(order, dimension, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection, 2, ordering);
   GridFunction state(&space);
   VectorFunctionCoefficient polynomial(2, ContinuousPolynomial);
   state.ProjectCoefficient(polynomial);

   std::set<int> dof_counts;
   std::set<int> geometries;
   for (int e = 0; e < space.GetNE(); ++e)
   {
      dof_counts.insert(space.GetFE(e)->GetDof());
      geometries.insert(static_cast<int>(space.GetFE(e)->GetGeomType()));
   }
   if (require_multiple_signatures)
   {
      Require(dof_counts.size() > 1,
              name + " does not exercise unequal element DOF counts");
      Require(geometries.size() > 1,
              name + " does not contain multiple element geometries");
   }

   VectorFunctionCoefficient velocity(dimension, ConstantVelocity);
   auto physics = std::make_shared<AdvectionFacePhysics>(&velocity);
   OFDG ofdg(&space, BasisType::GaussLobatto, physics);
   OEDG2024 oedg(&space, BasisType::GaussLobatto, physics);
   KXRCFIndicator kxrcf(&space, physics, 1e-10);

   Vector stabilization, decay, oedg_decay, indicators;
   Array<bool> active;
   ofdg.ComputeStabilization(state, stabilization);
   ofdg.CompDecay(state, decay, 0.03);
   oedg.CompDecay(state, oedg_decay, 0.03);
   kxrcf.Compute(state, active, &indicators);

   Require(stabilization.Normlinf() < 2e-8,
           name + " stabilized a continuous affine polynomial: " +
           std::to_string(stabilization.Normlinf()));
   decay -= state;
   oedg_decay -= state;
   Require(decay.Normlinf() < 2e-8,
           name + " OFDG changed a continuous affine polynomial");
   Require(oedg_decay.Normlinf() < 2e-8,
           name + " OEDG changed a continuous affine polynomial");
   Require(indicators.Normlinf() < 2e-8,
           name + " KXRCF marked a continuous affine polynomial");

   state.ProjectCoefficient(polynomial);
   Array<int> scalar_dofs;
   space.GetElementDofs(0, scalar_dofs);
   for (int i = 0; i < scalar_dofs.Size(); ++i)
   {
      state(space.DofToVDof(scalar_dofs[i], 0)) += 0.4;
      state(space.DofToVDof(scalar_dofs[i], 1)) -= 0.2;
   }

   kxrcf.Compute(state, active, &indicators);
   ofdg.ComputeStabilization(state, stabilization, &active);
   ofdg.CompDecay(state, decay, 0.03, &active);
   oedg.CompDecay(state, oedg_decay, 0.03);
   Require(std::isfinite(stabilization.Norml2()) &&
           std::isfinite(decay.Norml2()) &&
           std::isfinite(oedg_decay.Norml2()) &&
           std::isfinite(indicators.Norml2()),
           name + " produced a non-finite result");
   if (require_multiple_signatures)
   {
      Require(indicators.Normlinf() > 1e-8,
              name + " KXRCF missed a geometry-local discontinuity");
   }

   GridFunction ofdg_result(&space, decay);
   GridFunction oedg_result(&space, oedg_decay);
   for (int e = 0; e < space.GetNE(); ++e)
   {
      for (int c = 0; c < 2; ++c)
      {
         Require(std::abs(ElementIntegral(state, e, c) -
                          ElementIntegral(ofdg_result, e, c)) < 2e-10,
                 name + " OFDG changed an element mean");
         Require(std::abs(ElementIntegral(state, e, c) -
                          ElementIntegral(oedg_result, e, c)) < 2e-10,
                 name + " OEDG changed an element mean");
      }
   }
}

} // namespace

int main()
{
   try
   {
      for (int order = 1; order <= 3; ++order)
      {
         Mesh mixed_2d = MakeMixedTwoDimensionalMesh();
         ValidateMixedMesh(mixed_2d, order, Ordering::byVDIM,
                           "mixed triangle-quadrilateral P" +
                           std::to_string(order));

         Mesh mixed_3d(std::string(OFDG_MFEM_DATA_DIR) +
                       "/fichera-mixed.mesh");
         ValidateMixedMesh(mixed_3d, order, Ordering::byVDIM,
                           "mixed tetrahedron-hexahedron-prism P" +
                           std::to_string(order));
      }

      Mesh mixed_2d_by_nodes = MakeMixedTwoDimensionalMesh();
      ValidateMixedMesh(mixed_2d_by_nodes, 2, Ordering::byNODES,
                        "mixed triangle-quadrilateral P2 byNODES");
      Mesh mixed_3d_by_nodes(std::string(OFDG_MFEM_DATA_DIR) +
                             "/fichera-mixed.mesh");
      ValidateMixedMesh(mixed_3d_by_nodes, 2, Ordering::byNODES,
                        "mixed tetrahedron-hexahedron-prism P2 byNODES");

      Mesh tetrahedra = Mesh::MakeCartesian3D(
         1, 1, 1, Element::TETRAHEDRON);
      ValidateMixedMesh(tetrahedra, 2, Ordering::byVDIM,
                        "homogeneous tetrahedra P2", false);
      Mesh hexahedra = Mesh::MakeCartesian3D(
         2, 1, 1, Element::HEXAHEDRON);
      ValidateMixedMesh(hexahedra, 2, Ordering::byVDIM,
                        "homogeneous hexahedra P2", false);
      Mesh prisms = MakePrismMesh();
      ValidateMixedMesh(prisms, 2, Ordering::byVDIM,
                        "homogeneous prisms P2", false);
   }
   catch (const std::exception &error)
   {
      std::cerr << "Mixed geometry tests failed: " << error.what() << '\n';
      return 1;
   }

   std::cout << "Mixed affine 2D and 3D filter tests passed.\n";
   return 0;
}
