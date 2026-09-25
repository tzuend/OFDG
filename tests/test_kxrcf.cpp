#include "mfem.hpp"

#include "../src/kxrcf.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace mfem;
using namespace ofdg;

namespace
{

void Require(bool condition, const std::string &message)
{
   if (!condition)
   {
      throw std::runtime_error(message);
   }
}

bool Near(real_t actual, real_t expected,
          real_t relative_tolerance = 1e-12,
          real_t absolute_tolerance = 1e-13)
{
   const real_t scale = std::max(real_t(1.0),
                                 std::max(std::abs(actual), std::abs(expected)));
   return std::abs(actual - expected) <=
          absolute_tolerance + relative_tolerance * scale;
}

int CountActive(const Array<bool> &active)
{
   int count = 0;
   for (int e = 0; e < active.Size(); ++e)
   {
      count += active[e] ? 1 : 0;
   }
   return count;
}

Mesh MakeMixedMesh()
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

void SetGeometryState(FiniteElementSpace &space, Vector &state)
{
   state = 0.0;
   Array<int> dofs;
   for (int e = 0; e < space.GetNE(); ++e)
   {
      space.GetElementDofs(e, dofs);
      const real_t label = space.GetFE(e)->GetGeomType() ==
                           Geometry::TRIANGLE ? 2.0 : 1.0;
      for (int j = 0; j < dofs.Size(); ++j)
      {
         state(space.DofToVDof(dofs[j], 0)) = 1.0;
         state(space.DofToVDof(dofs[j], 1)) = label;
      }
   }
}

VectorFunctionCoefficient ConstantVelocity(int dim, real_t x_velocity)
{
   return VectorFunctionCoefficient(
      dim, [dim, x_velocity](const Vector &, Vector &velocity)
      {
         velocity.SetSize(dim);
         velocity = 0.0;
         velocity(0) = x_velocity;
      });
}

void TestZeroAndConstantStates()
{
   Mesh mesh = Mesh::MakeCartesian2D(4, 3, Element::QUADRILATERAL,
                                     true, 1.0, 1.0);
   L2_FECollection fec(2, 2);
   FiniteElementSpace fes(&mesh, &fec, 2, Ordering::byNODES);
   GridFunction state(&fes);
   auto velocity = ConstantVelocity(2, 1.0);
   KXRCFIndicator indicator(&fes, &velocity);
   Array<bool> active;
   Vector values;

   state = 0.0;
   indicator.Compute(state, active, &values);
   Require(CountActive(active) == 0, "zero state activated cells");
   Require(values.Normlinf() == 0.0, "zero state produced an indicator");

   state = 2.5;
   indicator.Compute(state, active, &values);
   Require(CountActive(active) == 0, "constant state activated cells");
   Require(values.Normlinf() < 1e-12,
           "constant state produced a nonzero indicator");
}

void TestOneDimensionalFingerprint()
{
   Mesh mesh = Mesh::MakeCartesian1D(4, 1.0);
   L2_FECollection fec(1, 1);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction state(&fes);
   FunctionCoefficient jump([](const Vector &x)
   {
      return x(0) < 0.5 ? 1.0 : 2.0;
   });
   state.ProjectCoefficient(jump);

   for (const real_t direction : {real_t(1.0), real_t(-1.0)})
   {
      auto velocity = ConstantVelocity(1, direction);
      KXRCFIndicator indicator(&fes, &velocity);
      Array<bool> active;
      Vector values;
      DenseMatrix component_values;
      indicator.Compute(state, active, &values, &component_values);

      Require(CountActive(active) == 1,
              "1D discontinuity must activate exactly one inflow cell");

      const int expected_element = direction > 0.0 ? 2 : 1;
      const real_t expected_indicator = direction > 0.0 ? 4.0 : 8.0;

      for (int e = 0; e < mesh.GetNE(); ++e)
      {
         Require(active[e] == (e == expected_element),
                 "1D active-cell fingerprint changed");
         const real_t expected = e == expected_element ? expected_indicator : 0.0;
         Require(Near(values(e), expected),
                 "1D pooled indicator fingerprint changed");
         Require(Near(component_values(0, e), expected),
                 "1D component indicator fingerprint changed");
      }
   }
}

void TestTwoDimensionalInflowLayer()
{
   constexpr int nx = 8;
   constexpr int ny = 6;
   Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL,
                                     true, 1.0, 1.0);
   L2_FECollection fec(1, 2);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction state(&fes);
   FunctionCoefficient jump([](const Vector &x)
   {
      return x(0) < 0.5 ? 1.0 : 2.0;
   });
   state.ProjectCoefficient(jump);

   for (const real_t direction : {real_t(1.0), real_t(-1.0)})
   {
      auto velocity = ConstantVelocity(2, direction);
      KXRCFIndicator indicator(&fes, &velocity);
      Array<bool> active;
      indicator.Compute(state, active);
      Require(CountActive(active) == ny,
              "2D discontinuity must activate one complete inflow layer");

      const real_t expected_x = direction > 0.0
                                ? 0.5 + 0.5 / nx
                                : 0.5 - 0.5 / nx;
      for (int e = 0; e < mesh.GetNE(); ++e)
      {
         Vector center;
         ElementTransformation *transformation =
            mesh.GetElementTransformation(e);
         transformation->Transform(
            Geometries.GetCenter(mesh.GetElementGeometry(e)), center);
         Require(active[e] == (std::abs(center(0) - expected_x) < 1e-12),
                 "2D inflow-layer fingerprint changed");
      }
   }
}

void TestMultiplicativeScalingAndPooling()
{
   Mesh mesh = Mesh::MakeCartesian2D(6, 4, Element::TRIANGLE,
                                     true, 1.0, 1.0);
   L2_FECollection fec(1, 2);
   FiniteElementSpace fes(&mesh, &fec, 3, Ordering::byNODES);
   GridFunction state(&fes);
   VectorFunctionCoefficient initial(3, [](const Vector &x, Vector &value)
   {
      value.SetSize(3);
      value(0) = 1.0;
      value(1) = x(0) < 0.5 ? 1.0 : 3.0;
      value(2) = 2.0;
   });
   state.ProjectCoefficient(initial);
   auto velocity = ConstantVelocity(2, 1.0);
   KXRCFIndicator indicator(&fes, &velocity);

   Array<bool> active;
   Vector values;
   DenseMatrix components;
   indicator.Compute(state, active, &values, &components);
   Require(CountActive(active) > 0,
           "a jump in one component must activate pooled cells");

   Array<bool> scaled_active;
   Vector scaled_values;
   DenseMatrix scaled_components;
   state *= 17.0;
   indicator.Compute(state, scaled_active, &scaled_values, &scaled_components);

   Require(active.Size() == scaled_active.Size(), "scaled mask size changed");
   for (int e = 0; e < active.Size(); ++e)
   {
      Require(active[e] == scaled_active[e], "multiplicative scaling changed mask");
      Require(Near(values(e), scaled_values(e)),
              "multiplicative scaling changed pooled indicator");
      for (int c = 0; c < components.Height(); ++c)
      {
         Require(Near(components(c, e), scaled_components(c, e)),
                 "multiplicative scaling changed component indicator");
      }
   }
}

void TestVectorOrderingInvariance()
{
   Mesh mesh = Mesh::MakeCartesian2D(5, 4, Element::QUADRILATERAL,
                                     true, 1.0, 1.0);
   L2_FECollection fec(2, 2);
   FiniteElementSpace nodes_fes(&mesh, &fec, 2, Ordering::byNODES);
   FiniteElementSpace vdim_fes(&mesh, &fec, 2, Ordering::byVDIM);
   GridFunction nodes_state(&nodes_fes);
   GridFunction vdim_state(&vdim_fes);
   VectorFunctionCoefficient initial(2, [](const Vector &x, Vector &value)
   {
      value.SetSize(2);
      value(0) = 1.0 + 0.2 * std::sin(2.0 * M_PI * x(0));
      value(1) = x(1) < 0.5 ? 1.0 : 2.0;
   });
   nodes_state.ProjectCoefficient(initial);
   vdim_state.ProjectCoefficient(initial);
   auto velocity = ConstantVelocity(2, 0.4);

   KXRCFIndicator nodes_indicator(&nodes_fes, &velocity);
   KXRCFIndicator vdim_indicator(&vdim_fes, &velocity);
   Array<bool> nodes_active, vdim_active;
   Vector nodes_values, vdim_values;
   nodes_indicator.Compute(nodes_state, nodes_active, &nodes_values);
   vdim_indicator.Compute(vdim_state, vdim_active, &vdim_values);

   Require(nodes_active.Size() == vdim_active.Size(), "ordering changed mask size");
   for (int e = 0; e < nodes_active.Size(); ++e)
   {
      Require(nodes_active[e] == vdim_active[e], "ordering changed active mask");
      Require(Near(nodes_values(e), vdim_values(e), 2e-12),
              "ordering changed indicator values");
   }
}

void TestMixedTriangleQuadrilateralMesh()
{
   Mesh mesh = MakeMixedMesh();
   L2_FECollection fec(2, 2);
   FiniteElementSpace space(&mesh, &fec, 2, Ordering::byNODES);
   GridFunction state(&space);
   auto velocity = ConstantVelocity(2, 1.0);
   KXRCFIndicator indicator(&space, &velocity, 1e-10);
   Array<bool> active;
   Vector values;
   DenseMatrix components;

   bool found_triangle = false;
   bool found_quadrilateral = false;
   for (int e = 0; e < space.GetNE(); ++e)
   {
      found_triangle |= space.GetFE(e)->GetDof() == 6;
      found_quadrilateral |= space.GetFE(e)->GetDof() == 9;
   }
   Require(found_triangle && found_quadrilateral,
           "mixed P2 mesh does not contain both element DOF counts");

   state = 0.0;
   indicator.Compute(state, active, &values);
   Require(CountActive(active) == 0 && values.Normlinf() == 0.0,
           "mixed mesh zero state produced an indicator");

   state = 2.5;
   indicator.Compute(state, active, &values);
   Require(CountActive(active) == 0 && values.Normlinf() < 1e-12,
           "mixed mesh constant state produced an indicator");

   SetGeometryState(space, state);
   indicator.Compute(state, active, &values, &components);
   Require(CountActive(active) > 0,
           "mixed-face discontinuity did not activate any cells");
   for (int e = 0; e < space.GetNE(); ++e)
   {
      Require(std::isfinite(values(e)),
              "mixed mesh produced a non-finite indicator");
      Require(Near(components(0, e), 0.0),
              "constant component produced a mixed-mesh indicator");
      Require(Near(values(e), components(1, e)),
              "mixed-mesh component pooling changed");
   }

   const Vector reference_values(values);
   const Array<bool> reference_active(active);
   state *= 19.0;
   indicator.Compute(state, active, &values);
   for (int e = 0; e < space.GetNE(); ++e)
   {
      Require(active[e] == reference_active[e],
              "scaling changed the mixed-mesh active mask");
      Require(Near(values(e), reference_values(e), 2e-12),
              "scaling changed a mixed-mesh indicator");
   }
}

// Cached geometry must not freeze flow physics or depend on MFEM's shared
// transformation scratch, which other operators overwrite between calls.
void TestCachedGeometryWithChangingVelocity()
{
   Mesh mesh = Mesh::MakeCartesian2D(3, 3, Element::QUADRILATERAL, true);
   mesh.SetCurvature(3);
   mesh.Transform([](const Vector &x, Vector &y)
   {
      y = x;
      y(0) += 0.03 * std::sin(3.0 * x(0)) * std::sin(3.0 * x(1));
   });
   L2_FECollection fec(2, 2);
   FiniteElementSpace space(&mesh, &fec);
   GridFunction state(&space);
   FunctionCoefficient initial([](const Vector &x)
   {
      return x(0) < 0.5 ? 1.0 : 2.0;
   });
   state.ProjectCoefficient(initial);
   real_t direction = 1.0;
   VectorFunctionCoefficient velocity(2, [&direction](const Vector &x, Vector &v)
   {
      v.SetSize(2);
      v(0) = direction * (1.0 + x(1));
      v(1) = 0.2 * direction * x(0);
   });
   KXRCFIndicator reused(&space, &velocity);
   Vector positive_values;
   for (real_t sign : {1.0, -1.0})
   {
      direction = sign;
      for (int f = 0; f < mesh.GetNumFaces(); ++f)
      {
         mesh.GetFaceElementTransformations(f);
      }
      Array<bool> active, expected_active;
      Vector values, expected;
      reused.Compute(state, active, &values);
      KXRCFIndicator fresh(&space, &velocity);
      fresh.Compute(state, expected_active, &expected);
      for (int e = 0; e < space.GetNE(); ++e)
      {
         Require(active[e] == expected_active[e] && Near(values(e), expected(e)),
                 "cached face geometry changed the indicator after a flow update");
      }
      if (sign > 0.0) { positive_values = values; }
      else
      {
         values -= positive_values;
         Require(values.Normlinf() > 1e-8,
                 "changing velocity did not change the inflow indicator");
      }
   }
}

} // namespace

int main()
{
   try
   {
      TestZeroAndConstantStates();
      TestOneDimensionalFingerprint();
      TestTwoDimensionalInflowLayer();
      TestMultiplicativeScalingAndPooling();
      TestVectorOrderingInvariance();
      TestMixedTriangleQuadrilateralMesh();
      TestCachedGeometryWithChangingVelocity();
   }
   catch (const std::exception &error)
   {
      std::cerr << "KXRCF tests failed: " << error.what() << '\n';
      return 1;
   }

   std::cout << "KXRCF baseline tests passed.\n";
   return 0;
}
