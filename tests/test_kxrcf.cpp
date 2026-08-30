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
   }
   catch (const std::exception &error)
   {
      std::cerr << "KXRCF tests failed: " << error.what() << '\n';
      return 1;
   }

   std::cout << "KXRCF baseline tests passed.\n";
   return 0;
}
