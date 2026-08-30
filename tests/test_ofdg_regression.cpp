#include "mfem.hpp"

#include "../src/ofdg_serial_optimized.hpp"

#include <algorithm>
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
   if (!condition) { throw std::runtime_error(message); }
}

bool Near(real_t actual, real_t expected,
          real_t relative_tolerance = 2e-12,
          real_t absolute_tolerance = 2e-13)
{
   const real_t scale = std::max(real_t(1.0),
                                 std::max(std::abs(actual), std::abs(expected)));
   return std::abs(actual - expected) <=
          absolute_tolerance + relative_tolerance * scale;
}

void RequireVectorNear(const Vector &actual, const real_t *expected,
                       int size, const std::string &name)
{
   Require(actual.Size() == size, name + " size changed");
   for (int i = 0; i < size; ++i)
   {
      Require(Near(actual(i), expected[i]),
              name + " changed at degree of freedom " + std::to_string(i));
   }
}

real_t ElementIntegral(const GridFunction &state, int element)
{
   const FiniteElementSpace *fes = state.FESpace();
   const FiniteElement *fe = fes->GetFE(element);
   ElementTransformation *transformation =
      fes->GetMesh()->GetElementTransformation(element);
   const IntegrationRule &rule =
      IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 1);

   real_t integral = 0.0;
   for (int q = 0; q < rule.GetNPoints(); ++q)
   {
      const IntegrationPoint &ip = rule.IntPoint(q);
      transformation->SetIntPoint(&ip);
      integral += ip.weight * transformation->Weight() *
                  state.GetValue(element, ip);
   }
   return integral;
}

void TestFrozenOneDimensionalFingerprint()
{
   Mesh mesh = Mesh::MakeCartesian1D(4, 1.0);
   DG_FECollection fec(2, 1, BasisType::GaussLobatto);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction state(&fes);
   FunctionCoefficient jump([](const Vector &x)
   {
      return x(0) < 0.5 ? 1.0 : 2.0;
   });
   state.ProjectCoefficient(jump);

   VectorFunctionCoefficient velocity(1, [](const Vector &, Vector &value)
   {
      value.SetSize(1);
      value(0) = 1.0;
   });
   OFDG ofdg(&fes, BasisType::GaussLobatto, velocity);

   Vector stabilization, decay;
   ofdg.ComputeStabilization(state, stabilization);
   ofdg.CompDecay(state, decay, 0.0375);
   GridFunction decayed_state(&fes, decay);

   const real_t expected_stabilization[] = {
      2.9419311152784784e-15, -3.2361242268063259e-16,
      -2.9419311152784784e-15, 4.9463978198199943,
      -5.6530260797942731, 17.665706499357086,
      5.8838622305569568e-15, -3.1184469821951871e-15,
      -5.8838622305569568e-15, 0.0, 0.0, 0.0
   };
   const real_t expected_decay[] = {
      0.99999999999999989, 1.0, 1.0000000000000002,
      0.94976732237602679, 1.1199516839719805, 1.5704259417360513,
      1.9999999999999998, 2.0, 2.0, 2.0, 2.0, 2.0
   };

   RequireVectorNear(stabilization, expected_stabilization, 12,
                     "OFDG stabilization fingerprint");
   RequireVectorNear(decay, expected_decay, 12,
                     "OFDG decay fingerprint");

   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      Require(Near(ElementIntegral(state, e),
                   ElementIntegral(decayed_state, e)),
              "OFDG decay no longer preserves the element mean");
   }
}

void TestActiveMaskLeavesInactiveElementsUntouched()
{
   Mesh mesh = Mesh::MakeCartesian1D(4, 1.0);
   DG_FECollection fec(2, 1, BasisType::GaussLobatto);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction state(&fes);
   FunctionCoefficient wave([](const Vector &x)
   {
      return std::sin(3.0 * x(0)) + (x(0) > 0.5 ? 0.4 : 0.0);
   });
   state.ProjectCoefficient(wave);
   OFDG ofdg(&fes, BasisType::GaussLobatto);

   Array<bool> active(mesh.GetNE());
   active = false;
   active[1] = true;

   Vector stabilization, decay;
   ofdg.ComputeStabilization(state, stabilization, &active);
   ofdg.CompDecay(state, decay, 0.08, &active);

   Array<int> vdofs;
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      if (active[e]) { continue; }
      fes.GetElementVDofs(e, vdofs);
      for (int i = 0; i < vdofs.Size(); ++i)
      {
         Require(std::abs(stabilization(vdofs[i])) < 1e-13,
                 "inactive element received stabilization");
         Require(Near(decay(vdofs[i]), state(vdofs[i])),
                 "inactive element was changed by decay");
      }
   }
}

} // namespace

int main()
{
   try
   {
      TestFrozenOneDimensionalFingerprint();
      TestActiveMaskLeavesInactiveElementsUntouched();
   }
   catch (const std::exception &error)
   {
      std::cerr << "OFDG regression tests failed: " << error.what() << '\n';
      return 1;
   }

   std::cout << "OFDG frozen regression tests passed.\n";
   return 0;
}
