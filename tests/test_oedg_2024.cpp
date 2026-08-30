#include "mfem.hpp"

#include "../src/oedg_2024.hpp"

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

real_t ElementIntegral(const GridFunction &state, int element)
{
   const FiniteElement *fe = state.FESpace()->GetFE(element);
   ElementTransformation *transformation =
      state.FESpace()->GetMesh()->GetElementTransformation(element);
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

void TestConstantAndConservation()
{
   Mesh mesh = Mesh::MakeCartesian1D(4, 1.0);
   DG_FECollection collection(2, 1, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection);
   GridFunction state(&space);
   OEDG2024 oedg(&space, BasisType::GaussLobatto,
                 std::make_shared<UnitFacePhysics>());

   state = 3.25;
   Vector decay;
   oedg.CompDecay(state, decay, 0.07);
   decay -= state;
   Require(decay.Normlinf() < 2e-12,
           "OEDG changed a constant state");

   FunctionCoefficient discontinuity([](const Vector &x)
   {
      return std::sin(3.0 * x(0)) + (x(0) > 0.5 ? 0.8 : 0.0);
   });
   state.ProjectCoefficient(discontinuity);
   oedg.CompDecay(state, decay, 0.03);
   GridFunction filtered(&space, decay);
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      Require(std::abs(ElementIntegral(state, e) -
                       ElementIntegral(filtered, e)) < 3e-12,
              "OEDG changed an element mean");
   }
}

void TestScaleAndShiftInvariance()
{
   Mesh mesh = Mesh::MakeCartesian1D(5, 1.0);
   DG_FECollection collection(3, 1, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection);
   GridFunction state(&space);
   FunctionCoefficient wave([](const Vector &x)
   {
      return std::sin(2.0 * M_PI * x(0)) + (x(0) > 0.4 ? 0.3 : 0.0);
   });
   state.ProjectCoefficient(wave);

   OEDG2024 oedg(&space, BasisType::GaussLobatto,
                 std::make_shared<UnitFacePhysics>());
   Vector reference;
   oedg.CompDecay(state, reference, 0.02);

   const real_t scale = 17.0;
   const real_t shift = -4.5;
   Vector transformed(state);
   transformed *= scale;
   transformed += shift;
   Vector transformed_decay;
   oedg.CompDecay(transformed, transformed_decay, 0.02);
   transformed_decay -= shift;
   transformed_decay /= scale;
   transformed_decay -= reference;
   Require(transformed_decay.Normlinf() < 3e-11,
           "OEDG lost scale/shift invariance");
}

} // namespace

int main()
{
   try
   {
      TestConstantAndConservation();
      TestScaleAndShiftInvariance();
   }
   catch (const std::exception &error)
   {
      std::cerr << "OEDG 2024 tests failed: " << error.what() << '\n';
      return 1;
   }
   std::cout << "OEDG 2024 tests passed.\n";
   return 0;
}
