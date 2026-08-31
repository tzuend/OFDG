#include "mfem.hpp"

#include "../src/oedg_2024.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

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

class StateMagnitudePhysics : public FacePhysics
{
public:
   FacePhysicsSample Evaluate(
      const Vector &state1, const Vector &state2, const Vector &,
      FaceElementTransformations &) const override
   {
      return {std::abs(state1(0)), std::abs(state2(0)), 0.0};
   }
};

void SetAffineElementState(FiniteElementSpace &space, GridFunction &state,
                           const std::vector<std::vector<real_t>> &offsets,
                           real_t slope)
{
   Array<int> scalar_dofs;
   Vector point;
   for (int e = 0; e < space.GetNE(); ++e)
   {
      const FiniteElement *element = space.GetFE(e);
      const IntegrationRule &nodes = element->GetNodes();
      ElementTransformation *transformation =
         space.GetMesh()->GetElementTransformation(e);
      space.GetElementDofs(e, scalar_dofs);
      for (int i = 0; i < scalar_dofs.Size(); ++i)
      {
         transformation->Transform(nodes.IntPoint(i), point);
         for (int c = 0; c < space.GetVDim(); ++c)
         {
            state(space.DofToVDof(scalar_dofs[i], c)) =
               slope * point(0) + offsets[c][e];
         }
      }
   }
}

void TestCellMeanWaveSpeed()
{
   Mesh mesh = Mesh::MakeCartesian1D(2, 1.0);
   DG_FECollection collection(1, 1, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection);
   GridFunction state(&space);
   SetAffineElementState(space, state, {{-0.5, 1.0}}, 2.0);

   OEDG2024 oedg(&space, BasisType::GaussLobatto,
                 std::make_shared<StateMagnitudePhysics>());
   Vector decay;
   oedg.CompDecay(state, decay, 0.1);

   Array<int> dofs;
   space.GetElementVDofs(0, dofs);
   for (int i = 0; i < dofs.Size(); ++i)
   {
      Require(std::abs(decay(dofs[i]) - state(dofs[i])) < 2e-12,
              "OEDG used a face trace instead of the zero cell mean for beta");
   }
}

void TestRelativeConstantTolerance()
{
   Mesh mesh = Mesh::MakeCartesian1D(2, 1.0);
   DG_FECollection collection(1, 1, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection);
   GridFunction state(&space);
   SetAffineElementState(space, state, {{1.0e12, 1.0e12 + 1.0e-3}}, 0.0);

   OEDG2024 oedg(&space, BasisType::GaussLobatto,
                 std::make_shared<UnitFacePhysics>());
   Vector decay;
   oedg.CompDecay(state, decay, 0.1);
   decay -= state;
   Require(decay.Normlinf() < 2e-12,
           "OEDG did not apply the paper's relative constant-state tolerance");
}

real_t ComponentAmplitude(const GridFunction &state, int component)
{
   const FiniteElementSpace *space = state.FESpace();
   real_t integral = 0.0;
   real_t volume = 0.0;
   real_t minimum = std::numeric_limits<real_t>::infinity();
   real_t maximum = -minimum;
   Vector value;
   for (int e = 0; e < space->GetNE(); ++e)
   {
      const FiniteElement *element = space->GetFE(e);
      ElementTransformation *transformation =
         space->GetMesh()->GetElementTransformation(e);
      const IntegrationRule &rule =
         IntRules.Get(element->GetGeomType(), 2 * element->GetOrder() + 1);
      for (int q = 0; q < rule.GetNPoints(); ++q)
      {
         const IntegrationPoint &point = rule.IntPoint(q);
         transformation->SetIntPoint(&point);
         state.GetVectorValue(e, point, value);
         const real_t weight = point.weight * transformation->Weight();
         integral += weight * value(component);
         volume += weight;
         minimum = std::min(minimum, value(component));
         maximum = std::max(maximum, value(component));
      }
   }
   const real_t mean = integral / volume;
   return std::max(mean - minimum, maximum - mean);
}

void TestPoolingOccursPerFace()
{
   Mesh mesh = Mesh::MakeCartesian1D(3, 1.0);
   DG_FECollection collection(1, 1, BasisType::GaussLobatto);
   FiniteElementSpace scalar_space(&mesh, &collection);
   FiniteElementSpace system_space(&mesh, &collection, 2,
                                   Ordering::byVDIM);
   GridFunction scalar_state(&scalar_space);
   GridFunction system_state(&system_space);
   SetAffineElementState(scalar_space, scalar_state,
                         {{1.0, 0.0, 0.0}}, 0.2);
   SetAffineElementState(system_space, system_state,
                         {{1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}, 0.2);

   auto physics = std::make_shared<UnitFacePhysics>();
   OEDG2024 scalar_oedg(&scalar_space, BasisType::GaussLobatto, physics);
   OEDG2024 system_oedg(&system_space, BasisType::GaussLobatto, physics);
   Vector scalar_stabilization, system_stabilization;
   scalar_oedg.ComputeStabilization(scalar_state, scalar_stabilization);
   system_oedg.ComputeStabilization(system_state, system_stabilization);
   const real_t amplitude0 = ComponentAmplitude(system_state, 0);
   const real_t amplitude1 = ComponentAmplitude(system_state, 1);
   const real_t expected_factor = 1.0 + amplitude0 / amplitude1;

   Array<int> scalar_dofs;
   Array<int> system_scalar_dofs;
   scalar_space.GetElementDofs(1, scalar_dofs);
   system_space.GetElementDofs(1, system_scalar_dofs);
   for (int i = 0; i < scalar_dofs.Size(); ++i)
   {
      const real_t scalar_value = scalar_stabilization(scalar_dofs[i]);
      const int system_dof = system_space.DofToVDof(system_scalar_dofs[i], 0);
      Require(std::abs(system_stabilization(system_dof) -
                       expected_factor * scalar_value) < 3e-11,
              "OEDG pooled components after summing different faces: system=" +
              std::to_string(system_stabilization(system_dof)) +
              ", scalar=" + std::to_string(scalar_value));
   }
}

} // namespace

int main()
{
   try
   {
      TestConstantAndConservation();
      TestScaleAndShiftInvariance();
      TestCellMeanWaveSpeed();
      TestRelativeConstantTolerance();
      TestPoolingOccursPerFace();
   }
   catch (const std::exception &error)
   {
      std::cerr << "OEDG 2024 tests failed: " << error.what() << '\n';
      return 1;
   }
   std::cout << "OEDG 2024 tests passed.\n";
   return 0;
}
