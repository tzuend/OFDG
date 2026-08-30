#include "mfem.hpp"

#include "../src/euler_positivity.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace mfem;

namespace
{

void Require(bool condition, const std::string &message)
{
   if (!condition) { throw std::runtime_error(message); }
}

Vector ElementMean(const FiniteElementSpace &space, const Vector &state,
                   int element)
{
   const FiniteElement *fe = space.GetFE(element);
   const int ndof = fe->GetDof();
   Array<int> dofs;
   space.GetElementDofs(element, dofs);
   Vector mean(space.GetVDim());
   mean = 0.0;
   Vector shape(ndof);
   ElementTransformation *transformation =
      space.GetMesh()->GetElementTransformation(element);
   const IntegrationRule &rule =
      IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 1);
   real_t volume = 0.0;
   for (int q = 0; q < rule.GetNPoints(); ++q)
   {
      const IntegrationPoint &ip = rule.IntPoint(q);
      transformation->SetIntPoint(&ip);
      fe->CalcShape(ip, shape);
      const real_t weight = ip.weight * transformation->Weight();
      volume += weight;
      for (int c = 0; c < space.GetVDim(); ++c)
      {
         for (int i = 0; i < ndof; ++i)
         {
            mean(c) += weight * shape(i) *
                       state(space.DofToVDof(dofs[i], c));
         }
      }
   }
   mean /= volume;
   return mean;
}

void SetNodeState(const FiniteElementSpace &space, Vector &state,
                  int scalar_dof, real_t density, real_t momentum,
                  real_t energy)
{
   state(space.DofToVDof(scalar_dof, 0)) = density;
   state(space.DofToVDof(scalar_dof, 1)) = momentum;
   state(space.DofToVDof(scalar_dof, 2)) = energy;
}

void TestLimiterPreservesMeanAndPositivity()
{
   Mesh mesh = Mesh::MakeCartesian1D(1, 1.0);
   DG_FECollection collection(1, 1, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection, 3, Ordering::byNODES);
   Array<int> dofs;
   space.GetElementDofs(0, dofs);
   Vector state(space.GetVSize());
   SetNodeState(space, state, dofs[0], -0.2, 0.0, 2.5);
   SetNodeState(space, state, dofs[1], 2.2, 0.0, 2.5);
   const Vector before = ElementMean(space, state, 0);

   EulerPositivityLimiter limiter(&space, 1.4, 1e-10, 1e-10);
   const EulerPositivityDiagnostics diagnostics = limiter.Apply(state);
   const Vector after = ElementMean(space, state, 0);
   Vector mean_difference(after);
   mean_difference -= before;
   Require(diagnostics.limited_elements == 1,
           "negative nodal density did not activate limiter");
   Require(mean_difference.Normlinf() < 2e-12,
           "positivity limiter changed the cell mean");
   Require(state(space.DofToVDof(dofs[0], 0)) >= 0.9e-10,
           "positivity limiter left a negative nodal density");
}

void TestAdmissibleAndImpossibleMeans()
{
   Mesh mesh = Mesh::MakeCartesian1D(1, 1.0);
   DG_FECollection collection(1, 1, BasisType::GaussLobatto);
   FiniteElementSpace space(&mesh, &collection, 3, Ordering::byNODES);
   Array<int> dofs;
   space.GetElementDofs(0, dofs);
   Vector state(space.GetVSize());
   for (int i = 0; i < dofs.Size(); ++i)
   {
      SetNodeState(space, state, dofs[i], 1.0, 0.1, 2.6);
   }
   Vector original(state);
   EulerPositivityLimiter limiter(&space, 1.4);
   EulerPositivityDiagnostics diagnostics = limiter.Apply(state);
   state -= original;
   Require(diagnostics.limited_elements == 0 && state.Normlinf() < 1e-14,
           "limiter changed a constant admissible state");

   for (int i = 0; i < dofs.Size(); ++i)
   {
      SetNodeState(space, state, dofs[i], -1.0, 0.0, 2.5);
   }
   diagnostics = limiter.Apply(state);
   Require(diagnostics.inadmissible_means == 1,
           "inadmissible cell mean was not reported");
}

} // namespace

int main()
{
   try
   {
      TestLimiterPreservesMeanAndPositivity();
      TestAdmissibleAndImpossibleMeans();
   }
   catch (const std::exception &error)
   {
      std::cerr << "Euler positivity tests failed: " << error.what() << '\n';
      return 1;
   }
   std::cout << "Euler positivity tests passed.\n";
   return 0;
}
