#include "mfem.hpp"

#include "../src/kxrcf.hpp"
#include "../src/oedg_2024.hpp"
#include "../src/euler_positivity.hpp"
#include "../src/ofdg.hpp"

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

VectorFunctionCoefficient PositiveVelocity()
{
   return VectorFunctionCoefficient(1, [](const Vector &, Vector &velocity)
   {
      velocity.SetSize(1);
      velocity(0) = 1.0;
   });
}

Mesh MakeMixedMesh()
{
   Mesh mesh(2, 6, 4);
   mesh.AddVertex(0.0, 0.0);
   mesh.AddVertex(2.0, 0.0);
   mesh.AddVertex(2.0, 1.0);
   mesh.AddVertex(0.0, 1.0);
   mesh.AddVertex(0.5, 0.5);
   mesh.AddVertex(1.5, 0.5);
   mesh.AddQuad(0, 1, 5, 4);
   mesh.AddTriangle(1, 2, 5);
   mesh.AddQuad(2, 3, 4, 5);
   mesh.AddTriangle(3, 0, 4);
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

real_t ElementIntegral(const GridFunction &state, int element)
{
   const FiniteElementSpace *space = state.FESpace();
   const FiniteElement *fe = space->GetFE(element);
   ElementTransformation *transformation =
      space->GetMesh()->GetElementTransformation(element);
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

void TestMixedMeshConsistency()
{
   Mesh serial_mesh = MakeMixedMesh();
   L2_FECollection serial_collection(2, 2);
   FiniteElementSpace serial_space(&serial_mesh, &serial_collection, 2,
                                   Ordering::byNODES);
   GridFunction serial_state(&serial_space);
   SetGeometryState(serial_space, serial_state);
   VectorFunctionCoefficient serial_velocity(
      2, [](const Vector &, Vector &velocity)
      {
         velocity.SetSize(2);
         velocity = 0.0;
         velocity(0) = 1.0;
      });
   KXRCFIndicator serial_indicator(&serial_space, &serial_velocity, 1e-10);
   Array<bool> serial_active;
   Vector serial_values;
   serial_indicator.Compute(serial_state, serial_active, &serial_values);

   int expected_active = 0;
   real_t expected_sum = 0.0;
   real_t expected_maximum = 0.0;
   for (int e = 0; e < serial_space.GetNE(); ++e)
   {
      expected_active += serial_active[e] ? 1 : 0;
      expected_sum += serial_values(e);
      expected_maximum = std::max(expected_maximum, serial_values(e));
   }

   int partitioning[4] = {0, 1, 0, 1};
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh, partitioning);
   L2_FECollection collection(2, 2);
   ParFiniteElementSpace space(&mesh, &collection, 2, Ordering::byNODES);
   ParGridFunction state(&space);
   SetGeometryState(space, state);
   VectorFunctionCoefficient velocity(
      2, [](const Vector &, Vector &value)
      {
         value.SetSize(2);
         value = 0.0;
         value(0) = 1.0;
      });
   KXRCFIndicator indicator(&space, &velocity, 1e-10);
   Array<bool> active;
   Vector values;
   indicator.Compute(state, active, &values);

   int local_active = 0;
   real_t local_sum = 0.0;
   real_t local_maximum = 0.0;
   for (int e = 0; e < space.GetNE(); ++e)
   {
      local_active += active[e] ? 1 : 0;
      local_sum += values(e);
      local_maximum = std::max(local_maximum, values(e));
   }

   int global_active = 0;
   real_t global_sum = 0.0;
   real_t global_maximum = 0.0;
   MPI_Allreduce(&local_active, &global_active, 1, MPI_INT, MPI_SUM,
                 mesh.GetComm());
   MPI_Allreduce(&local_sum, &global_sum, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, mesh.GetComm());
   MPI_Allreduce(&local_maximum, &global_maximum, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, mesh.GetComm());

   Require(global_active == expected_active,
           "distributed mixed-mesh active count differs from serial");
   Require(std::abs(global_sum - expected_sum) < 2e-12,
           "distributed mixed-mesh indicator sum differs from serial");
   Require(std::abs(global_maximum - expected_maximum) < 2e-12,
           "distributed mixed-mesh indicator maximum differs from serial");

   state = 3.0;
   indicator.Compute(state, active, &values);
   local_maximum = values.Size() == 0 ? 0.0 : values.Normlinf();
   MPI_Allreduce(&local_maximum, &global_maximum, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, mesh.GetComm());
   Require(global_maximum < 1e-12,
           "distributed mixed-mesh constant state produced an indicator");
}

void RunParallelChecks()
{
   Mesh serial_mesh = Mesh::MakeCartesian1D(4, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   DG_FECollection collection(2, 1, BasisType::GaussLobatto);
   ParFiniteElementSpace space(&mesh, &collection);
   ParGridFunction state(&space);
   Array<int> element_vdofs;
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      Vector center;
      mesh.GetElementTransformation(e)->Transform(
         Geometries.GetCenter(mesh.GetElementGeometry(e)), center);
      space.GetElementVDofs(e, element_vdofs);
      state.SetSubVector(element_vdofs, center(0) < 0.5 ? 1.0 : 2.0);
   }

   Require(mesh.GetNSharedFaces() > 0,
           "parallel test partition has no shared faces");

   // The positivity limiter must be safe as the first object that requests
   // shared-face transformations. This catches missing MFEM face-neighbor
   // initialization instead of letting an earlier filter hide it.
   ParFiniteElementSpace euler_space(
      &mesh, &collection, 3, Ordering::byVDIM);
   ParGridFunction euler_state(&euler_space);
   Array<int> scalar_dofs;
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      euler_space.GetElementDofs(e, scalar_dofs);
      for (int i = 0; i < scalar_dofs.Size(); ++i)
      {
         euler_state(euler_space.DofToVDof(scalar_dofs[i], 0)) = 1.0;
         euler_state(euler_space.DofToVDof(scalar_dofs[i], 1)) = 0.2;
         euler_state(euler_space.DofToVDof(scalar_dofs[i], 2)) = 2.52;
      }
   }
   Vector before = euler_state;
   EulerPositivityLimiter positivity(&euler_space, 1.4);
   const EulerPositivityDiagnostics diagnostics = positivity.Apply(euler_state);
   Require(diagnostics.limited_elements == 0 &&
           diagnostics.inadmissible_means == 0,
           "parallel positivity limiter changed an admissible constant state");
   euler_state -= before;
   Require(euler_state.Normlinf() < 2e-14,
           "parallel positivity limiter is inconsistent across shared faces");

   auto velocity = PositiveVelocity();
   KXRCFIndicator indicator(&space, &velocity);
   Array<bool> active;
   Vector indicator_values;
   indicator.Compute(state, active, &indicator_values);
   int local_active = 0;
   for (int e = 0; e < active.Size(); ++e)
   {
      local_active += active[e] ? 1 : 0;
   }
   int global_active = 0;
   MPI_Allreduce(&local_active, &global_active, 1, MPI_INT, MPI_SUM,
                 mesh.GetComm());
   if (global_active != 1)
   {
      std::cerr << "rank " << Mpi::WorldRank() << " indicators:";
      for (int e = 0; e < indicator_values.Size(); ++e)
      {
         std::cerr << ' ' << indicator_values(e);
      }
      std::cerr << '\n';
   }
   Require(global_active == 1,
           "KXRCF did not see the partition-boundary discontinuity exactly once");

   // Lobatto projection puts the interface value on both traces, but leaves a
   // derivative jump in the left cell. This specifically exercises OFDG's
   // shared-face derivative halo rather than only its solution-jump term.
   FunctionCoefficient projected_jump([](const Vector &x)
   {
      return x(0) < 0.5 ? 1.0 : 2.0;
   });
   state.ProjectCoefficient(projected_jump);

   OFDG ofdg(&space, BasisType::GaussLobatto, velocity);
   Vector stabilization, decay;
   ofdg.ComputeStabilization(state, stabilization);
   ofdg.CompDecay(state, decay, 0.02);

   real_t local_norm_squared = stabilization * stabilization;
   real_t global_norm_squared = 0.0;
   MPI_Allreduce(&local_norm_squared, &global_norm_squared, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, mesh.GetComm());
   Require(global_norm_squared > 1e-12,
           "OFDG missed the partition-boundary discontinuity");

   auto face_physics = std::make_shared<AdvectionFacePhysics>(&velocity);
   OEDG2024 oedg(&space, BasisType::GaussLobatto, face_physics);
   Vector oedg_decay;
   oedg.CompDecay(state, oedg_decay, 0.02);
   ParGridFunction oedg_decayed(&space, oedg_decay);
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      Require(std::abs(ElementIntegral(state, e) -
                       ElementIntegral(oedg_decayed, e)) < 2e-12,
              "parallel OEDG decay changed an element mean");
   }

   ParGridFunction decayed(&space, decay);
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      Require(std::abs(ElementIntegral(state, e) -
                       ElementIntegral(decayed, e)) < 2e-12,
              "parallel OFDG decay changed an element mean");
   }

   state = 3.0;
   indicator.Compute(state, active);
   for (int e = 0; e < active.Size(); ++e)
   {
      Require(!active[e], "constant parallel state was marked troubled");
   }
   ofdg.CompDecay(state, decay, 0.02);
   decay -= state;
   Require(decay.Normlinf() < 2e-12,
           "constant parallel state was changed by OFDG");

}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   try
   {
      Require(Mpi::WorldSize() == 2,
              "parallel consistency test must run with two MPI ranks");
      TestMixedMeshConsistency();
      RunParallelChecks();
   }
   catch (const std::exception &error)
   {
      std::cerr << "Parallel consistency tests failed on rank "
                << Mpi::WorldRank() << ": " << error.what() << '\n';
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
   }

   if (Mpi::Root())
   {
      std::cout << "Two-rank OFDG/OEDG/KXRCF/positivity tests passed.\n";
   }
   return 0;
}
