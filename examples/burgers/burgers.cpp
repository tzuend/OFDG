// Burgers accuracy and post-shock robustness examples using MFEM's
// hyperbolic DG integrators. OFDG is applied once after each complete MFEM
// Runge--Kutta step; the custom stage-by-stage RK4 remains a separate test.

#include "mfem.hpp"

#include "../euler/euler.hpp"
#include "../../src/experiment_rk.hpp"
#include "../../src/conservation.hpp"
#include "../../src/glvis_output.hpp"
#include "../../src/profile_output.hpp"
#include "../../src/study_filter.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

using namespace mfem;
using namespace std;

namespace
{

real_t InitialValue(const Vector &x)
{
   return 0.5 + std::sin(x(0));
}

real_t SmoothExactValue(const Vector &x, real_t time)
{
   // Before the first shock (t=1), characteristics give
   // u(x,t)=u_0(x-u t). Newton's method is well conditioned for t<1.
   real_t value = InitialValue(x);
   for (int iteration = 0; iteration < 30; ++iteration)
   {
      const real_t foot = x(0) - time * value;
      const real_t residual = value - 0.5 - std::sin(foot);
      const real_t derivative = 1.0 + time * std::cos(foot);
      const real_t update = residual / derivative;
      value -= update;
      if (std::abs(update) < 1e-14) { break; }
   }
   return value;
}

Mesh PeriodicMesh(int elements)
{
   Mesh mesh = Mesh::MakeCartesian1D(elements, 2.0 * M_PI);
   std::vector<int> vertex_map(mesh.GetNV());
   for (int v = 0; v < mesh.GetNV(); ++v) { vertex_map[v] = v; }
   vertex_map.back() = 0;
   return Mesh::MakePeriodic(mesh, vertex_map);
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   int problem = 1;
   int elements = 80;
   int order = 3;
   int ode_solver_type = 4;
   real_t final_time = -1.0;
   real_t dt = -1.0;
   real_t cfl = 0.15;
   bool stabilization = true;
   bool use_kxrcf = true;
   real_t kxrcf_threshold = 1.0;
   string method_name;
   string cadence_name = "auto";
   string profile_prefix;
   bool visualization = false;
   int vis_steps = 50;

   OptionsParser args(argc, argv);
   args.AddOption(&problem, "-p", "--problem",
                  "1: smooth accuracy (t<1), 2: post-shock robustness.");
   args.AddOption(&elements, "-n", "--elements", "Number of mesh elements.");
   args.AddOption(&order, "-o", "--order", "DG polynomial degree.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  ODESolver::ExplicitTypes.c_str());
   args.AddOption(&final_time, "-tf", "--t-final", "Final time.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Fixed time step; a non-positive value uses the CFL rule.");
   args.AddOption(&cfl, "-c", "--cfl-number", "CFL number.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable live GLVis visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&stabilization, "-stab", "--stabilization", "-no-stab",
                  "--no-stabilization", "Enable OFDG stabilization.");
   args.AddOption(&use_kxrcf, "-kxrcf", "--use-kxrcf", "-no-kxrcf",
                  "--no-use-kxrcf", "Restrict OFDG using KXRCF.");
   args.AddOption(&kxrcf_threshold, "-kt", "--kxrcf-threshold",
                  "KXRCF troubled-cell threshold.");
   args.AddOption(&method_name, "-method", "--method",
                  "Study method: dg, ofdg, ofdg-kxrcf, or oedg.");
   args.AddOption(&cadence_name, "-cadence", "--filter-cadence",
                  "Filter cadence: auto, step, or stage.");
   args.AddOption(&profile_prefix, "-profile", "--profile-prefix",
                  "Write rank-local CSV samples using this file prefix.");
   args.ParseCheck();

   MFEM_VERIFY(problem == 1 || problem == 2,
               "Burgers problem must be 1 or 2.");
   MFEM_VERIFY(elements > 0 && order >= 1,
               "Burgers requires positive elements and order >= 1.");
   MFEM_VERIFY(vis_steps > 0,
               "The GLVis visualization interval must be positive.");

   if (final_time < 0.0) { final_time = problem == 1 ? 0.6 : 1.5; }

   Mesh mesh = PeriodicMesh(elements);
   ParMesh parallel_mesh(MPI_COMM_WORLD, mesh);
   DG_FECollection collection(order, 1);
   ParFiniteElementSpace space(&parallel_mesh, &collection);
   const HYPRE_BigInt global_unknowns = space.GlobalTrueVSize();
   ParGridFunction solution(&space);
   FunctionCoefficient initial(InitialValue);
   solution.ProjectCoefficient(initial);
   const real_t initial_integral = GlobalScalarIntegral(solution);

   BurgersFlux flux(1);
   RusanovFlux numerical_flux(flux);
   DGHyperbolicConservationLaws evolution(
      space, std::make_unique<HyperbolicFormIntegrator>(numerical_flux, 1),
      true);

   StudyMethod method = method_name.empty()
                           ? (!stabilization ? StudyMethod::DG
                              : (use_kxrcf ? StudyMethod::OFDGKXRCF
                                           : StudyMethod::OFDG))
                           : ParseStudyMethod(method_name);
   const FilterCadence requested_cadence = ParseFilterCadence(cadence_name);
   auto face_physics = std::make_shared<BurgersFacePhysics>();
   StudyFilter filter(&space, BasisType::GaussLegendre, face_physics,
                      method, requested_cadence, kxrcf_threshold);

   GLVisOutput glvis(parallel_mesh, visualization);
   if (glvis.Enabled())
   {
      glvis.Send(parallel_mesh, solution, "Burgers, t = 0", "pause\n");
      if (Mpi::Root())
      {
         cout << "GLVis visualization paused. Press space in the GLVis "
              << "window to resume it.\n";
      }
   }

   std::unique_ptr<ODESolver> ode_solver =
      ODESolver::SelectExplicit(ode_solver_type);
   real_t time = 0.0;
   evolution.SetTime(time);
   ode_solver->Init(evolution);
   ExperimentRungeKutta experiment_solver(
      ExperimentRungeKutta::Scheme::ClassicalRK4,
      [&filter](Vector &candidate, real_t step, bool final_stage)
      {
         return filter.Apply(candidate, step, final_stage);
      });
   experiment_solver.Init(evolution);

   const real_t h = 2.0 * M_PI / elements;
   const bool use_cfl_time_step = dt <= 0.0;
   Vector residual(solution.Size());
   evolution.Mult(solution, residual);
   real_t maximum_speed = evolution.GetMaxCharSpeed();
   MPI_Allreduce(MPI_IN_PLACE, &maximum_speed, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX,
                 parallel_mesh.GetComm());
   if (use_cfl_time_step)
   {
      dt = cfl * h / ((2 * order + 1) * maximum_speed);
   }

   tic_toc.Clear();
   tic_toc.Start();
   int steps = 0;
   while (time < final_time)
   {
      real_t step_size = std::min(dt, final_time - time);
      if (filter.Cadence() == FilterCadence::Stage)
      {
         MFEM_VERIFY(experiment_solver.TryStep(solution, time, step_size),
                     "Burgers experiment step was unexpectedly rejected.");
      }
      else
      {
         ode_solver->Step(solution, time, step_size);
         filter.Apply(solution, step_size, true);
      }

      int finite_state = std::isfinite(solution.Normlinf()) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &finite_state, 1, MPI_INT, MPI_MIN,
                    parallel_mesh.GetComm());
      if (!finite_state)
      {
         cerr << "Burgers solution became non-finite at t=" << time << '\n';
         return 2;
      }

      if (use_cfl_time_step)
      {
         maximum_speed = evolution.GetMaxCharSpeed();
         MPI_Allreduce(MPI_IN_PLACE, &maximum_speed, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX,
                       parallel_mesh.GetComm());
         dt = cfl * h / ((2 * order + 1) * maximum_speed);
      }
      ++steps;
      if (glvis.Enabled() &&
          (steps % vis_steps == 0 || time >= final_time))
      {
         ostringstream title;
         title << "Burgers, t = " << time;
         glvis.Send(parallel_mesh, solution, title.str());
      }
   }
   tic_toc.Stop();

   real_t runtime = tic_toc.RealTime();
   MPI_Allreduce(MPI_IN_PLACE, &runtime, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MAX, parallel_mesh.GetComm());
   WriteScalarSamples(solution, profile_prefix, Mpi::WorldRank());

   real_t l1_error = std::numeric_limits<real_t>::quiet_NaN();
   real_t l2_error = std::numeric_limits<real_t>::quiet_NaN();
   real_t linf_error = std::numeric_limits<real_t>::quiet_NaN();
   if (final_time < 1.0)
   {
      FunctionCoefficient exact([final_time](const Vector &x)
      {
         return SmoothExactValue(x, final_time);
      });
      l1_error = solution.ComputeL1Error(exact);
      l2_error = solution.ComputeL2Error(exact);
      linf_error = solution.ComputeMaxError(exact);
   }

   real_t minimum = solution.Min();
   real_t maximum = solution.Max();
   MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MIN, parallel_mesh.GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MAX, parallel_mesh.GetComm());
   const real_t conservation_drift = RelativeConservationDrift(
      initial_integral, GlobalScalarIntegral(solution));
   const StudyFilterStatistics filter_stats = filter.Statistics();

   if (Mpi::Root())
   {
      cout.precision(16);
      cout << "method=" << StudyMethodName(method)
           << " cadence=" << FilterCadenceName(filter.Cadence())
           << " problem=" << problem << " elements=" << elements
           << " order=" << order << " steps=" << steps
           << " dofs=" << global_unknowns
           << " time=" << final_time
           << " l1_error=" << l1_error
           << " l2_error=" << l2_error
           << " linf_error=" << linf_error
           << " conservation_drift=" << conservation_drift
           << " min=" << minimum << " max=" << maximum
           << " filter_applications=" << filter_stats.applications
           << " active_elements=" << filter_stats.active_elements
           << " runtime_seconds=" << runtime << '\n';
   }
   return 0;
}
