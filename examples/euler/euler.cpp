// Derived from MFEM Example 18 (BSD 3-Clause). See
// ../../THIRD_PARTY_NOTICES.md for copyright and licence terms.
//                                MFEM Example 18
//
// Compile with: make ex18
//
// Sample runs:
//
//       ex18 -p 1 -r 2 -o 1 -s 3
//       ex18 -p 1 -r 1 -o 3 -s 4
//       ex18 -p 1 -r 0 -o 5 -s 6
//       ex18 -p 2 -r 1 -o 1 -s 3 -mf
//       ex18 -p 2 -r 0 -o 3 -s 3 -mf
//
// Description:  This example code solves the compressible Euler system of
//               equations, a model nonlinear hyperbolic PDE, with a
//               discontinuous Galerkin (DG) formulation.
//
//                (u_t, v)_T - (F(u), ∇ v)_T + <F̂(u,n), [[v]]>_F = 0
//
//               where (⋅,⋅)_T is volume integration, and <⋅,⋅>_F is face
//               integration, F is the Euler flux function, and F̂ is the
//               numerical flux.
//
//               Specifically, it solves for an exact solution of the equations
//               whereby a vortex is transported by a uniform flow. Since all
//               boundaries are periodic here, the method's accuracy can be
//               assessed by measuring the difference between the solution and
//               the initial condition at a later time when the vortex returns
//               to its initial location.
//
//               Note that as the order of the spatial discretization increases,
//               the timestep must become smaller. This example currently uses a
//               simple estimate derived by Cockburn and Shu for the 1D RKDG
//               method. An additional factor can be tuned by passing the --cfl
//               (or -c shorter) flag.
//
//               The example demonstrates usage of DGHyperbolicConservationLaws
//               that wraps NonlinearFormIntegrators containing element and face
//               integration schemes. In this case the system also involves an
//               external approximate Riemann solver for the DG interface flux.
//               By default, weak-divergence is pre-assembled in element-wise
//               manner, which corresponds to (I_h(F(u_h)), ∇ v). This yields
//               better performance and similar accuracy for the included test
//               problems. This can be turned off and use nonlinear assembly
//               similar to matrix-free assembly when -mf flag is provided.
//               It also demonstrates how to use GLVis for in-situ visualization
//               of vector grid function and how to set top-view.
//
//               We recommend viewing examples 9, 14 and 17 before viewing this
//               example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include "euler.hpp"
#include "../../src/profile_output.hpp"

#include "../../src/euler_positivity.hpp"
#include "../../src/conservation.hpp"
#include "../../src/experiment_rk.hpp"
#include "../../src/glvis_output.hpp"
#include "../../src/study_filter.hpp"

using namespace std;
using namespace mfem;
using namespace ofdg;

namespace
{

void WriteEulerSamples(const ParGridFunction &solution,
                       const ParFiniteElementSpace &scalar_space,
                       const ParMesh &mesh, real_t gamma,
                       const string &prefix, int rank)
{
   if (prefix.empty()) { return; }
   ostringstream name;
   name << prefix << ".rank" << setfill('0') << setw(6) << rank << ".csv";
   ofstream output(name.str());
   MFEM_VERIFY(output, "Unable to open Euler profile output " << name.str());
   output << "rank,element,sample,x,y,density,velocity_x,velocity_y,pressure\n";
   output << setprecision(17);
   const int dim = mesh.Dimension();
   const int scalar_dofs = scalar_space.GetNDofs();
   for (int element = 0; element < mesh.GetNE(); ++element)
   {
      const FiniteElement *fe = scalar_space.GetFE(element);
      ElementTransformation *transformation =
         scalar_space.GetMesh()->GetElementTransformation(element);
      const IntegrationRule &nodes = fe->GetNodes();
      auto evaluate = [&](const IntegrationPoint &ip, Vector &conservative)
      {
         for (int component = 0; component < dim + 2; ++component)
         {
            GridFunction field(const_cast<ParFiniteElementSpace *>(&scalar_space),
                               const_cast<real_t *>(solution.GetData()) +
                                  component * scalar_dofs);
            conservative(component) = field.GetValue(element, ip);
         }
      };
      auto write_state = [&](const IntegrationPoint &ip, const char *sample,
                             const Vector &conservative) {
         Vector point;
         transformation->Transform(ip, point);
         const real_t density = conservative(0);
         real_t momentum_squared = 0.0;
         for (int d = 0; d < dim; ++d)
         {
            momentum_squared += conservative(1 + d) * conservative(1 + d);
         }
         const real_t pressure = (gamma - 1.0) *
            (conservative(dim + 1) - 0.5 * momentum_squared / density);
         output << rank << ',' << element << ',' << sample << ','
                << point(0) << ',' << (dim > 1 ? point(1) : 0.0) << ','
                << density << ',' << conservative(1) / density << ','
                << (dim > 1 ? conservative(2) / density : 0.0) << ','
                << pressure << '\n';
      };
      auto write_point = [&](const IntegrationPoint &ip, const char *sample) {
         Vector conservative(dim + 2);
         evaluate(ip, conservative);
         write_state(ip, sample, conservative);
      };
      if (dim == 1) {
         for (int q = 0; q <= 20; ++q) {
            IntegrationPoint ip;
            ip.Set1w(q / 20.0, 1.0);
            write_point(ip, "polynomial");
         }
      } else {
         for (int q = 0; q < nodes.GetNPoints(); ++q) {
            write_point(nodes.IntPoint(q), "polynomial");
         }
      }
      const auto &center = Geometries.GetCenter(fe->GetGeomType());
      write_point(center, "center");
      const auto average = PhysicalCellAverage(*transformation, fe->GetOrder(), dim + 2, evaluate);
      write_state(center, "cell_average", average);
   }
}

bool CheckPhysicalState(const GridFunction &solution, int dim, real_t gamma,
                        real_t time, const char *stage,
                        real_t *reported_density = nullptr,
                        real_t *reported_pressure = nullptr)
{
   const int scalar_dofs = solution.FESpace()->GetNDofs();
   real_t minimum_density = std::numeric_limits<real_t>::infinity();
   real_t minimum_pressure = std::numeric_limits<real_t>::infinity();

   for (int i = 0; i < scalar_dofs; ++i)
   {
      const real_t density = solution(i);
      real_t momentum_squared = 0.0;
      for (int d = 0; d < dim; ++d)
      {
         const real_t momentum = solution((1 + d) * scalar_dofs + i);
         momentum_squared += momentum * momentum;
      }
      const real_t energy = solution((dim + 1) * scalar_dofs + i);
      const real_t pressure = (gamma - 1.0) *
                              (energy - 0.5 * momentum_squared / density);
      minimum_density = std::min(minimum_density, density);
      minimum_pressure = std::min(minimum_pressure, pressure);
   }

#ifdef MFEM_USE_MPI
   if (const auto *parallel_space = dynamic_cast<const ParFiniteElementSpace *>(
          solution.FESpace()))
   {
      MPI_Allreduce(MPI_IN_PLACE, &minimum_density, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    parallel_space->GetComm());
      MPI_Allreduce(MPI_IN_PLACE, &minimum_pressure, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    parallel_space->GetComm());
   }
#endif

   if (reported_density) { *reported_density = minimum_density; }
   if (reported_pressure) { *reported_pressure = minimum_pressure; }

   if (minimum_density > 0.0 && minimum_pressure > 0.0 &&
       std::isfinite(minimum_density) && std::isfinite(minimum_pressure))
   {
      return true;
   }

   cerr << "Nonphysical Euler state " << stage << " at t=" << time
        << ": min(rho)=" << minimum_density
        << ", min(p)=" << minimum_pressure << '\n';
   return false;
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int rank = Mpi::WorldRank();
   const int process_count = Mpi::WorldSize();

   // 1. Parse command-line options.
   int problem = 1;
   const real_t specific_heat_ratio = 1.4;
   const real_t gas_constant = 1.0;

   string mesh_file = "";
   int IntOrderOffset = 1;
   int ref_levels = 1;
   int elements_override = 0;
   int order = 3;
   int ode_solver_type = 4;
   real_t t_final = -1.0;
   real_t dt = -0.01;
   real_t cfl = 0.3;
   bool visualization = false;
   bool preassembleWeakDiv = true;
   bool save_output = false;
   int vis_steps = 50;

   bool stabilization = true;
   bool use_kxrcf = true;
   real_t kxrcf_threshold = 1.0;
   string method_name;
   string cadence_name = "auto";
   bool use_positivity = false;
   bool outflow_boundary = false;
   real_t state_scale = 1.0;
   int maximum_step_retries = 8;
   string profile_prefix;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use. If not provided, then a periodic square"
                  " mesh will be used.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem setup to use. See EulerInitialCondition().");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&elements_override, "-n", "--elements",
                  "Override the default elements per coordinate direction.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  ODESolver::ExplicitTypes.c_str());
   args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step. Positive number skips CFL timestep calculation.");
   args.AddOption(&cfl, "-c", "--cfl-number",
                  "CFL number for timestep calculation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&preassembleWeakDiv, "-ea", "--element-assembly-divergence",
                  "-mf", "--matrix-free-divergence",
                  "Weak divergence assembly level\n"
                  "    ea - Element assembly with interpolated F\n"
                  "    mf - Nonlinear assembly in matrix-free manner");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&save_output, "-save", "--save-output", "-no-save",
                  "--no-save-output", "Write rank-local MFEM mesh/solution files.");
   args.AddOption(&stabilization, "-stab", "--stabilization", "-no-stab",
                  "--no-stabilization",
                  "Enable or disable stabilization.");
   args.AddOption(&use_kxrcf, "-kxrcf", "--use-kxrcf", "-no-kxrcf",
                  "--no-use-kxrcf",
                  "Restrict OFDG to cells selected by KXRCF.");
   args.AddOption(&kxrcf_threshold, "-kt", "--kxrcf-threshold",
                  "KXRCF troubled-cell threshold.");
   args.AddOption(&method_name, "-method", "--method",
                  "Study method: dg, ofdg, ofdg-kxrcf, or oedg.");
   args.AddOption(&cadence_name, "-cadence", "--filter-cadence",
                  "Filter cadence: auto, step, or stage.");
   args.AddOption(&use_positivity, "-positivity", "--positivity-limiter",
                  "-no-positivity", "--no-positivity-limiter",
                  "Apply the common conservative Euler positivity limiter.");
   args.AddOption(&state_scale, "-scale", "--state-scale",
                  "Multiply the complete initial conservative state.");
   args.AddOption(&maximum_step_retries, "-retries", "--step-retries",
                  "Maximum rejected-step halvings for inadmissible means.");
   args.AddOption(&outflow_boundary, "-outflow", "--outflow-boundary",
                  "-fixed-boundary", "--fixed-boundary", "Use transmissive exterior states (non-wall problems).");
   args.AddOption(&profile_prefix, "-profile", "--profile-prefix",
                  "Write rank-local CSV samples using this file prefix.");

   args.ParseCheck();

   if (t_final < 0.0)
   {
      t_final = problem == 5 ? 0.038
                : (problem == 6 ? 0.25
                   : (problem == 7 ? 1.3
                      : (problem == 8 ? 1.8 :
                         (problem == 3 ? 2.0 : 1.1))));
   }
   const bool use_cfl_time_step = dt <= 0.0;

   // 2. Read the mesh from the given mesh file. When the user does not provide
   //    mesh file, use the default mesh file for the problem.
   MFEM_VERIFY(elements_override >= 0,
               "The element override must be nonnegative.");
   MFEM_VERIFY(mesh_file.empty() || elements_override == 0,
               "Do not combine a mesh file with an element override.");
   Mesh mesh = mesh_file.empty() ? EulerMesh(problem, elements_override)
                                 : Mesh(mesh_file);
   const int dim = mesh.Dimension();
   const int num_equations = dim + 2;

   // Refine the mesh to increase the resolution. In this example we do
   // 'ref_levels' of uniform refinement, where 'ref_levels' is a command-line
   // parameter.
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh.UniformRefinement();
   }
   ParMesh parallel_mesh(MPI_COMM_WORLD, mesh);

   // 3. Define the ODE solver used for time integration. Several explicit
   //    Runge-Kutta methods are available.
   unique_ptr<ODESolver> ode_solver = ODESolver::SelectExplicit(ode_solver_type);

   // 4. Define the discontinuous DG finite element space of the given
   //    polynomial order on the refined mesh.
   DG_FECollection fec(order, dim);
   // Finite element space for a scalar (thermodynamic quantity)
   ParFiniteElementSpace fes(&parallel_mesh, &fec);
   // Finite element space for a mesh-dim vector quantity (momentum)
   ParFiniteElementSpace dfes(&parallel_mesh, &fec, dim, Ordering::byNODES);
   // Finite element space for all variables together (total thermodynamic state)
   ParFiniteElementSpace vfes(&parallel_mesh, &fec, num_equations,
                              Ordering::byNODES);

   // if (ode_solver_type != 4)
   // {
   //    cout << "Error: Only RK4 is supported for this example." << endl;
   //    return 1;
   // }
   // // Note: GaussLegendre is default in fec
   // auto ode_solver_stab_nd =
   //       std::make_unique<OEDG_RK4Solver<OFDG>>(fes, BasisType::GaussLegendre, velocity_function(problem));

   // This example depends on this ordering of the space.
   MFEM_ASSERT(fes.GetOrdering() == Ordering::byNODES, "");

   const HYPRE_BigInt global_unknowns = vfes.GlobalTrueVSize();
   if (Mpi::Root())
   {
      cout << "Number of unknowns: " << global_unknowns << endl;
   }

   // 5. Define the initial conditions, save the corresponding mesh and grid
   //    functions to files. These can be opened with GLVis using:
   //    "glvis -m euler-mesh.mesh -g euler-1-init.gf" (for x-momentum).

   // Initialize the state.
   VectorFunctionCoefficient unscaled_u0 = EulerInitialCondition(
      problem, specific_heat_ratio, gas_constant);
   ScalarVectorProductCoefficient u0(state_scale, unscaled_u0);
   ParGridFunction sol(&vfes);
   sol.ProjectCoefficient(u0);
   Vector initial_integrals(num_equations);
   for (int component = 0; component < num_equations; ++component)
   {
      ParGridFunction component_view(
         &fes, sol.GetData() + component * fes.GetNDofs());
      initial_integrals(component) = GlobalScalarIntegral(component_view);
   }
   ParGridFunction mom(&dfes, sol.GetData() + fes.GetNDofs());
   // Output the initial solution.
   if (save_output)
   {
      ostringstream mesh_name;
      mesh_name << "euler-mesh." << setfill('0') << setw(6) << rank;
      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(precision);
      mesh_ofs << parallel_mesh;

      for (int k = 0; k < num_equations; k++)
      {
         ParGridFunction uk(&fes, sol.GetData() + k * fes.GetNDofs());
         ostringstream sol_name;
         sol_name << "euler-" << k << "-init." << setfill('0')
                  << setw(6) << rank;
         ofstream sol_ofs(sol_name.str().c_str());
         sol_ofs.precision(precision);
         sol_ofs << uk;
      }
   }

   // 6. Set up the nonlinear form with euler flux and numerical flux
   EulerFlux flux(dim, specific_heat_ratio);
   RusanovFlux numericalFlux(flux);
   std::unique_ptr<EulerBoundaryIntegrator> boundary_integrator;
   if (problem >= 5)
   {
      boundary_integrator =
         std::make_unique<EulerBoundaryIntegrator>(
            numericalFlux, &u0, problem == 5, IntOrderOffset, outflow_boundary);
   }
   DGHyperbolicConservationLaws euler(
      vfes, std::unique_ptr<HyperbolicFormIntegrator>(
         new HyperbolicFormIntegrator(numericalFlux, IntOrderOffset)),
      preassembleWeakDiv, std::move(boundary_integrator));

   // 7. Visualize momentum with its magnitude.
   GLVisOutput glvis(parallel_mesh, visualization, precision);
   if (glvis.Enabled())
   {
      glvis.Send(parallel_mesh, mom, "momentum, t = 0",
                 "view 0 0\nkeys jlm\npause\n");
      if (Mpi::Root())
      {
         cout << "GLVis visualization paused."
              << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   // 8. Time integration

   // When dt is not specified, use CFL condition.
   // Compute h_min and initial maximum characteristic speed
   real_t hmin = infinity();
   if (use_cfl_time_step)
   {
      for (int i = 0; i < parallel_mesh.GetNE(); i++)
      {
         hmin = min(parallel_mesh.GetElementSize(i, 1), hmin);
      }
      MPI_Allreduce(MPI_IN_PLACE, &hmin, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    parallel_mesh.GetComm());
      // Find a safe dt, using a temporary vector. Calling Mult() computes the
      // maximum char speed at all quadrature points on all faces (and all
      // elements with -mf).
      Vector z(sol.Size());
      euler.Mult(sol, z);

      real_t max_char_speed = euler.GetMaxCharSpeed();
      MPI_Allreduce(MPI_IN_PLACE, &max_char_speed, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    parallel_mesh.GetComm());
      dt = cfl * hmin / max_char_speed / (2 * order + 1);
   }

   // Start the timer.
   tic_toc.Clear();
   tic_toc.Start();

   // Init time integration
   real_t t = 0.0;
   euler.SetTime(t);
   ode_solver->Init(euler);

   StudyMethod method = method_name.empty()
                           ? (!stabilization ? StudyMethod::DG
                              : (use_kxrcf ? StudyMethod::OFDGKXRCF
                                           : StudyMethod::OFDG))
                           : ParseStudyMethod(method_name);
   const FilterCadence requested_cadence = ParseFilterCadence(cadence_name);
   auto euler_face_physics =
      std::make_shared<EulerFacePhysics>(dim, specific_heat_ratio);
   StudyFilter filter(&vfes, BasisType::GaussLegendre, euler_face_physics,
                      method, requested_cadence, kxrcf_threshold);
   unique_ptr<EulerPositivityLimiter> positivity_limiter;
   if (use_positivity)
   {
      positivity_limiter =
         make_unique<EulerPositivityLimiter>(&vfes, specific_heat_ratio);
   }
   long limited_elements = 0;
   long rejected_steps = 0;
   long positivity_applications = 0;

   auto postprocess = [&](Vector &candidate, real_t step_size,
                          bool final_stage)
   {
      if (use_positivity)
      {
         const EulerPositivityDiagnostics before =
            positivity_limiter->Apply(candidate);
         limited_elements += before.limited_elements;
         ++positivity_applications;
         if (before.inadmissible_means > 0) { return false; }
      }
      else
      {
         GridFunction candidate_view(&vfes, candidate.GetData());
         if (!CheckPhysicalState(candidate_view, dim, specific_heat_ratio,
                                 t + step_size, "before filtering"))
         {
            return false;
         }
      }

      filter.Apply(candidate, step_size, final_stage);

      if (use_positivity && filter.Enabled() &&
          (filter.Cadence() == FilterCadence::Stage || final_stage))
      {
         const EulerPositivityDiagnostics after =
            positivity_limiter->Apply(candidate);
         limited_elements += after.limited_elements;
         ++positivity_applications;
         if (after.inadmissible_means > 0) { return false; }
      }
      return true;
   };

   const bool discontinuous_problem = problem >= 5;
   ExperimentRungeKutta experiment_solver(
      discontinuous_problem ? ExperimentRungeKutta::Scheme::SSPRK3
                            : ExperimentRungeKutta::Scheme::ClassicalRK4,
      postprocess);
   experiment_solver.Init(euler);
   const bool use_experiment_solver =
      filter.Cadence() == FilterCadence::Stage || use_positivity;

   // Integrate in time.
   bool done = false;
   int ti = 0;
   for (; !done;)
   {
      real_t dt_real = min(dt, t_final - t);
      if (use_experiment_solver)
      {
         bool accepted = false;
         const int retry_limit = use_positivity ? maximum_step_retries : 0;
         for (int retry = 0; retry <= retry_limit; ++retry)
         {
            if (experiment_solver.TryStep(sol, t, dt_real))
            {
               accepted = true;
               break;
            }
            ++rejected_steps;
            dt_real *= 0.5;
         }
         if (!accepted)
         {
            cerr << "Euler step remained inadmissible after "
                 << retry_limit << " retries at t=" << t << '\n';
            return 2;
         }
      }
      else
      {
         ode_solver->Step(sol, t, dt_real);
         if (filter.Enabled() &&
             !CheckPhysicalState(sol, dim, specific_heat_ratio, t,
                                 "before filtering"))
         {
            return 2;
         }
         filter.Apply(sol, dt_real, true);
      }

      if (!CheckPhysicalState(sol, dim, specific_heat_ratio, t,
                              "after accepted step"))
      {
         return 2;
      }

      if (use_cfl_time_step) // update time step size with CFL
      {
         real_t max_char_speed = euler.GetMaxCharSpeed();
         MPI_Allreduce(MPI_IN_PLACE, &max_char_speed, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX,
                       parallel_mesh.GetComm());
         dt = cfl * hmin / max_char_speed / (2 * order + 1);
      }
      ++ti;

      done = (t >= t_final - 1e-8 * dt);
      if (done || ti % vis_steps == 0)
      {
         if (Mpi::Root())
         {
            cout << "time step: " << ti << ", time: " << t << endl;
         }
         if (glvis.Enabled())
         {
            ostringstream title;
            title << "momentum, t = " << t;
            glvis.Send(parallel_mesh, mom, title.str());
         }
      }
   }

   tic_toc.Stop();
   real_t runtime = tic_toc.RealTime();
   MPI_Allreduce(MPI_IN_PLACE, &runtime, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MAX, parallel_mesh.GetComm());
   if (Mpi::Root()) { cout << " done, " << runtime << "s." << endl; }

   // 9. Save the final solution. This output can be viewed later using GLVis:
   //    "glvis -m euler-mesh-final.mesh -g euler-1-final.gf" (for x-momentum).
   if (save_output)
   {
      ostringstream mesh_name;
      mesh_name << "euler-mesh-final." << setfill('0') << setw(6) << rank;
      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(precision);
      mesh_ofs << parallel_mesh;

      for (int k = 0; k < num_equations; k++)
      {
         ParGridFunction uk(&fes, sol.GetData() + k * fes.GetNDofs());
         ostringstream sol_name;
         sol_name << "euler-" << k << "-final." << setfill('0')
                  << setw(6) << rank;
         ofstream sol_ofs(sol_name.str().c_str());
         sol_ofs.precision(precision);
         sol_ofs << uk;
      }
   }

   WriteEulerSamples(sol, fes, parallel_mesh, specific_heat_ratio,
                     profile_prefix, rank);

   // 10. Compute exact smooth-solution errors where available.
   real_t l1_error = std::numeric_limits<real_t>::quiet_NaN();
   real_t l2_error = std::numeric_limits<real_t>::quiet_NaN();
   real_t linf_error = std::numeric_limits<real_t>::quiet_NaN();
   real_t state_l2_error = std::numeric_limits<real_t>::quiet_NaN();
   if (problem >= 1 && problem <= 4)
   {
      VectorFunctionCoefficient exact = EulerVortexExactCondition(
         problem, t_final, specific_heat_ratio, gas_constant);
      FunctionCoefficient density_exact = EulerExactDensityCondition(
         problem, t_final, specific_heat_ratio, gas_constant);
      ParGridFunction density(&fes, sol.GetData());
      l1_error = density.ComputeL1Error(density_exact);
      l2_error = density.ComputeL2Error(density_exact);
      linf_error = density.ComputeMaxError(density_exact);
      state_l2_error = sol.ComputeLpError(2, exact);
   }
   real_t minimum_density = 0.0;
   real_t minimum_pressure = 0.0;
   CheckPhysicalState(sol, dim, specific_heat_ratio, t, "final state",
                      &minimum_density, &minimum_pressure);
   real_t conservation_drift = 0.0;
   for (int component = 0; component < num_equations; ++component)
   {
      ParGridFunction component_view(
         &fes, sol.GetData() + component * fes.GetNDofs());
      conservation_drift = std::max(
         conservation_drift,
         RelativeConservationDrift(initial_integrals(component),
                                   GlobalScalarIntegral(component_view)));
   }
   const StudyFilterStatistics filter_stats = filter.Statistics();
   if (Mpi::Root())
   {
      cout.precision(16);
      cout << "method=" << StudyMethodName(method)
           << " cadence=" << FilterCadenceName(filter.Cadence())
           << " problem=" << problem
           << " order=" << order
           << " dofs=" << global_unknowns
           << " ranks=" << process_count
           << " steps=" << ti
           << " time=" << t
           << " l1_error=" << l1_error
           << " l2_error=" << l2_error
           << " linf_error=" << linf_error
           << " state_l2_error=" << state_l2_error
           << " conservation_drift=" << conservation_drift
           << " min_density=" << minimum_density
           << " min_pressure=" << minimum_pressure
           << " filter_applications=" << filter_stats.applications
           << " active_elements=" << filter_stats.active_elements
           << " positivity_applications=" << positivity_applications
           << " limited_elements=" << limited_elements
           << " rejected_steps=" << rejected_steps
           << " runtime_seconds=" << runtime << '\n';
      cout << "Solution error: " << l2_error << endl;
   }

   return 0;
}
