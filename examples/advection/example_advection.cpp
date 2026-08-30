// Periodic linear-advection accuracy example. The spatial operator and time
// integration are standard MFEM components; OFDG is applied once after each
// complete Runge--Kutta step, optionally restricted by KXRCF.

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
#include <utility>
#include <vector>

using namespace mfem;
using namespace std;

namespace
{

Mesh PeriodicMesh(int dimension, int elements, bool triangular)
{
   if (dimension == 1)
   {
      Mesh mesh = Mesh::MakeCartesian1D(elements, 1.0);
      std::vector<int> vertex_map(mesh.GetNV());
      for (int v = 0; v < mesh.GetNV(); ++v) { vertex_map[v] = v; }
      vertex_map.back() = 0;
      return Mesh::MakePeriodic(mesh, vertex_map);
   }

   Mesh mesh = Mesh::MakeCartesian2D(elements, elements,
                                      triangular ? Element::TRIANGLE
                                                 : Element::QUADRILATERAL,
                                      true, 1.0, 1.0);
   Vector x_translation({1.0, 0.0});
   Vector y_translation({0.0, 1.0});
   std::vector<Vector> translations = {x_translation, y_translation};
   return Mesh::MakePeriodic(mesh,
      mesh.CreatePeriodicVertexMapping(translations));
}

real_t InitialValue(const Vector &x)
{
   if (x.Size() == 1)
   {
      const real_t wave = std::sin(2.0 * M_PI * x(0));
      return wave * wave;
   }
   const real_t wave = std::sin(M_PI * (x(0) + x(1)));
   return wave * wave;
}

void ProjectGaussRadau(ParGridFunction &solution, Coefficient &coefficient)
{
   const FiniteElementSpace *space = solution.FESpace();
   MFEM_VERIFY(space->GetMesh()->Dimension() == 1,
               "Gauss--Radau initialization is implemented for 1D only.");
   Array<int> dofs;
   Vector shape;
   for (int element = 0; element < space->GetNE(); ++element)
   {
      const FiniteElement *fe = space->GetFE(element);
      const int degree = fe->GetOrder();
      const int ndof = fe->GetDof();
      DenseMatrix system(ndof);
      Vector right_hand_side(ndof);
      system = 0.0;
      right_hand_side = 0.0;
      ElementTransformation *transformation =
         space->GetMesh()->GetElementTransformation(element);
      const IntegrationRule &rule = IntRules.Get(
         fe->GetGeomType(), 2 * degree + 6);
      shape.SetSize(ndof);
      for (int q = 0; q < rule.GetNPoints(); ++q)
      {
         const IntegrationPoint &point = rule.IntPoint(q);
         transformation->SetIntPoint(&point);
         fe->CalcShape(point, shape);
         const real_t weight = point.weight * transformation->Weight();
         for (int moment = 0; moment < degree; ++moment)
         {
            const real_t test = std::pow(point.x, moment);
            for (int i = 0; i < ndof; ++i)
            {
               system(moment, i) += weight * test * shape(i);
            }
            right_hand_side(moment) +=
               weight * test * coefficient.Eval(*transformation, point);
         }
      }
      IntegrationPoint downwind;
      downwind.x = 1.0;
      transformation->SetIntPoint(&downwind);
      fe->CalcShape(downwind, shape);
      for (int i = 0; i < ndof; ++i) { system(degree, i) = shape(i); }
      right_hand_side(degree) = coefficient.Eval(*transformation, downwind);

      DenseMatrixInverse inverse(system);
      Vector coefficients(ndof);
      inverse.Mult(right_hand_side, coefficients);
      space->GetElementDofs(element, dofs);
      solution.SetSubVector(dofs, coefficients);
   }
}

std::pair<real_t, real_t> SuperconvergenceErrors(
   const ParGridFunction &solution, Coefficient &exact)
{
   const FiniteElementSpace *space = solution.FESpace();
   if (space->GetMesh()->Dimension() != 1)
   {
      const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
      return {nan, nan};
   }
   real_t mean_error_squared = 0.0;
   real_t downwind_error_squared = 0.0;
   for (int element = 0; element < space->GetNE(); ++element)
   {
      const FiniteElement *fe = space->GetFE(element);
      ElementTransformation *transformation =
         space->GetMesh()->GetElementTransformation(element);
      const IntegrationRule &rule = IntRules.Get(
         fe->GetGeomType(), 2 * fe->GetOrder() + 6);
      real_t numerical_integral = 0.0;
      real_t exact_integral = 0.0;
      real_t volume = 0.0;
      for (int q = 0; q < rule.GetNPoints(); ++q)
      {
         const IntegrationPoint &point = rule.IntPoint(q);
         transformation->SetIntPoint(&point);
         const real_t weight = point.weight * transformation->Weight();
         volume += weight;
         numerical_integral += weight * solution.GetValue(element, point);
         exact_integral += weight * exact.Eval(*transformation, point);
      }
      const real_t mean_difference =
         (numerical_integral - exact_integral) / volume;
      mean_error_squared += volume * mean_difference * mean_difference;

      IntegrationPoint downwind;
      downwind.x = 1.0;
      transformation->SetIntPoint(&downwind);
      const real_t point_difference = solution.GetValue(element, downwind) -
                                      exact.Eval(*transformation, downwind);
      downwind_error_squared += volume * point_difference * point_difference;
   }
   real_t errors[2] = {mean_error_squared, downwind_error_squared};
   MPI_Allreduce(MPI_IN_PLACE, errors, 2,
                 MPITypeMap<real_t>::mpi_type, MPI_SUM,
                 dynamic_cast<const ParFiniteElementSpace *>(space)->GetComm());
   return {std::sqrt(errors[0]), std::sqrt(errors[1])};
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   int dimension = 1;
   int elements = 40;
   int order = 3;
   int ode_solver_type = 4;
   real_t final_time = 1.0;
   real_t dt = -1.0;
   real_t cfl = 0.2;
   bool stabilization = true;
   bool use_kxrcf = true;
   real_t kxrcf_threshold = 1.0;
   string method_name;
   string cadence_name = "auto";
   string profile_prefix;
   bool triangular = false;
   bool radau_initialization = false;
   bool visualization = false;
   int vis_steps = 50;

   OptionsParser args(argc, argv);
   args.AddOption(&dimension, "-d", "--dimension", "Spatial dimension: 1 or 2.");
   args.AddOption(&elements, "-n", "--elements",
                  "Elements per coordinate direction.");
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
   args.AddOption(&triangular, "-triangular", "--triangular-mesh",
                  "-quadrilateral", "--quadrilateral-mesh",
                  "Use triangles instead of quadrilaterals in 2D.");
   args.AddOption(&radau_initialization, "-radau-init",
                  "--gauss-radau-initialization", "-l2-init",
                  "--l2-initialization",
                  "Use the 1D right Gauss--Radau projection.");
   args.ParseCheck();

   MFEM_VERIFY(dimension == 1 || dimension == 2,
               "Advection dimension must be 1 or 2.");
   MFEM_VERIFY(elements > 0 && order >= 1,
               "Advection requires positive elements and order >= 1.");
   MFEM_VERIFY(vis_steps > 0,
               "The GLVis visualization interval must be positive.");

   MFEM_VERIFY(!radau_initialization || dimension == 1,
               "Gauss--Radau initialization is a 1D experiment.");
   Mesh serial_mesh = PeriodicMesh(dimension, elements, triangular);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   DG_FECollection collection(order, dimension);
   ParFiniteElementSpace space(&mesh, &collection);
   const HYPRE_BigInt global_unknowns = space.GlobalTrueVSize();
   ParGridFunction solution(&space);
   FunctionCoefficient initial(InitialValue);
   if (radau_initialization) { ProjectGaussRadau(solution, initial); }
   else { solution.ProjectCoefficient(initial); }
   const real_t initial_integral = GlobalScalarIntegral(solution);

   VectorFunctionCoefficient velocity(dimension,
      [dimension](const Vector &, Vector &value)
   {
      value.SetSize(dimension);
      value(0) = dimension == 1 ? 1.0 : 0.7;
      if (dimension == 2) { value(1) = 0.3; }
   });
   AdvectionFlux flux(velocity);
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
   auto face_physics = std::make_shared<AdvectionFacePhysics>(&velocity);
   StudyFilter filter(&space, BasisType::GaussLegendre, face_physics,
                      method, requested_cadence, kxrcf_threshold);

   GLVisOutput glvis(mesh, visualization);
   if (glvis.Enabled())
   {
      glvis.Send(mesh, solution, "advection, t = 0", "pause\n");
      if (Mpi::Root())
      {
         cout << "GLVis visualization paused. Press space in the GLVis "
              << "window to resume it.\n";
      }
   }

   std::unique_ptr<ODESolver> solver =
      ODESolver::SelectExplicit(ode_solver_type);
   real_t time = 0.0;
   evolution.SetTime(time);
   solver->Init(evolution);
   ExperimentRungeKutta experiment_solver(
      ExperimentRungeKutta::Scheme::ClassicalRK4,
      [&filter](Vector &candidate, real_t step, bool final_stage)
      {
         return filter.Apply(candidate, step, final_stage);
      });
   experiment_solver.Init(evolution);

   Vector residual(solution.Size());
   evolution.Mult(solution, residual);
   const bool use_cfl_time_step = dt <= 0.0;
   real_t maximum_speed = evolution.GetMaxCharSpeed();
   MPI_Allreduce(MPI_IN_PLACE, &maximum_speed, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, mesh.GetComm());
   const real_t h = 1.0 / elements;
   if (use_cfl_time_step)
   {
      const real_t temporal_scale = radau_initialization
         ? std::pow(h, std::max(1.0, (order + 2.0) / 4.0)) : h;
      dt = cfl * temporal_scale /
           ((2 * order + 1) * maximum_speed);
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
                     "Scalar experiment step was unexpectedly rejected.");
      }
      else
      {
         solver->Step(solution, time, step_size);
         filter.Apply(solution, step_size, true);
      }
      if (use_cfl_time_step)
      {
         maximum_speed = evolution.GetMaxCharSpeed();
         MPI_Allreduce(MPI_IN_PLACE, &maximum_speed, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX,
                       mesh.GetComm());
         const real_t temporal_scale = radau_initialization
            ? std::pow(h, std::max(1.0, (order + 2.0) / 4.0)) : h;
         dt = cfl * temporal_scale /
              ((2 * order + 1) * maximum_speed);
      }
      ++steps;
      if (glvis.Enabled() &&
          (steps % vis_steps == 0 || time >= final_time))
      {
         ostringstream title;
         title << "advection, t = " << time;
         glvis.Send(mesh, solution, title.str());
      }
   }
   tic_toc.Stop();

   FunctionCoefficient exact([dimension, final_time](const Vector &x)
   {
      Vector foot(x);
      foot(0) -= (dimension == 1 ? 1.0 : 0.7) * final_time;
      if (dimension == 2) { foot(1) -= 0.3 * final_time; }
      return InitialValue(foot);
   });
   const real_t l1_error = solution.ComputeL1Error(exact);
   const real_t l2_error = solution.ComputeL2Error(exact);
   const real_t linf_error = solution.ComputeMaxError(exact);
   const auto [cell_average_error, downwind_error] =
      SuperconvergenceErrors(solution, exact);
   const real_t conservation_drift = RelativeConservationDrift(
      initial_integral, GlobalScalarIntegral(solution));
   real_t minimum = solution.Min();
   real_t maximum = solution.Max();
   MPI_Allreduce(MPI_IN_PLACE, &minimum, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MIN, mesh.GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &maximum, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, mesh.GetComm());
   real_t runtime = tic_toc.RealTime();
   MPI_Allreduce(MPI_IN_PLACE, &runtime, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, mesh.GetComm());
   WriteScalarSamples(solution, profile_prefix, Mpi::WorldRank());
   const StudyFilterStatistics filter_stats = filter.Statistics();

   if (Mpi::Root())
   {
      cout.precision(16);
      cout << "method=" << StudyMethodName(method)
           << " cadence=" << FilterCadenceName(filter.Cadence())
           << " dimension=" << dimension << " elements=" << elements
           << " order=" << order << " steps=" << steps
           << " dofs=" << global_unknowns
           << " time=" << final_time
           << " l1_error=" << l1_error
           << " l2_error=" << l2_error
           << " linf_error=" << linf_error
           << " cell_average_l2_error=" << cell_average_error
           << " downwind_l2_error=" << downwind_error
           << " conservation_drift=" << conservation_drift
           << " minimum=" << minimum
           << " maximum=" << maximum
           << " filter_applications=" << filter_stats.applications
           << " active_elements=" << filter_stats.active_elements
           << " runtime_seconds=" << runtime << '\n';
   }
   return 0;
}
