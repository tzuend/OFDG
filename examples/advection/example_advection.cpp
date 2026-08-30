// In this file it is tested that the stabilization doesn't degrade the solution over time.
// i.e. on a traveling wave over periodic mesh after traversing the domain 100 times
// should look very similar to its intial condition.
// additionally simple traveling wave has a analytical solution to which we can compare to get the error.

// importantly both standard DG and stabilized DG should be compared to eachother to make sense of the result.

// This is a strictly 1D test for now

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

#include "../../src/ofdg_serial_optimized.hpp"
#include "../../src/rk4.hpp"
#include "../../src/kxrcf.hpp"

using namespace std;
using namespace mfem;

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
real_t u0_function(const Vector &x);

// Inflow boundary condition
real_t inflow_function(const Vector &x);

// Mesh bounding box
Vector bb_min, bb_max;

class DG_Solver : public Solver
{
private:
   SparseMatrix &M, &K, A;
   GMRESSolver linear_solver;
   BlockILU prec;
   real_t dt;
public:
   DG_Solver(SparseMatrix &M_, SparseMatrix &K_, const FiniteElementSpace &fes)
      : M(M_),
        K(K_),
        prec(fes.GetTypicalFE()->GetDof(),
             BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
        dt(-1.0)
   {
      linear_solver.iterative_mode = false;
      linear_solver.SetRelTol(1e-9);
      linear_solver.SetAbsTol(0.0);
      linear_solver.SetMaxIter(100);
      linear_solver.SetPrintLevel(0);
      linear_solver.SetPreconditioner(prec);
   }

   void SetTimeStep(real_t dt_)
   {
      if (dt_ != dt)
      {
         dt = dt_;
         // Form operator A = M - dt*K
         A = K;
         A *= -dt;
         A += M;

         // this will also call SetOperator on the preconditioner
         linear_solver.SetOperator(A);
      }
   }

   void SetOperator(const Operator &op) override
   {
      linear_solver.SetOperator(op);
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      linear_solver.Mult(x, y);
   }
};

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
template <typename T>
class FE_Evolution : public TimeDependentOperator
{
private:
   BilinearForm &M, &K;
   const Vector &b;
   const FiniteElementSpace *fes; // NEW: For stabilization S(u)
   Solver *M_prec;
   CGSolver M_solver;
   DG_Solver *dg_solver;

   mutable Vector z;
   mutable Vector S; // NEW: Stabilization vector

   bool use_stabilization = true; // NEW: Flag to enable/disable stabilization
   T ofdg; // NEW: Object to compute stabilization S(u)

public:
   FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_, bool use_stabilization_, const int btype_, VectorFunctionCoefficient &velocity_);

   void Mult(const Vector &x, Vector &y) const override;
   void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   ~FE_Evolution() override;

   T& GetStabilizer() { return ofdg; }
};


int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   problem = 0;
   const char *mesh_file = "../../../mfem/data/periodic-hexagon.mesh";
   // const char *mesh_file = "";
   int ref_levels = 0;
   int order = 3;
   bool pa = false;
   bool ea = false;
   bool fa = false;
   const char *device_config = "cpu";
   int ode_solver_type = 4; // should be RK4
   real_t t_final = 1;
   real_t dt = 0.0001;
   bool visualization = true;
   bool visit = false;  
   bool paraview = false;
   bool binary = false;
   int vis_steps = 5;
   bool solve_implicit_state = false;
   bool use_stabilization = true;
   int mesh_res = 50; // number of elements in the periodic mesh if no mesh file is provided
   const char *output_prefix = "solution";
   const char *output_dir = "out";
   bool only_final = false; // only save final solution, not initial condition
   double sigma_factor = 1.0; // factor to multiply the stabilization sigma values with, for testing
   int view_mode = 0; // 0: ofdg1d, 1: ofdg, 2: difference between ofdg and stabilization
   bool use_kxrcf = true; // whether to use KXRCF indicator for stabilization

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem setup to use. See options in velocity_function().");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&ea, "-ea", "--element-assembly", "-no-ea",
                  "--no-element-assembly", "Enable Element Assembly.");
   args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                  "--no-full-assembly", "Enable Full Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  ODESolver::Types.c_str());
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&solve_implicit_state, "-imp-state", "--implicit-state",
                  "-imp-slope", "--implicit-slope",
                  "Implicitly solve for stage state or slope.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                  "--no-paraview-datafiles",
                  "Save data files for ParaView (paraview.org) visualization.");
   args.AddOption(&binary, "-binary", "--binary-datafiles", "-ascii",
                  "--ascii-datafiles",
                  "Use binary (Sidre) or ascii format for VisIt data files.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&use_stabilization, "-stab", "--use-stabilization", "-no-stab", "--no-stabilization",
                  "Enable or disable the new stabilization term S(u).");
   args.AddOption(&mesh_res, "-mr", "--mesh-resolution",
                  "Number of elements in the periodic mesh if no mesh file is provided.");
   args.AddOption(&output_prefix, "-op", "--output-prefix",
                  "Prefix for output files (mesh and solution).");
   args.AddOption(&output_dir, "-od", "--output-dir",
                  "Directory for output files.");
   args.AddOption(&only_final, "-of", "--only-final", "-no-of", "--no-only-final",
                  "Only save final solution, not initial condition.");
   args.AddOption(&sigma_factor, "-sf", "--sigma-factor",
                  "Factor to multiply the stabilization sigma values with, for testing.");
   args.AddOption(&view_mode, "-vm", "--view-mode",
                  "Mode for viewing results: 0: ofdg1d, 1: ofdg, 2: difference between ofdg and stabilization.");
   args.AddOption(&use_kxrcf, "-kxrcf", "--use-kxrcf", "-no-kxrcf", "--no-use-kxrcf",
                  "Enable or disable the use of KXRCF indicator for stabilization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   Device device(device_config);
   device.Print();

   // 2. Read the mesh from the given mesh file. We can handle geometrically
   //    periodic meshes in this code.
   Mesh mesh;
   if (strlen(mesh_file) == 0 || (mesh_file[0] == '1' && mesh_file[1] == 'D'))
   {
      cout << "used 1D mesh" << endl;
      int ne = mesh_res; // number of elements
      Mesh mesh_init = Mesh::MakeCartesian1D(ne);// Make a mesh of the unit interval with 10 elements
      // Create the vertex mapping. To begin, create the identity mapping.
      std::vector<int> v2v(mesh_init.GetNV());
      for (int i = 0; i < mesh_init.GetNV(); ++i)
      {
         v2v[i] = i;
      }
      // Modify the mapping so that the last vertex gets mapped to the first vertex.
      v2v.back() = 0;
      mesh = Mesh::MakePeriodic(mesh_init, v2v); // Create the periodic mesh

      real_t h = 1.0 / ne;
      dt = 0.1 * h / (2*order+1);
   } else if (mesh_file[0] == '2' && mesh_file[1] == 'D') {
      // 2D mesh
      int ne = mesh_res; // number of elements

      // Create a 10x10 quad mesh of the unit square;
      Mesh mesh_init = Mesh::MakeCartesian2D(ne, ne, Element::QUADRILATERAL);
      // Create translation vectors defining the periodicity
      Vector x_translation({1.0, 0.0});
      Vector y_translation({0.0, 1.0});
      std::vector<Vector> translations = {x_translation, y_translation};
      // Create the periodic mesh using the vertex mapping defined by the translation vectors
      mesh = Mesh::MakePeriodic(mesh_init, mesh_init.CreatePeriodicVertexMapping(translations));
   }
   else 
   {
      mesh = Mesh(mesh_file, 1, 1);
   }
   int dim = mesh.Dimension();

   // 3. Define the ODE solver used for time integration. Several explicit
   //    Runge-Kutta methods are available.
   // unique_ptr<ODESolver> ode_solver = ODESolver::Select(ode_solver_type);
   // unique_ptr<ODESolver> ode_solver_nd = ODESolver::Select(ode_solver_type);
   // unique_ptr<ODESolver> ode_solver_stab_nd = ODESolver::Select(ode_solver_type);

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement, where 'ref_levels' is a
   //    command-line parameter. If the mesh is of NURBS type, we convert it to
   //    a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh.UniformRefinement();
   }
   if (mesh.NURBSext)
   {
      mesh.SetCurvature(max(order, 1));
   }
   mesh.GetBoundingBox(bb_min, bb_max, max(order, 1));

   // 5. Define the discontinuous DG finite element space of the given
   //    polynomial order on the refined mesh.
   DG_FECollection fec(order, dim, BasisType::GaussLobatto);
   FiniteElementSpace fes(&mesh, &fec);

   cout << "Number of unknowns: " << fes.GetVSize() << endl;

   // 6. Set up and assemble the bilinear and linear forms corresponding to the
   //    DG discretization. The DGTraceIntegrator involves integrals over mesh
   //    interior faces.
   VectorFunctionCoefficient velocity(dim, velocity_function);
   FunctionCoefficient inflow(inflow_function);
   FunctionCoefficient u0(u0_function);

   BilinearForm m(&fes);
   BilinearForm k(&fes);
   if (pa)
   {
      m.SetAssemblyLevel(AssemblyLevel::PARTIAL);
      k.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   }
   else if (ea)
   {
      m.SetAssemblyLevel(AssemblyLevel::ELEMENT);
      k.SetAssemblyLevel(AssemblyLevel::ELEMENT);
   }
   else if (fa)
   {
      m.SetAssemblyLevel(AssemblyLevel::FULL);
      k.SetAssemblyLevel(AssemblyLevel::FULL);
   }
   m.AddDomainIntegrator(new MassIntegrator);
   constexpr real_t alpha = -1.0;
   k.AddDomainIntegrator(new ConvectionIntegrator(velocity, alpha));
   k.AddInteriorFaceIntegrator(
      new NonconservativeDGTraceIntegrator(velocity, alpha));
   k.AddBdrFaceIntegrator(
      new NonconservativeDGTraceIntegrator(velocity, alpha));

   LinearForm b(&fes);
   b.AddBdrFaceIntegrator(
      new BoundaryFlowIntegrator(inflow, velocity, alpha));

   m.Assemble();
   int skip_zeros = 0;
   k.Assemble(skip_zeros);
   b.Assemble();
   m.Finalize();
   k.Finalize(skip_zeros);

   // 7. Define the initial conditions, save the corresponding grid function to
   //    a file and (optionally) save data in the VisIt format and initialize
   //    GLVis visualization.
   GridFunction u(&fes);
   u.ProjectCoefficient(u0);

   if (!only_final)
   {
      ofstream omesh(std::string(output_dir) + "/" + output_prefix + std::string(".mesh"));
      omesh.precision(precision);
      mesh.Print(omesh);
      ofstream osol(std::string(output_dir) + "/" + output_prefix + std::string("_init.gf"));
      osol.precision(precision);
      u.Save(osol);
   }

   // Create data collection for solution output: either VisItDataCollection for
   // ascii data files, or SidreDataCollection for binary data files.
   DataCollection *dc = NULL;
   if (visit)
   {
      if (binary)
      {
#ifdef MFEM_USE_SIDRE
         dc = new SidreDataCollection("Example9", &mesh);
#else
         MFEM_ABORT("Must build with MFEM_USE_SIDRE=YES for binary output.");
#endif
      }
      else
      {
         dc = new VisItDataCollection("Example9", &mesh);
         dc->SetPrecision(precision);
      }
      dc->RegisterField("solution", &u);
      dc->SetCycle(0);
      dc->SetTime(0.0);
      dc->Save();
   }

   ParaViewDataCollection *pd = NULL;
   if (paraview)
   {
      pd = new ParaViewDataCollection("Example9", &mesh);
      pd->SetPrefixPath("ParaView");
      pd->RegisterField("solution", &u);
      pd->SetLevelsOfDetail(order);
      pd->SetDataFormat(VTKFormat::BINARY);
      pd->SetHighOrderOutput(true);
      pd->SetCycle(0);
      pd->SetTime(0.0);
      pd->Save();
   }

   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout)
      {
         cout << "Unable to connect to GLVis server at "
              << vishost << ':' << visport << endl;
         visualization = false;
         cout << "GLVis visualization disabled.\n";
      }
      else
      {
         sout.precision(precision);
         sout << "solution\n" << mesh << u;
         sout << "pause\n";
         sout << flush;
         cout << "GLVis visualization paused."
              << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   // 8. Define the time-dependent evolution operator describing the ODE
   //    right-hand side, and perform time-integration (looping over the time
   //    iterations, ti, with a time-step dt).
   
   FE_Evolution<OFDG> adv_nd(m, k, b, false, fec.GetBasisType(), velocity);
   using ImplicitVariableType = FE_Evolution<OFDG>::ImplicitVariableType;
   ImplicitVariableType imp_var = solve_implicit_state ?
                                  ImplicitVariableType::STATE
                                  : ImplicitVariableType::SLOPE;
   adv_nd.SetImplicitVariableType(imp_var);

   FE_Evolution<OFDG> adv_stab_nd(m, k, b, false, fec.GetBasisType(), velocity);
   using ImplicitVariableType = FE_Evolution<OFDG>::ImplicitVariableType;
   ImplicitVariableType imp_var_stab = solve_implicit_state ?
                                  ImplicitVariableType::STATE
                                  : ImplicitVariableType::SLOPE;                
   adv_stab_nd.SetImplicitVariableType(imp_var_stab);

   unique_ptr<ODESolver> ode_solver_nd = ODESolver::Select(ode_solver_type);
   // auto ode_solver_stab_nd =
   //       std::make_unique<OEDG_RK4Solver<OFDG>>(adv_stab_nd.GetStabilizer());
   unique_ptr<ODESolver> ode_solver_stab_nd = ODESolver::Select(ode_solver_type);
   OFDG& stabilizer = adv_stab_nd.GetStabilizer();
   KXRCFIndicator kxrcf(&fes, &velocity, 1.0);
   Vector u_stab_tmp(fes.GetVSize());
   Array<bool> active_elements;

   real_t t = 0.0;

   adv_nd.SetTime(t);
   ode_solver_nd->Init(adv_nd);

   adv_stab_nd.SetTime(t);
   ode_solver_stab_nd->Init(adv_stab_nd);

   GridFunction u_stab(&fes); 
   u_stab = u;  // copy initial condition to another grid function for stabilized case

   double area = 0.0;
   int counter = 0;

   double max_diff = 0.0;

   bool done = false;
   for (int ti = 0; !done; )
   {
      real_t dt_real = min(dt, t_final - t);
      real_t t_copy = t; // copy of dt for stabilized case, in case we want to use a different time step for stabilitys
      ode_solver_nd->Step(u, t, dt_real);
      ode_solver_stab_nd->Step(u_stab, t_copy, dt_real);

      if (use_kxrcf)
      {
         kxrcf.Compute(u_stab, active_elements);
         stabilizer.CompDecay(u_stab, u_stab_tmp, dt_real, &active_elements);
      }
      else
      {
         stabilizer.CompDecay(u_stab, u_stab_tmp, dt_real);
      }
      u_stab = u_stab_tmp;
      ti++;
      
      done = (t >= t_final - 1e-8*dt);

      if (done || ti % vis_steps == 0)
      {
         cout << "time step: " << ti << ", time: " << t << endl;

         if (visualization)
         {
            if (view_mode == 0)
            {
               sout << "solution\n" << mesh << u << flush;
            }
            else if (view_mode == 1)
            {
               sout << "solution\n" << mesh << u_stab << flush;
            }
            else if (view_mode == 2)
            {
               GridFunction diff(&fes);
               diff = u_stab;
               diff -= u;

               if (max_diff < diff.Max())
               {
                  max_diff = diff.Max();
               } else if (max_diff < -diff.Min())
               {
                  max_diff = -diff.Min();
               }

               sout << "solution\n" << mesh << diff << flush;
            }
         }

         if (visit)
         {
            dc->SetCycle(ti);
            dc->SetTime(t);
            dc->Save();
         }

         if (paraview)
         {
            pd->SetCycle(ti);
            pd->SetTime(t);
            pd->Save();
         }
      }

      // if (counter >= 0) {
      //    done = true; // only a few steps for testing
      // }
      // counter++;
   }

   if (view_mode == 2)
   {
      std::cout << "max difference between stabilized and non-stabilized solution: " << max_diff << std::endl;
   }

   // 9. Save the final solution. This output can be viewed later using GLVis:
   //    "glvis -m ex9.mesh -g ex9-final.gf".
   {
      ofstream osol_nd(std::string(output_dir) + "/" + output_prefix + std::string("_final.gf"));
      osol_nd.precision(precision);
      u.Save(osol_nd);

      ofstream osol_stab_nd(std::string(output_dir) + "/" + output_prefix + std::string("_final_stab.gf"));
      osol_stab_nd.precision(precision);
      u_stab.Save(osol_stab_nd);
   }

   // 10. Free the used memory.
   delete pd;
   delete dc;

   return 0;
}


// Implementation of class FE_Evolution
template <typename T>
FE_Evolution<T>::FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_, bool use_stabilization_, const int btype_, VectorFunctionCoefficient &velocity_)
   : TimeDependentOperator(M_.FESpace()->GetTrueVSize()),
     M(M_), K(K_), b(b_), z(height), S(height), fes(M_.FESpace()), use_stabilization(use_stabilization_), ofdg(fes, btype_, velocity_)
{
   Array<int> ess_tdof_list;
   if (M.GetAssemblyLevel() == AssemblyLevel::LEGACY)
   {
      M_prec = new DSmoother(M.SpMat());
      M_solver.SetOperator(M.SpMat());
      dg_solver = new DG_Solver(M.SpMat(), K.SpMat(), *M.FESpace());
   }
   else
   {
      M_prec = new OperatorJacobiSmoother(M, ess_tdof_list);
      M_solver.SetOperator(M);
      dg_solver = NULL;
   }
   M_solver.SetPreconditioner(*M_prec);
   M_solver.iterative_mode = false;
   M_solver.SetRelTol(1e-9);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);

   std::cout << "\nFE_Evolution: Using stabilization S(u) = sum( delta^l*h*( u-P^l(u) ) )\n" << std::endl;
}


// We add the new stabilization term to the rhs as S(u)
// The time derivative is then du/dt = M^{-1} (K u + b + S(u))
// S(u) is sum( delta^l*h*( u-P^l(u) ) )
// delta^l
// h is the element size
// P^l(u) is the projection of u to the space of degree l polynomials on each element
template <typename T>
void FE_Evolution<T>::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b) - S(x) 
   K.Mult(x, z);
   z += b;
   
   
   M_solver.Mult(z, y);

   if (use_stabilization)
   {
      // NEW: Compute stabilization S(x) and add to rhs
      ofdg.ComputeStabilization(x, S);
      y -= S;
   }
}

template <typename T>
void FE_Evolution<T>::ImplicitSolve(const real_t dt, const Vector &x, Vector &k)
{
   MFEM_VERIFY(dg_solver != NULL,
               "Implicit time integration is not supported with partial assembly");
   // Construct current right-hand side for stage state vs. slope solve
   if (ImplicitVarTypeIsState())
   {
      // k, on return, is the stage value u
      M.Mult(x, z);
   }
   else
   {
      // k, on return, is the stage slope du/dt
      K.Mult(x, z);
   }
   z += b;
   dg_solver->SetTimeStep(dt);
   dg_solver->Mult(z, k);
}

template <typename T>
FE_Evolution<T>::~FE_Evolution()
{
   delete M_prec;
   delete dg_solver;
}

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      real_t center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   switch (problem)
   {
      case 0:
      {
         // Translations in 1D, 2D, and 3D
         switch (dim)
         {
            case 1: v(0) = 2; break;
            // case 2: v(0) = sqrt(2./3.); v(1) = sqrt(1./3.); break;
            // case 2: v(0) = 1.; v(1) = 1.; break;
            case 2: v(0) = 1; v(1) = 1; break;
            case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
               break;
         }
         break;
      }
      case 1:
      case 2:
      {
         // Clockwise rotation in 2D around the origin
         const real_t w = M_PI/2;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = w*X(1); v(1) = -w*X(0); break;
            case 3: v(0) = w*X(1); v(1) = -w*X(0); v(2) = 0.0; break;
         }
         break;
      }
      case 3:
      {
         // Clockwise twisting rotation in 2D around the origin
         const real_t w = M_PI/2;
         real_t d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
         d = d*d;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = d*w*X(1); v(1) = -d*w*X(0); break;
            case 3: v(0) = d*w*X(1); v(1) = -d*w*X(0); v(2) = 0.0; break;
         }
         break;
      }
   }
}

// Initial condition
real_t u0_function(const Vector &x)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      real_t center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   // print domain bounds for debugging
   std::cout << "Domain bounds: ";
   for (int i = 0; i < dim; i++)
   {
      std::cout << "[" << bb_min[i] << ", " << bb_max[i] << "] ";
   }
   std::cout << std::endl;

   switch (problem)
   {
      case 0:
         switch (dim)
         {
            case 1:
               return sin(X(0)*M_PI); // smooth
               // if (X(0) >= -0.5 && X(0) <= 0.5) {
               //    return 1.0;
               // } else {
               //    return 0.0;
               // }
            case 2:
               return sin(M_PI*X(0))*sin(M_PI*X(1)); // smooth
               // return sin(M_PI*X(0));
            case 3:
               return sin(M_PI*X(0))*sin(M_PI*X(1))*sin(M_PI*X(2)); // smooth
         }
      case 1:
      {
         switch (dim)
         {
            case 1:
               return exp(-40.*pow(X(0)-0.5,2));
            case 2:
            case 3:
            {
               real_t rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
               if (dim == 3)
               {
                  const real_t s = (1. + 0.25*cos(2*M_PI*X(2)));
                  rx *= s;
                  ry *= s;
               }
               return ( std::erfc(w*(X(0)-cx-rx))*std::erfc(-w*(X(0)-cx+rx)) *
                        std::erfc(w*(X(1)-cy-ry))*std::erfc(-w*(X(1)-cy+ry)) )/16;
            }
         }
      }
      case 2:
      {
         real_t x_ = X(0), y_ = X(1), rho, phi;
         rho = std::hypot(x_, y_);
         phi = atan2(y_, x_);
         return pow(sin(M_PI*rho),2)*sin(3*phi);
      }
      case 3:
      {
         const real_t f = M_PI;
         return sin(f*X(0))*sin(f*X(1));
      }
   }
   return 0.0;
}

// Inflow boundary condition (zero for the problems considered in this example)
real_t inflow_function(const Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3: return 0.0;
   }
   return 0.0;
}
