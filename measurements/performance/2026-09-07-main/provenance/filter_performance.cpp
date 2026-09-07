// End-to-end scalar RK4 timing with production DG/OFDG/KXRCF operators.
// Single rank isolates mesh/degree growth; clocks are outside kernel calls.
#include "mfem.hpp"
#include "../examples/euler/euler.hpp"
#include "../src/ofdg.hpp"
#include "../src/kxrcf.hpp"
#include "../src/conservation.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace mfem;
using namespace ofdg;

int main(int argc, char **argv)
{
   Mpi::Init(argc, argv); Hypre::Init();
   MFEM_VERIFY(Mpi::WorldSize()==1, "This benchmark isolates single-rank cost.");
   int n=16, p=2, repeats=5;
   bool curved=false, jump=false;
   double tf=0.05, minimum_sample=0.25;
   OptionsParser args(argc,argv);
   args.AddOption(&n,"-n","--resolution","Cells per direction.");
   args.AddOption(&p,"-o","--order","Solution degree.");
   args.AddOption(&repeats,"-r","--repeats","Measured repeats, plus one warm-up.");
   args.AddOption(&tf,"-tf","--final-time","Physical time per reset trajectory.");
   args.AddOption(&minimum_sample,"-ms","--minimum-sample","Minimum timed seconds per sample.");
   args.AddOption(&curved,"-curved","--curved","-affine","--affine","Cubic deformed geometry.");
   args.AddOption(&jump,"-jump","--jump","-smooth","--smooth","Periodic discontinuous initial data.");
   args.ParseCheck();
   MFEM_VERIFY(n>0 && p>0 && repeats>0 && tf>0 && minimum_sample>=0,"Invalid settings.");
   double setup_start=MPI_Wtime();
   Mesh serial=Mesh::MakeCartesian2D(n,n,Element::QUADRILATERAL,true,1.,1.);
   if(curved) {
      serial.SetCurvature(3);
      serial.Transform([](const Vector &x,Vector &y) {
         y=x; const double d=.04*std::sin(2*M_PI*x(0))*std::sin(2*M_PI*x(1));
         y(0)+=d; y(1)+=d;
      });
   }
   std::vector<Vector> translations={Vector({1.,0.}),Vector({0.,1.})};
   Mesh periodic=Mesh::MakePeriodic(serial,serial.CreatePeriodicVertexMapping(translations));
   ParMesh mesh(MPI_COMM_WORLD,periodic);
   DG_FECollection collection(p,2);
   ParFiniteElementSpace space(&mesh,&collection);
   ParGridFunction initial(&space),state(&space);
   FunctionCoefficient initial_value([jump](const Vector &x) {
      return jump ? (x(0)<.5 ? 1. : 1.4) : std::pow(std::sin(M_PI*(x(0)+x(1))),2);
   });
   initial.ProjectCoefficient(initial_value);
   VectorFunctionCoefficient velocity(2,[](const Vector &,Vector &v) { v.SetSize(2);v(0)=.7;v(1)=.3; });
   AdvectionFlux flux(velocity); RusanovFlux numerical_flux(flux);
   DGHyperbolicConservationLaws evolution(space,std::make_unique<HyperbolicFormIntegrator>(numerical_flux,1),true);
   auto physics=std::make_shared<AdvectionFacePhysics>(&velocity);
   Vector residual(space.GetVSize()),filtered(space.GetVSize());
   evolution.Mult(initial,residual);
   const double dt=.15/(n*(2*p+1)*evolution.GetMaxCharSpeed());
   const double common_setup=MPI_Wtime()-setup_start;
   const double initial_mass=GlobalScalarIntegral(initial);
   FunctionCoefficient exact([jump,tf](const Vector &x) {
      double xx=x(0)-.7*tf, yy=x(1)-.3*tf;xx-=std::floor(xx);
      return jump ? (xx<.5 ? 1. : 1.4) : std::pow(std::sin(M_PI*(xx+yy)),2);
   });
   const char *names[]={"dg","ofdg","ofdg-kxrcf"};
   // Rotate method order between repetitions to limit systematic drift bias.
   for(int repeat=-1;repeat<repeats;++repeat) for(int slot=0;slot<3;++slot) {
      const int method=(slot+repeat+1)%3;
      const double filter_setup_start=MPI_Wtime();
      std::unique_ptr<OFDG> filter;
      std::unique_ptr<KXRCFIndicator> detector;
      if(method) filter=std::make_unique<OFDG>(&space,BasisType::GaussLegendre,physics);
      if(method==2) detector=std::make_unique<KXRCFIndicator>(&space,physics);
      const double filter_setup=MPI_Wtime()-filter_setup_start;
      Array<bool> mask;
      double dg_seconds=0,detector_seconds=0,decay_seconds=0,total_seconds=0;
      long steps=0,cycles=0,active_count=0;
      do {
         state=initial; RK4Solver solver; solver.Init(evolution);
         double t=0; evolution.SetTime(0);
         const double start=MPI_Wtime();
         while(t<tf) {
            double step=std::min(dt,tf-t);
            const double a=MPI_Wtime(); solver.Step(state,t,step);
            const double b=MPI_Wtime(); dg_seconds+=b-a;
            if(detector) {
               detector->Compute(state,mask);
               for(int e=0;e<mask.Size();++e) active_count+=mask[e] ? 1 : 0;
            }
            const double c=MPI_Wtime(); detector_seconds+=detector ? c-b : 0.;
            if(filter) { filter->CompDecay(state,filtered,step,detector ? &mask : nullptr);state=filtered; }
            const double d=MPI_Wtime(); decay_seconds+=filter ? d-c : 0.;
            ++steps;
         }
         total_seconds+=MPI_Wtime()-start; ++cycles;
      } while(total_seconds<minimum_sample);
      const double mass_drift=std::abs(GlobalScalarIntegral(state)-initial_mass);
      const double error=state.ComputeL2Error(exact);
      MFEM_VERIFY(std::isfinite(error) && mass_drift<1e-9,"Numerical sanity check failed.");
      std::cout<<std::setprecision(16)<<"BENCH n="<<n<<" order="<<p<<" curved="<<curved
         <<" jump="<<jump<<" repeat="<<repeat<<" method="<<names[method]<<" cells="<<n*n
         <<" dofs="<<space.GlobalTrueVSize()<<" cycles="<<cycles<<" steps="<<steps<<" tf="<<tf
         <<" dt="<<dt<<" common_setup="<<common_setup<<" filter_setup="<<filter_setup
         <<" dg_seconds="<<dg_seconds<<" detector_seconds="<<detector_seconds
         <<" decay_seconds="<<decay_seconds<<" total_seconds="<<total_seconds
         <<" active_count="<<active_count<<" mass_drift="<<mass_drift<<" l2_error="<<error<<std::endl;
   }
}
