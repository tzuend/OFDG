// Isolated 3D volume transport showcase. Uses production StudyFilter unchanged.
#include "mfem.hpp"
#include "../examples/euler/euler.hpp"
#include "../src/study_filter.hpp"
#include "../src/conservation.hpp"
#include <fstream>
#include <iomanip>
using namespace mfem;
using namespace ofdg;

real_t Initial(const Vector &x, int problem)
{
   if (problem == 2) { return 1.; }
   if (problem == 0) { return 1. + .2*sin(2*M_PI*x(0))*sin(2*M_PI*x(1))*sin(2*M_PI*x(2)); }
   real_t r2=0.;
   const double center[3]={.3,.35,.4};
   for(int d=0;d<3;++d) { double a=x(d)-center[d]; a-=std::round(a); r2+=a*a; }
   return r2 < .04 ? 1. : 0.;
}

void Export(ParGridFunction &u, const std::string &prefix, double time, int problem)
{
   std::ofstream file(prefix+".csv"); file<<std::setprecision(17)<<"element,x,y,z,u,exact\n";
   auto *space=u.ParFESpace();
   for(int e=0;e<space->GetNE();++e) {
      auto *tr=space->GetElementTransformation(e);
      for(int k=0;k<5;++k) for(int j=0;j<5;++j) for(int i=0;i<5;++i) {
         IntegrationPoint ip; ip.Set3((i+.5)/5.,(j+.5)/5.,(k+.5)/5.);
         Vector x; tr->Transform(ip,x); Vector foot(x);
         foot(0)-=.5*time; foot(1)-=.3*time; foot(2)-=.2*time;
         file<<e<<','<<x(0)<<','<<x(1)<<','<<x(2)<<','<<u.GetValue(e,ip)<<','<<Initial(foot,problem)<<'\n';
      }
   }
   std::ofstream meshfile(prefix+".mesh"), gridfile(prefix+".gf");
   meshfile<<std::setprecision(17); gridfile<<std::setprecision(17);
   space->GetParMesh()->Print(meshfile); u.Save(gridfile);
}

int main(int argc, char **argv)
{
   Mpi::Init(argc,argv); Hypre::Init();
   MFEM_VERIFY(argc==7 && Mpi::WorldSize()==1,"Usage: volume_showcase N method problem dt T prefix; one rank");
   int n=std::stoi(argv[1]), problem=std::stoi(argv[3]);
   double dt=std::stod(argv[4]), final=std::stod(argv[5]); std::string prefix=argv[6];
   Mesh base=Mesh::MakeCartesian3D(n,n,n,Element::HEXAHEDRON,1.,1.,1.);
   std::vector<Vector> translations={Vector({1.,0.,0.}),Vector({0.,1.,0.}),Vector({0.,0.,1.})};
   Mesh periodic=Mesh::MakePeriodic(base,base.CreatePeriodicVertexMapping(translations));
   ParMesh mesh(MPI_COMM_WORLD,periodic);
   DG_FECollection collection(2,3); ParFiniteElementSpace space(&mesh,&collection); ParGridFunction u(&space);
   FunctionCoefficient initial([problem](const Vector &x){return Initial(x,problem);});
   u.ProjectCoefficient(initial); // Same nodal initialization for both methods.
   double mass0=GlobalScalarIntegral(u);
   VectorFunctionCoefficient velocity(3,[](const Vector &,Vector &v){v.SetSize(3);v(0)=.5;v(1)=.3;v(2)=.2;});
   AdvectionFlux flux(velocity); RusanovFlux numerical_flux(flux);
   DGHyperbolicConservationLaws evolution(space,std::make_unique<HyperbolicFormIntegrator>(numerical_flux,1),true);
   auto physics=std::make_shared<AdvectionFacePhysics>(&velocity);
   StudyFilter filter(&space,BasisType::GaussLegendre,physics,ParseStudyMethod(argv[2]),FilterCadence::Step,1.);
   RK4Solver solver; solver.Init(evolution); double t=0.; int steps=0;
   if(prefix!="-") Export(u,prefix+"-initial",0.,problem);
   StopWatch timer; timer.Start();
   while(t<final) { double step=std::min(dt,final-t); solver.Step(u,t,step); filter.Apply(u,step,true); ++steps; }
   timer.Stop();
   double error2=0.,norm2=0.,minimum=1e100,maximum=-1e100;
   const auto &rule=IntRules.Get(Geometry::CUBE,16);
   for(int e=0;e<space.GetNE();++e) {
      auto *tr=space.GetElementTransformation(e);
      for(int q=0;q<rule.GetNPoints();++q) {
         const auto &ip=rule.IntPoint(q); tr->SetIntPoint(&ip); Vector x; tr->Transform(ip,x);
         x(0)-=.5*t; x(1)-=.3*t; x(2)-=.2*t;
         double truth=Initial(x,problem), val=u.GetValue(e,ip),w=ip.weight*tr->Weight();
         error2+=w*(val-truth)*(val-truth);norm2+=w*truth*truth;
         minimum=std::min(minimum,val);maximum=std::max(maximum,val);
      }
   }
   MFEM_VERIFY(std::isfinite(error2)&&std::isfinite(minimum)&&std::isfinite(maximum),"Nonfinite result");
   if(prefix!="-") Export(u,prefix+"-final",t,problem);
   std::cout<<std::setprecision(17)<<"n,method,problem,dofs,steps,time,dt,l2,relative_l2,mass_drift,min,max,filter_applications,active_elements,solver_seconds\n"
     <<n<<','<<argv[2]<<','<<problem<<','<<space.GlobalTrueVSize()<<','<<steps<<','<<t<<','<<dt<<','<<sqrt(error2)<<','<<sqrt(error2/norm2)<<','<<GlobalScalarIntegral(u)-mass0<<','<<minimum<<','<<maximum<<','<<filter.Statistics().applications<<','<<filter.Statistics().active_elements<<','<<timer.RealTime()<<'\n';
}
