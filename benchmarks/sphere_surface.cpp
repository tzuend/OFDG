// Isolated surface-FEM screened Poisson study, adapted from pinned MFEM ex7.
// A 2D surface embedded in R^3; no production DG/OFDG source is changed.
#include "mfem.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace mfem;

std::unique_ptr<Mesh> SphereMesh(int elem_type, int geometry_order, int levels)
{
   int Nvert = 8, Nelem = 6;
   if (elem_type == 0)
   {
      Nvert = 6;
      Nelem = 8;
   }
   Mesh *mesh = new Mesh(2, Nvert, Nelem, 0, 3);

   if (elem_type == 0) // inscribed octahedron
   {
      const real_t tri_v[6][3] =
      {
         { 1,  0,  0}, { 0,  1,  0}, {-1,  0,  0},
         { 0, -1,  0}, { 0,  0,  1}, { 0,  0, -1}
      };
      const int tri_e[8][3] =
      {
         {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4},
         {1, 0, 5}, {2, 1, 5}, {3, 2, 5}, {0, 3, 5}
      };

      for (int j = 0; j < Nvert; j++)
      {
         mesh->AddVertex(tri_v[j]);
      }
      for (int j = 0; j < Nelem; j++)
      {
         int attribute = j + 1;
         mesh->AddTriangle(tri_e[j], attribute);
      }
      mesh->FinalizeTriMesh(1, 1, true);
   }
   else // inscribed cube
   {
      const real_t quad_v[8][3] =
      {
         {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
         {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}
      };
      const int quad_e[6][4] =
      {
         {3, 2, 1, 0}, {0, 1, 5, 4}, {1, 2, 6, 5},
         {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}
      };

      for (int j = 0; j < Nvert; j++)
      {
         mesh->AddVertex(quad_v[j]);
      }
      for (int j = 0; j < Nelem; j++)
      {
         int attribute = j + 1;
         mesh->AddQuad(quad_e[j], attribute);
      }
      mesh->FinalizeQuadMesh(1, 1, true);
   }


   mesh->SetCurvature(geometry_order, false, 3, Ordering::byVDIM);
   for (int level=0; level<levels; ++level) { mesh->UniformRefinement(); }
   // Match ex7's default: snap the final high-order nodes once.
   auto &nodes=*mesh->GetNodes();
   auto *fes=nodes.FESpace();
   for (int i=0; i<fes->GetNDofs(); ++i)
   {
      double radius=0;
      for (int d=0; d<3; ++d) { radius+=std::pow(nodes(fes->DofToVDof(i,d)),2); }
      radius=std::sqrt(radius);
      for (int d=0; d<3; ++d) { nodes(fes->DofToVDof(i,d))/=radius; }
   }
   return std::unique_ptr<Mesh>(mesh);
}

double Exact(const Vector &x) { return x(0)*x(1)/(x*x); }

// Export a tessellated rendering of the actual approximate surface, not a
// perfect-sphere replacement. The original element boundary is retained by id.
void ExportSurface(Mesh &mesh, const GridFunction &u, const std::string &prefix)
{
   std::ofstream meshfile(prefix+".mesh"), gridfile(prefix+".gf"), samples(prefix+"-surface.csv");
   meshfile<<std::setprecision(17); mesh.Print(meshfile);
   gridfile<<std::setprecision(17); u.Save(gridfile);
   samples<<"element,triangle,vertex,x,y,z,u,exact,error\n"<<std::setprecision(17);
   int tid=0; const int subdivisions=6;
   for (int e=0; e<mesh.GetNE(); ++e)
   {
      auto *tr=mesh.GetElementTransformation(e);
      auto emit=[&](double a,double b,double c,double d,double f,double g)
      {
         const double points[3][2]={{a,b},{c,d},{f,g}};
         for (int j=0; j<3; ++j)
         {
            IntegrationPoint ip; ip.Set2(points[j][0],points[j][1]);
            Vector x; tr->Transform(ip,x); double value=u.GetValue(e,ip), truth=Exact(x);
            samples<<e<<','<<tid<<','<<j<<','<<x(0)<<','<<x(1)<<','<<x(2)<<','
                   <<value<<','<<truth<<','<<value-truth<<'\n';
         }
         ++tid;
      };
      const double h=1./subdivisions;
      if (mesh.GetElementBaseGeometry(e)==Geometry::TRIANGLE)
      {
         for(int j=0;j<subdivisions;++j) for(int i=0;i<subdivisions-j;++i)
         {
            emit(i*h,j*h,(i+1)*h,j*h,i*h,(j+1)*h);
            if(i+j<subdivisions-1) emit((i+1)*h,j*h,(i+1)*h,(j+1)*h,i*h,(j+1)*h);
         }
      }
      else
      {
         for(int j=0;j<subdivisions;++j) for(int i=0;i<subdivisions;++i)
         {
            emit(i*h,j*h,(i+1)*h,j*h,(i+1)*h,(j+1)*h);
            emit(i*h,j*h,(i+1)*h,(j+1)*h,i*h,(j+1)*h);
         }
      }
   }
}

int main(int argc,char **argv)
{
   Mpi::Init(argc,argv);
   if(argc!=8) { std::cerr<<"family levels solution_order geometry_order extra constant output_prefix_or_dash\n"; return 2; }
   const std::string family=argv[1],prefix=argv[7];
   const int levels=std::stoi(argv[2]),p=std::stoi(argv[3]),g=std::stoi(argv[4]);
   const int extra=std::stoi(argv[5]); const bool constant=std::stoi(argv[6]);
   const auto start=std::chrono::steady_clock::now();
   auto mesh=SphereMesh(family=="tri"?0:1,g,levels);
   H1_FECollection fec(p,2); FiniteElementSpace space(mesh.get(),&fec);
   FunctionCoefficient exact([&](const Vector &x){return constant?1.:Exact(x);});
   FunctionCoefficient forcing([&](const Vector &x){return constant?1.:7*Exact(x);});
   ConstantCoefficient one(1.);
   const auto geom=family=="tri"?Geometry::TRIANGLE:Geometry::SQUARE;
   const int order=2*p+2*g+6+extra;
   const auto &assembly=IntRules.Get(geom,order);
   LinearForm rhs(&space); auto *load=new DomainLFIntegrator(forcing);load->SetIntRule(&assembly);rhs.AddDomainIntegrator(load);rhs.Assemble();
   BilinearForm form(&space);
   auto *diffusion=new DiffusionIntegrator(one);diffusion->SetIntRule(&assembly);
   auto *mass=new MassIntegrator(one);mass->SetIntRule(&assembly);
   form.AddDomainIntegrator(diffusion);form.AddDomainIntegrator(mass);form.Assemble();
   GridFunction u(&space);u=0.; Array<int> essential;SparseMatrix A;Vector X,B;
   form.FormLinearSystem(essential,u,rhs,A,X,B);
   GSSmoother preconditioner(A);CGSolver solver;
   solver.SetOperator(A);solver.SetPreconditioner(preconditioner);
   solver.SetRelTol(1e-12);solver.SetAbsTol(1e-14);solver.SetMaxIter(1000);solver.SetPrintLevel(0);
   solver.Mult(B,X);form.RecoverFEMSolution(X,rhs,u);
   Vector residual(B.Size()); A.Mult(X,residual);residual-=B;
   double l2=0,truthnorm=0,gradienterror=0,area=0,radiuserror=0,radialmax=0,mean=0,loadmean=0,maximumerror=0,minweight=1e100;
   const auto &evaluation=IntRules.Get(geom,order+8);
   for(int e=0;e<mesh->GetNE();++e)
   {
      auto *tr=mesh->GetElementTransformation(e);
      for(int q=0;q<evaluation.GetNPoints();++q)
      {
         const auto &ip=evaluation.IntPoint(q);tr->SetIntPoint(&ip);Vector x,gradient;
         tr->Transform(ip,x);const double weight=ip.weight*tr->Weight(),r2=x*x,r=std::sqrt(r2);
         const double truth=constant?1.:Exact(x),value=u.GetValue(e,ip),error=value-truth;
         u.GetGradient(*tr,gradient);
         Vector exact_gradient(3);exact_gradient=0.;
         if(!constant)
         {
            exact_gradient(0)=x(1)/r2;exact_gradient(1)=x(0)/r2;
            exact_gradient.Add(-2*x(0)*x(1)/(r2*r2),x);
         }
         // Project radial-extension gradient onto the discrete surface tangent.
         auto &J=tr->Jacobian();Vector normal(3);
         normal(0)=J(1,0)*J(2,1)-J(2,0)*J(1,1);
         normal(1)=J(2,0)*J(0,1)-J(0,0)*J(2,1);
         normal(2)=J(0,0)*J(1,1)-J(1,0)*J(0,1);normal/=normal.Norml2();
         exact_gradient.Add(-(exact_gradient*normal),normal);gradient-=exact_gradient;
         l2+=weight*error*error;truthnorm+=weight*truth*truth;gradienterror+=weight*(gradient*gradient);
         area+=weight;radiuserror+=weight*(r-1)*(r-1);radialmax=std::max(radialmax,std::abs(r-1));
         maximumerror=std::max(maximumerror,std::abs(error));mean+=weight*value;loadmean+=weight*(constant?1.:7*truth);minweight=std::min(minweight,tr->Weight());
      }
   }
   double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
   const double relative_residual=residual.Norml2()/B.Norml2();
   if(!solver.GetConverged()||!std::isfinite(l2)||minweight<=0||relative_residual>1e-9) {std::cerr<<"Solver/geometry check failed\n";return 3;}
   if(prefix!="-")ExportSurface(*mesh,u,prefix);
   std::cout<<"family,levels,p,g,extra,constant,elements,dofs,l2,relative_l2,h1_seminorm,energy_error,sampled_linf,area,relative_area_error,radial_rms,radial_max,mean_balance_defect,relative_residual,iterations,solve_and_measure_seconds\n"<<std::setprecision(17)
            <<family<<','<<levels<<','<<p<<','<<g<<','<<extra<<','<<constant<<','<<mesh->GetNE()<<','<<space.GetTrueVSize()<<','<<std::sqrt(l2)<<','<<std::sqrt(l2/truthnorm)<<','<<std::sqrt(gradienterror)<<','<<std::sqrt(l2+gradienterror)<<','<<maximumerror<<','<<area<<','<<std::abs(area-4*std::acos(-1.))/(4*std::acos(-1.))<<','<<std::sqrt(radiuserror/area)<<','<<radialmax<<','<<std::abs(mean-loadmean)<<','<<relative_residual<<','<<solver.GetNumIterations()<<','<<seconds<<'\n';
   return 0;
}
