// Isolated P1/Q1 surface-DG and OFDG adaptation on the exact unit sphere.
// Production volume OFDG is deliberately not modified or silently reused.
#include "mfem.hpp"
#include "../src/curved_geometry.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
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


struct SurfacePoint {
    Vector x=Vector(3), velocity=Vector(3);
    DenseMatrix J= DenseMatrix(3,2), inverse= DenseMatrix(2,3);
    double area;
};
SurfacePoint GeometryAt(ElementTransformation &tr,const IntegrationPoint &ip)
{
    tr.SetIntPoint(&ip);Vector y;tr.Transform(ip,y);double radius=y.Norml2();
    SurfacePoint s;s.x=y;s.x/=radius;const auto &base=tr.Jacobian();
    for(int a=0;a<2;++a){double normalpart=0;for(int i=0;i<3;++i)normalpart+=s.x(i)*base(i,a);for(int i=0;i<3;++i)s.J(i,a)=(base(i,a)-s.x(i)*normalpart)/radius;}
    double g00=0,g01=0,g11=0;
    for(int i=0;i<3;++i){g00+=s.J(i,0)*s.J(i,0);g01+=s.J(i,0)*s.J(i,1);g11+=s.J(i,1)*s.J(i,1);}
    double determinant=g00*g11-g01*g01;s.area=std::sqrt(determinant);
    for(int i=0;i<3;++i){s.inverse(0,i)=(g11*s.J(i,0)-g01*s.J(i,1))/determinant;s.inverse(1,i)=(-g01*s.J(i,0)+g00*s.J(i,1))/determinant;}
    s.velocity(0)=-s.x(1);s.velocity(1)=s.x(0);s.velocity(2)=0;
    return s;
}
double Exact(const Vector &x,double time,int problem)
{
    if(problem==2)return 1.;
    const double initial_x=std::cos(time)*x(0)+std::sin(time)*x(1);
    return problem==0?1.+.5*initial_x:(initial_x>.5?1.:0.);
}
struct Sample { Vector shape;SurfacePoint geometry;double weight; };
struct ElementData {
    Array<int> dofs;DenseMatrix mass,inverse;std::array<DenseMatrix,3> derivatives;
    Vector mean;double area=0,h=0;std::vector<Sample> samples;
};
struct FaceSample { Vector left,right,x;double weight; };
struct FaceData {int left,right;double length=0,speed=0;std::vector<FaceSample> samples;};
class SurfaceOperator : public TimeDependentOperator {
public:
    Mesh &mesh;FiniteElementSpace &space;std::vector<ElementData> elements;std::vector<FaceData> faces;
    SparseMatrix transport;mutable Vector work;
    double totalarea=0,minimum_h=1e100,flux_mismatch=0,constant_residual=0,filter_mean_defect=0;
    SurfaceOperator(Mesh &m,FiniteElementSpace &sp,int quadrature):TimeDependentOperator(sp.GetVSize()),mesh(m),space(sp),elements(m.GetNE()),transport(sp.GetVSize()),work(sp.GetVSize())
    {
        for(int e=0;e<mesh.GetNE();++e){auto &ed=elements[e];auto *fe=space.GetFE(e);auto *tr=space.GetElementTransformation(e);int nd=fe->GetDof();space.GetElementDofs(e,ed.dofs);ed.mass.SetSize(nd);ed.mass=0.;ed.mean.SetSize(nd);ed.mean=0.;std::array<DenseMatrix,3> rhs;for(auto &a:rhs){a.SetSize(nd);a=0.;}
            const auto &rule=ofdg::detail::CurvedRule(fe->GetGeomType(),quadrature);
            for(int q=0;q<rule.GetNPoints();++q){const auto &ip=rule.IntPoint(q);auto g=GeometryAt(*tr,ip);Vector shape(nd);fe->CalcShape(ip,shape);DenseMatrix ds(nd,2),grad(nd,3);fe->CalcDShape(ip,ds);mfem::Mult(ds,g.inverse,grad);double weight=ip.weight*g.area;ed.area+=weight;ed.mean.Add(weight,shape);AddMult_a_VVt(weight,shape,ed.mass);
                for(int i=0;i<nd;++i)for(int j=0;j<nd;++j){double convection=0;for(int a=0;a<3;++a){convection+=grad(i,a)*g.velocity(a);rhs[a](i,j)+=weight*shape(i)*grad(j,a);}transport.Add(ed.dofs[i],ed.dofs[j],weight*convection*shape(j));}
                ed.samples.push_back({shape,g,weight});
            }
            ed.inverse=ed.mass;ed.inverse.Invert();ed.mean/=ed.area;ed.h=std::sqrt(ed.area/Geometry::Volume[fe->GetGeomType()]);minimum_h=std::min(minimum_h,ed.h);totalarea+=ed.area;
            for(int a=0;a<3;++a){ed.derivatives[a].SetSize(nd);mfem::Mult(ed.inverse,rhs[a],ed.derivatives[a]);}
        }
        for(int face=0;face<mesh.GetNumFaces();++face){auto *tr=mesh.GetFaceElementTransformations(face);MFEM_VERIFY(tr&&tr->Elem2No>=0,"Sphere must have no boundary");FaceData fd;fd.left=tr->Elem1No;fd.right=tr->Elem2No;const auto &rule=IntRules.Get(Geometry::SEGMENT,quadrature);
            for(int q=0;q<rule.GetNPoints();++q){const auto &ip=rule.IntPoint(q);tr->SetAllIntPoints(&ip);std::array<Vector,2> shape;std::array<double,2> flux,length;Vector point;
                for(int side=0;side<2;++side){int e=side?fd.right:fd.left;auto &local=side?tr->GetElement2IntPoint():tr->GetElement1IntPoint();auto &loc=side?tr->Loc2.Transf:tr->Loc1.Transf;loc.SetIntPoint(&ip);const auto &edge=loc.Jacobian();auto g=GeometryAt(*(side?tr->Elem2:tr->Elem1),local);shape[side].SetSize(space.GetFE(e)->GetDof());space.GetFE(e)->CalcShape(local,shape[side]);
                    double nx=edge(1,0),ny=-edge(0,0);const auto &center=Geometries.GetCenter(space.GetFE(e)->GetGeomType());if(nx*(local.x-center.x)+ny*(local.y-center.y)<0){nx=-nx;ny=-ny;}
                    Vector contrav(2);g.inverse.Mult(g.velocity,contrav);flux[side]=g.area*(contrav(0)*nx+contrav(1)*ny);double len2=0;for(int d=0;d<3;++d)len2+=std::pow(g.J(d,0)*edge(0,0)+g.J(d,1)*edge(1,0),2);length[side]=std::sqrt(len2);if(!side)point=g.x;
                }
                flux_mismatch=std::max(flux_mismatch,std::abs(flux[0]+flux[1]));double beta=.5*(flux[0]-flux[1]),measure=.5*(length[0]+length[1]);fd.length+=ip.weight*measure;fd.speed=std::max(fd.speed,std::abs(beta)/measure);fd.samples.push_back({shape[0],shape[1],point,ip.weight*measure});
                int upwind=beta>=0?fd.left:fd.right;const Vector &up=shape[beta>=0?0:1];for(int side=0;side<2;++side){int e=side?fd.right:fd.left;double sign=side?1.:-1.;for(int i=0;i<shape[side].Size();++i)for(int j=0;j<up.Size();++j)transport.Add(elements[e].dofs[i],elements[upwind].dofs[j],sign*ip.weight*beta*shape[side](i)*up(j));}
            }
            for(auto &s:fd.samples)s.weight/=fd.length;faces.push_back(std::move(fd));
        }
        transport.Finalize();Vector one(height),zero(height);one=1.;Mult(one,zero);constant_residual=zero.Normlinf();
    }
    void Mult(const Vector &u,Vector &du) const override
    {
        transport.Mult(u,work);du.SetSize(height);Vector local,result;
        for(const auto &e:elements){work.GetSubVector(e.dofs,local);result.SetSize(local.Size());e.inverse.Mult(local,result);du.SetSubVector(e.dofs,result);}
    }
    void Project(Vector &u,int problem) const
    {
        u.SetSize(height);for(const auto &e:elements){Vector rhs(e.dofs.Size()),coeff(e.dofs.Size());rhs=0.;for(const auto &s:e.samples)rhs.Add(s.weight*Exact(s.geometry.x,0,problem),s.shape);e.inverse.Mult(rhs,coeff);u.SetSubVector(e.dofs,coeff);}
    }
    double Integral(const Vector &u) const
    {
        double value=0;Vector coeff;for(const auto &e:elements){u.GetSubVector(e.dofs,coeff);value+=e.area*(e.mean*coeff);}return value;
    }
    void Filter(Vector &u,double dt)
    {
        const double global_mean=Integral(u)/totalarea;double amplitude=0,scale=1;
        std::vector<Vector> local(elements.size());std::vector<std::array<Vector,3>> gradient(elements.size());std::vector<double> rates(elements.size(),0.);
        for(size_t e=0;e<elements.size();++e){const auto &ed=elements[e];u.GetSubVector(ed.dofs,local[e]);for(int a=0;a<3;++a){gradient[e][a].SetSize(local[e].Size());ed.derivatives[a].Mult(local[e],gradient[e][a]);}for(const auto &s:ed.samples){double v=s.shape*local[e];amplitude=std::max(amplitude,std::abs(v-global_mean));scale=std::max(scale,std::abs(v));}}
        for(const auto &f:faces)for(const auto &s:f.samples){amplitude=std::max(amplitude,std::abs(s.left*local[f.left]-global_mean));amplitude=std::max(amplitude,std::abs(s.right*local[f.right]-global_mean));}
        if(amplitude<=1e-12*scale)return;
        for(const auto &f:faces){double jump0=0,jump1=0;
            for(const auto &s:f.samples){jump0+=s.weight*std::abs(s.left*local[f.left]-s.right*local[f.right]);double norm1=0;
                // Ambient Cartesian components of the projected tangential
                // gradient: the same componentwise absolute sum as volume OFDG.
                for(int a=0;a<3;++a)norm1+=std::abs(s.left*gradient[f.left][a]-s.right*gradient[f.right][a]);jump1+=s.weight*norm1;
            }
            for(int e:{f.left,f.right})rates[e]+=f.speed*(.5*jump0/elements[e].h+1.5*jump1)/amplitude;
        }
        for(size_t e=0;e<elements.size();++e){double mean=elements[e].mean*local[e],decay=std::exp(-dt*rates[e]);for(int i=0;i<local[e].Size();++i)local[e](i)=mean+decay*(local[e](i)-mean);filter_mean_defect=std::max(filter_mean_defect,std::abs(elements[e].mean*local[e]-mean));u.SetSubVector(elements[e].dofs,local[e]);}
    }
};
void SaveSamples(SurfaceOperator &op,const Vector &u,double time,int problem,const std::string &path)
{
    std::ofstream out(path);out<<"element,triangle,vertex,x,y,z,u,exact,error\n"<<std::setprecision(17);int tid=0;
    for(int e=0;e<op.mesh.GetNE();++e){auto *fe=op.space.GetFE(e);auto *tr=op.space.GetElementTransformation(e);Vector coeff;u.GetSubVector(op.elements[e].dofs,coeff);int n=5;double h=1./n;
        auto emit=[&](double a,double b,double c,double d,double f,double g){double points[3][2]={{a,b},{c,d},{f,g}};for(int j=0;j<3;++j){IntegrationPoint ip;ip.Set2(points[j][0],points[j][1]);auto geo=GeometryAt(*tr,ip);Vector shape(coeff.Size());fe->CalcShape(ip,shape);double v=shape*coeff,t=Exact(geo.x,time,problem);out<<e<<','<<tid<<','<<j<<','<<geo.x(0)<<','<<geo.x(1)<<','<<geo.x(2)<<','<<v<<','<<t<<','<<v-t<<'\n';}++tid;};
        for(int j=0;j<n;++j)for(int i=0;i<n;++i){if(fe->GetGeomType()==Geometry::TRIANGLE){if(i+j>=n)continue;emit(i*h,j*h,(i+1)*h,j*h,i*h,(j+1)*h);if(i+j<n-1)emit((i+1)*h,j*h,(i+1)*h,(j+1)*h,i*h,(j+1)*h);}else{emit(i*h,j*h,(i+1)*h,j*h,(i+1)*h,(j+1)*h);emit(i*h,j*h,(i+1)*h,(j+1)*h,i*h,(j+1)*h);}}
    }
}
int main(int argc,char **argv)
{
    Mpi::Init(argc,argv);if(argc!=8){std::cerr<<"family levels method problem cfl quadrature output_prefix_or_dash\n";return 2;}
    const std::string family=argv[1],method=argv[3],prefix=argv[7];int level=std::stoi(argv[2]),problem=std::stoi(argv[4]),quadrature=std::stoi(argv[6]);double cfl=std::stod(argv[5]);
    auto mesh=SphereMesh(family=="tri"?0:1,1,level);DG_FECollection fec(1,2,BasisType::GaussLobatto);FiniteElementSpace space(mesh.get(),&fec);SurfaceOperator op(*mesh,space,quadrature);Vector u;op.Project(u,problem);double initialmass=op.Integral(u),initialmin=1e100,initialmax=-1e100;
    for(const auto &e:op.elements){Vector c;u.GetSubVector(e.dofs,c);for(const auto &s:e.samples){double v=s.shape*c;initialmin=std::min(initialmin,v);initialmax=std::max(initialmax,v);}}
    if(prefix!="-")SaveSamples(op,u,0,problem,prefix+"-initial.csv");
    const double finaltime=.5*std::acos(-1.);const int steps=std::ceil(finaltime/(cfl*op.minimum_h/3));double dt=finaltime/steps;RK4Solver solver;solver.Init(op);double time=0;
    for(int i=0;i<steps;++i){solver.Step(u,time,dt);if(method=="ofdg")op.Filter(u,dt);}
    double error=0,l1=0,truth=0,minimum=1e100,maximum=-1e100;
    for(int e=0;e<mesh->GetNE();++e){Vector coeff;u.GetSubVector(op.elements[e].dofs,coeff);auto *fe=space.GetFE(e);auto *tr=space.GetElementTransformation(e);const auto &rule=ofdg::detail::CurvedRule(fe->GetGeomType(),problem==1?std::max(60,quadrature+8):quadrature+8);
        for(int q=0;q<rule.GetNPoints();++q){const auto &ip=rule.IntPoint(q);auto g=GeometryAt(*tr,ip);Vector shape(coeff.Size());fe->CalcShape(ip,shape);double v=shape*coeff,t=Exact(g.x,finaltime,problem),w=ip.weight*g.area;error+=w*(v-t)*(v-t);l1+=w*std::abs(v-t);truth+=w*t*t;minimum=std::min(minimum,v);maximum=std::max(maximum,v);}
    }
    // Include both one-sided edge samples in the reported extrema.
    for(const auto &f:op.faces){Vector a,b;u.GetSubVector(op.elements[f.left].dofs,a);u.GetSubVector(op.elements[f.right].dofs,b);for(const auto &s:f.samples){double v1=s.left*a,v2=s.right*b;minimum=std::min({minimum,v1,v2});maximum=std::max({maximum,v1,v2});}}
    const double drift=std::abs(op.Integral(u)-initialmass);if(!std::isfinite(error)||error<0||drift>1e-8||op.constant_residual>1e-8||op.flux_mismatch>1e-9){std::cerr<<"Control failure "<<error<<' '<<drift<<' '<<op.constant_residual<<' '<<op.flux_mismatch<<'\n';return 3;}
    if(prefix!="-")SaveSamples(op,u,finaltime,problem,prefix+"-final.csv");
    std::cout<<"family,level,method,problem,cfl,quadrature,elements,dofs,steps,final_time,l2,relative_l2,l1,minimum,maximum,initial_minimum,initial_maximum,mass_drift,filter_mean_defect,constant_operator_residual,flux_mismatch,area\n"<<std::setprecision(17)<<family<<','<<level<<','<<method<<','<<problem<<','<<cfl<<','<<quadrature<<','<<mesh->GetNE()<<','<<space.GetVSize()<<','<<steps<<','<<time<<','<<std::sqrt(error)<<','<<std::sqrt(error/truth)<<','<<l1<<','<<minimum<<','<<maximum<<','<<initialmin<<','<<initialmax<<','<<drift<<','<<op.filter_mean_defect<<','<<op.constant_residual<<','<<op.flux_mismatch<<','<<op.totalarea<<'\n';return 0;
}
