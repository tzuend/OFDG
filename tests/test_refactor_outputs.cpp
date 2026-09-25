// Emit full coefficient vectors for before/after refactor comparisons.
#include "mfem.hpp"
#include "../src/ofdg.hpp"
#include "../src/oedg_2024.hpp"
#include "../src/kxrcf.hpp"
#include <fstream>
#include <iomanip>
#include <string>
using namespace mfem;

void Bend(const Vector &x, Vector &y)
{
    y=x;
    y(0)+=0.18*x(0)*(1-x(0))*x(1);
    y(1)+=0.12*x(1)*(1-x(1))*x(0);
    if(x.Size()==3)y(2)+=0.1*x(2)*(1-x(2))*x(0);
}
void Write(std::ostream &out,const std::string &name,const Vector &v)
{
    out<<name<<' '<<v.Size();
    for(int i=0;i<v.Size();++i)out<<' '<<std::hexfloat<<v(i);
    out<<'\n';
}
void Check(std::ostream &out,int dim,Element::Type type,int degree,int basis,bool curved,int ranks,int rank)
{
    Mesh mesh=dim==1?Mesh::MakeCartesian1D(4):dim==2?
        Mesh::MakeCartesian2D(2,2,type,true):Mesh::MakeCartesian3D(2,1,1,type);
    if(curved){mesh.SetCurvature(3);mesh.Transform(Bend);}
    DG_FECollection fec(degree,dim,basis);
    FiniteElementSpace serial_space(&mesh,&fec,2,Ordering::byNODES);
    GridFunction serial_state(&serial_space);
    VectorFunctionCoefficient function(2,[](const Vector &x,Vector &u){u(0)=2+.2*x(0)+.13*std::sin(3*x(0));u(1)=1+.15*x(0)*x(0);if(x.Size()>1){u(0)+=.2*x(1);u(1)-=.1*x(1);}});
    serial_state.ProjectCoefficient(function);
    Array<int>dofs;
    for(int e=0;e<mesh.GetNE();++e){serial_space.GetElementVDofs(e,dofs);for(int j=0;j<dofs.Size();++j)serial_state(dofs[j])+=.01*(e%3)*(j+1)/dofs.Size();}
    Array<int>partition(mesh.GetNE());for(int e=0;e<mesh.GetNE();++e)partition[e]=e%ranks;
    // A one-cell 3D mesh cannot be partitioned onto two nonempty ranks.
    if(ranks>mesh.GetNE())return;
    std::unique_ptr<ParMesh> pm;
    std::unique_ptr<ParFiniteElementSpace> ps;
    std::unique_ptr<ParGridFunction> pu;
    FiniteElementSpace *space=&serial_space;const Vector *state=&serial_state;
    if(ranks>1){pm=std::make_unique<ParMesh>(MPI_COMM_WORLD,mesh,partition.GetData());ps=std::make_unique<ParFiniteElementSpace>(pm.get(),&fec,2,Ordering::byNODES);pu=std::make_unique<ParGridFunction>(pm.get(),&serial_state,partition.GetData());space=ps.get();state=pu.get();}
    out<<"case "<<dim<<' '<<int(type)<<' '<<degree<<' '<<basis<<' '<<curved<<' '<<rank<<'\n';
    VectorFunctionCoefficient velocity(dim,[](const Vector &x,Vector &v){for(int d=0;d<v.Size();++d)v(d)=.7+.1*d+.2*x(d);});
    auto physics=std::make_shared<ofdg::AdvectionFacePhysics>(&velocity);
    ofdg::OFDG filter(space,basis,physics);
    ofdg::KXRCFIndicator indicator(space,physics,.01);
    Array<bool>active;Vector values,result;DenseMatrix components;
    indicator.Compute(*state,active,&values,&components);Write(out,"indicator",values);
    for(int c=0;c<components.Width();++c){Vector column(components.GetColumn(c),components.Height());Write(out,"component",column);}
    for(int masked=0;masked<2;++masked){
        filter.ComputeStabilization(*state,result,masked?&active:nullptr);Write(out,"residual",result);
        filter.CompDecay(*state,result,.037,masked?&active:nullptr);Write(out,"decay",result);
    }
    for(int e=0;e<active.Size();++e)active[e]=(e%2==0);
    filter.ComputeStabilization(*state,result,&active);Write(out,"partial-residual",result);
    filter.CompDecay(*state,result,.037,&active);Write(out,"partial-decay",result);
    // Exercise inactive cells even when this indicator selects every element.
    active=false;filter.CompDecay(*state,result,.037,&active);Write(out,"inactive",result);
    if(!curved){ofdg::OEDG2024 oedg(space,basis,physics);oedg.ComputeStabilization(*state,result);Write(out,"oedg-residual",result);oedg.CompDecay(*state,result,.037);Write(out,"oedg-decay",result);}
}
int main(int argc,char **argv)
{
    Mpi::Init(argc,argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if(argc!=2)return 2;
    std::ofstream out(std::string(argv[1])+"-rank"+std::to_string(rank)+".txt");
    for(int basis:{BasisType::GaussLobatto,BasisType::GaussLegendre}){
        for(int p:{1,2,3})Check(out,1,Element::SEGMENT,p,basis,false,ranks,rank);
        for(auto type:{Element::TRIANGLE,Element::QUADRILATERAL})for(int p:{1,2,3})for(bool curved:{false,true})Check(out,2,type,p,basis,curved,ranks,rank);
        for(auto type:{Element::TETRAHEDRON,Element::HEXAHEDRON,Element::WEDGE})for(bool curved:{false,true})Check(out,3,type,2,basis,curved,ranks,rank);
    }
}
