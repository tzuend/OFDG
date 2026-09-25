#include "reused_derivatives.hpp"
#include <iostream>
#include <iomanip>
using namespace mfem;
// A valid non-partition-of-unity basis: the constant coefficient vector is not 1.
struct ScaledBasis:ScalarFiniteElement {
 const FiniteElement &base;
 ScaledBasis(const FiniteElement &f):ScalarFiniteElement(f.GetDim(),f.GetGeomType(),f.GetDof(),f.GetOrder(),f.Space()),base(f){Nodes=f.GetNodes();}
 double Scale(int i)const{return .75+.03125*(i%11);}
 void CalcShape(const IntegrationPoint &ip,Vector &s)const override{base.CalcShape(ip,s);for(int i=0;i<GetDof();++i)s(i)*=Scale(i);}
 void CalcDShape(const IntegrationPoint &ip,DenseMatrix &s)const override{base.CalcDShape(ip,s);for(int i=0;i<GetDof();++i)for(int d=0;d<GetDim();++d)s(i,d)*=Scale(i);}
};
double fact(int n){double a=1;for(int i=2;i<=n;++i)a*=i;return a;}
int main(int argc,char**argv){Mpi::Init(argc,argv);std::cout<<"dimension,basis,scaled,centered,max_reference_error,max_constant_error,max_q_error\n"<<std::setprecision(17);
 for(int dim:{2,3})for(int basis:{BasisType::GaussLobatto,BasisType::GaussLegendre,BasisType::Positive}){
 DG_FECollection collection(5,dim,basis);const auto &base=*collection.FiniteElementForGeometry(dim==2?Geometry::SQUARE:Geometry::CUBE);ScaledBasis scaled(base);
 for(int scale=0;scale<2;++scale){const FiniteElement &fe=scale?static_cast<const FiniteElement&>(scaled):base;reused::Plan plan(dim,5);reused::Operators op(fe,plan);int nd=fe.GetDof();DenseMatrix v(nd),values(nd,2),coeff(nd,2);Vector shape(nd);
 for(int i=0;i<nd;++i){const auto &ip=fe.GetNodes().IntPoint(i);fe.CalcShape(ip,shape);for(int j=0;j<nd;++j)v(i,j)=shape(j);values(i,0)=1;values(i,1)=std::pow(ip.x,5)+std::pow(ip.y,5)+(dim==3?std::pow(ip.z,5):0)+ip.x*ip.x*ip.y*ip.y;}
 DenseMatrixInverse solve(v);solve.Mult(values,coeff);
 for(int center=0;center<2;++center){reused::Field f(op,coeff,center);double err=0,constant=0,qerr=0;
 for(double t:{.13,.51,.89}){IntegrationPoint ip;ip.Set3(t,.37,.61);const double x[3]={ip.x,ip.y,ip.z};auto val=f.Evaluate(ip);fe.CalcShape(ip,shape);qerr=std::max(qerr,std::abs(shape*op.constant-1));
 for(int i=1;i<plan.Size();++i){auto a=plan.indices[i];double exact=0;for(int d=0;d<dim;++d){bool pure=true;for(int k=0;k<dim;++k)if(k!=d&&a[k])pure=false;if(pure)exact+=fact(5)/fact(5-a[d])*std::pow(x[d],5-a[d]);}if(a[0]<=2&&a[1]<=2&&a[2]==0)exact+=fact(2)/fact(2-a[0])*fact(2)/fact(2-a[1])*std::pow(x[0],2-a[0])*std::pow(x[1],2-a[1]);err=std::max(err,std::abs(val(i,1)-exact)/std::max(1.,std::abs(exact)));constant=std::max(constant,std::abs(val(i,0)));}}
 std::cout<<dim<<','<<basis<<','<<scale<<','<<center<<','<<err<<','<<constant<<','<<qerr<<'\n';if(err>1e-7||constant>1e-7||qerr>1e-12)return 1;}
 }}return 0;}
