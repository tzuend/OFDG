// Small volume-only comparison. The FE evaluator sees only stored mesh fields.
#include "mfem.hpp"
#include "../src/curved_geometry.hpp"
#include "fe_mapping_jets.hpp"
#include "reused_derivatives.hpp"
#include <iostream>
#include <iomanip>
#include <string>
using namespace mfem;
using namespace fe_mapping;
namespace {
struct Fixture {
 Vector origin;
 DenseMatrix affine;
 double a,b;
 Jets Map(const Jets &r)const{
  Jets q;for(int d=0;d<3;++d){q[d]=Jet(origin(d));for(int k=0;k<3;++k)q[d]=q[d]+affine(d,k)*r[k];}
  return {q[0]+a*q[1]*(Jet(1)-q[1]),q[1]+b*q[2]*(Jet(1)-q[2]),q[2]};
 }
};
// Explicit global inverse supplies the exact represented field independently of
// local formal inversion. Constants and sine have direct analytic derivatives.
Jets InverseGlobal(const Jets& x,double a,double b){
 Jet r=x[2],t=x[1]-b*r*(Jet(1)-r),s=x[0]-a*t*(Jet(1)-t);return {s,t,r};
}
Jet Represented(const Jets&q){return q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[0]*q[1]*q[2];}
double Exact(int field,const Vector&x,int index,double a,double b){
 if(field==0)return index==0?1:0;
 if(field==2){auto p=powers[index];return std::pow(2.,p[0])*std::pow(3.,p[1])*std::sin(2*x(0)+3*x(1)+x(2)+Degree(index)*std::acos(-1.)/2);}
 Jets q={Jet(x(0)),Jet(x(1)),Jet(x(2))};for(int d=0;d<3;++d)q[d].c[d+1]=1;
 return Factor(index)*Represented(InverseGlobal(q,a,b)).c[index];
}
std::array<std::array<double,N>,3> ExactValues(const Vector &x,double a,double b){
 std::array<std::array<double,N>,3> v{};v[0][0]=1;
 Jets q={Jet(x(0)),Jet(x(1)),Jet(x(2))};for(int d=0;d<3;++d)q[d].c[d+1]=1;
 auto u=Represented(InverseGlobal(q,a,b));for(int i=0;i<N;++i){v[1][i]=Factor(i)*u.c[i];v[2][i]=Exact(2,x,i,a,b);}return v;
}
std::string Word(int i){std::string s;for(int d=0;d<3;++d)s+=std::string(powers[i][d],"xyz"[d]);return s;}
struct ElementData {
 Fixture fixture;
 std::unique_ptr<GeometryField> geometry;
 std::array<std::unique_ptr<reused::Field>,2> rg,ru;
 std::array<Vector,3> original,polynomial;
 std::array<std::array<Vector,N>,3> recursive;
};
struct Row {int field,index,mode;double error=0,truth=0,linf=0,face=0,facetruth=0,jump=0,jumpabs=0,op=0,directnorm=0;};
}
int main(int argc,char**argv){
 Mpi::Init(argc,argv);
 if(argc!=5){std::cerr<<"n p curved extra\n";return 2;}
 int n=std::stoi(argv[1]),p=std::stoi(argv[2]),curved=std::stoi(argv[3]),extra=std::stoi(argv[4]);double h=1./n,a=curved?.2:0,b=curved?.15:0;
 Mesh mesh=Mesh::MakeCartesian3D(n,n,n,Element::HEXAHEDRON);
 std::vector<ElementData> data(mesh.GetNE());IntegrationPoint zero;zero.Set3(0,0,0);
 for(int e=0;e<mesh.GetNE();++e){auto *map=mesh.GetElementTransformation(e);map->SetIntPoint(&zero);map->Transform(zero,data[e].fixture.origin);data[e].fixture.affine=map->Jacobian();data[e].fixture.a=a;data[e].fixture.b=b;}
 mesh.SetCurvature(2);mesh.Transform([&](const Vector&s,Vector&x){x=s;x(0)+=a*s(1)*(1-s(1));x(1)+=b*s(2)*(1-s(2));});
 reused::Plan reuse_plan(3,3);
 double setup=Seconds();auto *iso=dynamic_cast<IsoparametricTransformation*>(mesh.GetElementTransformation(0));
 reused::Operators geometry_ops(*iso->GetFE(),reuse_plan);
 auto gb=std::make_shared<PolynomialBasis>(*iso->GetFE());
 DG_FECollection fec(p,3,BasisType::GaussLobatto);FiniteElementSpace space(&mesh,&fec);PolynomialBasis ub(*space.GetFE(0));
 reused::Operators solution_ops(*space.GetFE(0),reuse_plan);
 double reused_inverse_residual=0,reused_first_residual=0,reused_geometry_residual=0;
 std::vector<Row> rows;for(int f=0;f<3;++f)for(int i=1;i<N;++i)for(int mode=0;mode<5;++mode)rows.push_back({f,i,mode});
 double initialization=Seconds()-setup,evaluation=0,input[3]={},geom=0,jac=0,gd=0,residual=0,fe_residual=0,agreement=0,first=0,representation=0,facegeom=0,minjac=1e99,faces=0,interior=0;
 int qb=0;
 for(int e=0;e<mesh.GetNE();++e){setup=Seconds();auto &ed=data[e];auto *map=space.GetElementTransformation(e);auto *fe=space.GetFE(e);int nd=fe->GetDof();
  ed.geometry=std::make_unique<GeometryField>(*dynamic_cast<IsoparametricTransformation*>(map),gb);
  qb=ofdg::detail::CurvedQuadratureOrder(p,*map);auto D=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,qb);
  DenseMatrix mass(nd);mass=0.;std::array<Vector,3> rhs={Vector(nd),Vector(nd),Vector(nd)};for(auto &v:rhs)v=0.;
  Vector shape(nd),x;auto &ir=ofdg::detail::CurvedRule(fe->GetGeomType(),qb+16);
  for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);fe->CalcShape(ip,shape);double w=ip.weight*map->Weight();AddMult_a_VVt(w,shape,mass);for(int f=0;f<3;++f)rhs[f].Add(w*Exact(f,x,0,a,b),shape);}
  DenseMatrixInverse inverse(mass);
  for(int f=0;f<3;++f){ed.original[f].SetSize(nd);ed.polynomial[f].SetSize(nd);inverse.Mult(rhs[f],ed.original[f]);ub.coefficients.Mult(ed.original[f],ed.polynomial[f]);
   for(int i=1;i<N;++i){Vector v(ed.original[f]),z(nd);for(char axis:Word(i)){D[axis-'x'].Mult(v,z);v=z;}ed.recursive[f][i]=v;}}
  auto coordinates=reused::GeometryCoefficients(*dynamic_cast<IsoparametricTransformation*>(map));
  DenseMatrix coefficients(nd,3);for(int f=0;f<3;++f)for(int j=0;j<nd;++j)coefficients(j,f)=ed.original[f](j);
  for(int c=0;c<2;++c){ed.rg[c]=std::make_unique<reused::Field>(geometry_ops,coordinates,c);ed.ru[c]=std::make_unique<reused::Field>(solution_ops,coefficients,c);}
  initialization+=Seconds()-setup;
 }
 auto evaluate=[&](int e,const IntegrationPoint&ip,std::array<std::array<double,N>,3>&direct,std::array<std::array<double,N>,3>&feder,Vector&shape,std::array<DenseMatrix,2>&rv){
  double start=Seconds();auto &ed=data[e];auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);map->SetIntPoint(&ip);shape.SetSize(fe->GetDof());fe->CalcShape(ip,shape);
  auto af=[&](const Jets&r){return ed.fixture.Map(r);};auto ff=[&](const Jets&r){return ed.geometry->Map(r);};
  auto ar=Inverse(af,ip,3,h,residual),fr=Inverse(ff,ip,3,h,fe_residual);auto av=ub.Evaluate(ar),fv=ub.Evaluate(fr);
  auto ag=af(Coordinates(ip,3,true)),fg=ff(Coordinates(ip,3,true));Vector x;map->Transform(ip,x);minjac=std::min(minjac,map->Jacobian().Det());
  for(int d=0;d<3;++d){geom=std::max(geom,std::abs(fg[d].c[0]-x(d))/h);for(int k=0;k<N;++k)gd=std::max(gd,Factor(k)*std::abs(fg[d].c[k]-ag[d].c[k])/h);for(int k=0;k<3;++k)jac=std::max(jac,std::abs(fg[d].c[k+1]-map->Jacobian()(d,k))/h);}
  DenseMatrix ds(fe->GetDof(),3);fe->CalcPhysDShape(*map,ds);
  for(int f=0;f<3;++f){auto u=Evaluate(av,ed.polynomial[f]),v=Evaluate(fv,ed.polynomial[f]);representation=std::max(representation,std::abs(v.c[0]-shape*ed.original[f]));
   for(int i=0;i<N;++i){double scale=Factor(i)/std::pow(h,Degree(i));direct[f][i]=scale*u.c[i];feder[f][i]=scale*v.c[i];agreement=std::max(agreement,Factor(i)*std::abs(v.c[i]-u.c[i])/std::max(1.,std::abs(Factor(i)*u.c[i])));}
   for(int d=0;d<3;++d){double val=0;for(int k=0;k<fe->GetDof();++k)val+=ds(k,d)*ed.original[f](k);first=std::max(first,std::abs(val-feder[f][d+1])/std::max(1.,std::abs(val)));}}
  for(int c=0;c<2;++c){reused::Geometry g(*ed.rg[c],ip,h);rv[c]=reused::Physical(g,*ed.ru[c],ip);reused_inverse_residual=std::max(reused_inverse_residual,g.inverse_residual);
   for(int d=0;d<3;++d)for(int k=0;k<3;++k)reused_geometry_residual=std::max(reused_geometry_residual,std::abs(g.K[0][d*3+k]-map->Jacobian()(d,k)/h));
   for(int f=0;f<3;++f)for(int d=0;d<3;++d){double expected=0;for(int k=0;k<fe->GetDof();++k)expected+=ds(k,d)*ed.original[f](k);reused_first_residual=std::max(reused_first_residual,std::abs(rv[c](d+1,f)-expected)/std::max(1.,std::abs(expected)));}}
  evaluation+=Seconds()-start;
 };
 for(int e=0;e<mesh.GetNE();++e){auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);const auto &ir=ofdg::detail::CurvedRule(fe->GetGeomType(),qb+16+extra);Vector shape,x;
  for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);std::array<std::array<double,N>,3> ad,fd;std::array<DenseMatrix,2> rv;evaluate(e,ip,ad,fd,shape,rv);map->Transform(ip,x);double w=ip.weight*map->Weight();
   for(int f=0;f<3;++f){double err=shape*data[e].original[f]-Exact(f,x,0,a,b);input[f]+=w*err*err;}
   auto refs=ExactValues(x,a,b);for(auto &r:rows){double ref=refs[r.field][r.index],dv=ad[r.field][r.index],v=r.mode>=3?rv[r.mode-3](reuse_plan.Index(powers[r.index]),r.field):r.mode==2?fd[r.field][r.index]:(r.mode==1?dv:shape*data[e].recursive[r.field][r.index]),err=v-ref;
    r.error+=w*err*err;r.truth+=w*ref*ref;r.linf=std::max(r.linf,std::abs(err));r.op+=w*(v-dv)*(v-dv);r.directnorm+=w*dv*dv;}}
 }
 for(int face=0;face<mesh.GetNumFaces();++face){auto *tr=mesh.GetFaceElementTransformations(face);if(!tr)continue;auto &ir=IntRules.Get(tr->GetGeometryType(),qb+16+extra);
  for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);tr->SetAllIntPoints(&ip);Vector x;tr->Face->Transform(ip,x);double w=ip.weight*tr->Face->Weight();std::vector<double> left(rows.size());auto refs=ExactValues(x,a,b);
   for(int side=0;side<2;++side){int e=side?tr->Elem2No:tr->Elem1No;if(e<0)continue;faces+=w;auto &local=side?tr->GetElement2IntPoint():tr->GetElement1IntPoint();Vector shape;std::array<std::array<double,N>,3> ad,fd;std::array<DenseMatrix,2> rv;evaluate(e,local,ad,fd,shape,rv);
    auto fx=data[e].geometry->Map(Coordinates(local,3));for(int d=0;d<3;++d)facegeom=std::max(facegeom,std::abs(fx[d].c[0]-x(d))/h);
    for(size_t i=0;i<rows.size();++i){auto &r=rows[i];double v=r.mode>=3?rv[r.mode-3](reuse_plan.Index(powers[r.index]),r.field):r.mode==2?fd[r.field][r.index]:(r.mode==1?ad[r.field][r.index]:shape*data[e].recursive[r.field][r.index]),ref=refs[r.field][r.index];r.face+=w*(v-ref)*(v-ref);r.facetruth+=w*ref*ref;if(side){r.jump+=w*(v-left[i])*(v-left[i]);r.jumpabs+=w*std::abs(v-left[i]);}else left[i]=v;}if(side)interior+=w;}
  }
 }
 if(reused_inverse_residual>1e-10||reused_geometry_residual>1e-10||reused_first_residual>1e-8){std::cerr<<"Reuse 3D control failed "<<reused_inverse_residual<<' '<<reused_geometry_residual<<' '<<reused_first_residual<<'\n';return 5;}
 if(geom>1e-10||jac>1e-10||gd>1e-10||facegeom>1e-10||fe_residual>1e-10||agreement>1e-8||first>1e-8||minjac<=0){std::cerr<<"3D controls failed "<<geom<<' '<<jac<<' '<<gd<<' '<<facegeom<<' '<<fe_residual<<' '<<agreement<<' '<<first<<'\n';return 4;}
 std::cout<<"family,n,p,amplitude,basis,map,extra,field,derivative,method,error_l2,relative_l2,sampled_linf,truth_l2,face_rms,face_relative_l2,jump_rms,jump_mean_abs,input_error_l2,operator_error_l2,operator_relative_l2,inverse_residual,representation_residual,first_derivative_residual,geometry_residual,min_jacobian,elements,geometry_derivative_residual,geometry_jacobian_residual,fe_inverse_residual,method_agreement,fe_first_residual,face_geometry_residual,initialization_seconds,evaluation_seconds,reused_inverse_residual,reused_first_residual,reused_geometry_residual\n"<<std::setprecision(17);
 const char *fields[]={"constant","represented_polynomial","sine"};const char *methods[]={"recursive","direct","fe_geometry","reused","reused_centered"};
 for(auto &r:rows)std::cout<<"hex,"<<n<<','<<p<<','<<a<<",gll,shear,"<<extra<<','<<fields[r.field]<<','<<Word(r.index)<<','<<methods[r.mode]<<','<<std::sqrt(r.error)<<','<<(r.truth>1e-24?std::sqrt(r.error/r.truth):NAN)<<','<<r.linf<<','<<std::sqrt(r.truth)<<','<<std::sqrt(r.face/faces)<<','<<(r.facetruth>1e-24?std::sqrt(r.face/r.facetruth):NAN)<<','<<std::sqrt(r.jump/interior)<<','<<r.jumpabs/interior<<','<<std::sqrt(input[r.field])<<','<<std::sqrt(r.op)<<','<<(r.directnorm>1e-24?std::sqrt(r.op/r.directnorm):NAN)<<','<<residual<<','<<representation<<','<<first<<','<<geom<<','<<minjac<<','<<mesh.GetNE()<<','<<gd<<','<<jac<<','<<fe_residual<<','<<agreement<<','<<first<<','<<facegeom<<','<<initialization<<','<<evaluation<<','<<reused_inverse_residual<<','<<reused_first_residual<<','<<reused_geometry_residual<<'\n';
 return 0;
}
