#include "mfem.hpp"
#include "../src/curved_geometry.hpp"
#include "fe_mapping_high_order.hpp"
#include "reused_derivatives.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
using namespace mfem;
template<int D,int O> int Run(int n,bool curved,const std::string &dump){
 using J=high_order::Jet<D,O>;using L=high_order::Layout<D,O>;using Js=high_order::Jets<D,O>;
 double h=1./n;const double origin[3]={.2,.3,.15};
 Mesh mesh=D==2?Mesh::MakeCartesian2D(1,1,Element::QUADRILATERAL,true):Mesh::MakeCartesian3D(1,1,1,Element::HEXAHEDRON);
 auto mapping=[&](const Js&r){Js s;for(int d=0;d<D;++d)s[d]=J(origin[d])+h*r[d];if(curved){if constexpr(D==2){s[0]=s[0]+.7*s[0]*(J(1)-s[0]);s[1]=s[1]+.49*s[1]*(J(1)-s[1]);}else{J x=s[0]+.2*s[1]*(J(1)-s[1]),y=s[1]+.15*s[2]*(J(1)-s[2]);s[0]=x;s[1]=y;}}return s;};
 mesh.SetCurvature(D==2?4:2);mesh.Transform([&](const Vector&r,Vector&x){Js s;for(int d=0;d<D;++d)s[d]=J(r(d));auto y=mapping(s);for(int d=0;d<D;++d)x(d)=y[d].c[0];});
 double start=fe_mapping::Seconds();auto *map=dynamic_cast<IsoparametricTransformation*>(mesh.GetElementTransformation(0));
 auto gb=std::make_shared<fe_mapping::PolynomialBasis>(*map->GetFE());fe_mapping::GeometryField gf(*map,gb);
 DG_FECollection fec(5,D,BasisType::GaussLobatto);FiniteElementSpace space(&mesh,&fec);auto *fe=space.GetFE(0);fe_mapping::PolynomialBasis ub(*fe);int nd=fe->GetDof();
 auto polynomial=[](const double *s){double u=std::pow(s[0],5)+std::pow(s[1],5)+s[0]*s[0]*s[1]*s[1];if constexpr(D==3)u+=std::pow(s[2],5)+s[0]*s[1]*s[2];return u;};
 std::array<Vector,2> original={Vector(nd),Vector(nd)},coeff={Vector(nd),Vector(nd)};
 for(int i=0;i<nd;++i){auto &ip=fe->GetNodes().IntPoint(i);double s[3]={origin[0]+h*ip.x,origin[1]+h*ip.y,origin[2]+h*ip.z};original[0](i)=1;original[1](i)=polynomial(s);}
 for(int f=0;f<2;++f)ub.coefficients.Mult(original[f],coeff[f]);
 auto deriv=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,ofdg::detail::CurvedQuadratureOrder(5,*map));
 std::array<std::vector<Vector>,2> recursive;
 for(int f=0;f<2;++f){recursive[f].resize(L::N);recursive[f][0]=original[f];for(int i=1;i<L::N;++i){auto p=L::Powers()[i];int axis=D-1;while(!p[axis])--axis;--p[axis];int prev=L::Index(p[0],p[1],p[2]);recursive[f][i].SetSize(nd);deriv[axis].Mult(recursive[f][prev],recursive[f][i]);}}
 reused::Plan plan(D,O);reused::Operators gop(*map->GetFE(),plan),uop(*fe,plan);
 DenseMatrix uc(nd,2);for(int i=0;i<nd;++i)for(int f=0;f<2;++f)uc(i,f)=original[f](i);
 std::array<std::unique_ptr<reused::Field>,2> geometry,solution;
 for(int c=0;c<2;++c){geometry[c]=std::make_unique<reused::Field>(gop,reused::GeometryCoefficients(*map),c);solution[c]=std::make_unique<reused::Field>(uop,uc,c);}
 if(!dump.empty()) {std::ofstream out(dump);out<<std::setprecision(17)<<D<<' '<<n<<' '<<curved<<'\n';
  for(int which=0;which<2;++which){const auto &element=which?*fe:*map->GetFE();DenseMatrix cs=which?uc:reused::GeometryCoefficients(*map);out<<element.GetDof()<<' '<<cs.Width()<<'\n';
   for(int i=0;i<element.GetDof();++i){const auto &ip=element.GetNodes().IntPoint(i);out<<ip.x<<' '<<ip.y<<' '<<ip.z;for(int c=0;c<cs.Width();++c)out<<' '<<cs(i,c);out<<'\n';}}
 }
 double initialization=fe_mapping::Seconds()-start;
 const double points[3][3]={{.23,.31,.37},{.61,.47,.53},{0.,.4,.7}};
 std::cout<<"dimension,n,curved,capacity,point,xi,eta,zeta,x,y,z,field,dx,dy,dz,order,method,value,input_error,geometry_residual,inverse_residual,initialization_seconds,evaluation_seconds,reused_inverse_residual\n"<<std::setprecision(17);
 for(int point=0;point<3;++point){IntegrationPoint ip;ip.Set3(points[point][0],points[point][1],points[point][2]);start=fe_mapping::Seconds();double residual=0,ar=0;
  std::array<DenseMatrix,2> rv;double reuse_residual=0;
  for(int c=0;c<2;++c){reused::Geometry g(*geometry[c],ip,h);rv[c]=reused::Physical(g,*solution[c],ip);reuse_residual=std::max(reuse_residual,g.inverse_residual);}
  auto fm=[&](const Js&r){return high_order::Map<D,O>(gf,r);};auto r=high_order::Inverse<D,O>(fm,ip,h,residual),s=high_order::Inverse<D,O>(mapping,ip,h,ar);
  auto fv=high_order::Basis<D,O>(ub,r),av=high_order::Basis<D,O>(ub,s);Vector x,shape(nd);map->Transform(ip,x);fe->CalcShape(ip,shape);auto gx=fm(high_order::Coordinates<D,O>(ip));double gr=0;for(int d=0;d<D;++d)gr=std::max(gr,std::abs(gx[d].c[0]-x(d))/h);
  double local[3]={origin[0]+h*ip.x,origin[1]+h*ip.y,origin[2]+h*ip.z};
  for(int f=0;f<2;++f){auto u=high_order::Evaluate<D,O>(fv,coeff[f]),v=high_order::Evaluate<D,O>(av,coeff[f]);double elapsed=fe_mapping::Seconds()-start,input=std::abs(shape*original[f]-(f?polynomial(local):1));
   for(int i=1;i<L::N;++i){auto p=L::Powers()[i];for(int method=0;method<5;++method){double val=method>=3?rv[method-3](plan.Index(p),f):method==0?shape*recursive[f][i]:L::Factor(i)*(method==1?v.c[i]:u.c[i])/std::pow(h,L::Degree(i));
    std::cout<<D<<','<<n<<','<<curved<<','<<O<<','<<point<<','<<ip.x<<','<<ip.y<<','<<ip.z<<','<<x(0)<<','<<x(1)<<','<<(D==3?x(2):0)<<','<<(f?"quintic":"constant")<<','<<p[0]<<','<<p[1]<<','<<p[2]<<','<<L::Degree(i)<<','<<(method==0?"recursive":method==1?"direct":method==2?"fe_geometry":method==3?"reused":"reused_centered")<<','<<val<<','<<input<<','<<gr<<','<<std::max(ar,residual)<<','<<initialization<<','<<elapsed<<','<<reuse_residual<<'\n';}
   }
  }
 }
 return 0;
}
int main(int argc,char**argv){Mpi::Init(argc,argv);if(argc!=5&&argc!=6){std::cerr<<"dimension n curved capacity\n";return 2;}int d=std::stoi(argv[1]),n=std::stoi(argv[2]),c=std::stoi(argv[3]),o=std::stoi(argv[4]);
 if(d==2&&o==5)return Run<2,5>(n,c,argc==6?argv[5]:"");if(d==2&&o==10)return Run<2,10>(n,c,argc==6?argv[5]:"");if(d==3&&o==5)return Run<3,5>(n,c,argc==6?argv[5]:"");if(d==3&&o==10)return Run<3,10>(n,c,argc==6?argv[5]:"");return 2;
}
