// Independent physical finite-difference references live in the analysis script.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#define main preserved_direct_study_main
#include "direct_derivatives.cpp"
#undef main
#pragma clang diagnostic pop
#include "fe_mapping_jets.hpp"
int main(){
 std::cout<<std::setprecision(17);
 H1_FECollection fec(4,2);const auto *fe=fec.FiniteElementForGeometry(Geometry::SQUARE);
 auto basis=std::make_shared<fe_mapping::PolynomialBasis>(*fe);
 for(bool coupled:{false,true}){
  direct::Geometry analytic{.12,.18,.20,.06,-.03,.22,coupled?1.2:.7,coupled};
  DenseMatrix points(2,fe->GetDof());for(int i=0;i<fe->GetDof();++i){auto ip=fe->GetNodes().IntPoint(i);auto x=analytic.Map(direct::Jet(ip.x),direct::Jet(ip.y));for(int d=0;d<2;++d)points(d,i)=x[d].c[0];}
  IsoparametricTransformation tr;tr.SetFE(fe);tr.SetPointMat(points);fe_mapping::GeometryField field(tr,basis);
  IntegrationPoint ip;ip.Set2(.23,.31);double residual=0,h=.25;
  auto inverse=fe_mapping::Inverse([&](const fe_mapping::Jets&r){return field.Map(r);},ip,2,h,residual);
  auto r=inverse[0],s=inverse[1];auto u=r*r*r*s+.4*r*s*s+.7*r*r-.2*s;
  for(int i=1;i<direct::count;++i){int x=direct::ax[i],y=direct::ay[i];std::cout<<(coupled?"coupled":"separable")<<','<<x<<','<<y<<','<<u.c[fe_mapping::Index(x,y)]*direct::Factorial(x)*direct::Factorial(y)/std::pow(h,x+y)<<'\n';}
 }
 return 0;
}
