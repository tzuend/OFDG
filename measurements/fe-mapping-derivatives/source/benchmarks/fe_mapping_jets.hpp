// Experiment-only degree-three Taylor arithmetic and polynomial FE fields.
// Geometry evaluation has no access to fixture formulas or physical derivatives.
#pragma once
#include "mfem.hpp"
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <chrono>

namespace fe_mapping {
using namespace mfem;
constexpr int N=20;
inline const std::array<std::array<int,3>,N> powers={{
 {0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},{1,1,0},{1,0,1},
 {0,2,0},{0,1,1},{0,0,2},{3,0,0},{2,1,0},{2,0,1},{1,2,0},
 {1,1,1},{1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}}};
inline int Index(int x,int y,int z=0) {
 for(int i=0;i<N;++i)if(powers[i]==std::array<int,3>{x,y,z})return i;
 return -1;
}
inline int Degree(int i){return powers[i][0]+powers[i][1]+powers[i][2];}
inline double Factorial(int k){return k<2?1:(k==2?2:6);}
inline double Factor(int i){double f=1;for(int k:powers[i])f*=Factorial(k);return f;}
struct Jet {
 std::array<double,N> c{};
 Jet(double v=0){c[0]=v;}
};
inline Jet operator+(Jet a,const Jet& b){for(int i=0;i<N;++i)a.c[i]+=b.c[i];return a;}
inline Jet operator-(Jet a,const Jet& b){for(int i=0;i<N;++i)a.c[i]-=b.c[i];return a;}
inline Jet operator*(double s,Jet a){for(double &v:a.c)v*=s;return a;}
inline Jet operator*(const Jet&a,const Jet&b){
 static const auto products=[](){std::vector<std::array<int,3>> v;
  for(int i=0;i<N;++i)for(int j=0;j<N;++j){int k=Index(powers[i][0]+powers[j][0],powers[i][1]+powers[j][1],powers[i][2]+powers[j][2]);if(k>=0)v.push_back({i,j,k});}return v;}();
 Jet r;for(auto ijk:products)r.c[ijk[2]]+=a.c[ijk[0]]*b.c[ijk[1]];return r;
}
using Jets=std::array<Jet,3>;
inline Jets Coordinates(const IntegrationPoint &ip,int dim,bool variables=false){
 Jets r={Jet(ip.x),Jet(ip.y),Jet(dim==3?ip.z:0)};
 if(variables)for(int d=0;d<dim;++d)r[d].c[1+d]=1;
 return r;
}
// Change of polynomial basis only, not a projection. Retain every tensor term,
// including terms above total degree three before Taylor evaluation at a point.
struct PolynomialBasis {
 int dim,degree;
 std::vector<std::array<int,3>> exponents;
 DenseMatrix coefficients;
 PolynomialBasis(const FiniteElement &fe):dim(fe.GetDim()),degree(fe.GetOrder()) {
  const bool triangle=fe.GetGeomType()==Geometry::TRIANGLE;
  MFEM_VERIFY(dim==2||dim==3,"Only 2D/3D volume polynomial FE are supported");
  for(int k=0;k<=(dim==3?degree:0);++k)for(int j=0;j<=degree;++j)for(int i=0;i<=degree;++i)
   if(!triangle||i+j<=degree)exponents.push_back({i,j,k});
  int nd=fe.GetDof();MFEM_VERIFY(int(exponents.size())==nd,"Unsupported polynomial space");
  DenseMatrix v(nd),samples(nd);Vector shape(nd);
  for(int row=0;row<nd;++row){auto ijk=exponents[row];IntegrationPoint ip;
   ip.Set3(double(ijk[0])/degree,double(ijk[1])/degree,double(ijk[2])/degree);
   fe.CalcShape(ip,shape);double x[3]={2*ip.x-1,2*ip.y-1,2*ip.z-1};
   for(int col=0;col<nd;++col){v(row,col)=1;for(int d=0;d<dim;++d)v(row,col)*=std::pow(x[d],exponents[col][d]);samples(row,col)=shape(col);}
  }
  DenseMatrixInverse inverse(v);coefficients.SetSize(nd);Vector rhs(nd),sol(nd);
  for(int j=0;j<nd;++j){for(int i=0;i<nd;++i)rhs(i)=samples(i,j);inverse.Mult(rhs,sol);for(int i=0;i<nd;++i)coefficients(i,j)=sol(i);}
 }
 std::vector<Jet> Evaluate(const Jets&r)const{
  std::array<std::vector<Jet>,3> xp;
  for(int d=0;d<dim;++d){xp[d].resize(degree+1);xp[d][0]=Jet(1);for(int k=1;k<=degree;++k)xp[d][k]=xp[d][k-1]*(2*r[d]-Jet(1));}
  std::vector<Jet> out;out.reserve(exponents.size());
  for(auto ijk:exponents){Jet v=xp[0][ijk[0]]*xp[1][ijk[1]];if(dim==3)v=v*xp[2][ijk[2]];out.push_back(v);}return out;
 }
};
inline Jet Evaluate(const std::vector<Jet>&basis,const Vector&coeff){Jet v;for(int i=0;i<coeff.Size();++i)v=v+coeff(i)*basis[i];return v;}
struct GeometryField {
 std::shared_ptr<const PolynomialBasis> basis;
 std::vector<Vector> coefficients;
 GeometryField(IsoparametricTransformation &map,std::shared_ptr<const PolynomialBasis> b):basis(std::move(b)){
  const auto &p=map.GetPointMat();MFEM_VERIFY(p.Height()==basis->dim,"Volume geometry required");
  for(int d=0;d<p.Height();++d){Vector original(p.Width()),converted(p.Width());for(int j=0;j<p.Width();++j)original(j)=p(d,j);basis->coefficients.Mult(original,converted);coefficients.push_back(converted);}
 }
 Jets Map(const Jets&r)const{auto values=basis->Evaluate(r);Jets result;for(int d=0;d<basis->dim;++d)result[d]=Evaluate(values,coefficients[d]);return result;}
};
// Formal local inverse, resolving a Taylor degree per correction. The physical
// increments are scaled by h to avoid conflating roundoff with element size.
template<class Map> Jets Inverse(const Map &map,const IntegrationPoint &ip,int dim,double h,double &residual){
 auto g=map(Coordinates(ip,dim,true));DenseMatrix j(dim),inverse(dim);
 for(int d=0;d<dim;++d)for(int k=0;k<dim;++k)j(d,k)=g[d].c[1+k];
 CalcInverse(j,inverse);Jets target=g;for(int d=0;d<dim;++d){target[d]=Jet(g[d].c[0]);target[d].c[d+1]=h;}
 auto r=Coordinates(ip,dim);
 for(int iteration=0;iteration<4;++iteration){auto value=map(r);Jets next=r;for(int d=0;d<dim;++d)for(int k=0;k<dim;++k)next[d]=next[d]-inverse(d,k)*(value[k]-target[k]);r=next;}
 auto value=map(r);for(int d=0;d<dim;++d)for(int k=0;k<N;++k)residual=std::max(residual,std::abs(value[d].c[k]-target[d].c[k])/h);
 return r;
}
inline double Seconds(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
}
