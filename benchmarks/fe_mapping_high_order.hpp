// Experiment-only arbitrary finite-order Taylor evaluator. Instantiate order 5
// for the main target and order 10 only for the bounded stress probes.
#pragma once
#include "fe_mapping_jets.hpp"
namespace high_order {
template<int D,int O> struct Layout {
 static constexpr int Count(){int n=1;for(int i=1;i<=D;++i)n=n*(O+i)/i;return n;}
 static constexpr int N=Count();
 static const auto& Powers(){static const auto v=[](){std::array<std::array<int,3>,N> a{};int i=0;
  for(int total=0;total<=O;++total)for(int x=total;x>=0;--x)for(int y=total-x;y>=0;--y){int z=total-x-y;if(D==3||z==0)a[i++]={x,y,z};}return a;}();return v;}
 static int Index(int x,int y,int z=0){
  if(x+y+z>O||x<0||y<0||z<0||(D==2&&z))return -1;
  static const auto indices=[](){std::array<int,(O+1)*(O+1)*(O+1)> a{};a.fill(-1);int i=0;for(auto p:Powers())a[(p[0]*(O+1)+p[1])*(O+1)+p[2]]=i++;return a;}();
  return indices[(x*(O+1)+y)*(O+1)+z];
 }
 static int Degree(int i){auto p=Powers()[i];return p[0]+p[1]+p[2];}
 static double Factor(int i){double v=1;for(int a:Powers()[i])for(int k=2;k<=a;++k)v*=k;return v;}
 static const auto& Products(){static const auto v=[](){std::vector<std::array<int,3>> p;for(int i=0;i<N;++i)for(int j=0;j<N;++j){auto a=Powers()[i],b=Powers()[j];int k=Index(a[0]+b[0],a[1]+b[1],a[2]+b[2]);if(k>=0)p.push_back({i,j,k});}return p;}();return v;}
};
template<int D,int O> struct Jet {
 using L=Layout<D,O>;std::array<double,L::N> c{};Jet(double x=0){c[0]=x;}
 friend Jet operator+(Jet a,const Jet&b){for(int i=0;i<L::N;++i)a.c[i]+=b.c[i];return a;}
 friend Jet operator-(Jet a,const Jet&b){for(int i=0;i<L::N;++i)a.c[i]-=b.c[i];return a;}
 friend Jet operator*(double x,Jet a){for(auto &v:a.c)v*=x;return a;}
 friend Jet operator*(const Jet&a,const Jet&b){Jet c;for(auto p:L::Products())c.c[p[2]]+=a.c[p[0]]*b.c[p[1]];return c;}
};
template<int D,int O> using Jets=std::array<Jet<D,O>,D>;
template<int D,int O> Jets<D,O> Coordinates(const mfem::IntegrationPoint &ip,bool vars=false){
 Jets<D,O> r;double x[3]={ip.x,ip.y,ip.z};for(int d=0;d<D;++d){r[d]=Jet<D,O>(x[d]);if(vars){int p[3]={};p[d]=1;r[d].c[Layout<D,O>::Index(p[0],p[1],p[2])]=1;}}return r;
}
template<int D,int O> std::vector<Jet<D,O>> Basis(const fe_mapping::PolynomialBasis &basis,const Jets<D,O>&r){
 using J=Jet<D,O>;std::array<std::vector<J>,D> xp;
 for(int d=0;d<D;++d){xp[d].resize(basis.degree+1);xp[d][0]=J(1);for(int k=1;k<=basis.degree;++k)xp[d][k]=xp[d][k-1]*(2*r[d]-J(1));}
 std::vector<J> out;for(auto p:basis.exponents){J v(1);for(int d=0;d<D;++d)v=v*xp[d][p[d]];out.push_back(v);}return out;
}
template<int D,int O> Jet<D,O> Evaluate(const std::vector<Jet<D,O>>&basis,const mfem::Vector&coeff){Jet<D,O> v;for(int i=0;i<coeff.Size();++i)v=v+coeff(i)*basis[i];return v;}
template<int D,int O> Jets<D,O> Map(const fe_mapping::GeometryField &g,const Jets<D,O>&r){auto values=Basis<D,O>(*g.basis,r);Jets<D,O> out;for(int d=0;d<D;++d)out[d]=Evaluate<D,O>(values,g.coefficients[d]);return out;}
template<int D,int O,class Mapping> Jets<D,O> Inverse(const Mapping &mapping,const mfem::IntegrationPoint&ip,double h,double &residual){
 using J=Jet<D,O>;using L=Layout<D,O>;auto g=mapping(Coordinates<D,O>(ip,true));mfem::DenseMatrix a(D),inverse(D);Jets<D,O> target;
 for(int d=0;d<D;++d){target[d]=J(g[d].c[0]);for(int k=0;k<D;++k){int p[3]={};p[k]=1;int i=L::Index(p[0],p[1],p[2]);a(d,k)=g[d].c[i];if(d==k)target[d].c[i]=h;}}
 mfem::CalcInverse(a,inverse);auto r=Coordinates<D,O>(ip);
 for(int iteration=0;iteration<O+1;++iteration){auto v=mapping(r);auto next=r;for(int d=0;d<D;++d)for(int k=0;k<D;++k)next[d]=next[d]-inverse(d,k)*(v[k]-target[k]);r=next;}
 auto v=mapping(r);for(int d=0;d<D;++d)for(int i=0;i<L::N;++i)residual=std::max(residual,std::abs(v[d].c[i]-target[d].c[i])/h);return r;
}
}
