// Experiment-only curved derivatives from the production affine reference
// matrices. No analytic mapping, monomial conversion, or inverse-map series.
#pragma once
#include "mfem.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

// Defined by the unchanged production object in libofdg.a. This declaration
// keeps the experiment separate from the production public interface.
namespace ofdg { namespace detail {
std::vector<mfem::DenseMatrix> AssembleReferenceDerivatives(
    const mfem::FiniteElement &fe);
}}

namespace reused {
using Multi = std::array<int,3>;
inline int Degree(const Multi &a) { return a[0]+a[1]+a[2]; }
inline double Choose(int n,int k) {
    double v=1; for(int i=1;i<=k;++i) v=v*(n-k+i)/i; return v;
}
struct Split { int left,right; double binomial; };

// Same canonical parent rule as the production derivative graph, extended to
// the requested order (including known polynomial-zero states).
struct Plan {
    int dim,order;
    std::vector<Multi> indices;
    std::vector<int> parent,axis;
    std::vector<std::vector<Split>> splits;
    std::vector<std::array<int,3>> raised;
    Plan(int d,int m):dim(d),order(m) {
        for(int total=0;total<=m;++total)
            for(int x=total;x>=0;--x)
                for(int y=total-x;y>=0;--y) {
                    int z=total-x-y;
                    if(d==3 || z==0) indices.push_back({x,y,z});
                }
        const int n=Size(); parent.resize(n,-1); axis.resize(n,-1);
        splits.resize(n); raised.resize(n);
        for(int i=0;i<n;++i) {
            auto a=indices[i];
            if(i) {
                int d=dim-1; while(a[d]==0) --d;
                --a[d]; parent[i]=Index(a); axis[i]=d;
            }
            for(int d=0;d<3;++d) {
                a=indices[i]; ++a[d]; raised[i][d]=Index(a);
            }
            for(int j=0;j<n;++j) {
                Multi b{}; double factor=1; bool inside=true;
                for(int d=0;d<3;++d) {
                    if(indices[j][d]>indices[i][d]) { inside=false; break; }
                    b[d]=indices[i][d]-indices[j][d];
                    factor*=Choose(indices[i][d],indices[j][d]);
                }
                if(inside) splits[i].push_back({j,Index(b),factor});
            }
        }
    }
    int Size() const { return static_cast<int>(indices.size()); }
    int Index(const Multi &a) const {
        for(int i=0;i<Size();++i) if(indices[i]==a) return i;
        return -1;
    }
};

struct Operators {
    const mfem::FiniteElement &fe;
    const Plan &plan;
    std::vector<mfem::DenseMatrix> D;
    mfem::Vector constant,center_shape;
    Operators(const mfem::FiniteElement &element,const Plan &p)
        :fe(element),plan(p),D(ofdg::detail::AssembleReferenceDerivatives(fe)),
         constant(fe.GetDof()),center_shape(fe.GetDof()) {
        const int n=fe.GetDof(); mfem::DenseMatrix V(n); mfem::Vector shape(n),one(n);
        one=1.0;
        for(int q=0;q<n;++q) {
            fe.CalcShape(fe.GetNodes().IntPoint(q),shape);
            for(int j=0;j<n;++j) V(q,j)=shape(j);
        }
        mfem::DenseMatrixInverse solve(V); solve.Mult(one,constant);
        fe.CalcShape(mfem::Geometries.GetCenter(fe.GetGeomType()),center_shape);
    }
    bool Zero(int i) const {
        const auto &a=plan.indices[i]; const int p=fe.GetOrder();
        if(fe.GetGeomType()==mfem::Geometry::TRIANGLE ||
           fe.GetGeomType()==mfem::Geometry::TETRAHEDRON) return Degree(a)>p;
        for(int d=0;d<plan.dim;++d) if(a[d]>p) return true;
        return false;
    }
};

struct Field {
    const Operators &op;
    mfem::DenseMatrix original;
    std::vector<mfem::DenseMatrix> states;
    Field(const Operators &ops,const mfem::DenseMatrix &coeff,bool center)
        :op(ops),original(coeff),states(ops.plan.Size()) {
        states[0]=coeff;
        if(center) {
            for(int c=0;c<coeff.Width();++c) {
                double value=0;
                for(int i=0;i<coeff.Height();++i) value+=op.center_shape(i)*coeff(i,c);
                for(int i=0;i<coeff.Height();++i)
                    states[0](i,c)-=value*op.constant(i);
            }
        }
        for(int i=1;i<op.plan.Size();++i) {
            states[i].SetSize(coeff.Height(),coeff.Width());
            if(op.Zero(i)) states[i]=0.0;
            else mfem::Mult(op.D[op.plan.axis[i]],states[op.plan.parent[i]],states[i]);
        }
    }
    mfem::DenseMatrix Evaluate(const mfem::IntegrationPoint &ip) const {
        mfem::Vector shape(op.fe.GetDof()); op.fe.CalcShape(ip,shape);
        mfem::DenseMatrix values(op.plan.Size(),original.Width());
        for(int i=0;i<op.plan.Size();++i) {
            const auto &state=i ? states[i] : original;
            for(int c=0;c<original.Width();++c) {
                double v=0; for(int k=0;k<shape.Size();++k) v+=shape(k)*state(k,c);
                values(i,c)=v;
            }
        }
        return values;
    }
};

inline mfem::DenseMatrix GeometryCoefficients(mfem::IsoparametricTransformation &map) {
    const auto &p=map.GetPointMat(); mfem::DenseMatrix c(p.Width(),p.Height());
    for(int i=0;i<p.Width();++i) for(int d=0;d<p.Height();++d) c(i,d)=p(d,i);
    return c;
}

using Matrix=std::array<double,9>;
inline Matrix Product(const Matrix &a,const Matrix &b,int dim) {
    Matrix c{};
    for(int i=0;i<dim;++i) for(int j=0;j<dim;++j)
        for(int k=0;k<dim;++k) c[i*dim+j]+=a[i*dim+k]*b[k*dim+j];
    return c;
}

struct Geometry {
    const Plan &plan;
    double h;
    mfem::DenseMatrix derivatives;
    std::vector<Matrix> K,B;
    double inverse_residual=0;
    Geometry(const Field &field,const mfem::IntegrationPoint &ip,double scale)
        :plan(field.op.plan),h(scale),derivatives(field.Evaluate(ip)),
         K(plan.Size()),B(plan.Size()) {
        const int dim=plan.dim;
        for(int b=0;b<plan.Size();++b) if(Degree(plan.indices[b])<plan.order)
            for(int i=0;i<dim;++i) for(int a=0;a<dim;++a)
                K[b][i*dim+a]=derivatives(plan.raised[b][a],i)/h;
        mfem::DenseMatrix k0(dim),b0(dim);
        for(int i=0;i<dim;++i) for(int j=0;j<dim;++j) k0(i,j)=K[0][i*dim+j];
        MFEM_VERIFY(std::abs(k0.Det())>1e-12,"Singular scaled geometry Jacobian");
        mfem::CalcInverse(k0,b0);
        for(int i=0;i<dim;++i) for(int j=0;j<dim;++j) B[0][i*dim+j]=b0(i,j);
        for(int b=1;b<plan.Size();++b) if(Degree(plan.indices[b])<plan.order) {
            Matrix rhs{};
            for(const auto &s:plan.splits[b]) if(s.left) {
                auto p=Product(K[s.left],B[s.right],dim);
                for(int k=0;k<dim*dim;++k) rhs[k]+=s.binomial*p[k];
            }
            B[b]=Product(B[0],rhs,dim); for(double &v:B[b]) v=-v;
        }
        // Check every computed derivative of K B = I, scaled by the sum of
        // term magnitudes, so high derivative magnitudes do not hide failures.
        for(int b=0;b<plan.Size();++b) if(Degree(plan.indices[b])<plan.order) {
            Matrix sum{},magnitude{};
            for(const auto &s:plan.splits[b]) {
                auto p=Product(K[s.left],B[s.right],dim);
                for(int k=0;k<dim*dim;++k) {
                    double v=s.binomial*p[k]; sum[k]+=v; magnitude[k]+=std::abs(v);
                }
            }
            for(int i=0;i<dim;++i) for(int j=0;j<dim;++j) {
                const int k=i*dim+j; double target=(b==0 && i==j)?1:0;
                inverse_residual=std::max(inverse_residual,
                    std::abs(sum[k]-target)/std::max(1.0,magnitude[k]));
            }
        }
    }
};

// G(alpha,beta) = reference derivative beta of the scaled physical derivative
// alpha. Only |alpha|+|beta| <= m is used. These are point values, not FE states.
inline mfem::DenseMatrix Physical(const Geometry &geometry,const Field &solution,
                                 const mfem::IntegrationPoint &ip) {
    const auto &p=geometry.plan; const int n=p.Size(),dim=p.dim;
    auto reference=solution.Evaluate(ip); const int components=reference.Width();
    std::vector<double> g(static_cast<size_t>(n)*n*components,0.0);
    auto at=[&](int a,int b,int c)->double& {return g[(a*n+b)*components+c];};
    for(int b=0;b<n;++b) for(int c=0;c<components;++c) at(0,b,c)=reference(b,c);
    for(int a=1;a<n;++a) {
        int parent=p.parent[a],direction=p.axis[a],remaining=p.order-Degree(p.indices[a]);
        for(int b=0;b<n && Degree(p.indices[b])<=remaining;++b)
            for(const auto &split:p.splits[b]) for(int k=0;k<dim;++k) {
                const int raised=p.raised[split.right][k];
                const double factor=split.binomial*geometry.B[split.left][k*dim+direction];
                for(int c=0;c<components;++c) at(a,b,c)+=factor*at(parent,raised,c);
            }
    }
    mfem::DenseMatrix result(n,components);
    for(int a=0;a<n;++a) for(int c=0;c<components;++c)
        result(a,c)=at(a,0,c)/std::pow(geometry.h,Degree(p.indices[a]));
    return result;
}
}
