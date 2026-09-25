// Isolated direct differentiation experiment; no production code is modified.
// A degree-three Taylor jet in physical increments is pulled back through the
// known polynomial geometry by formal inversion. No derivative is projected.
#include "mfem.hpp"
#include "../src/curved_geometry.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace mfem;
Vector Apply(const std::vector<DenseMatrix>& D,const Vector& u,const std::string& word){Vector v(u),z(u.Size());for(char c:word){D[c=='x'?0:1].Mult(v,z);v=z;}return v;}
double PowerDerivative(double s,int degree,int order){if(order>degree)return 0;double v=1;for(int j=0;j<order;++j)v*=degree-j;return v*std::pow(s,degree-order);}
// Exact derivatives of [T^{-1}(x)]^degree for T(s)=s+a*s*(1-s).
double InversePower(double x,double a,int degree,int order){double s=2*x/(1+a+std::sqrt((1+a)*(1+a)-4*a*x));double j=1+a-2*a*s;double d1=1/j,d2=2*a/std::pow(j,3),d3=12*a*a/std::pow(j,5);if(order==0)return PowerDerivative(s,degree,0);if(order==1)return PowerDerivative(s,degree,1)*d1;if(order==2)return PowerDerivative(s,degree,2)*d1*d1+PowerDerivative(s,degree,1)*d2;return PowerDerivative(s,degree,3)*d1*d1*d1+3*PowerDerivative(s,degree,2)*d1*d2+PowerDerivative(s,degree,1)*d3;}
double Exact(int field,const Vector& x,int nx,int ny,double a){double pi=std::acos(-1.);if(field==5)return PowerDerivative(x(0),2,nx)*PowerDerivative(x(1),0,ny)+PowerDerivative(x(0),1,nx)*PowerDerivative(x(1),1,ny)+PowerDerivative(x(0),0,nx)*PowerDerivative(x(1),2,ny);if(field==0)return std::pow(2.,nx)*std::pow(3.,ny)*std::sin(2*x(0)+3*x(1)+(nx+ny)*pi/2);if(field==1)return std::pow(8.,nx)*std::pow(6.,ny)*std::sin(8*x(0)+6*x(1)+(nx+ny)*pi/2);if(field==2)return std::pow(.5,ny)*std::exp(x(0)+.5*x(1));if(field==3)return (nx+ny==0)?1.:0.;return InversePower(x(0),a,2,nx)*InversePower(x(1),.7*a,0,ny)+InversePower(x(0),a,1,nx)*InversePower(x(1),.7*a,1,ny)+InversePower(x(0),a,0,nx)*InversePower(x(1),.7*a,2,ny);}


namespace direct {
constexpr int count=10;
const int ax[count]={0,1,0,2,1,0,3,2,1,0};
const int ay[count]={0,0,1,0,1,2,0,1,2,3};
int Index(int x,int y) { for(int i=0;i<count;++i) if(ax[i]==x && ay[i]==y)return i; return -1; }
struct Jet {
    std::array<double,count> c{};
    Jet(double value=0) { c[0]=value; }
};
Jet operator+(Jet a,const Jet &b) { for(int i=0;i<count;++i)a.c[i]+=b.c[i];return a; }
Jet operator-(Jet a,const Jet &b) { for(int i=0;i<count;++i)a.c[i]-=b.c[i];return a; }
Jet operator*(double a,Jet b) { for(double &v:b.c)v*=a;return b; }
Jet operator*(const Jet &a,const Jet &b) {
    Jet z;
    for(int i=0;i<count;++i)for(int j=0;j<count;++j) {
        int x=ax[i]+ax[j],y=ay[i]+ay[j];
        if(x+y<=3)z.c[Index(x,y)]+=a.c[i]*b.c[j];
    }
    return z;
}
struct Geometry {
    // Original Cartesian element expressed in its local reference coordinates.
    double s0,t0,sx,sy,tx,ty,amplitude;
    bool coupled;
    std::array<Jet,2> Map(const Jet &xi,const Jet &eta) const {
        Jet s=Jet(s0)+sx*xi+sy*eta,t=Jet(t0)+tx*xi+ty*eta;
        if(coupled) { Jet b=amplitude*s*(Jet(1)-s)*t*(Jet(1)-t);return {s+b,t+.7*b}; }
        return {s+amplitude*s*(Jet(1)-s),t+(.7*amplitude)*t*(Jet(1)-t)};
    }
    std::array<Jet,2> InverseJet(double xi,double eta,double h,double &residual) const {
        Jet r(xi),s(eta);r.c[1]=1;s.c[2]=1;
        auto g=Map(r,s);
        const double a=g[0].c[1],b=g[0].c[2],c=g[1].c[1],d=g[1].c[2],det=a*d-b*c;
        Jet targetx(g[0].c[0]),targety(g[1].c[0]);targetx.c[1]=h;targety.c[2]=h;
        r=Jet(xi);s=Jet(eta);
        // Fixed base-point Jacobian: each correction resolves one Taylor order.
        // Four passes include a redundant residual correction for roundoff.
        for(int it=0;it<4;++it) {
            auto value=Map(r,s);Jet ex=value[0]-targetx,ey=value[1]-targety;
            r=r-(d/det)*ex+(b/det)*ey;s=s+(c/det)*ex-(a/det)*ey;
        }
        auto value=Map(r,s);
        for(int i=0;i<count;++i)residual=std::max(residual,std::max(std::abs(value[0].c[i]-targetx.c[i]),std::abs(value[1].c[i]-targety.c[i]))/h);
        return {r,s};
    }
};
struct PolynomialBasis {
    std::vector<std::array<int,2>> powers;
    DenseMatrix coefficients; // monomial coefficients of each original FE basis function
    PolynomialBasis(const FiniteElement &fe,bool triangle) {
        const int p=fe.GetOrder(),nd=fe.GetDof();
        for(int j=0;j<=p;++j)for(int i=0;i<=p;++i)if(!triangle||i+j<=p)powers.push_back({i,j});
        MFEM_VERIFY(int(powers.size())==nd,"Polynomial dimension mismatch");
        DenseMatrix v(nd),samples(nd);Vector shape(nd);int row=0;
        for(int j=0;j<=p;++j)for(int i=0;i<=p;++i)if(!triangle||i+j<=p) {
            IntegrationPoint ip;ip.Set2(double(i)/p,double(j)/p);fe.CalcShape(ip,shape);
            for(int k=0;k<nd;++k) {v(row,k)=std::pow(2*ip.x-1,powers[k][0])*std::pow(2*ip.y-1,powers[k][1]);samples(row,k)=shape(k);}
            ++row;
        }
        DenseMatrixInverse inverse(v);coefficients.SetSize(nd);Vector rhs(nd),sol(nd);
        for(int j=0;j<nd;++j){for(int i=0;i<nd;++i)rhs(i)=samples(i,j);inverse.Mult(rhs,sol);for(int i=0;i<nd;++i)coefficients(i,j)=sol(i);}
    }
    std::vector<Jet> Evaluate(const Jet &r,const Jet &s,int p) const {
        std::vector<Jet> xp(p+1),yp(p+1),out;xp[0]=Jet(1);yp[0]=Jet(1);
        for(int i=1;i<=p;++i){xp[i]=xp[i-1]*(2*r-Jet(1));yp[i]=yp[i-1]*(2*s-Jet(1));}
        for(auto ij:powers)out.push_back(xp[ij[0]]*yp[ij[1]]);return out;
    }
};
double Factorial(int i){return i<2?1:(i==2?2:6);}
struct Row {int field,nx,ny,mode;std::string word;double error=0,truth=0,linf=0,face=0,facetruth=0,jump=0,jumpabs=0,operatorerror=0,directnorm=0;};
struct ElementData { Geometry geometry;std::vector<Vector> original,polynomial,recursive; };
}
int main(int argc,char **argv) {
    Mpi::Init(argc,argv);
    if(argc!=8){std::cerr<<"family n p amplitude basis map extra\n";return 2;}
    std::string family=argv[1],basis=argv[5],mapping=argv[6];int n=std::stoi(argv[2]),p=std::stoi(argv[3]),extra=std::stoi(argv[7]);double amplitude=std::stod(argv[4]),h=1./n;
    bool triangle=family=="tri";
    Mesh mesh=Mesh::MakeCartesian2D(n,n,triangle?Element::TRIANGLE:Element::QUADRILATERAL,true);
    std::vector<direct::ElementData> data(mesh.GetNE());
    for(int e=0;e<mesh.GetNE();++e){auto *tr=mesh.GetElementTransformation(e);IntegrationPoint ip;ip.Set2(0,0);Vector x;tr->Transform(ip,x);tr->SetIntPoint(&ip);auto &j=tr->Jacobian();data[e].geometry={x(0),x(1),j(0,0),j(0,1),j(1,0),j(1,1),amplitude,mapping=="coupled"};}
    mesh.SetCurvature(4);mesh.Transform([&](const Vector&s,Vector&x){x=s;if(mapping=="coupled"){double b=amplitude*s(0)*(1-s(0))*s(1)*(1-s(1));x(0)+=b;x(1)+=.7*b;}else{x(0)+=amplitude*s(0)*(1-s(0));x(1)+=.7*amplitude*s(1)*(1-s(1));}});
    int bt=basis=="gll"?BasisType::GaussLobatto:(basis=="gl"?BasisType::GaussLegendre:BasisType::Positive);
    DG_FECollection fec(p,2,bt);FiniteElementSpace space(&mesh,&fec);direct::PolynomialBasis poly(*space.GetFE(0),triangle);
    std::vector<int> fields=mapping=="coupled"?std::vector<int>{0,1,2,3}:std::vector<int>{4,3};if(amplitude==0)fields.push_back(5);
    std::vector<std::string> words={"x","y","xx","xy","yy"};if(p>=3)for(auto s:{"xxx","xxy","xyy","yyy"})words.push_back(s);
    std::vector<direct::Row> rows;
    for(int f:fields)for(auto word:words)for(int mode=0;mode<2;++mode){int nx=std::count(word.begin(),word.end(),'x');rows.push_back({f,nx,int(word.size())-nx,mode,word});}
    double input[6]={},inverse_residual=0,representation_residual=0,first_derivative_residual=0,geometry_residual=0,minjac=1e100,face_measure=0,interior_measure=0;int qb=0;
    auto evaluate=[&](int e,const IntegrationPoint &ip,Vector &shape,std::array<direct::Jet,6>&values) {
        const auto &ed=data[e];auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);fe->CalcShape(ip,shape);map->SetIntPoint(&ip);
        auto rs=ed.geometry.InverseJet(ip.x,ip.y,h,inverse_residual);auto monomials=poly.Evaluate(rs[0],rs[1],p);
        DenseMatrix ds(fe->GetDof(),2);fe->CalcPhysDShape(*map,ds);Vector x;map->Transform(ip,x);
        auto xy=ed.geometry.Map(direct::Jet(ip.x),direct::Jet(ip.y));
        geometry_residual=std::max(geometry_residual,std::max(std::abs(xy[0].c[0]-x(0)),std::abs(xy[1].c[0]-x(1)))/h);
        for(int f:fields){values[f]=direct::Jet();for(int j=0;j<int(monomials.size());++j)values[f]=values[f]+ed.polynomial[f](j)*monomials[j];
            representation_residual=std::max(representation_residual,std::abs(values[f].c[0]-shape*ed.original[f]));
            for(int axis=0;axis<2;++axis){double v=0;for(int j=0;j<fe->GetDof();++j)v+=ds(j,axis)*ed.original[f](j);first_derivative_residual=std::max(first_derivative_residual,std::abs(values[f].c[axis+1]/h-v)/std::max(1.,std::abs(v)));}
        }
    };
    auto derivative=[&](const direct::Jet &value,const direct::Row &r){return value.c[direct::Index(r.nx,r.ny)]*direct::Factorial(r.nx)*direct::Factorial(r.ny)/std::pow(h,r.nx+r.ny);};
    for(int e=0;e<space.GetNE();++e){auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);int nd=fe->GetDof();auto &ed=data[e];qb=ofdg::detail::CurvedQuadratureOrder(p,*map);
        // Exactly the existing input projection and production derivative matrices.
        auto D=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,qb);const auto &init=ofdg::detail::CurvedRule(fe->GetGeomType(),qb+16);DenseMatrix mass(nd);mass=0.;std::vector<Vector> rhs(6,Vector(nd));for(auto &v:rhs)v=0.;Vector shape(nd),x;
        for(int q=0;q<init.GetNPoints();++q){auto &ip=init.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);fe->CalcShape(ip,shape);double w=ip.weight*map->Weight();AddMult_a_VVt(w,shape,mass);for(int f:fields)rhs[f].Add(w*Exact(f,x,0,0,amplitude),shape);}
        DenseMatrixInverse solve(mass);ed.original.resize(6,Vector(nd));ed.polynomial.resize(6,Vector(nd));ed.recursive.resize(rows.size());
        for(int f:fields){solve.Mult(rhs[f],ed.original[f]);poly.coefficients.Mult(ed.original[f],ed.polynomial[f]);}
        for(size_t i=0;i<rows.size();++i)if(rows[i].mode==0)ed.recursive[i]=Apply(D,ed.original[rows[i].field],rows[i].word);
        const auto &ir=ofdg::detail::CurvedRule(fe->GetGeomType(),qb+16+extra);
        for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);std::array<direct::Jet,6> values;evaluate(e,ip,shape,values);map->Transform(ip,x);double w=ip.weight*map->Weight();minjac=std::min(minjac,map->Jacobian().Det());
            for(int f:fields){double d=shape*ed.original[f]-Exact(f,x,0,0,amplitude);input[f]+=w*d*d;}
            for(size_t i=0;i<rows.size();++i){auto &r=rows[i];double exact=Exact(r.field,x,r.nx,r.ny,amplitude),dv=derivative(values[r.field],r),v=r.mode?dv:shape*ed.recursive[i],err=v-exact;r.error+=w*err*err;r.truth+=w*exact*exact;r.linf=std::max(r.linf,std::abs(err));r.operatorerror+=w*(v-dv)*(v-dv);r.directnorm+=w*dv*dv;}
        }
    }
    for(int face=0;face<mesh.GetNumFaces();++face){auto *tr=mesh.GetFaceElementTransformations(face);if(!tr)continue;const auto &ir=IntRules.Get(tr->GetGeometryType(),qb+16+extra);
        for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);tr->SetAllIntPoints(&ip);Vector x;tr->Face->Transform(ip,x);double w=ip.weight*tr->Face->Weight();std::vector<double> left(rows.size());
            for(int side=0;side<2;++side){int e=side?tr->Elem2No:tr->Elem1No;if(e<0)continue;face_measure+=w;const auto &local=side?tr->GetElement2IntPoint():tr->GetElement1IntPoint();Vector shape(space.GetFE(e)->GetDof());std::array<direct::Jet,6> values;evaluate(e,local,shape,values);
                for(size_t i=0;i<rows.size();++i){auto &r=rows[i];double v=r.mode?derivative(values[r.field],r):shape*data[e].recursive[i],t=Exact(r.field,x,r.nx,r.ny,amplitude);r.face+=w*(v-t)*(v-t);r.facetruth+=w*t*t;if(side){r.jump+=w*(v-left[i])*(v-left[i]);r.jumpabs+=w*std::abs(v-left[i]);}else left[i]=v;}
                if(side)interior_measure+=w;
            }
        }
    }
    if(minjac<=0||inverse_residual>1e-9||representation_residual>1e-8||first_derivative_residual>1e-7||geometry_residual>1e-10){std::cerr<<"Control failed: "<<minjac<<' '<<inverse_residual<<' '<<representation_residual<<' '<<first_derivative_residual<<' '<<geometry_residual<<'\n';return 3;}
    std::cout<<"family,n,p,amplitude,basis,map,extra,field,derivative,method,error_l2,relative_l2,sampled_linf,truth_l2,face_rms,face_relative_l2,jump_rms,jump_mean_abs,input_error_l2,operator_error_l2,operator_relative_l2,inverse_residual,representation_residual,first_derivative_residual,geometry_residual,min_jacobian,elements\n"<<std::setprecision(17);
    const char *names[]={"sine","shortwave","exponential","constant","represented_quadratic","physical_quadratic"};
    for(const auto &r:rows)std::cout<<family<<','<<n<<','<<p<<','<<amplitude<<','<<basis<<','<<mapping<<','<<extra<<','<<names[r.field]<<','<<r.word<<','<<(r.mode?"direct":"recursive")<<','<<std::sqrt(r.error)<<','<<(r.truth>1e-24?std::sqrt(r.error/r.truth):NAN)<<','<<r.linf<<','<<std::sqrt(r.truth)<<','<<std::sqrt(r.face/face_measure)<<','<<(r.facetruth>1e-24?std::sqrt(r.face/r.facetruth):NAN)<<','<<std::sqrt(r.jump/interior_measure)<<','<<r.jumpabs/interior_measure<<','<<std::sqrt(input[r.field])<<','<<std::sqrt(r.operatorerror)<<','<<(r.directnorm>1e-24?std::sqrt(r.operatorerror/r.directnorm):NAN)<<','<<inverse_residual<<','<<representation_residual<<','<<first_derivative_residual<<','<<geometry_residual<<','<<minjac<<','<<space.GetNE()<<'\n';
}
