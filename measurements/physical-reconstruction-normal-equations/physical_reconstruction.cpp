// Isolated sensor experiment; does not change production OFDG.
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

// Derivatives of Legendre polynomials with respect to their argument.
std::vector<std::array<double,4>> Legendre(double x, int degree)
{
    std::vector<std::array<double,4>> v(degree+1);
    v[0] = {1,0,0,0};
    if (degree) { v[1] = {x,1,0,0}; }
    for (int n=2; n<=degree; ++n) {
        for (int d=0; d<=3; ++d) {
            v[n][d] = ((2*n-1)*(x*v[n-1][d] + (d ? d*v[n-1][d-1] : 0))
                       -(n-1)*v[n-2][d])/n;
        }
    }
    return v;
}

struct PhysicalBasis {
    int degree;
    std::array<double,2> center, scale;
    int Size() const { return (degree+1)*(degree+2)/2; }
    Vector Evaluate(const Vector &x, int nx=0, int ny=0) const
    {
        auto px=Legendre((x(0)-center[0])/scale[0],degree);
        auto py=Legendre((x(1)-center[1])/scale[1],degree);
        Vector values(Size()); int i=0;
        const double factor=std::pow(scale[0],-nx)*std::pow(scale[1],-ny);
        for (int total=0; total<=degree; ++total) {
            for (int b=0; b<=total; ++b) {
                values(i++)=factor*px[total-b][nx]*py[b][ny];
            }
        }
        return values;
    }
};
struct Sample { Vector x, shape; double weight; };
std::vector<Sample> Samples(const FiniteElement &fe, ElementTransformation &map, int order)
{
    std::vector<Sample> samples;
    const auto &rule=ofdg::detail::CurvedRule(fe.GetGeomType(),order);
    for (int q=0; q<rule.GetNPoints(); ++q) {
        const auto &ip=rule.IntPoint(q);
        map.SetIntPoint(&ip);
        Sample s; s.shape.SetSize(fe.GetDof());
        map.Transform(ip,s.x); fe.CalcShape(ip,s.shape);
        s.weight=ip.weight*map.Weight(); samples.push_back(std::move(s));
    }
    return samples;
}
struct Row {
    int field, mode, nx, ny; std::string word;
    double error=0, truth=0, maximum=0, faceerror=0, facetruth=0;
    double jump2=0, jump1=0;
};
struct ElementData {
    std::vector<Vector> derivatives;
    std::vector<std::vector<Vector>> reconstructed;
    std::vector<PhysicalBasis> basis;
};
int main(int argc, char **argv)
{
    Mpi::Init(argc,argv);
    if (argc!=8) { std::cerr<<"family n p amplitude basis map extra\n"; return 2; }
    const std::string family=argv[1], basis_name=argv[5], mapping=argv[6];
    const int n=std::stoi(argv[2]), p=std::stoi(argv[3]), extra=std::stoi(argv[7]);
    const double amplitude=std::stod(argv[4]);
    Mesh mesh=Mesh::MakeCartesian2D(n,n,family=="tri"?Element::TRIANGLE:Element::QUADRILATERAL,true);
    mesh.SetCurvature(4);
    mesh.Transform([&](const Vector &s,Vector &x) {
        x=s;
        if (mapping=="coupled") {
            const double b=amplitude*s(0)*(1-s(0))*s(1)*(1-s(1));
            x(0)+=b; x(1)+=.7*b;
        } else {
            x(0)+=amplitude*s(0)*(1-s(0));
            x(1)+=.7*amplitude*s(1)*(1-s(1));
        }
    });
    int bt=basis_name=="gll"?BasisType::GaussLobatto:
           (basis_name=="gl"?BasisType::GaussLegendre:BasisType::Positive);
    DG_FECollection collection(p,2,bt); FiniteElementSpace space(&mesh,&collection);
    std::vector<int> fields=mapping=="coupled"?std::vector<int>{0,1,2,3}:std::vector<int>{4,3};
    if (amplitude==0) { fields.push_back(5); }
    std::vector<std::string> words={"x","y","xx","xy","yy"};
    if (p>=3) { for (auto w:{"xxx","xxy","xyy","yyy"}) { words.push_back(w); } }
    std::vector<Row> rows;
    for (int f:fields) {
        for (const auto &w:words) {
            int nx=std::count(w.begin(),w.end(),'x');
            for (int mode=0; mode<5; ++mode) { rows.push_back({f,mode,nx,int(w.size())-nx,w}); }
        }
    }
    std::vector<ElementData> data(space.GetNE());
    double input_error[6]={}, reconstruction_error[6][5]={};
    double face_measure=0, interior_measure=0, mean_defect=0, solve_residual=0;
    // A degree-q physical polynomial composed with degree-4 geometry has
    // degree <=4q. The mass integrand has degree <=8q plus the determinant.
    const int assembly_order=8*(p+2)+8+extra, evaluation_order=assembly_order+8;
    for (int e=0; e<space.GetNE(); ++e) {
        const auto &fe=*space.GetFE(e); auto &map=*space.GetElementTransformation(e);
        const int nd=fe.GetDof(); auto &ed=data[e];
        const int production_order=ofdg::detail::CurvedQuadratureOrder(p,map);
        auto D=ofdg::detail::ProjectedPhysicalDerivatives(fe,map,production_order);
        auto input_samples=Samples(fe,map,production_order+16); // unchanged initialization
        DenseMatrix mass(nd); mass=0.; std::vector<Vector> rhs(6,Vector(nd)), U(6,Vector(nd));
        for (auto &v:rhs) { v=0.; }
        for (const auto &s:input_samples) {
            AddMult_a_VVt(s.weight,s.shape,mass);
            for (int f:fields) { rhs[f].Add(s.weight*Exact(f,s.x,0,0,amplitude),s.shape); }
        }
        DenseMatrixInverse input_solve(mass);
        for (int f:fields) { input_solve.Mult(rhs[f],U[f]); }
        ed.derivatives.resize(rows.size());
        for (size_t i=0; i<rows.size(); ++i) {
            if (rows[i].mode==0) { ed.derivatives[i]=Apply(D,U[rows[i].field],rows[i].word); }
        }
        auto assembly=Samples(fe,map,assembly_order);
        auto evaluation=Samples(fe,map,evaluation_order);
        std::array<double,2> low={1e100,1e100}, high={-1e100,-1e100};
        for (const auto &s:assembly) {
            for (int j=0;j<2;++j) { low[j]=std::min(low[j],s.x(j)); high[j]=std::max(high[j],s.x(j)); }
        }
        ed.reconstructed.resize(3,std::vector<Vector>(6));
        std::vector<Vector> oracle(6);
        for (int enrichment=0; enrichment<3; ++enrichment) {
            PhysicalBasis pb{p+enrichment,{}, {}};
            for (int j=0;j<2;++j) { pb.center[j]=(low[j]+high[j])/2; pb.scale[j]=(high[j]-low[j])/2; }
            ed.basis.push_back(pb); int nb=pb.Size();
            DenseMatrix M(nb); M=0.; std::vector<Vector> b(6,Vector(nb)), exact_rhs(6,Vector(nb));
            for (auto &v:b) { v=0.; } for (auto &v:exact_rhs) { v=0.; }
            for (const auto &s:assembly) {
                Vector psi=pb.Evaluate(s.x);
                AddMult_a_VVt(s.weight,psi,M);
                for (int f:fields) {
                    b[f].Add(s.weight*(s.shape*U[f]),psi);
                    if (enrichment==2) { exact_rhs[f].Add(s.weight*Exact(f,s.x,0,0,amplitude),psi); }
                }
            }
            DenseMatrixInverse inverse(M);
            for (int f:fields) {
                auto &c=ed.reconstructed[enrichment][f]; c.SetSize(nb); inverse.Mult(b[f],c);
                Vector residual(nb); M.Mult(c,residual); residual-=b[f];
                solve_residual=std::max(solve_residual,residual.Normlinf()/std::max(1.,b[f].Normlinf()));
                if (enrichment==2) { oracle[f].SetSize(nb); inverse.Mult(exact_rhs[f],oracle[f]); }
            }
        }
        // The oracle is diagnostic only: it receives the analytic field rather
        // than U, so it must never enter the feasible-method win statistics.
        ed.reconstructed.push_back(oracle);
        double element_mean[6][5]={};
        for (const auto &s:evaluation) {
            for (int f:fields) {
                double u=s.shape*U[f], truth=Exact(f,s.x,0,0,amplitude);
                input_error[f]+=s.weight*(u-truth)*(u-truth);
                for (int mode=1;mode<5;++mode) {
                    int b=std::min(mode-1,2);
                    double difference=ed.basis[b].Evaluate(s.x)*ed.reconstructed[mode-1][f]-u;
                    reconstruction_error[f][mode]+=s.weight*difference*difference;
                    element_mean[f][mode]+=s.weight*difference;
                }
            }
            // Cache one physical evaluation vector per degree and derivative.
            for (const auto &word:words) {
                int nx=std::count(word.begin(),word.end(),'x'), ny=word.size()-nx;
                std::vector<Vector> psi;
                for (const auto &b:ed.basis) { psi.push_back(b.Evaluate(s.x,nx,ny)); }
                for (size_t i=0;i<rows.size();++i) {
                    auto &r=rows[i]; if (r.word!=word) { continue; }
                    double v=r.mode==0?s.shape*ed.derivatives[i]:psi[std::min(r.mode-1,2)]*ed.reconstructed[r.mode-1][r.field];
                    double truth=Exact(r.field,s.x,nx,ny,amplitude), error=v-truth;
                    r.error+=s.weight*error*error; r.truth+=s.weight*truth*truth;
                    r.maximum=std::max(r.maximum,std::abs(error));
                }
            }
        }
        for (int f:fields) { for (int m=1;m<4;++m) { mean_defect=std::max(mean_defect,std::abs(element_mean[f][m])); } }
    }
    for (int f=0;f<mesh.GetNumFaces();++f) {
        auto *tr=mesh.GetFaceElementTransformations(f); if (!tr) { continue; }
        const auto &rule=IntRules.Get(tr->GetGeometryType(),evaluation_order);
        for (int q=0;q<rule.GetNPoints();++q) {
            const auto &ip=rule.IntPoint(q); tr->SetAllIntPoints(&ip);
            Vector x; tr->Face->Transform(ip,x); double weight=ip.weight*tr->Face->Weight();
            std::vector<double> first(rows.size());
            for (int side=0;side<2;++side) {
                int e=side?tr->Elem2No:tr->Elem1No; if (e<0) { continue; }
                auto &ed=data[e]; Vector shape(space.GetFE(e)->GetDof());
                space.GetFE(e)->CalcShape(side?tr->GetElement2IntPoint():tr->GetElement1IntPoint(),shape);
                face_measure+=weight;
                for (const auto &word:words) {
                    int nx=std::count(word.begin(),word.end(),'x'),ny=word.size()-nx;
                    std::vector<Vector> psi;
                    for (const auto &b:ed.basis) { psi.push_back(b.Evaluate(x,nx,ny)); }
                    for (size_t i=0;i<rows.size();++i) {
                        auto &r=rows[i]; if (r.word!=word) { continue; }
                        double v=r.mode==0?shape*ed.derivatives[i]:psi[std::min(r.mode-1,2)]*ed.reconstructed[r.mode-1][r.field];
                        double truth=Exact(r.field,x,nx,ny,amplitude);
                        r.faceerror+=weight*(v-truth)*(v-truth); r.facetruth+=weight*truth*truth;
                        if (!side) { first[i]=v; }
                        else { double jump=v-first[i]; r.jump2+=weight*jump*jump; r.jump1+=weight*std::abs(jump); }
                    }
                }
            }
            if (tr->Elem2No>=0) { interior_measure+=weight; }
        }
    }
    if (!std::isfinite(mean_defect)||mean_defect>1e-9||solve_residual>1e-10) {
        std::cerr<<"Projection control failure "<<mean_defect<<" "<<solve_residual<<'\n'; return 3;
    }
    const char *names[]={"sine","shortwave","exponential","constant","represented_quadratic","physical_quadratic"};
    const char *methods[]={"recursive","physical_k","physical_k1","physical_k2","oracle_k2"};
    std::cout<<"family,n,p,amplitude,basis,map,extra,field,derivative,method,q,error_l2,relative_l2,sampled_linf,truth_l2,face_rms,face_relative_l2,jump_rms,jump_mean_abs,input_error_l2,reconstruction_l2,mean_defect,solve_residual,assembly_order,elements\n"<<std::setprecision(17);
    for (const auto &r:rows) {
        std::cout<<family<<','<<n<<','<<p<<','<<amplitude<<','<<basis_name<<','<<mapping<<','<<extra<<','<<names[r.field]<<','<<r.word<<','<<methods[r.mode]<<','<<(r.mode==0?p:p+std::min(r.mode-1,2))<<','<<std::sqrt(r.error)<<','<<(r.truth>1e-24?std::sqrt(r.error/r.truth):NAN)<<','<<r.maximum<<','<<std::sqrt(r.truth)<<','<<std::sqrt(r.faceerror/face_measure)<<','<<(r.facetruth>1e-24?std::sqrt(r.faceerror/r.facetruth):NAN)<<','<<std::sqrt(r.jump2/interior_measure)<<','<<r.jump1/interior_measure<<','<<std::sqrt(input_error[r.field])<<','<<std::sqrt(reconstruction_error[r.field][r.mode])<<','<<mean_defect<<','<<solve_residual<<','<<assembly_order<<','<<space.GetNE()<<'\n';
    }
}
