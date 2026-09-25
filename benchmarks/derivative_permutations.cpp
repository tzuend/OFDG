// Isolated manufactured-derivative experiment. No stabilization or PDE evolution.
#include "mfem.hpp"
#include "../src/curved_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace mfem;

// Words list chronological application order, e.g. xxy = Dy Dx Dx U.
std::vector<std::string> Permutations(std::string word) {
    std::vector<std::string> result;
    do { result.push_back(word); } while(std::next_permutation(word.begin(),word.end()));
    return result;
}
Vector Apply(const std::vector<DenseMatrix> &D, const Vector &u, const std::string &word) {
    Vector v(u), next(u.Size());
    for(char c:word) { D[c=='x'?0:1].Mult(v,next); v=next; }
    return v;
}
double Exact(int field,const Vector &x,int nx,int ny) {
    if(field==0) return std::pow(2.,nx)*std::pow(3.,ny)*std::sin(2*x(0)+3*x(1)+(nx+ny)*std::acos(-1.)/2);
    if(field==1) return std::pow(.5,ny)*std::exp(x(0)+.5*x(1));
    // x^2 y + x y^2: nonzero xxy and xyy derivatives, both exactly 2.
    auto monomial=[&](int a,int b) { if(nx>a||ny>b)return 0.; double r=1.;for(int j=0;j<nx;++j)r*=a-j;for(int j=0;j<ny;++j)r*=b-j;return r*std::pow(x(0),a-nx)*std::pow(x(1),b-ny); };
    return monomial(2,1)+monomial(1,2);
}
struct Row {int field;std::string multi,path; double error=0,shift=0,exact=0,floor=0,input=0;};
int main(int argc,char **argv) {
    Mpi::Init(argc,argv);
    if(argc!=6) {std::cerr<<"family n p amplitude qextra\n";return 2;}
    std::string family=argv[1];int n=std::stoi(argv[2]),p=std::stoi(argv[3]),extra=std::stoi(argv[5]);double amplitude=std::stod(argv[4]);
    Mesh mesh=Mesh::MakeCartesian2D(n,n,family=="tri"?Element::TRIANGLE:Element::QUADRILATERAL,true);
    mesh.SetCurvature(4);
    mesh.Transform([amplitude](const Vector &s,Vector &x){x=s;double b=amplitude*s(0)*(1-s(0))*s(1)*(1-s(1));x(0)+=b;x(1)+=.7*b;});
    DG_FECollection fec(p,2,BasisType::GaussLobatto);FiniteElementSpace space(&mesh,&fec);
    std::vector<Row> rows;
    for(int f=0;f<3;++f)for(std::string multi:{"xy","xxy","xyy"}) {auto paths=Permutations(multi);paths.push_back("mean");for(auto path:paths)rows.push_back({f,multi,path});}
    double weak=0,constant=0,min_weight=1e100;int qbase=0;
    for(int e=0;e<space.GetNE();++e) {
        auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);int nd=fe->GetDof();
        qbase=ofdg::detail::CurvedQuadratureOrder(p,*map);
        auto D=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,qbase+extra);
        // Independent physical mass and RHS assembly at a higher quadrature order.
        const auto &rule=ofdg::detail::CurvedRule(fe->GetGeomType(),qbase+16);
        DenseMatrix M(nd);M=0.;std::vector<Vector> rhs(3,Vector(nd));for(auto &v:rhs)v=0.;
        Vector sh(nd),x;DenseMatrix grad(nd,2);
        for(int q=0;q<rule.GetNPoints();++q) {auto &ip=rule.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);fe->CalcShape(ip,sh);double w=ip.weight*map->Weight();min_weight=std::min(min_weight,map->Jacobian().Det());for(int i=0;i<nd;++i)for(int j=0;j<nd;++j)M(i,j)+=w*sh(i)*sh(j);for(int f=0;f<3;++f)rhs[f].Add(w*Exact(f,x,0,0),sh);}
        DenseMatrixInverse solve(M);std::vector<Vector> U(3,Vector(nd));for(int f=0;f<3;++f)solve.Mult(rhs[f],U[f]);
        Vector one(nd),z(nd);one=1.;for(auto &d:D){d.Mult(one,z);constant=std::max(constant,z.Normlinf());}
        for(int f=0;f<3;++f) {
            // Independent weak first-derivative check on the projected input U.
            for(int j=0;j<2;++j) {
                Vector b(nd),mz(nd);b=0.;D[j].Mult(U[f],z);M.Mult(z,mz);
                for(int q=0;q<rule.GetNPoints();++q){auto &ip=rule.IntPoint(q);map->SetIntPoint(&ip);fe->CalcShape(ip,sh);fe->CalcPhysDShape(*map,grad);double direct=0.;for(int i=0;i<nd;++i)direct+=grad(i,j)*U[f](i);b.Add(ip.weight*map->Weight()*direct,sh);}
                mz-=b;weak=std::max(weak,mz.Normlinf()/std::max(1.,b.Normlinf()));
            }
            for(std::string multi:{"xy","xxy","xyy"}) {
                auto paths=Permutations(multi);std::vector<Vector> values;Vector mean(nd);mean=0.;for(auto path:paths){values.push_back(Apply(D,U[f],path));mean.Add(1./paths.size(),values.back());}values.push_back(mean);paths.push_back("mean");
                int nx=std::count(multi.begin(),multi.end(),'x'),ny=multi.size()-nx;Vector ref_rhs(nd),ref(nd);ref_rhs=0.;
                for(int q=0;q<rule.GetNPoints();++q){auto &ip=rule.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);fe->CalcShape(ip,sh);ref_rhs.Add(ip.weight*map->Weight()*Exact(f,x,nx,ny),sh);}solve.Mult(ref_rhs,ref);
                for(int q=0;q<rule.GetNPoints();++q) {
                    auto &ip=rule.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);fe->CalcShape(ip,sh);double w=ip.weight*map->Weight(),truth=Exact(f,x,nx,ny),avg=sh*mean,fl=sh*ref-truth,initial=sh*U[f]-Exact(f,x,0,0);
                    for(size_t r=0;r<paths.size();++r)for(auto &row:rows)if(row.field==f&&row.multi==multi&&row.path==paths[r]){double v=sh*values[r];row.error+=w*(v-truth)*(v-truth);row.shift+=w*(v-avg)*(v-avg);row.exact+=w*truth*truth;row.floor+=w*fl*fl;row.input+=w*initial*initial;}
                }
            }
        }
    }
    if(!std::isfinite(weak)||weak>1e-9||constant>1e-8||min_weight<=0) {std::cerr<<"failed weak/constant/geometry control "<<weak<<" "<<constant<<" "<<min_weight<<"\n";return 3;}
    std::cout<<"family,n,p,amplitude,qextra,field,multi,path,error_l2,path_minus_mean_l2,exact_l2,best_projection_error_l2,input_error_l2,first_weak_relative,constant_max,min_jacobian_det,quadrature_order,elements\n"<<std::setprecision(17);
    const char *names[]={"sine","exponential","cubic"};
    for(auto &r:rows) std::cout<<family<<','<<n<<','<<p<<','<<amplitude<<','<<extra<<','<<names[r.field]<<','<<r.multi<<','<<r.path<<','<<std::sqrt(r.error)<<','<<std::sqrt(r.shift)<<','<<std::sqrt(r.exact)<<','<<std::sqrt(r.floor)<<','<<std::sqrt(r.input)<<','<<weak<<','<<constant<<','<<min_weight<<','<<qbase+extra<<','<<space.GetNE()<<'\n';
}
