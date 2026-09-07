#include "mfem.hpp"
#include "../src/ofdg.hpp"
#include "../src/kxrcf.hpp"
#include "../src/curved_geometry.hpp"
#include "../src/euler_positivity.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace mfem;

void Require(bool condition, const std::string &message)
{ if (!condition) { throw std::runtime_error(message); } }

void Bend(const Vector &x, Vector &y)
{
    y = x;
    y(0) += 0.18 * x(0) * (1.0-x(0)) * x(1);
    y(1) += 0.12 * x(1) * (1.0-x(1)) * x(0);
    if (x.Size() == 3) { y(2) += 0.1*x(2)*(1.0-x(2))*x(0); }
}

Mesh PolynomialMesh(int dim, Element::Type type, int n = 2)
{
    Mesh mesh = dim == 2 ? Mesh::MakeCartesian2D(n,n,type,true) :
                          Mesh::MakeCartesian3D(n,n,n,type);
    mesh.SetCurvature(3);
    mesh.Transform(Bend);
    return mesh;
}

Vector Integrals(FiniteElementSpace &space, const Vector &values)
{
    GridFunction state(&space); state = values;
    Vector integrals(space.GetNE()*space.GetVDim()); integrals = 0.0;
    for (int e=0; e<space.GetNE(); ++e) {
        auto *map = space.GetElementTransformation(e);
        const auto &rule = ofdg::detail::CurvedRule(space.GetFE(e)->GetGeomType(), 30);
        for (int q=0; q<rule.GetNPoints(); ++q) {
            const auto &ip=rule.IntPoint(q); map->SetIntPoint(&ip);
            for (int c=0; c<space.GetVDim(); ++c) {
                integrals(e*space.GetVDim()+c) += ip.weight*map->Weight()*state.GetValue(e,ip,c+1);
            }
        }
    }
    return integrals;
}

void SetState(FiniteElementSpace &space, GridFunction &state)
{
    VectorFunctionCoefficient function(space.GetVDim(), [](const Vector &x, Vector &u) {
        for (int c=0; c<u.Size(); ++c) { u(c)=2.0+c+0.2*x(0)+0.1*x(1)*x(1); }
    });
    state.ProjectCoefficient(function);
    Array<int> dofs;
    for (int e=0; e<space.GetNE(); ++e) {
        space.GetElementDofs(e,dofs);
        for (int c=0; c<space.GetVDim(); ++c) {
            for (int i=0; i<dofs.Size(); ++i) {
                state(space.DofToVDof(dofs[i],c)) += 0.03*(e%3)*(i+1)/dofs.Size();
            }
        }
    }
}

void ValidateFilter(Mesh &mesh, int degree, const std::string &name)
{
    DG_FECollection collection(degree,mesh.Dimension(),BasisType::GaussLobatto);
    FiniteElementSpace space(&mesh,&collection,2,Ordering::byNODES);
    GridFunction state(&space); SetState(space,state);
    ofdg::OFDG filter(&space,BasisType::GaussLobatto);
    Vector result, residual;
    filter.CompDecay(state,result,0.05);
    filter.ComputeStabilization(state,residual);
    Require(std::isfinite(result.Norml2()),name+" nonfinite decay");
    Vector difference=Integrals(space,result); difference-=Integrals(space,state);
    Require(difference.Normlinf()<2e-10,name+" mean drift");
    Require(Integrals(space,residual).Normlinf()<2e-10,name+" residual mean drift");
    Array<bool> active(space.GetNE()); active=false;
    filter.CompDecay(state,result,0.05,&active); result-=state;
    Require(result.Normlinf()==0.0,name+" inactive cells changed");
    state=2.0; filter.CompDecay(state,result,0.05); result-=state;
    if (!(result.Normlinf()<1e-12)) { std::cerr << "constant defect=" << result.Normlinf() << " degree=" << degree << "\n"; }
    Require(result.Normlinf()<1e-12,name+" constant changed");

    // Same polynomial field represented in two coefficient bases.
    DG_FECollection other_collection(degree,mesh.Dimension(),BasisType::GaussLegendre);
    FiniteElementSpace other_space(&mesh,&other_collection,2,Ordering::byNODES);
    SetState(space,state);
    GridFunction other_state(&other_space); other_state.ProjectGridFunction(state);
    ofdg::OFDG other_filter(&other_space,BasisType::GaussLegendre);
    Vector other_result; other_filter.CompDecay(other_state,other_result,0.05);
    filter.CompDecay(state,result,0.05);
    GridFunction filtered(&space); filtered=result;
    GridFunction expected(&other_space); expected.ProjectGridFunction(filtered);
    other_result-=expected;
    Require(other_result.Normlinf()<3e-10,name+" coefficient basis dependence");
    std::cout << name << " P" << degree << " conservation/basis/constants passed\n";
}

void ValidateMPI(Mesh &serial_mesh, const std::string &name)
{
    DG_FECollection collection(2,serial_mesh.Dimension(),BasisType::GaussLobatto);
    FiniteElementSpace serial_space(&serial_mesh,&collection,2,Ordering::byNODES);
    GridFunction serial_state(&serial_space); SetState(serial_space,serial_state);
    int ranks; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    Array<int> partition(serial_mesh.GetNE());
    for (int e=0;e<partition.Size();++e) { partition[e]=e%ranks; }
    ParMesh mesh(MPI_COMM_WORLD,serial_mesh,partition.GetData());
    ParFiniteElementSpace space(&mesh,&collection,2,Ordering::byNODES);
    ParGridFunction state(&mesh,&serial_state,partition.GetData());
    VectorFunctionCoefficient velocity(mesh.Dimension(), [](const Vector &,Vector &v) { v=1.0; });
    auto physics=std::make_shared<ofdg::AdvectionFacePhysics>(&velocity);
    ofdg::OFDG serial_filter(&serial_space,BasisType::GaussLobatto,physics);
    ofdg::OFDG filter(&space,BasisType::GaussLobatto,physics);
    ofdg::KXRCFIndicator serial_indicator(&serial_space,physics,0.001);
    ofdg::KXRCFIndicator indicator(&space,physics,0.001);
    Array<bool> serial_active,active; Vector serial_values,values;
    serial_indicator.Compute(serial_state,serial_active,&serial_values);
    indicator.Compute(state,active,&values);
    // MFEM global element numbers are rank-contiguous, not the original
    // serial element numbering. Transfer labels through the same partition.
    GridFunction serial_labels(&serial_space); Array<int> label_dofs;
    for(int e=0;e<serial_space.GetNE();++e) {
        serial_space.GetElementVDofs(e,label_dofs);
        for(int i=0;i<label_dofs.Size();++i){serial_labels(label_dofs[i])=e;}
    }
    ParGridFunction labels(&mesh,&serial_labels,partition.GetData());
    for(int e=0;e<space.GetNE();++e) {
        int global=static_cast<int>(std::lround(labels.GetValue(e,Geometries.GetCenter(space.GetFE(e)->GetGeomType()))));
        Require(active[e]==serial_active[global],name+" MPI mask mismatch");
        if (std::abs(values(e)-serial_values(global))>=3e-10) {
            std::cerr << name << " e=" << global << " local=" << values(e)
                      << " serial=" << serial_values(global) << " difference="
                      << values(e)-serial_values(global) << '\n';
        }
        Require(std::abs(values(e)-serial_values(global))<3e-10,name+" MPI indicator mismatch");
    }
    for (int masked=0;masked<2;++masked) {
        Vector result,serial_result;
        filter.CompDecay(state,result,0.05,masked?&active:nullptr);
        serial_filter.CompDecay(serial_state,serial_result,0.05,masked?&serial_active:nullptr);
        GridFunction serial_filtered(&serial_space);serial_filtered=serial_result;
        ParGridFunction expected(&mesh,&serial_filtered,partition.GetData());
        result-=expected;
        double local=result.Normlinf(),global;
        MPI_Allreduce(&local,&global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        Require(global<3e-10,name+" MPI decay mismatch");
    }
    std::cout<<name<<" MPI OFDG/KXRCF passed\n";
}

void ValidateDerivatives()
{
    // Fixed smooth geometry and u=sin(x): exact derivatives are analytic.
    double previous[3]={0,0,0};
    for(int n : {2,4,8}) {
        Mesh mesh=PolynomialMesh(2,Element::QUADRILATERAL,n);
        DG_FECollection collection(4,2,BasisType::GaussLobatto);
        FiniteElementSpace space(&mesh,&collection);
        GridFunction state(&space); FunctionCoefficient function([](const Vector &x){return std::sin(x(0));});
        state.ProjectCoefficient(function);
        double errors[3]={0,0,0};
        for(int e=0;e<space.GetNE();++e) {
            auto *fe=space.GetFE(e); auto *map=space.GetElementTransformation(e);
            int qorder=ofdg::detail::CurvedQuadratureOrder(4,*map);
            auto derivative=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,qorder);
            auto finer=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,qorder+8);
            finer[0]-=derivative[0];
            Require(finer[0].MaxMaxNorm()<1e-9,"derivative quadrature sensitivity");
            Array<int>dofs;space.GetElementDofs(e,dofs);Vector value;state.GetSubVector(dofs,value);
            Vector repeated[3];
            for(int d=0;d<3;++d) { repeated[d].SetSize(fe->GetDof()); derivative[0].Mult(d?repeated[d-1]:value,repeated[d]); }
            // Independently assemble the weak residual of the first derivative.
            Vector weak(fe->GetDof());weak=0.0;
            Vector shape(fe->GetDof()),x;DenseMatrix grad(fe->GetDof(),2);
            const auto &rule=ofdg::detail::CurvedRule(fe->GetGeomType(),qorder+8);
            for(int q=0;q<rule.GetNPoints();++q) {
                const auto &ip=rule.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);
                fe->CalcShape(ip,shape);fe->CalcPhysDShape(*map,grad);
                double direct=0.0;for(int j=0;j<fe->GetDof();++j){direct+=grad(j,0)*value(j);}
                double weight=ip.weight*map->Weight();
                weak.Add(weight*(shape*repeated[0]-direct),shape);
                double exact[3]={std::cos(x(0)),-std::sin(x(0)),-std::cos(x(0))};
                for(int d=0;d<3;++d){double error=shape*repeated[d]-exact[d];errors[d]+=weight*error*error;}
            }
            Require(weak.Normlinf()<2e-12,"first derivative weak projection residual");
        }
        std::cout<<"derivatives n="<<n;
        for(int d=0;d<3;++d){errors[d]=std::sqrt(errors[d]);std::cout<<" D"<<d+1<<"="<<errors[d];
            if(n>2){Require(errors[d]<previous[d]/2.5,"projected derivative convergence");}
            previous[d]=errors[d];}
        std::cout<<'\n';
    }
}

void ValidatePositivity(Mesh &mesh)
{
    const int dim=mesh.Dimension();
    DG_FECollection collection(2,dim,BasisType::GaussLobatto);
    FiniteElementSpace space(&mesh,&collection,dim+2,Ordering::byNODES);
    GridFunction state(&space);state=0.0;
    Array<int>dofs;
    for(int e=0;e<space.GetNE();++e){space.GetElementDofs(e,dofs);
        for(int i=0;i<dofs.Size();++i){state(space.DofToVDof(dofs[i],0))=1.0;
            state(space.DofToVDof(dofs[i],dim+1))=2.5;}}
    space.GetElementDofs(0,dofs);
    state(space.DofToVDof(dofs[0],0))=-0.2;
    Vector before=Integrals(space,state);
    EulerPositivityLimiter limiter(&space,1.4,1e-10,1e-10);
    auto diagnostics=limiter.Apply(state);
    Vector defect=Integrals(space,state);defect-=before;
    Require(diagnostics.inadmissible_means==0 && diagnostics.limited_elements>0,
            "curved positivity activation");
    Require(defect.Normlinf()<2e-10,"curved positivity mean drift");
    Require(limiter.Apply(state).inadmissible_means==0,"curved positivity inadmissible result");
}

void ValidateContinuousState(Mesh &mesh)
{
    DG_FECollection collection(3,mesh.Dimension(),BasisType::GaussLobatto);
    FiniteElementSpace space(&mesh,&collection);
    GridFunction state(&space);
    FunctionCoefficient linear([](const Vector &x){return 1+x(0)+0.5*x(1);});
    state.ProjectCoefficient(linear);
    ofdg::OFDG filter(&space,BasisType::GaussLobatto);
    Vector residual;filter.ComputeStabilization(state,residual);
    Require(residual.Normlinf()<2e-9,"continuous physical linear state has derivative jumps");
}

void ValidateFaceMeasure()
{
    Mesh mesh=PolynomialMesh(2,Element::QUADRILATERAL);
    bool distinguishes_reference_weights=false;
    for(int f=0;f<mesh.GetNumFaces();++f) {
        int e1,e2;mesh.GetFaceElements(f,&e1,&e2);if(e2<0){continue;}
        auto *face=mesh.GetFaceElementTransformations(f);
        const auto &rule=ofdg::detail::FaceRule(2,*face);
        Vector weights=ofdg::detail::NormalizedPhysicalFaceWeights(*face,rule);
        Require(std::abs(weights.Sum()-1)<2e-13,"face weights do not normalize");
        double moment=0,reference=0,reference_measure=0;Vector x;
        for(int q=0;q<rule.GetNPoints();++q){
            face->Face->Transform(rule.IntPoint(q),x);
            double value=x(0)*x(0)+x(1)*x(1);
            moment+=weights(q)*value;
            reference+=rule.IntPoint(q).weight*value;
            reference_measure+=rule.IntPoint(q).weight;
        }
        const auto &fine=IntRules.Get(Geometry::SEGMENT,50);
        double oracle=0,area=0;
        for(int q=0;q<fine.GetNPoints();++q){
            face->SetAllIntPoints(&fine.IntPoint(q));face->Face->Transform(fine.IntPoint(q),x);
            double weight=fine.IntPoint(q).weight*face->Face->Weight();
            area+=weight;oracle+=weight*(x(0)*x(0)+x(1)*x(1));
        }
        Require(std::abs(moment-oracle/area)<2e-12,"physical face moment mismatch");
        distinguishes_reference_weights |= std::abs(moment-reference/reference_measure)>1e-6;
    }
    Require(distinguishes_reference_weights,"face fixture does not distinguish physical weights");
}

void ValidateParameterizations()
{
    // The physical cells are identical squares, but one map has an interior
    // bubble. The physical linear field belongs to both mapped spaces.
    Mesh straight=Mesh::MakeCartesian2D(2,2,Element::QUADRILATERAL,true);
    Mesh reparameterized(straight);
    reparameterized.SetCurvature(3,true);
    auto *nodes=reparameterized.GetNodes();auto *space=nodes->FESpace();
    Array<int>dofs;
    for(int e=0;e<space->GetNE();++e){space->GetElementDofs(e,dofs);
        const auto &points=space->GetFE(e)->GetNodes();
        for(int i=0;i<dofs.Size();++i){const auto &ip=points.IntPoint(i);
            (*nodes)(space->DofToVDof(dofs[i],0)) += 0.1*ip.x*(1-ip.x)*ip.y*(1-ip.y);}}
    reparameterized.NodesUpdated();
    ValidateContinuousState(straight);ValidateContinuousState(reparameterized);
    ValidateFilter(reparameterized,3,"equivalent square parameterization");
}

int main(int argc,char **argv)
{
    MPI_Init(&argc,&argv);
    std::cout << std::unitbuf;
    try {
        int ranks;MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        for(auto type:{Element::TRIANGLE,Element::QUADRILATERAL,Element::TETRAHEDRON,Element::HEXAHEDRON,Element::WEDGE}) {
            int dim=(type==Element::TRIANGLE||type==Element::QUADRILATERAL)?2:3;
            Mesh mesh=PolynomialMesh(dim,type);
            std::string name="curved geometry "+std::to_string(type);
            if(ranks==1){for(int p:{1,2,3}){ValidateFilter(mesh,p,name);}
                ValidateContinuousState(mesh); ValidatePositivity(mesh);}
            ValidateMPI(mesh,name);
        }
        for(const char *file:{"disc-nurbs.mesh","ball-nurbs.mesh"}) {
            Mesh mesh((std::string(OFDG_MFEM_DATA_DIR)+"/"+file).c_str(),1,1);
            {
                DG_FECollection collection(2,mesh.Dimension());
                FiniteElementSpace space(&mesh,&collection);
                bool rejected=false;
                try { ofdg::OFDG filter(&space,BasisType::GaussLegendre); }
                catch(const std::invalid_argument &) { rejected=true; }
                Require(rejected,"native NURBS must report its support boundary");
            }
            mesh.SetCurvature(3); // Explicit polynomial approximation, not exact NURBS.
            if(ranks==1){ValidateFilter(mesh,2,std::string("projected ")+file); ValidatePositivity(mesh);}
            ValidateMPI(mesh,std::string("projected ")+file);
        }
        if(ranks==1){ValidateDerivatives();ValidateParameterizations();ValidateFaceMeasure();}
    } catch(const std::exception &e) {
        std::cerr<<"FAILED: "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();
}
