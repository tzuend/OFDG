#include "mfem.hpp"
#include "../examples/euler/euler.hpp"
#include "../src/ofdg.hpp"
#include "../src/kxrcf.hpp"
#include "../src/euler_positivity.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace mfem;

Vector Totals(ParFiniteElementSpace &space, const Vector &state)
{
    GridFunction field(&space);field=state;
    Vector total(4);total=0.0;
    for(int e=0;e<space.GetNE();++e){auto *map=space.GetElementTransformation(e);
        const auto &rule=ofdg::detail::VolumeRule(2,*map);
        for(int q=0;q<rule.GetNPoints();++q){const auto &ip=rule.IntPoint(q);map->SetIntPoint(&ip);
            for(int c=0;c<4;++c){total(c)+=ip.weight*map->Weight()*field.GetValue(e,ip,c+1);}}}
    MPI_Allreduce(MPI_IN_PLACE,total.GetData(),4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return total;
}

int main(int argc,char **argv)
{
    MPI_Init(&argc,&argv);
    try {
        Mesh original=Mesh::MakeCartesian2D(4,4,Element::QUADRILATERAL,true);
        original.SetCurvature(3);
        original.Transform([](const Vector &x,Vector &y){y=x;double b=.04*std::sin(2*M_PI*x(0))*std::sin(2*M_PI*x(1));y(0)+=b;y(1)+=b;});
        std::vector<Vector> translations={Vector({1.,0.}),Vector({0.,1.})};
        Mesh serial=Mesh::MakePeriodic(original,original.CreatePeriodicVertexMapping(translations));
        ParMesh mesh(MPI_COMM_WORLD,serial);
        DG_FECollection collection(2,2,BasisType::GaussLegendre);
        ParFiniteElementSpace space(&mesh,&collection,4,Ordering::byNODES);
        ParGridFunction state(&space);
        VectorFunctionCoefficient initial(4,[](const Vector &x,Vector &u){
            double rho=1+.1*std::sin(2*M_PI*x(0));u(0)=rho;u(1)=.3*rho;u(2)=-.2*rho;u(3)=2.5+.065*rho;});
        state.ProjectCoefficient(initial);
        Vector before=Totals(space,state);
        EulerFlux flux(2,1.4);RusanovFlux numerical_flux(flux);
        DGHyperbolicConservationLaws evolution(space,std::make_unique<HyperbolicFormIntegrator>(numerical_flux,1),true);
        auto physics=std::make_shared<ofdg::EulerFacePhysics>(2,1.4);
        ofdg::OFDG filter(&space,BasisType::GaussLegendre,physics);
        ofdg::KXRCFIndicator indicator(&space,physics);
        EulerPositivityLimiter limiter(&space,1.4);
        RK4Solver solver;solver.Init(evolution);
        double time=0.0;Array<bool>active;Vector filtered;
        for(int step=0;step<20;++step){double dt=.0005;solver.Step(state,time,dt);
            indicator.Compute(state,active);filter.CompDecay(state,filtered,dt,&active);state=filtered;
            if(limiter.Apply(state).inadmissible_means){throw std::runtime_error("Euler inadmissible mean");}}
        Vector drift=Totals(space,state);drift-=before;
        if(!std::isfinite(state.Norml2()) || drift.Normlinf()>1e-9){throw std::runtime_error("Euler curved conservation");}
        if(Mpi::Root()){std::cout<<"curved Euler: 20 steps, conservation defect="<<drift.Normlinf()<<'\n';}
    } catch(const std::exception &e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
