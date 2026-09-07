#include "../src/profile_output.hpp"
#include "../examples/euler/euler.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

void Require(bool condition, const char *message) {
   if (!condition) { throw std::runtime_error(message); }
}
int main(int argc, char **argv) {
   mfem::Mpi::Init(argc, argv);
   mfem::Hypre::Init();
   // x(xi) is quadratic, hence x^2 is represented exactly at P4.
   auto mesh = mfem::Mesh::MakeCartesian1D(1);
   mesh.SetCurvature(2);
   mfem::VectorFunctionCoefficient warp(1, [](const mfem::Vector &x, mfem::Vector &y) {
      y.SetSize(1); y(0)=x(0)+.2*x(0)*(1-x(0));
   });
   mesh.Transform(warp);
   mfem::DG_FECollection fec(4,1);
   mfem::FiniteElementSpace fes(&mesh,&fec);
   mfem::GridFunction field(&fes);
   mfem::FunctionCoefficient square([](const mfem::Vector &x) {return x(0)*x(0);});
   field.ProjectCoefficient(square);
   auto *map=mesh.GetElementTransformation(0);
   auto average=PhysicalCellAverage(*map,4,1,[&](const mfem::IntegrationPoint &ip,mfem::Vector &v) {v(0)=field.GetValue(0,ip);});
   Require(std::abs(average(0)-1./3)<1e-13,"Incorrect physical mean on curved element");
   Require(std::abs(field.GetValue(0,mfem::Geometries.GetCenter(mfem::Geometry::SEGMENT))-average(0))>.01,"Fixture fails to distinguish center from average");
   WriteScalarSamples(field,"preview_export_check",0);
   // A transmissive face gives exactly the physical flux for the interior state,
   // regardless of an intentionally different supplied exterior coefficient.
   auto line=mfem::Mesh::MakeCartesian1D(1);
   mfem::DG_FECollection linear(1,1);
   mfem::FiniteElementSpace space(&line,&linear);
   mfem::EulerFlux flux(1,1.4);
   mfem::RusanovFlux numerical(flux);
   mfem::Vector state(3);state(0)=1.;state(1)=.2;state(2)=2.52;
   mfem::VectorConstantCoefficient same(state);
   mfem::Vector wrong(state);wrong(0)=2.;
   mfem::VectorConstantCoefficient different(wrong);
   mfem::EulerBoundaryIntegrator outflow(numerical,&different,false,1,true);
   mfem::EulerBoundaryIntegrator matching(numerical,&same,false,1);
   mfem::EulerBoundaryIntegrator fixed(numerical,&different,false,1);
   const auto &fe=*space.GetFE(0);
   mfem::Vector values(fe.GetDof()*3),a,b,c;
   for(int component=0;component<3;++component)
      for(int j=0;j<fe.GetDof();++j)values(component*fe.GetDof()+j)=state(component);
   for(int face=0;face<line.GetNBE();++face) {
      auto &trans=*line.GetBdrFaceTransformations(face);
      outflow.AssembleFaceVector(fe,fe,trans,values,a);
      matching.AssembleFaceVector(fe,fe,trans,values,b);
      fixed.AssembleFaceVector(fe,fe,trans,values,c);
      a-=b; c-=b;
      Require(a.Norml2()<1e-13,"Outflow does not match physical boundary flux");
      Require(c.Norml2()>.01,"Fixed boundary behavior was lost");
   }
   std::cout<<"Preview physical-average and outflow checks passed\n";
}
