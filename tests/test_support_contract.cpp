#include "mfem.hpp"

#include "../src/ofdg.hpp"

#include <iostream>
#include <string>

using namespace mfem;

namespace
{

void ConstructFilter(FiniteElementSpace &space)
{
   ofdg::OFDG filter(&space, BasisType::GaussLegendre);
}

} // namespace

int main(int argc, char *argv[])
{
   MFEM_VERIFY(argc >= 2, "Expected a support-contract test mode.");
   const std::string mode = argv[1];

   if (mode == "supported")
   {
      Mesh mesh = Mesh::MakeCartesian2D(
         2, 2, Element::QUADRILATERAL, true, 1.0, 1.0);
      DG_FECollection collection(2, 2);
      FiniteElementSpace space(&mesh, &collection);
      ConstructFilter(space);
      return 0;
   }

   if (mode == "mixed" || mode == "curved")
   {
      MFEM_VERIFY(argc == 3, "Expected the MFEM mesh path.");
      Mesh mesh(argv[2], 1, 1);
      DG_FECollection collection(2, mesh.Dimension());
      FiniteElementSpace space(&mesh, &collection);
      ConstructFilter(space);
      return 0;
   }

   if (mode == "variable-order")
   {
      Mesh mesh = Mesh::MakeCartesian1D(3, 1.0);
      DG_FECollection collection(2, 1);
      FiniteElementSpace space(&mesh, &collection);
      space.SetElementOrder(1, 1);
      ConstructFilter(space);
      return 0;
   }

   if (mode == "modified")
   {
      Mesh mesh = Mesh::MakeCartesian1D(3, 1.0);
      DG_FECollection collection(2, 1);
      FiniteElementSpace space(&mesh, &collection);
      ofdg::OFDG filter(&space, BasisType::GaussLegendre);
      GridFunction state(&space);
      state = 1.0;
      mesh.UniformRefinement();
      Vector means;
      filter.ComputeMean(state, means);
      return 0;
   }

   std::cerr << "Unknown support-contract mode: " << mode << '\n';
   return 2;
}
