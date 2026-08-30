#pragma once

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace ofdg::detail
{

/** Snapshot and validation of the geometry contract supported by the filters. */
class SupportedSpaceContract
{
private:
   const mfem::FiniteElementSpace *space = nullptr;
   long mesh_sequence = -1;
   long nodes_sequence = -1;
   long space_sequence = -1;
   int elements = 0;
   int faces = 0;
   int vector_size = 0;

   static bool IsAffine(mfem::ElementTransformation &transformation,
                        mfem::Geometry::Type geometry)
   {
      const int quadrature_order =
         std::max(4, 2 * transformation.OrderJ() + 2);
      const mfem::IntegrationRule &rule =
         mfem::IntRules.Get(geometry, quadrature_order);
      const mfem::IntegrationPoint &center =
         mfem::Geometries.GetCenter(geometry);
      transformation.SetIntPoint(&center);
      mfem::DenseMatrix reference(transformation.Jacobian());
      const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
         std::max(1.0, reference.FNorm());

      for (int q = 0; q < rule.GetNPoints(); ++q)
      {
         transformation.SetIntPoint(&rule.IntPoint(q));
         mfem::DenseMatrix difference(transformation.Jacobian());
         difference -= reference;
         if (difference.FNorm() > tolerance) { return false; }
      }
      return true;
   }

public:
   SupportedSpaceContract() = default;

   SupportedSpaceContract(const mfem::FiniteElementSpace *space_,
                          const char *consumer)
      : space(space_)
   {
      MFEM_VERIFY(space != nullptr,
                  consumer << " requires a valid finite element space.");
      mfem::Mesh *mesh = space->GetMesh();
      MFEM_VERIFY(mesh != nullptr,
                  consumer << " finite element space has no mesh.");
      MFEM_VERIFY(space->GetNE() > 0,
                  consumer << " requires at least one mesh element.");
      MFEM_VERIFY(mesh->Dimension() >= 1 && mesh->Dimension() <= 3,
                  consumer << " supports dimensions 1, 2, and 3.");
      MFEM_VERIFY(mesh->SpaceDimension() == mesh->Dimension(),
                  consumer << " requires a full-dimensional mesh.");
      MFEM_VERIFY(space->GetVDim() >= 1,
                  consumer << " requires at least one solution component.");
      MFEM_VERIFY(!space->IsVariableOrder(),
                  consumer << " does not yet support variable polynomial order.");

      const mfem::FiniteElement *first = space->GetFE(0);
      MFEM_VERIFY(first != nullptr,
                  consumer << " encountered a null finite element.");
      MFEM_VERIFY(first->GetOrder() >= 1,
                  consumer << " requires polynomial degree at least one.");
      MFEM_VERIFY(first->GetRangeType() == mfem::FiniteElement::SCALAR,
                  consumer << " requires scalar finite elements repeated over components.");
      MFEM_VERIFY(first->GetMapType() == mfem::FiniteElement::VALUE,
                  consumer << " requires VALUE-mapped scalar finite elements.");

      for (int element = 0; element < space->GetNE(); ++element)
      {
         const mfem::FiniteElement *finite_element = space->GetFE(element);
         MFEM_VERIFY(finite_element != nullptr,
                     consumer << " encountered a null finite element at element "
                              << element << '.');
         MFEM_VERIFY(finite_element->GetGeomType() == first->GetGeomType(),
                     consumer << " does not yet support meshes containing multiple "
                              << "element geometries (element " << element << ").");
         MFEM_VERIFY(finite_element->GetOrder() == first->GetOrder(),
                     consumer << " does not yet support variable polynomial order "
                              << "(element " << element << ").");
         MFEM_VERIFY(finite_element->GetDof() == first->GetDof() &&
                        finite_element->GetDim() == first->GetDim(),
                     consumer << " requires a uniform element-local DOF layout.");
         MFEM_VERIFY(finite_element->GetRangeType() == mfem::FiniteElement::SCALAR &&
                        finite_element->GetMapType() == mfem::FiniteElement::VALUE,
                     consumer << " requires VALUE-mapped scalar finite elements.");
         mfem::ElementTransformation *transformation =
            mesh->GetElementTransformation(element);
         MFEM_VERIFY(transformation != nullptr,
                     consumer << " encountered a null element transformation.");
         MFEM_VERIFY(IsAffine(*transformation, finite_element->GetGeomType()),
                     consumer << " does not yet support curvilinear or non-affine "
                              << "element mappings (element " << element << ").");
      }

      mesh_sequence = mesh->GetSequence();
      nodes_sequence = mesh->GetNodesSequence();
      space_sequence = space->GetSequence();
      elements = space->GetNE();
      faces = mesh->GetNumFaces();
      vector_size = space->GetVSize();
   }

   void VerifyUnchanged(const char *consumer) const
   {
      MFEM_VERIFY(space != nullptr,
                  consumer << " has no validated finite element space.");
      const mfem::Mesh *mesh = space->GetMesh();
      MFEM_VERIFY(mesh != nullptr,
                  consumer << " finite element space has no mesh.");
      MFEM_VERIFY(mesh->GetSequence() == mesh_sequence &&
                     mesh->GetNodesSequence() == nodes_sequence &&
                     space->GetSequence() == space_sequence &&
                     space->GetNE() == elements &&
                     mesh->GetNumFaces() == faces &&
                     space->GetVSize() == vector_size,
                  consumer << " caches are invalid because the mesh, geometry, or "
                           << "finite element space changed after construction. "
                           << "Reconstruct the filter after modifying the space.");
   }
};

} // namespace ofdg::detail
