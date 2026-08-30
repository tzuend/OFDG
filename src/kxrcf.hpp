#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "face_physics.hpp"
#include "support_contract.hpp"

namespace ofdg
{

using mfem::Array;
using mfem::CalcOrtho;
using mfem::DenseMatrix;
using mfem::ElementTransformation;
using mfem::FaceElementTransformations;
using mfem::FiniteElement;
using mfem::FiniteElementSpace;
using mfem::Geometries;
using mfem::Geometry;
using mfem::IntegrationPoint;
using mfem::IntegrationRule;
using mfem::IntRules;
using mfem::Mesh;
using mfem::MPITypeMap;
using mfem::ParFiniteElementSpace;
using mfem::ParGridFunction;
using mfem::ParMesh;
using mfem::real_t;
using mfem::Vector;
using mfem::VectorCoefficient;

namespace detail
{

struct KXRCFInternalTiming
{
   double element_scales = 0.0;
   double face_integration = 0.0;
   double normalization = 0.0;
   double total = 0.0;
};

#ifdef KXRCF_INTERNAL_TIMING
using KXRCFTimingClock = std::chrono::steady_clock;
inline double KXRCFSecondsSince(const KXRCFTimingClock::time_point &begin)
{
   return std::chrono::duration<double>(KXRCFTimingClock::now() - begin).count();
}
#endif

} // namespace detail


/**
 * Multidimensional / multicomponent KXRCF troubled-cell indicator.
 *
 * For each element K and component c,
 *
 *                | int_{Gamma_K^-} (u_K,c - u_N,c) ds |
 * I_{K,c} = --------------------------------------------------
 *             h_K^((k+1)/2) |Gamma_K^-| ||u_K,c||
 *
 * where
 *
 *   Gamma_K^- = { x in dK : v(x) . n_K(x) < 0 }.
 *
 * The final element indicator is
 *
 *   I_K = max_c I_{K,c},
 *
 * and K is marked troubled if
 *
 *   I_K > threshold.
 *
 *
 * Normalization:
 *
 *   dim == 1:
 *       ||u_K,c|| = |cell average|
 *
 *   dim >= 2:
 *       ||u_K,c|| = max_q |u_K,c(x_q)|
 *
 *
 * Element radius h:
 *
 *   dim == 1:
 *       h = element_length / 2
 *
 *   dim >= 2:
 *       h is generalized to the largest distance from the mapped
 *       reference-element center to any mapped reference vertex.
 *
 * This gives the exact old 1D definition and a geometry-independent
 * circumscribing radius for arbitrary MFEM element types.
 *
 *
 * Boundary treatment:
 *
 *   Only faces with two neighboring elements are considered.
 *   Physical boundary faces are ignored.
 *
 *
 * Requirements:
 *
 *   - Scalar finite elements, possibly repeated through vdim.
 *   - VALUE mapping.
 *   - Mesh dimensions 1, 2, or 3.
 *   - Physical space dimension equal to mesh dimension.
 *   - A common VectorCoefficient velocity.
 */
class KXRCFIndicator
{
private:
   FiniteElementSpace *fes;
   detail::SupportedSpaceContract contract;
   std::shared_ptr<const FacePhysics> face_physics;

   int dim;
   int ncomp;

   real_t threshold;
   real_t relative_scale_floor;
   struct ElementCache
   {
      Array<int> vdofs;
      DenseMatrix evaluation;
      Vector physical_weights;
      real_t volume = 0.0;
      real_t radius = 0.0;
   };
   std::vector<ElementCache> element_cache;
   struct LocalFaceCache
   {
      int face = -1;
      int element1 = -1;
      int element2 = -1;
      Array<int> vdofs1;
      Array<int> vdofs2;
      DenseMatrix evaluation1;
      DenseMatrix evaluation2;
      DenseMatrix unit_normals;
      Vector physical_weights;
      const IntegrationRule *rule = nullptr;
   };
   std::vector<LocalFaceCache> local_face_cache;
   mutable detail::KXRCFInternalTiming internal_timing;

   struct FaceScratch
   {
      Vector shape1;
      Vector shape2;
      Vector state1;
      Vector state2;
      Vector normal;
   };
   mutable FaceScratch face_scratch;


   /**
    * Evaluate all components of an element-local state at one
    * integration point.
    *
    * GetElementVDofs() returns the local vector in byNODES ordering:
    *
    *   dof 0: component 0, component 1, ...
    *   dof 1: component 0, component 1, ...
    *
    * Therefore
    *
    *   local_index = j * ncomp + c.
    */
   void EvaluateElementState(const Vector &local_state,
                             const Vector &shape,
                             Vector &value) const;


   /**
    * Compute the normalization scale S_{K,c}.
    *
    * 1D:
    *
    *   S_{K,c} = | average_K u_c |
    *
    * 2D / 3D:
    *
    *   S_{K,c} = max_q |u_c(x_q)|.
    *
    * Output is indexed as
    *
    *   scales(c, e).
    */
   void ComputeElementScales(const Vector &x,
                             DenseMatrix &scales) const;


   /**
    * Compute the generalized KXRCF element radius.
    *
    * In 1D this is exactly half the element length, matching the
    * original implementation.
    *
    * In 2D/3D we map the center and vertices of the reference
    * element and take
    *
    *   h_K = max_vertex |x_vertex - x_center|.
    *
    * This is well defined for triangles, quads, tetrahedra,
    * hexahedra, prisms, pyramids, etc.
    */
   real_t ComputeElementRadius(int e) const;
   void BuildElementCache();
   void BuildLocalFaceCache();

   void AccumulateCachedFace(const LocalFaceCache &cache,
                             FaceElementTransformations &transformations,
                             const Vector &local_state1,
                             const Vector &local_state2,
                             DenseMatrix &jump_integral,
                             Vector &inflow_measure) const;

   void AccumulateFace(FaceElementTransformations &transformations,
                       const FiniteElement &element1,
                       const FiniteElement &element2,
                       const Vector &local_state1,
                       const Vector &local_state2,
                       int element1_number,
                       int element2_number,
                       bool accumulate_element2,
                       DenseMatrix &jump_integral,
                       Vector &inflow_measure) const;


public:
   KXRCFIndicator(FiniteElementSpace *fes_,
                  VectorCoefficient *velocity_,
                  real_t threshold_ = 1.0,
                  real_t relative_scale_floor_ = 1e-14)
      : KXRCFIndicator(
           fes_,
           std::make_shared<AdvectionFacePhysics>(velocity_),
           threshold_,
           relative_scale_floor_)
   {
      MFEM_VERIFY(velocity_->GetVDim() == dim,
                  "KXRCF velocity dimension must match mesh dimension.");
   }

   KXRCFIndicator(FiniteElementSpace *fes_,
                  std::shared_ptr<const FacePhysics> face_physics_,
                  real_t threshold_ = 1.0,
                  real_t relative_scale_floor_ = 1e-14)
      : fes(fes_),
        contract(fes_, "KXRCF"),
        face_physics(std::move(face_physics_)),
        dim(0),
        ncomp(0),
        threshold(threshold_),
        relative_scale_floor(relative_scale_floor_)
   {
      MFEM_VERIFY(fes != NULL,
                  "KXRCF requires a valid finite element space.");

      MFEM_VERIFY(face_physics != nullptr,
                  "KXRCF requires valid face physics.");

      Mesh *mesh = fes->GetMesh();

      MFEM_VERIFY(mesh != NULL,
                  "KXRCF finite element space has no mesh.");

      dim = mesh->Dimension();
      ncomp = fes->GetVDim();

      MFEM_VERIFY(dim >= 1 && dim <= 3,
                  "KXRCF only supports dimensions 1, 2, and 3.");

      /*
       * CalcOrtho() below assumes a standard full-dimensional
       * computational domain.
       */
      MFEM_VERIFY(mesh->SpaceDimension() == dim,
                  "KXRCF currently requires mesh dimension and "
                  "physical space dimension to agree.");

      MFEM_VERIFY(ncomp >= 1,
                  "KXRCF requires at least one solution component.");

      MFEM_VERIFY(threshold >= 0.0,
                  "KXRCF threshold must be non-negative.");

      MFEM_VERIFY(relative_scale_floor >= 0.0,
                  "KXRCF scale floor must be non-negative.");

      /*
       * This implementation represents a system as vdim copies
       * of a scalar finite element.
       *
       * This is exactly the layout normally used by DG systems
       * such as Euler in MFEM.
       */
      for (int e = 0; e < fes->GetNE(); ++e)
      {
         const FiniteElement *fe = fes->GetFE(e);

         MFEM_VERIFY(fe != NULL,
                     "KXRCF encountered a null finite element.");

         MFEM_VERIFY(fe->GetRangeType() == FiniteElement::SCALAR,
                     "KXRCF expects scalar finite elements repeated "
                     "over the vector dimension.");

         MFEM_VERIFY(fe->GetMapType() == FiniteElement::VALUE,
                     "KXRCF currently expects VALUE-mapped scalar "
                     "finite elements.");
      }

      BuildElementCache();
      BuildLocalFaceCache();
   }

   void ResetInternalTimings() const
   {
      internal_timing = detail::KXRCFInternalTiming{};
   }

   void PrintInternalTimings(std::ostream &out = std::cout) const
   {
#ifdef KXRCF_INTERNAL_TIMING
      out << "KXRCF internal timings\n"
          << "  element scales:   " << internal_timing.element_scales << " s\n"
          << "  face integration: " << internal_timing.face_integration << " s\n"
          << "  normalization:    " << internal_timing.normalization << " s\n"
          << "  total:            " << internal_timing.total << " s\n";
#else
      out << "KXRCF internal timing is disabled. "
          << "Compile with -DKXRCF_INTERNAL_TIMING to enable it.\n";
#endif
   }


   /**
    * Compute the KXRCF troubled-cell mask.
    *
    * x:
    *     DG state vector.
    *
    * active:
    *     One bool per element. True iff at least one component has
    *     indicator > threshold.
    *
    * indicator_values:
    *     Optional pooled indicator
    *
    *         I_K = max_c I_{K,c}.
    *
    * component_indicator_values:
    *     Optional matrix containing all component-wise indicators:
    *
    *         (*component_indicator_values)(c, e) = I_{K,c}.
    */
   void Compute(
      const Vector &x,
      Array<bool> &active,
      Vector *indicator_values = NULL,
      DenseMatrix *component_indicator_values = NULL) const;
};


// --------------------------------------------------------------------------
// Local state evaluation
// --------------------------------------------------------------------------

void KXRCFIndicator::EvaluateElementState(
   const Vector &local_state,
   const Vector &shape,
   Vector &value) const
{
   const int ndof = shape.Size();

   MFEM_VERIFY(
      local_state.Size() == ndof * ncomp,
      "Unexpected local vector size in KXRCF.");

   value.SetSize(ncomp);
   value = 0.0;

   /*
    * GetElementVDofs() always returns the element-local vector
    * DOFs in MFEM's byNODES ordering, independently of the
    * global ordering of the finite element space.
    *
    * MFEM byNODES layout is:
    *
    *   component 0: dof 0, dof 1, ..., dof ndof-1
    *   component 1: dof 0, dof 1, ..., dof ndof-1
    *   ...
    *
    * Therefore:
    *
    *   local_index = c * ndof + j.
    */
   for (int c = 0; c < ncomp; ++c)
   {
      for (int j = 0; j < ndof; ++j)
      {
         value(c) +=
            shape(j) *
            local_state(c * ndof + j);
      }
   }
}

// --------------------------------------------------------------------------
// Element solution scales
// --------------------------------------------------------------------------

void KXRCFIndicator::ComputeElementScales(
   const Vector &x,
   DenseMatrix &scales) const
{
   const int ne = fes->GetNE();

   scales.SetSize(ncomp, ne);
   scales = 0.0;

   Vector local_state;
   DenseMatrix local_matrix;
   DenseMatrix values;

   for (int e = 0; e < ne; ++e)
   {
      const ElementCache &cache = element_cache[e];
      const int ndof = cache.evaluation.Width();
      local_state.SetSize(cache.vdofs.Size());
      x.GetSubVector(cache.vdofs, local_state);
      local_matrix.UseExternalData(local_state.GetData(), ndof, ncomp);
      values.SetSize(cache.evaluation.Height(), ncomp);
      Mult(cache.evaluation, local_matrix, values);

      if (dim == 1)
      {
         /*
          * Original 1D KXRCF normalization:
          *
          *     | element average |
          */
         for (int c = 0; c < ncomp; ++c)
         {
            real_t integral = 0.0;
            for (int q = 0; q < values.Height(); ++q)
            {
               integral += cache.physical_weights(q) * values(q, c);
            }
            scales(c, e) = std::abs(integral / cache.volume);
         }
      }
      else
      {
         /*
          * Multidimensional KXRCF normalization:
          *
          *     max_q |u(x_q)|
          */
         for (int q = 0; q < values.Height(); ++q)
         {
            for (int c = 0; c < ncomp; ++c)
            {
               scales(c, e) =
                  std::max(scales(c, e),
                           std::abs(values(q, c)));
            }
         }
      }
   }
}


void KXRCFIndicator::BuildElementCache()
{
   Mesh *mesh = fes->GetMesh();
   element_cache.resize(fes->GetNE());

   for (int e = 0; e < fes->GetNE(); ++e)
   {
      ElementCache &cache = element_cache[e];
      const FiniteElement *fe = fes->GetFE(e);
      ElementTransformation *transformation =
         mesh->GetElementTransformation(e);
      fes->GetElementVDofs(e, cache.vdofs);

      const IntegrationRule &rule =
         IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 2);
      cache.evaluation.SetSize(rule.GetNPoints(), fe->GetDof());
      cache.physical_weights.SetSize(rule.GetNPoints());
      cache.volume = 0.0;
      Vector shape(fe->GetDof());

      for (int q = 0; q < rule.GetNPoints(); ++q)
      {
         const IntegrationPoint &ip = rule.IntPoint(q);
         fe->CalcShape(ip, shape);
         for (int i = 0; i < shape.Size(); ++i)
         {
            cache.evaluation(q, i) = shape(i);
         }
         transformation->SetIntPoint(&ip);
         cache.physical_weights(q) = ip.weight * transformation->Weight();
         cache.volume += cache.physical_weights(q);
      }
      MFEM_VERIFY(cache.volume > 0.0,
                  "KXRCF encountered an element with non-positive volume.");
      cache.radius = ComputeElementRadius(e);
   }
}


void KXRCFIndicator::BuildLocalFaceCache()
{
   Mesh *mesh = fes->GetMesh();
   local_face_cache.clear();
   local_face_cache.reserve(mesh->GetNumFaces());

   for (int f = 0; f < mesh->GetNumFaces(); ++f)
   {
      FaceElementTransformations *transformations =
         mesh->GetFaceElementTransformations(f);
      if (!transformations || transformations->Elem1No < 0 ||
          transformations->Elem2No < 0)
      {
         continue;
      }

      local_face_cache.emplace_back();
      LocalFaceCache &cache = local_face_cache.back();
      cache.face = f;
      cache.element1 = transformations->Elem1No;
      cache.element2 = transformations->Elem2No;
      fes->GetElementVDofs(cache.element1, cache.vdofs1);
      fes->GetElementVDofs(cache.element2, cache.vdofs2);

      const FiniteElement *element1 = fes->GetFE(cache.element1);
      const FiniteElement *element2 = fes->GetFE(cache.element2);
      const int integration_order = dim == 1 ? 0 :
         2 * std::max(element1->GetOrder(), element2->GetOrder()) + 2;
      cache.rule = &IntRules.Get(transformations->FaceGeom,
                                 integration_order);
      const int point_count = cache.rule->GetNPoints();
      cache.evaluation1.SetSize(point_count, element1->GetDof());
      cache.evaluation2.SetSize(point_count, element2->GetDof());
      cache.unit_normals.SetSize(point_count, dim);
      cache.physical_weights.SetSize(point_count);
      Vector shape1(element1->GetDof());
      Vector shape2(element2->GetDof());
      Vector normal(dim);

      for (int q = 0; q < point_count; ++q)
      {
         const IntegrationPoint &ip = cache.rule->IntPoint(q);
         transformations->SetAllIntPoints(&ip);
         const IntegrationPoint &ip1 =
            transformations->GetElement1IntPoint();
         element1->CalcShape(ip1, shape1);
         element2->CalcShape(transformations->GetElement2IntPoint(), shape2);
         for (int i = 0; i < shape1.Size(); ++i)
         {
            cache.evaluation1(q, i) = shape1(i);
         }
         for (int i = 0; i < shape2.Size(); ++i)
         {
            cache.evaluation2(q, i) = shape2(i);
         }

         if (dim == 1)
         {
            transformations->Elem1->SetIntPoint(&ip1);
            const real_t orientation =
               transformations->Elem1->Jacobian()(0, 0) >= 0.0 ? 1.0 : -1.0;
            normal(0) = (ip1.x < 0.5 ? -1.0 : 1.0) * orientation;
            cache.physical_weights(q) = 1.0;
         }
         else
         {
            CalcOrtho(transformations->Jacobian(), normal);
            const real_t normal_norm = normal.Norml2();
            MFEM_VERIFY(normal_norm > 0.0,
                        "KXRCF encountered a degenerate face.");
            cache.physical_weights(q) = ip.weight * normal_norm;
            normal /= normal_norm;
         }
         for (int d = 0; d < dim; ++d)
         {
            cache.unit_normals(q, d) = normal(d);
         }
      }
   }
}


void KXRCFIndicator::AccumulateCachedFace(
   const LocalFaceCache &cache,
   FaceElementTransformations &transformations,
   const Vector &local_state1,
   const Vector &local_state2,
   DenseMatrix &jump_integral,
   Vector &inflow_measure) const
{
   Vector &state1 = face_scratch.state1;
   Vector &state2 = face_scratch.state2;
   Vector &normal = face_scratch.normal;
   state1.SetSize(ncomp);
   state2.SetSize(ncomp);
   normal.SetSize(dim);

   for (int q = 0; q < cache.rule->GetNPoints(); ++q)
   {
      state1 = 0.0;
      state2 = 0.0;
      for (int c = 0; c < ncomp; ++c)
      {
         for (int j = 0; j < cache.evaluation1.Width(); ++j)
         {
            state1(c) += cache.evaluation1(q, j) *
                         local_state1(c * cache.evaluation1.Width() + j);
         }
         for (int j = 0; j < cache.evaluation2.Width(); ++j)
         {
            state2(c) += cache.evaluation2(q, j) *
                         local_state2(c * cache.evaluation2.Width() + j);
         }
      }
      for (int d = 0; d < dim; ++d)
      {
         normal(d) = cache.unit_normals(q, d);
      }

      transformations.SetAllIntPoints(&cache.rule->IntPoint(q));
      const real_t transport =
         face_physics->Evaluate(state1, state2, normal,
                                transformations).normal_transport;
      const real_t weight = cache.physical_weights(q);

      if (transport < 0.0)
      {
         inflow_measure(cache.element1) += weight;
         for (int c = 0; c < ncomp; ++c)
         {
            jump_integral(c, cache.element1) +=
               weight * (state1(c) - state2(c));
         }
      }
      else if (transport > 0.0)
      {
         inflow_measure(cache.element2) += weight;
         for (int c = 0; c < ncomp; ++c)
         {
            jump_integral(c, cache.element2) +=
               weight * (state2(c) - state1(c));
         }
      }
   }
}


// --------------------------------------------------------------------------
// Generalized element radius
// --------------------------------------------------------------------------

real_t KXRCFIndicator::ComputeElementRadius(
   int e) const
{
   Mesh *mesh = fes->GetMesh();

   /*
    * Preserve the exact old 1D definition.
    */
   if (dim == 1)
   {
      const real_t h =
         0.5 * mesh->GetElementVolume(e);

      MFEM_VERIFY(h > 0.0,
                  "KXRCF encountered non-positive "
                  "1D element radius.");

      return h;
   }

   const Geometry::Type geom =
      mesh->GetElementGeometry(e);

   ElementTransformation *Tr =
      mesh->GetElementTransformation(e);

   /*
    * Reference element center and vertices.
    */
   const IntegrationPoint &center_ip =
      Geometries.GetCenter(geom);

   const IntegrationRule *vertex_rule =
      Geometries.GetVertices(geom);

   MFEM_VERIFY(vertex_rule != NULL,
               "KXRCF could not obtain reference vertices.");

   MFEM_VERIFY(vertex_rule->GetNPoints() > 0,
               "KXRCF element has no reference vertices.");

   Vector center(dim);
   Vector point(dim);

   Tr->Transform(center_ip, center);

   real_t radius = 0.0;

   for (int i = 0;
        i < vertex_rule->GetNPoints();
        ++i)
   {
      const IntegrationPoint &vip =
         vertex_rule->IntPoint(i);

      Tr->Transform(vip, point);

      real_t distance2 = 0.0;

      for (int d = 0; d < dim; ++d)
      {
         const real_t diff =
            point(d) - center(d);

         distance2 += diff * diff;
      }

      radius =
         std::max(radius,
                  std::sqrt(distance2));
   }

   MFEM_VERIFY(radius > 0.0,
               "KXRCF encountered non-positive "
               "element radius.");

   return radius;
}


void KXRCFIndicator::AccumulateFace(
   FaceElementTransformations &transformations,
   const FiniteElement &element1,
   const FiniteElement &element2,
   const Vector &local_state1,
   const Vector &local_state2,
   int element1_number,
   int element2_number,
   bool accumulate_element2,
   DenseMatrix &jump_integral,
   Vector &inflow_measure) const
{
   Vector &shape1 = face_scratch.shape1;
   Vector &shape2 = face_scratch.shape2;
   Vector &state1 = face_scratch.state1;
   Vector &state2 = face_scratch.state2;
   Vector &normal = face_scratch.normal;
   shape1.SetSize(element1.GetDof());
   shape2.SetSize(element2.GetDof());
   state1.SetSize(ncomp);
   state2.SetSize(ncomp);
   normal.SetSize(dim);

   const int integration_order = dim == 1 ? 0 :
      2 * std::max(element1.GetOrder(), element2.GetOrder()) + 2;
   const IntegrationRule &rule =
      IntRules.Get(transformations.FaceGeom, integration_order);

   for (int q = 0; q < rule.GetNPoints(); ++q)
   {
      const IntegrationPoint &ip = rule.IntPoint(q);
      transformations.SetAllIntPoints(&ip);
      const IntegrationPoint &ip1 = transformations.GetElement1IntPoint();
      const IntegrationPoint &ip2 = transformations.GetElement2IntPoint();

      element1.CalcShape(ip1, shape1);
      element2.CalcShape(ip2, shape2);
      EvaluateElementState(local_state1, shape1, state1);
      EvaluateElementState(local_state2, shape2, state2);

      real_t physical_weight;
      if (dim == 1)
      {
         transformations.Elem1->SetIntPoint(&ip1);
         const DenseMatrix &jacobian = transformations.Elem1->Jacobian();
         const real_t orientation = jacobian(0, 0) >= 0.0 ? 1.0 : -1.0;
         normal(0) = (ip1.x < 0.5 ? -1.0 : 1.0) * orientation;
         physical_weight = 1.0;
      }
      else
      {
         CalcOrtho(transformations.Jacobian(), normal);
         const real_t normal_norm = normal.Norml2();
         MFEM_VERIFY(normal_norm > 0.0,
                     "KXRCF encountered a degenerate face.");
         physical_weight = ip.weight * normal_norm;
         normal /= normal_norm;
      }

      const real_t normal_transport =
         face_physics->Evaluate(state1, state2, normal,
                                transformations).normal_transport;

      if (normal_transport < 0.0)
      {
         inflow_measure(element1_number) += physical_weight;
         for (int c = 0; c < ncomp; ++c)
         {
            jump_integral(c, element1_number) +=
               physical_weight * (state1(c) - state2(c));
         }
      }
      else if (normal_transport > 0.0 && accumulate_element2)
      {
         inflow_measure(element2_number) += physical_weight;
         for (int c = 0; c < ncomp; ++c)
         {
            jump_integral(c, element2_number) +=
               physical_weight * (state2(c) - state1(c));
         }
      }
   }
}


// --------------------------------------------------------------------------
// Main KXRCF computation
// --------------------------------------------------------------------------

void KXRCFIndicator::Compute(
   const Vector &x,
   Array<bool> &active,
   Vector *indicator_values,
DenseMatrix *component_indicator_values) const
{
   contract.VerifyUnchanged("KXRCF");
#ifdef KXRCF_INTERNAL_TIMING
   const auto total_timer_begin = detail::KXRCFTimingClock::now();
#endif
   MFEM_VERIFY(x.Size() == fes->GetVSize(),
               "KXRCF input size does not match "
               "the finite element space.");

   Mesh *mesh = fes->GetMesh();
   const int ne = fes->GetNE();

   // -----------------------------------------------------------------------
   // Output initialization
   // -----------------------------------------------------------------------

   active.SetSize(ne);
   active = false;

   if (indicator_values)
   {
      indicator_values->SetSize(ne);
      *indicator_values = 0.0;
   }

   if (component_indicator_values)
   {
      component_indicator_values->SetSize(ncomp, ne);
      *component_indicator_values = 0.0;
   }


   // -----------------------------------------------------------------------
   // 1. Compute solution normalization scales S_{K,c}
   // -----------------------------------------------------------------------

   DenseMatrix solution_scales;
#ifdef KXRCF_INTERNAL_TIMING
   const auto scales_timer_begin = detail::KXRCFTimingClock::now();
#endif
   ComputeElementScales(x, solution_scales);
#ifdef KXRCF_INTERNAL_TIMING
   internal_timing.element_scales +=
      detail::KXRCFSecondsSince(scales_timer_begin);
#endif


   /*
    * Numerical floors are computed separately for each component.
    *
    * This is important for systems whose components have very
    * different physical magnitudes.
    */
   Vector global_scales(ncomp);
   global_scales = 0.0;

   for (int c = 0; c < ncomp; ++c)
   {
      for (int e = 0; e < ne; ++e)
      {
         global_scales(c) =
            std::max(global_scales(c),
                     solution_scales(c, e));
      }
   }

#ifdef MFEM_USE_MPI
   if (const auto *parallel_space =
          dynamic_cast<const ParFiniteElementSpace *>(fes))
   {
      MPI_Allreduce(MPI_IN_PLACE, global_scales.GetData(), ncomp,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    parallel_space->GetComm());
   }
#endif

   bool any_nonzero_component = false;

   for (int c = 0; c < ncomp; ++c)
   {
      if (global_scales(c) > 0.0)
      {
         any_nonzero_component = true;
         break;
      }
   }

   /*
    * Every component is numerically zero according to the KXRCF
    * normalization. Nothing can be marked.
    *
    * This reproduces the behavior of the old 1D implementation.
    */
   if (!any_nonzero_component)
   {
      return;
   }

   Vector scale_floors(ncomp);

   for (int c = 0; c < ncomp; ++c)
   {
      scale_floors(c) =
         relative_scale_floor *
         global_scales(c);
   }


   // -----------------------------------------------------------------------
   // 2. Integrate signed jumps over inflow boundaries
   // -----------------------------------------------------------------------

   /*
    * jump_integral(c,e)
    *
    *     = integral_{Gamma_e^-}
    *         (u_e,c - u_neighbor,c) ds
    */
   DenseMatrix jump_integral(ncomp, ne);
   jump_integral = 0.0;

   /*
    * inflow_measure(e)
    *
    *     = |Gamma_e^-|
    */
   Vector inflow_measure(ne);
   inflow_measure = 0.0;


   Array<int> vdofs1;
   Array<int> vdofs2;

#ifdef KXRCF_INTERNAL_TIMING
   const auto faces_timer_begin = detail::KXRCFTimingClock::now();
#endif

   Vector local_state1;
   Vector local_state2;
   for (const LocalFaceCache &cache : local_face_cache)
   {
      FaceElementTransformations *FTr =
         mesh->GetFaceElementTransformations(cache.face);
      local_state1.SetSize(cache.vdofs1.Size());
      local_state2.SetSize(cache.vdofs2.Size());
      x.GetSubVector(cache.vdofs1, local_state1);
      x.GetSubVector(cache.vdofs2, local_state2);
      AccumulateCachedFace(cache, *FTr, local_state1, local_state2,
                           jump_integral, inflow_measure);
   }

#ifdef MFEM_USE_MPI
   if (auto *parallel_space = dynamic_cast<ParFiniteElementSpace *>(fes))
   {
      ParMesh *parallel_mesh = parallel_space->GetParMesh();
      ParGridFunction parallel_state(parallel_space);
      parallel_state = x;
      parallel_state.ExchangeFaceNbrData();

      for (int shared_face = 0;
           shared_face < parallel_mesh->GetNSharedFaces(); ++shared_face)
      {
         FaceElementTransformations *transformations =
            parallel_mesh->GetSharedFaceTransformations(shared_face, true);
         const int local_element = transformations->Elem1No;
         const int neighbor_element =
            transformations->Elem2No - parallel_mesh->GetNE();

         parallel_space->GetElementVDofs(local_element, vdofs1);
         parallel_space->GetFaceNbrElementVDofs(neighbor_element, vdofs2);
         Vector local_state(vdofs1.Size());
         Vector neighbor_state(vdofs2.Size());
         x.GetSubVector(vdofs1, local_state);
         parallel_state.FaceNbrData().GetSubVector(vdofs2, neighbor_state);

         AccumulateFace(*transformations,
                        *parallel_space->GetFE(local_element),
                        *parallel_space->GetFaceNbrFE(neighbor_element),
                        local_state, neighbor_state,
                        local_element, -1, false,
                        jump_integral, inflow_measure);
      }
   }
#endif

#ifdef KXRCF_INTERNAL_TIMING
   internal_timing.face_integration +=
      detail::KXRCFSecondsSince(faces_timer_begin);
   const auto normalization_timer_begin = detail::KXRCFTimingClock::now();
#endif


   // -----------------------------------------------------------------------
   // 3. Construct normalized indicators
   // -----------------------------------------------------------------------

   for (int e = 0; e < ne; ++e)
   {
      /*
       * No inflow boundary -> indicator remains zero.
       */
      if (inflow_measure(e) == 0.0)
      {
         continue;
      }

      const FiniteElement *fe =
         fes->GetFE(e);

      const int k =
         fe->GetOrder();

      const real_t h = element_cache[e].radius;

      const real_t h_factor =
         std::pow(
            h,
            0.5 * static_cast<real_t>(k + 1));

      MFEM_VERIFY(h_factor > 0.0,
                  "Invalid KXRCF h scaling.");

      real_t pooled_indicator = 0.0;

      for (int c = 0;
           c < ncomp;
           ++c)
      {
         /*
          * This component is identically zero according to the
          * global KXRCF scale.
          */
         if (global_scales(c) == 0.0)
         {
            continue;
         }

         const real_t solution_scale =
            std::max(solution_scales(c, e),
                     scale_floors(c));

         MFEM_VERIFY(solution_scale > 0.0,
                     "Invalid KXRCF solution scale.");

         const real_t denominator =
            h_factor *
            inflow_measure(e) *
            solution_scale;

         MFEM_VERIFY(denominator > 0.0,
                     "Invalid KXRCF denominator.");

         /*
          * Important:
          *
          * take the absolute value AFTER integrating the
          * signed jump:
          *
          *       | integral (u_K - u_N) ds |
          *
          * not
          *
          *       integral |u_K - u_N| ds.
          */
         const real_t component_indicator =
            std::abs(jump_integral(c, e))
            / denominator;

         if (component_indicator_values)
         {
            (*component_indicator_values)(c, e) =
               component_indicator;
         }

         pooled_indicator =
            std::max(pooled_indicator,
                     component_indicator);
      }

      if (indicator_values)
      {
         (*indicator_values)(e) =
            pooled_indicator;
      }

      /*
       * Any troubled component activates stabilization for the
       * whole state vector in this element.
       */
      active[e] =
         pooled_indicator > threshold;
   }

#ifdef KXRCF_INTERNAL_TIMING
   internal_timing.normalization +=
      detail::KXRCFSecondsSince(normalization_timer_begin);
   internal_timing.total += detail::KXRCFSecondsSince(total_timer_begin);
#endif
}

} // namespace ofdg
