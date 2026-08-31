#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "face_physics.hpp"

namespace ofdg
{

/** KXRCF troubled-cell indicator for scalar fields and systems.
 *
 * For element K and component c,
 *
 *   I_K,c = | integral_(Gamma_K^-) (u_K,c-u_N,c) ds |
 *           / (h_K^((k+1)/2) |Gamma_K^-| ||u_K,c||).
 *
 * The component maximum activates the element. Physical boundaries are
 * ignored. Inputs are preconditioned to be static, affine, full-dimensional
 * VALUE-mapped scalar DG elements, repeated over the vector dimension. A
 * uniform order is required; 2D triangles and quadrilaterals may be mixed.
 */
class KXRCFIndicator
{
private:
   mfem::FiniteElementSpace *fes;
   std::shared_ptr<const FacePhysics> face_physics;

   int dim;
   int ncomp;

   mfem::real_t threshold;
   mfem::real_t relative_scale_floor;
   struct ElementCache
   {
      mfem::Array<int> vdofs;
      mfem::DenseMatrix evaluation;
      mfem::Vector physical_weights;
      mfem::real_t volume = 0.0;
      mfem::real_t radius = 0.0;
   };
   std::vector<ElementCache> element_cache;
   struct LocalFaceCache
   {
      int face = -1;
      int element1 = -1;
      int element2 = -1;
      mfem::Array<int> vdofs1;
      mfem::Array<int> vdofs2;
      mfem::DenseMatrix evaluation1;
      mfem::DenseMatrix evaluation2;
      mfem::DenseMatrix unit_normals;
      mfem::Vector physical_weights;
      const mfem::IntegrationRule *rule = nullptr;
   };
   std::vector<LocalFaceCache> local_face_cache;

   struct FaceScratch
   {
      mfem::Vector shape1;
      mfem::Vector shape2;
      mfem::Vector state1;
      mfem::Vector state2;
      mfem::Vector normal;
   };
   mutable FaceScratch face_scratch;

   void EvaluateElementState(const mfem::Vector &local_state,
                             const mfem::Vector &shape,
                             mfem::Vector &value) const;
   void ComputeElementScales(const mfem::Vector &state,
                             mfem::DenseMatrix &scales) const;
   mfem::real_t ComputeElementRadius(int e) const;
   void BuildElementCache();
   void BuildLocalFaceCache();

   void AccumulateCachedFace(const LocalFaceCache &cache,
                             mfem::FaceElementTransformations &transformations,
                             const mfem::Vector &local_state1,
                             const mfem::Vector &local_state2,
                             mfem::DenseMatrix &jump_integral,
                             mfem::Vector &inflow_measure) const;

   void AccumulateFace(mfem::FaceElementTransformations &transformations,
                       const mfem::FiniteElement &element1,
                       const mfem::FiniteElement &element2,
                       const mfem::Vector &local_state1,
                       const mfem::Vector &local_state2,
                       int element1_number,
                       int element2_number,
                       bool accumulate_element2,
                       mfem::DenseMatrix &jump_integral,
                       mfem::Vector &inflow_measure) const;

   bool ComputeScales(const mfem::Vector &state,
                      mfem::DenseMatrix &solution_scales,
                      mfem::Vector &global_scales,
                      mfem::Vector &scale_floors) const;

   void AccumulateFaces(const mfem::Vector &state,
                        mfem::DenseMatrix &jump_integrals,
                        mfem::Vector &inflow_measures) const;

   void FinalizeIndicators(const mfem::DenseMatrix &solution_scales,
                           const mfem::Vector &global_scales,
                           const mfem::Vector &scale_floors,
                           const mfem::DenseMatrix &jump_integrals,
                           const mfem::Vector &inflow_measures,
                           mfem::Array<bool> &active,
                           mfem::Vector *indicator_values,
                           mfem::DenseMatrix *component_values) const;

public:
   KXRCFIndicator(mfem::FiniteElementSpace *fes_,
                  mfem::VectorCoefficient *velocity_,
                  mfem::real_t threshold_ = 1.0,
                  mfem::real_t relative_scale_floor_ = 1e-14)
      : KXRCFIndicator(
           fes_,
           std::make_shared<AdvectionFacePhysics>(velocity_),
           threshold_,
           relative_scale_floor_)
   {
   }

   KXRCFIndicator(mfem::FiniteElementSpace *fes_,
                  std::shared_ptr<const FacePhysics> face_physics_,
                  mfem::real_t threshold_ = 1.0,
                  mfem::real_t relative_scale_floor_ = 1e-14)
      : fes(fes_),
        face_physics(std::move(face_physics_)),
        dim(0),
        ncomp(0),
        threshold(threshold_),
        relative_scale_floor(relative_scale_floor_)
   {
      dim = fes->GetMesh()->Dimension();
      ncomp = fes->GetVDim();
      BuildElementCache();
      BuildLocalFaceCache();
   }

   /** Compute the pooled troubled-cell mask and optional indicator values. */
   void Compute(
      const mfem::Vector &x,
      mfem::Array<bool> &active,
      mfem::Vector *indicator_values = NULL,
      mfem::DenseMatrix *component_indicator_values = NULL) const;
};

// Local state evaluation

void KXRCFIndicator::EvaluateElementState(
   const mfem::Vector &local_state,
   const mfem::Vector &shape,
   mfem::Vector &value) const
{
   const int ndof = shape.Size();

   value.SetSize(ncomp);
   value = 0.0;

   // Element-local vector DOFs are component blocks, regardless of the
   // finite-element space global ordering: local_index = c * ndof + j.
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

// Element solution scales

void KXRCFIndicator::ComputeElementScales(
   const mfem::Vector &x,
   mfem::DenseMatrix &scales) const
{
   const int ne = fes->GetNE();

   scales.SetSize(ncomp, ne);
   scales = 0.0;

   mfem::Vector local_state;
   mfem::DenseMatrix local_matrix;
   mfem::DenseMatrix values;

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
         for (int c = 0; c < ncomp; ++c)
         {
            mfem::real_t integral = 0.0;
            for (int q = 0; q < values.Height(); ++q)
            {
               integral += cache.physical_weights(q) * values(q, c);
            }
            scales(c, e) = std::abs(integral / cache.volume);
         }
      }
      else
      {
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
   mfem::Mesh *mesh = fes->GetMesh();
   element_cache.resize(fes->GetNE());

   for (int e = 0; e < fes->GetNE(); ++e)
   {
      ElementCache &cache = element_cache[e];
      const mfem::FiniteElement *fe = fes->GetFE(e);
      mfem::ElementTransformation *transformation =
         mesh->GetElementTransformation(e);
      fes->GetElementVDofs(e, cache.vdofs);

      const mfem::IntegrationRule &rule =
         mfem::IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 2);
      cache.evaluation.SetSize(rule.GetNPoints(), fe->GetDof());
      cache.physical_weights.SetSize(rule.GetNPoints());
      cache.volume = 0.0;
      mfem::Vector shape(fe->GetDof());

      for (int q = 0; q < rule.GetNPoints(); ++q)
      {
         const mfem::IntegrationPoint &ip = rule.IntPoint(q);
         fe->CalcShape(ip, shape);
         for (int i = 0; i < shape.Size(); ++i)
         {
            cache.evaluation(q, i) = shape(i);
         }
         transformation->SetIntPoint(&ip);
         cache.physical_weights(q) = ip.weight * transformation->Weight();
         cache.volume += cache.physical_weights(q);
      }
      cache.radius = ComputeElementRadius(e);
   }
}

void KXRCFIndicator::BuildLocalFaceCache()
{
   mfem::Mesh *mesh = fes->GetMesh();
   local_face_cache.clear();
   local_face_cache.reserve(mesh->GetNumFaces());

   for (int f = 0; f < mesh->GetNumFaces(); ++f)
   {
      mfem::FaceElementTransformations *transformations =
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

      const mfem::FiniteElement *element1 = fes->GetFE(cache.element1);
      const mfem::FiniteElement *element2 = fes->GetFE(cache.element2);
      const int integration_order = dim == 1 ? 0 :
         2 * std::max(element1->GetOrder(), element2->GetOrder()) + 2;
      cache.rule = &mfem::IntRules.Get(transformations->FaceGeom,
                                 integration_order);
      const int point_count = cache.rule->GetNPoints();
      cache.evaluation1.SetSize(point_count, element1->GetDof());
      cache.evaluation2.SetSize(point_count, element2->GetDof());
      cache.unit_normals.SetSize(point_count, dim);
      cache.physical_weights.SetSize(point_count);
      mfem::Vector shape1(element1->GetDof());
      mfem::Vector shape2(element2->GetDof());
      mfem::Vector normal(dim);

      for (int q = 0; q < point_count; ++q)
      {
         const mfem::IntegrationPoint &ip = cache.rule->IntPoint(q);
         transformations->SetAllIntPoints(&ip);
         const mfem::IntegrationPoint &ip1 =
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
            const mfem::real_t orientation =
               transformations->Elem1->Jacobian()(0, 0) >= 0.0 ? 1.0 : -1.0;
            normal(0) = (ip1.x < 0.5 ? -1.0 : 1.0) * orientation;
            cache.physical_weights(q) = 1.0;
         }
         else
         {
            mfem::CalcOrtho(transformations->Jacobian(), normal);
            const mfem::real_t normal_norm = normal.Norml2();
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
   mfem::FaceElementTransformations &transformations,
   const mfem::Vector &local_state1,
   const mfem::Vector &local_state2,
   mfem::DenseMatrix &jump_integral,
   mfem::Vector &inflow_measure) const
{
   mfem::Vector &state1 = face_scratch.state1;
   mfem::Vector &state2 = face_scratch.state2;
   mfem::Vector &normal = face_scratch.normal;
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
      const mfem::real_t transport =
         face_physics->Evaluate(state1, state2, normal,
                                transformations).normal_transport;
      const mfem::real_t weight = cache.physical_weights(q);

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

// Generalized element radius

mfem::real_t KXRCFIndicator::ComputeElementRadius(
   int e) const
{
   mfem::Mesh *mesh = fes->GetMesh();

   if (dim == 1)
   {
      const mfem::real_t h =
         0.5 * mesh->GetElementVolume(e);

      return h;
   }

   const mfem::Geometry::Type geom =
      mesh->GetElementGeometry(e);

   mfem::ElementTransformation *Tr =
      mesh->GetElementTransformation(e);

   const mfem::IntegrationPoint &center_ip =
      mfem::Geometries.GetCenter(geom);

   const mfem::IntegrationRule *vertex_rule =
      mfem::Geometries.GetVertices(geom);

   mfem::Vector center(dim);
   mfem::Vector point(dim);

   Tr->Transform(center_ip, center);

   mfem::real_t radius = 0.0;

   for (int i = 0;
        i < vertex_rule->GetNPoints();
        ++i)
   {
      const mfem::IntegrationPoint &vip =
         vertex_rule->IntPoint(i);

      Tr->Transform(vip, point);

      mfem::real_t distance2 = 0.0;

      for (int d = 0; d < dim; ++d)
      {
         const mfem::real_t diff =
            point(d) - center(d);

         distance2 += diff * diff;
      }

      radius =
         std::max(radius,
                  std::sqrt(distance2));
   }

   return radius;
}

void KXRCFIndicator::AccumulateFace(
   mfem::FaceElementTransformations &transformations,
   const mfem::FiniteElement &element1,
   const mfem::FiniteElement &element2,
   const mfem::Vector &local_state1,
   const mfem::Vector &local_state2,
   int element1_number,
   int element2_number,
   bool accumulate_element2,
   mfem::DenseMatrix &jump_integral,
   mfem::Vector &inflow_measure) const
{
   mfem::Vector &shape1 = face_scratch.shape1;
   mfem::Vector &shape2 = face_scratch.shape2;
   mfem::Vector &state1 = face_scratch.state1;
   mfem::Vector &state2 = face_scratch.state2;
   mfem::Vector &normal = face_scratch.normal;
   shape1.SetSize(element1.GetDof());
   shape2.SetSize(element2.GetDof());
   state1.SetSize(ncomp);
   state2.SetSize(ncomp);
   normal.SetSize(dim);

   const int integration_order = dim == 1 ? 0 :
      2 * std::max(element1.GetOrder(), element2.GetOrder()) + 2;
   const mfem::IntegrationRule &rule =
      mfem::IntRules.Get(transformations.FaceGeom, integration_order);

   for (int q = 0; q < rule.GetNPoints(); ++q)
   {
      const mfem::IntegrationPoint &ip = rule.IntPoint(q);
      transformations.SetAllIntPoints(&ip);
      const mfem::IntegrationPoint &ip1 = transformations.GetElement1IntPoint();
      const mfem::IntegrationPoint &ip2 = transformations.GetElement2IntPoint();

      element1.CalcShape(ip1, shape1);
      element2.CalcShape(ip2, shape2);
      EvaluateElementState(local_state1, shape1, state1);
      EvaluateElementState(local_state2, shape2, state2);

      mfem::real_t physical_weight;
      if (dim == 1)
      {
         transformations.Elem1->SetIntPoint(&ip1);
         const mfem::DenseMatrix &jacobian = transformations.Elem1->Jacobian();
         const mfem::real_t orientation = jacobian(0, 0) >= 0.0 ? 1.0 : -1.0;
         normal(0) = (ip1.x < 0.5 ? -1.0 : 1.0) * orientation;
         physical_weight = 1.0;
      }
      else
      {
         mfem::CalcOrtho(transformations.Jacobian(), normal);
         const mfem::real_t normal_norm = normal.Norml2();
         physical_weight = ip.weight * normal_norm;
         normal /= normal_norm;
      }

      const mfem::real_t normal_transport =
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

// Main KXRCF computation

bool KXRCFIndicator::ComputeScales(
   const mfem::Vector &state,
   mfem::DenseMatrix &solution_scales,
   mfem::Vector &global_scales,
   mfem::Vector &scale_floors) const
{
   ComputeElementScales(state, solution_scales);
   global_scales.SetSize(ncomp);
   global_scales = 0.0;

   for (int c = 0; c < ncomp; ++c)
   {
      for (int e = 0; e < fes->GetNE(); ++e)
      {
         global_scales(c) =
            std::max(global_scales(c), solution_scales(c, e));
      }
   }

#ifdef MFEM_USE_MPI
   if (const auto *parallel_space =
          dynamic_cast<const mfem::ParFiniteElementSpace *>(fes))
   {
      MPI_Allreduce(MPI_IN_PLACE, global_scales.GetData(), ncomp,
                    mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_MAX,
                    parallel_space->GetComm());
   }
#endif

   scale_floors.SetSize(ncomp);
   bool nonzero = false;
   for (int c = 0; c < ncomp; ++c)
   {
      nonzero |= global_scales(c) > 0.0;
      scale_floors(c) = relative_scale_floor * global_scales(c);
   }
   return nonzero;
}

void KXRCFIndicator::AccumulateFaces(
   const mfem::Vector &state,
   mfem::DenseMatrix &jump_integrals,
   mfem::Vector &inflow_measures) const
{
   jump_integrals.SetSize(ncomp, fes->GetNE());
   jump_integrals = 0.0;
   inflow_measures.SetSize(fes->GetNE());
   inflow_measures = 0.0;

   mfem::Vector local_state1;
   mfem::Vector local_state2;
   mfem::Mesh *mesh = fes->GetMesh();
   for (const LocalFaceCache &cache : local_face_cache)
   {
      mfem::FaceElementTransformations *transformations =
         mesh->GetFaceElementTransformations(cache.face);
      local_state1.SetSize(cache.vdofs1.Size());
      local_state2.SetSize(cache.vdofs2.Size());
      state.GetSubVector(cache.vdofs1, local_state1);
      state.GetSubVector(cache.vdofs2, local_state2);
      AccumulateCachedFace(cache, *transformations, local_state1, local_state2,
                           jump_integrals, inflow_measures);
   }

#ifdef MFEM_USE_MPI
   if (auto *parallel_space =
          dynamic_cast<mfem::ParFiniteElementSpace *>(fes))
   {
      mfem::ParMesh *parallel_mesh = parallel_space->GetParMesh();
      mfem::ParGridFunction parallel_state(parallel_space);
      parallel_state = state;
      parallel_state.ExchangeFaceNbrData();

      mfem::Array<int> local_vdofs;
      mfem::Array<int> neighbor_vdofs;
      for (int face = 0; face < parallel_mesh->GetNSharedFaces(); ++face)
      {
         mfem::FaceElementTransformations *transformations =
            parallel_mesh->GetSharedFaceTransformations(face, true);
         const int local_element = transformations->Elem1No;
         const int neighbor_element =
            transformations->Elem2No - parallel_mesh->GetNE();

         parallel_space->GetElementVDofs(local_element, local_vdofs);
         parallel_space->GetFaceNbrElementVDofs(neighbor_element,
                                                neighbor_vdofs);
         mfem::Vector local_state(local_vdofs.Size());
         mfem::Vector neighbor_state(neighbor_vdofs.Size());
         state.GetSubVector(local_vdofs, local_state);
         parallel_state.FaceNbrData().GetSubVector(neighbor_vdofs,
                                                   neighbor_state);

         AccumulateFace(*transformations,
                        *parallel_space->GetFE(local_element),
                        *parallel_space->GetFaceNbrFE(neighbor_element),
                        local_state, neighbor_state, local_element, -1, false,
                        jump_integrals, inflow_measures);
      }
   }
#endif
}

void KXRCFIndicator::FinalizeIndicators(
   const mfem::DenseMatrix &solution_scales,
   const mfem::Vector &global_scales,
   const mfem::Vector &scale_floors,
   const mfem::DenseMatrix &jump_integrals,
   const mfem::Vector &inflow_measures,
   mfem::Array<bool> &active,
   mfem::Vector *indicator_values,
   mfem::DenseMatrix *component_values) const
{
   for (int e = 0; e < fes->GetNE(); ++e)
   {
      if (inflow_measures(e) == 0.0) { continue; }

      const int order = fes->GetFE(e)->GetOrder();
      const mfem::real_t h_factor =
         std::pow(element_cache[e].radius,
                  0.5 * static_cast<mfem::real_t>(order + 1));
      mfem::real_t pooled = 0.0;

      for (int c = 0; c < ncomp; ++c)
      {
         if (global_scales(c) == 0.0) { continue; }

         const mfem::real_t solution_scale =
            std::max(solution_scales(c, e), scale_floors(c));
         const mfem::real_t denominator =
            h_factor * inflow_measures(e) * solution_scale;

         // The KXRCF numerator is the absolute value of the signed
         // face integral, not the integral of the absolute jump.
         const mfem::real_t component_indicator =
            std::abs(jump_integrals(c, e)) / denominator;
         if (component_values)
         {
            (*component_values)(c, e) = component_indicator;
         }
         pooled = std::max(pooled, component_indicator);
      }

      if (indicator_values) { (*indicator_values)(e) = pooled; }
      active[e] = pooled > threshold;
   }
}

void KXRCFIndicator::Compute(
   const mfem::Vector &state,
   mfem::Array<bool> &active,
   mfem::Vector *indicator_values,
   mfem::DenseMatrix *component_values) const
{
   const int element_count = fes->GetNE();
   active.SetSize(element_count);
   active = false;

   if (indicator_values)
   {
      indicator_values->SetSize(element_count);
      *indicator_values = 0.0;
   }
   if (component_values)
   {
      component_values->SetSize(ncomp, element_count);
      *component_values = 0.0;
   }

   mfem::DenseMatrix solution_scales;
   mfem::Vector global_scales;
   mfem::Vector scale_floors;
   if (!ComputeScales(state, solution_scales, global_scales, scale_floors))
   {
      return;
   }

   mfem::DenseMatrix jump_integrals;
   mfem::Vector inflow_measures;
   AccumulateFaces(state, jump_integrals, inflow_measures);
   FinalizeIndicators(solution_scales, global_scales, scale_floors,
                      jump_integrals, inflow_measures, active,
                      indicator_values, component_values);
}

} // namespace ofdg
