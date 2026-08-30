#include <algorithm>
#include <cmath>

#include "mfem.hpp"

using namespace mfem;


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
   VectorCoefficient *velocity;

   int dim;
   int ncomp;

   real_t threshold;
   real_t relative_scale_floor;


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


public:
   KXRCFIndicator(FiniteElementSpace *fes_,
                  VectorCoefficient *velocity_,
                  real_t threshold_ = 1.0,
                  real_t relative_scale_floor_ = 1e-14)
      : fes(fes_),
        velocity(velocity_),
        dim(0),
        ncomp(0),
        threshold(threshold_),
        relative_scale_floor(relative_scale_floor_)
   {
      MFEM_VERIFY(fes != NULL,
                  "KXRCF requires a valid finite element space.");

      MFEM_VERIFY(velocity != NULL,
                  "KXRCF requires a valid velocity coefficient.");

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

      MFEM_VERIFY(velocity->GetVDim() == dim,
                  "KXRCF velocity dimension must match mesh dimension.");

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
   Mesh *mesh = fes->GetMesh();
   const int ne = fes->GetNE();

   scales.SetSize(ncomp, ne);
   scales = 0.0;

   Array<int> vdofs;

   for (int e = 0; e < ne; ++e)
   {
      const FiniteElement *fe = fes->GetFE(e);
      ElementTransformation *Tr =
         mesh->GetElementTransformation(e);

      const int ndof = fe->GetDof();

      fes->GetElementVDofs(e, vdofs);

      MFEM_VERIFY(vdofs.Size() == ndof * ncomp,
                  "Unexpected number of element VDofs in KXRCF.");

      Vector local_state(vdofs.Size());
      x.GetSubVector(vdofs, local_state);

      Vector shape(ndof);
      Vector state(ncomp);

      /*
       * More than sufficient for evaluating / integrating a
       * degree-k DG polynomial.
       */
      const IntegrationRule &ir =
         IntRules.Get(fe->GetGeomType(),
                      2 * fe->GetOrder() + 2);

      if (dim == 1)
      {
         /*
          * Original 1D KXRCF normalization:
          *
          *     | element average |
          */
         Vector integral(ncomp);
         integral = 0.0;

         real_t volume = 0.0;

         for (int q = 0; q < ir.GetNPoints(); ++q)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);

            Tr->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            EvaluateElementState(local_state,
                                 shape,
                                 state);

            const real_t weight =
               ip.weight * Tr->Weight();

            for (int c = 0; c < ncomp; ++c)
            {
               integral(c) +=
                  weight * state(c);
            }

            volume += weight;
         }

         MFEM_VERIFY(volume > 0.0,
                     "KXRCF encountered an element with "
                     "non-positive volume.");

         for (int c = 0; c < ncomp; ++c)
         {
            scales(c, e) =
               std::abs(integral(c) / volume);
         }
      }
      else
      {
         /*
          * Multidimensional KXRCF normalization:
          *
          *     max_q |u(x_q)|
          */
         for (int q = 0; q < ir.GetNPoints(); ++q)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);

            fe->CalcShape(ip, shape);

            EvaluateElementState(local_state,
                                 shape,
                                 state);

            for (int c = 0; c < ncomp; ++c)
            {
               scales(c, e) =
                  std::max(scales(c, e),
                           std::abs(state(c)));
            }
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


// --------------------------------------------------------------------------
// Main KXRCF computation
// --------------------------------------------------------------------------

void KXRCFIndicator::Compute(
   const Vector &x,
   Array<bool> &active,
   Vector *indicator_values,
   DenseMatrix *component_indicator_values) const
{
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
   ComputeElementScales(x, solution_scales);


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

   Vector velocity_value(dim);
   Vector normal(dim);

   for (int f = 0;
        f < mesh->GetNumFaces();
        ++f)
   {
      FaceElementTransformations *FTr =
         mesh->GetFaceElementTransformations(f);

      if (!FTr)
      {
         continue;
      }

      const int e1 = FTr->Elem1No;
      const int e2 = FTr->Elem2No;

      /*
       * First multidimensional version:
       *
       * only interfaces with a valid state on both sides.
       *
       * Periodically connected faces are included provided MFEM
       * exposes both neighboring elements.
       */
      if (e1 < 0 || e2 < 0)
      {
         continue;
      }

      const FiniteElement *fe1 =
         fes->GetFE(e1);

      const FiniteElement *fe2 =
         fes->GetFE(e2);

      const int ndof1 = fe1->GetDof();
      const int ndof2 = fe2->GetDof();

      fes->GetElementVDofs(e1, vdofs1);
      fes->GetElementVDofs(e2, vdofs2);

      MFEM_VERIFY(vdofs1.Size() == ndof1 * ncomp,
                  "Unexpected element-1 VDof count.");

      MFEM_VERIFY(vdofs2.Size() == ndof2 * ncomp,
                  "Unexpected element-2 VDof count.");

      Vector u_e1(vdofs1.Size());
      Vector u_e2(vdofs2.Size());

      x.GetSubVector(vdofs1, u_e1);
      x.GetSubVector(vdofs2, u_e2);

      Vector shape1(ndof1);
      Vector shape2(ndof2);

      Vector state1(ncomp);
      Vector state2(ncomp);


      // --------------------------------------------------------------------
      // Face quadrature
      // --------------------------------------------------------------------

      int face_integration_order = 0;

      if (dim > 1)
      {
         const int k_face =
            std::max(fe1->GetOrder(),
                     fe2->GetOrder());

         face_integration_order =
            2 * k_face + 2;
      }

      /*
       * In 1D the face is a point, so order zero gives the
       * single required integration point and reproduces the
       * original implementation.
       */
      const IntegrationRule &ir =
         IntRules.Get(FTr->FaceGeom,
                      face_integration_order);


      for (int q = 0;
           q < ir.GetNPoints();
           ++q)
      {
         const IntegrationPoint &ip =
            ir.IntPoint(q);

         /*
          * Updates the face point and its corresponding points
          * in both neighboring reference elements.
          */
         FTr->SetAllIntPoints(&ip);

         const IntegrationPoint &eip1 =
            FTr->GetElement1IntPoint();

         const IntegrationPoint &eip2 =
            FTr->GetElement2IntPoint();


         // -----------------------------------------------------------------
         // Evaluate traces from both elements
         // -----------------------------------------------------------------

         fe1->CalcShape(eip1, shape1);
         fe2->CalcShape(eip2, shape2);

         EvaluateElementState(u_e1,
                              shape1,
                              state1);

         EvaluateElementState(u_e2,
                              shape2,
                              state2);


         // -----------------------------------------------------------------
         // Determine outward normal of element 1 and physical ds
         // -----------------------------------------------------------------

         real_t vn1 = 0.0;
         real_t physical_face_weight = 0.0;

         if (dim == 1)
         {
            /*
             * Preserve exactly the old 1D normal convention.
             *
             * Reference segment:
             *
             *       x=0               x=1
             *        |-----------------|
             *       n=-1              n=+1
             */
            FTr->Elem1->SetIntPoint(&eip1);

            const DenseMatrix &J =
               FTr->Elem1->Jacobian();

            MFEM_VERIFY(
               J.Height() == 1 &&
               J.Width() == 1,
               "Unexpected 1D Jacobian dimensions.");

            const real_t orientation =
               (J(0, 0) >= 0.0)
               ? 1.0
               : -1.0;

            const real_t reference_normal =
               (eip1.x < 0.5)
               ? -1.0
               : 1.0;

            const real_t normal1 =
               reference_normal *
               orientation;

            velocity->Eval(velocity_value,
                           *FTr->Elem1,
                           eip1);

            vn1 =
               velocity_value(0) *
               normal1;

            /*
             * A zero-dimensional face has no ordinary length
             * measure. The 1D KXRCF formula counts inflow
             * endpoints, exactly as the previous implementation.
             */
            physical_face_weight = 1.0;
         }
         else
         {
            /*
             * MFEM's CalcOrtho gives the scaled face normal:
             *
             *     n_scaled = n_unit * J_face.
             *
             * Hence
             *
             *     |n_scaled| = physical face Jacobian.
             */
            CalcOrtho(FTr->Jacobian(),
                      normal);

            const real_t normal_norm =
               normal.Norml2();

            MFEM_VERIFY(normal_norm > 0.0,
                        "KXRCF encountered a degenerate face.");

            physical_face_weight =
               ip.weight *
               normal_norm;

            /*
             * Velocity magnitude does NOT enter the KXRCF
             * numerator.
             *
             * It is used only to determine which element is
             * receiving flow at this point.
             *
             * Since normal is merely scaled by a positive
             * surface Jacobian, its scaling does not affect
             * the sign.
             */
            velocity->Eval(velocity_value,
                           *FTr->Elem1,
                           eip1);

            vn1 =
               velocity_value *
               normal;
         }


         // -----------------------------------------------------------------
         // Accumulate into the element for which this quadrature
         // point lies on the inflow boundary.
         // -----------------------------------------------------------------

         if (vn1 < 0.0)
         {
            /*
             * Flow enters element 1:
             *
             *     jump = u_1 - u_2
             */
            inflow_measure(e1) +=
               physical_face_weight;

            for (int c = 0; c < ncomp; ++c)
            {
               jump_integral(c, e1) +=
                  physical_face_weight *
                  (state1(c) - state2(c));
            }
         }
         else if (vn1 > 0.0)
         {
            /*
             * Since n_2 = -n_1, positive v.n_1 means the flow
             * enters element 2:
             *
             *     jump = u_2 - u_1
             */
            inflow_measure(e2) +=
               physical_face_weight;

            for (int c = 0; c < ncomp; ++c)
            {
               jump_integral(c, e2) +=
                  physical_face_weight *
                  (state2(c) - state1(c));
            }
         }

         /*
          * vn1 == 0:
          *
          * characteristic/tangential point.
          * It belongs to neither inflow boundary.
          */
      }
   }


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

      const real_t h =
         ComputeElementRadius(e);

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
}