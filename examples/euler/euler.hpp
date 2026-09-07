#pragma once

// Derived from MFEM Example 18 (BSD 3-Clause). See
// ../../THIRD_PARTY_NOTICES.md for copyright and licence terms.
//                  MFEM Example 18 - Serial/Parallel Shared Code
//                      (Implementation of Time-dependent DG Operator)
//
// This code provide example problems for the Euler equations and implements
// the time-dependent DG operator given by the equation:
//
//            (u_t, v)_T - (F(u), ∇ v)_T + <F̂(u, n), [[v]]>_F = 0.
//
// This operator is designed for explicit time stepping methods. Specifically,
// the function DGHyperbolicConservationLaws::Mult implements the following
// transformation:
//
//                             u ↦ M⁻¹(-DF(u) + NF(u))
//
// where M is the mass matrix, DF is the weak divergence of flux, and NF is the
// interface flux. The inverse of the mass matrix is computed element-wise by
// leveraging the block-diagonal structure of the DG mass matrix. Additionally,
// the flux-related terms are computed using the HyperbolicFormIntegrator.
//
// The maximum characteristic speed is determined for each time step. For more
// details, refer to the documentation of DGHyperbolicConservationLaws::Mult.
//

#include <functional>
#include "mfem.hpp"
#include "../../src/face_physics.hpp"

namespace mfem
{

/** Boundary flux for fixed Euler states or slip/reflecting walls. */
class EulerBoundaryIntegrator : public NonlinearFormIntegrator
{
private:
   const NumericalFlux &numerical_flux;
   VectorCoefficient *fixed_state;
   bool reflecting;
   bool outflow;
   int integration_order_offset;
   int dimension;
   int equations;
   real_t max_char_speed = 0.0;

public:
   EulerBoundaryIntegrator(const NumericalFlux &flux,
                           VectorCoefficient *state,
                           bool reflecting_wall,
                           int order_offset, bool outflow_boundary = false)
      : numerical_flux(flux), fixed_state(state), reflecting(reflecting_wall), outflow(outflow_boundary),
        integration_order_offset(order_offset),
        dimension(flux.GetFluxFunction().dim),
        equations(flux.GetFluxFunction().num_equations)
   {
      MFEM_VERIFY(reflecting || outflow || fixed_state != nullptr,
                  "A fixed Euler boundary needs a state coefficient.");
   }

   void ResetMaxCharSpeed() { max_char_speed = 0.0; }
   real_t GetMaxCharSpeed() const { return max_char_speed; }

   void AssembleFaceVector(const FiniteElement &element,
                           const FiniteElement &,
                           FaceElementTransformations &transformations,
                           const Vector &element_state,
                           Vector &element_vector) override
   {
      MFEM_ASSERT(transformations.Elem2No < 0, "Expected a boundary face.");
      const int dofs = element.GetDof();
      Vector shape(dofs);
      Vector state_in(equations);
      Vector state_out(equations);
      Vector normal(transformations.GetSpaceDim());
      Vector numerical_flux_normal(equations);
      const DenseMatrix state_matrix(element_state.GetData(), dofs, equations);
      element_vector.SetSize(dofs * equations);
      element_vector = 0.0;
      DenseMatrix vector_matrix(element_vector.GetData(), dofs, equations);

      const IntegrationRule *rule = IntRule;
      if (!rule)
      {
         rule = &IntRules.Get(transformations.GetGeometryType(),
                             2 * element.GetOrder() +
                             integration_order_offset);
      }
      for (int q = 0; q < rule->GetNPoints(); ++q)
      {
         const IntegrationPoint &point = rule->IntPoint(q);
         transformations.SetAllIntPoints(&point);
         element.CalcShape(transformations.GetElement1IntPoint(), shape);
         state_matrix.MultTranspose(shape, state_in);
         if (normal.Size() == 1)
         {
            normal(0) = 2.0 * transformations.GetElement1IntPoint().x - 1.0;
         }
         else
         {
            CalcOrtho(transformations.Jacobian(), normal);
         }

         if (reflecting)
         {
            Vector unit_normal(normal);
            unit_normal /= unit_normal.Norml2();
            ofdg::ReflectEulerState(state_in, unit_normal, dimension, state_out);
         }
         else if (outflow)
         {
            state_out = state_in;
         }
         else
         {
            fixed_state->Eval(state_out, transformations, point);
         }

         max_char_speed = std::max(
            max_char_speed,
            numerical_flux.Eval(state_in, state_out, normal,
                                transformations, numerical_flux_normal));
         AddMult_a_VWt(-point.weight, shape, numerical_flux_normal,
                       vector_matrix);
      }
   }
};

/// @brief Time dependent DG operator for hyperbolic conservation laws
class DGHyperbolicConservationLaws : public TimeDependentOperator
{
private:
   const int num_equations; // the number of equations
   const int dim;
   FiniteElementSpace &vfes; // vector finite element space
   // Element integration form. Should contain ComputeFlux
   std::unique_ptr<HyperbolicFormIntegrator> formIntegrator;
   std::unique_ptr<EulerBoundaryIntegrator> boundaryIntegrator;
   // Base Nonlinear Form
   std::unique_ptr<NonlinearForm> nonlinearForm;
   // element-wise inverse mass matrix
   std::vector<DenseMatrix> invmass; // local scalar inverse mass
   std::vector<DenseMatrix> weakdiv; // local weak divergence (trial space ByDim)
   // global maximum characteristic speed. Updated by form integrators
   mutable real_t max_char_speed;
   // auxiliary variable used in Mult
   mutable Vector z;

   // Compute element-wise inverse mass matrix
   void ComputeInvMass();
   // Compute element-wise weak-divergence matrix
   void ComputeWeakDivergence();

public:
   /**
    * @brief Construct a new DGHyperbolicConservationLaws object
    *
    * @param vfes_ vector finite element space. Only tested for DG [Pₚ]ⁿ
    * @param formIntegrator_ integrator (F(u,x), grad v)
    * @param preassembleWeakDivergence preassemble weak divergence for faster
    *                                  assembly
    */
   DGHyperbolicConservationLaws(
      FiniteElementSpace &vfes_,
      std::unique_ptr<HyperbolicFormIntegrator> formIntegrator_,
      bool preassembleWeakDivergence=true,
      std::unique_ptr<EulerBoundaryIntegrator>
         boundaryIntegrator_=nullptr);
   /**
    * @brief Apply nonlinear form to obtain M⁻¹(DIVF + JUMP HAT(F))
    *
    * @param x current solution vector
    * @param y resulting dual vector to be used in an EXPLICIT solver
    */
   void Mult(const Vector &x, Vector &y) const override;
   // get global maximum characteristic speed to be used in CFL condition
   // where max_char_speed is updated during Mult.
   real_t GetMaxCharSpeed() { return max_char_speed; }
   void Update();

};

//////////////////////////////////////////////////////////////////
///        HYPERBOLIC CONSERVATION LAWS IMPLEMENTATION         ///
//////////////////////////////////////////////////////////////////

// Implementation of class DGHyperbolicConservationLaws
DGHyperbolicConservationLaws::DGHyperbolicConservationLaws(
   FiniteElementSpace &vfes_,
   std::unique_ptr<HyperbolicFormIntegrator> formIntegrator_,
   bool preassembleWeakDivergence,
   std::unique_ptr<EulerBoundaryIntegrator> boundaryIntegrator_)
   : TimeDependentOperator(vfes_.GetTrueVSize()),
     num_equations(formIntegrator_->num_equations),
     dim(vfes_.GetMesh()->SpaceDimension()),
     vfes(vfes_),
     formIntegrator(std::move(formIntegrator_)),
     boundaryIntegrator(std::move(boundaryIntegrator_)),
     z(vfes_.GetTrueVSize())
{
   // Standard local assembly and inversion for energy mass matrices.
   ComputeInvMass();
#ifndef MFEM_USE_MPI
   nonlinearForm.reset(new NonlinearForm(&vfes));
#else
   ParFiniteElementSpace *pvfes = dynamic_cast<ParFiniteElementSpace *>(&vfes);
   if (pvfes)
   {
      nonlinearForm.reset(new ParNonlinearForm(pvfes));
   }
   else
   {
      nonlinearForm.reset(new NonlinearForm(&vfes));
   }
#endif
   if (preassembleWeakDivergence)
   {
      ComputeWeakDivergence();
   }
   else
   {
      nonlinearForm->AddDomainIntegrator(formIntegrator.get());
   }
   nonlinearForm->AddInteriorFaceIntegrator(formIntegrator.get());
   if (boundaryIntegrator)
   {
      nonlinearForm->AddBdrFaceIntegrator(boundaryIntegrator.get());
   }
   nonlinearForm->UseExternalIntegrators();

}

void DGHyperbolicConservationLaws::ComputeInvMass()
{
   InverseIntegrator inv_mass(new MassIntegrator());

   invmass.resize(vfes.GetNE());
   for (int i=0; i<vfes.GetNE(); i++)
   {
      int dof = vfes.GetFE(i)->GetDof();
      invmass[i].SetSize(dof);
      inv_mass.AssembleElementMatrix(*vfes.GetFE(i),
                                     *vfes.GetElementTransformation(i),
                                     invmass[i]);
   }
}

void DGHyperbolicConservationLaws::ComputeWeakDivergence()
{
   TransposeIntegrator weak_div(new GradientIntegrator());
   DenseMatrix weakdiv_bynodes;

   weakdiv.resize(vfes.GetNE());
   for (int i=0; i<vfes.GetNE(); i++)
   {
      int dof = vfes.GetFE(i)->GetDof();
      weakdiv_bynodes.SetSize(dof, dof*dim);
      weak_div.AssembleElementMatrix2(*vfes.GetFE(i), *vfes.GetFE(i),
                                      *vfes.GetElementTransformation(i),
                                      weakdiv_bynodes);
      weakdiv[i].SetSize(dof, dof*dim);
      // Reorder so that trial space is ByDim.
      // This makes applying weak divergence to flux value simpler.
      for (int j=0; j<dof; j++)
      {
         for (int d=0; d<dim; d++)
         {
            weakdiv[i].SetCol(j*dim + d, weakdiv_bynodes.GetColumn(d*dof + j));
         }
      }

   }
}


void DGHyperbolicConservationLaws::Mult(const Vector &x, Vector &y) const
{
   // 0. Reset wavespeed computation before operator application.
   formIntegrator->ResetMaxCharSpeed();
   if (boundaryIntegrator) { boundaryIntegrator->ResetMaxCharSpeed(); }
   // 1. Apply Nonlinear form to obtain an auxiliary result
   //         z = - <F̂(u_h,n), [[v]]>_e
   //    If weak-divergence is not preassembled, we also have weak-divergence
   //         z = - <F̂(u_h,n), [[v]]>_e + (F(u_h), ∇v)
   nonlinearForm->Mult(x, z);
   if (!weakdiv.empty()) // if weak divergence is pre-assembled
   {
      // Apply weak divergence to F(u_h), and inverse mass to z_loc + weakdiv_loc
      Vector current_state; // view of current state at a node
      DenseMatrix current_flux(num_equations, dim); // flux of current state
      DenseMatrix flux; // element flux value. Whose column is ordered by dim.
      DenseMatrix current_xmat; // view of current states in an element, dof x num_eq
      DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
      DenseMatrix current_ymat; // view of element result, dof x num_eq
      const FluxFunction &fluxFunction = formIntegrator->GetFluxFunction();
      Array<int> vdofs;
      Vector xval, zval;
      for (int i=0; i<vfes.GetNE(); i++)
      {
         ElementTransformation* Tr = vfes.GetElementTransformation(i);
         int dof = vfes.GetFE(i)->GetDof();
         vfes.GetElementVDofs(i, vdofs);
         x.GetSubVector(vdofs, xval);
         current_xmat.UseExternalData(xval.GetData(), dof, num_equations);
         flux.SetSize(num_equations, dim*dof);
         for (int j=0; j<dof; j++) // compute flux for all nodes in the element
         {
            current_xmat.GetRow(j, current_state);
            fluxFunction.ComputeFlux(current_state, *Tr, current_flux);
            for (int d = 0; d < dim; ++d)
            {
               for (int equation = 0; equation < num_equations; ++equation)
               {
                  flux(equation, j * dim + d) = current_flux(equation, d);
               }
            }
         }
         // Compute weak-divergence and add it to auxiliary result, z
         // Recalling that weakdiv is reordered by dim, we can apply
         // weak-divergence to the transpose of flux.
         z.GetSubVector(vdofs, zval);
         current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
         mfem::AddMult_a_ABt(1.0, weakdiv[i], flux, current_zmat);
         // Apply inverse mass to auxiliary result to obtain the final result
         current_ymat.SetSize(dof, num_equations);
         mfem::Mult(invmass[i], current_zmat, current_ymat);
         y.SetSubVector(vdofs, current_ymat.GetData());
      }
   }
   else
   {
      // Apply block inverse mass
      Vector zval; // z_loc, dof*num_eq

      DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
      DenseMatrix current_ymat; // view of element result, dof x num_eq
      Array<int> vdofs;
      for (int i=0; i<vfes.GetNE(); i++)
      {
         int dof = vfes.GetFE(i)->GetDof();
         vfes.GetElementVDofs(i, vdofs);
         z.GetSubVector(vdofs, zval);
         current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
         current_ymat.SetSize(dof, num_equations);
         mfem::Mult(invmass[i], current_zmat, current_ymat);
         y.SetSubVector(vdofs, current_ymat.GetData());
      }
   }
   max_char_speed = formIntegrator->GetMaxCharSpeed();
   if (boundaryIntegrator)
   {
      max_char_speed = std::max(max_char_speed,
                                boundaryIntegrator->GetMaxCharSpeed());
   }
}

void DGHyperbolicConservationLaws::Update()
{
   nonlinearForm->Update();
   height = nonlinearForm->Height();
   width = height;
   z.SetSize(height);

   ComputeInvMass();
   if (!weakdiv.empty()) {ComputeWeakDivergence();}
}

std::function<void(const Vector&, Vector&)> GetMovingVortexInit(
   const real_t radius, const real_t Minf, const real_t beta,
   const real_t gas_constant, const real_t specific_heat_ratio)
{
   return [specific_heat_ratio,
           gas_constant, Minf, radius, beta](const Vector &x, Vector &y)
   {
      MFEM_ASSERT(x.Size() == 2, "");

      const real_t xc = 0.5, yc = 0.5;

      // Nice units
      const real_t vel_inf = 1.;
      const real_t den_inf = 1.;

      // Derive remainder of background state from this and Minf
      const real_t pres_inf = (den_inf / specific_heat_ratio) *
                              (vel_inf / Minf) * (vel_inf / Minf);
      const real_t temp_inf = pres_inf / (den_inf * gas_constant);

      real_t r2rad = 0.0;
      r2rad += (x(0) - xc) * (x(0) - xc);
      r2rad += (x(1) - yc) * (x(1) - yc);
      r2rad /= (radius * radius);

      const real_t shrinv1 = 1.0 / (specific_heat_ratio - 1.);

      const real_t velX =
         vel_inf * (1 - beta * (x(1) - yc) / radius * std::exp(-0.5 * r2rad));
      const real_t velY =
         vel_inf * beta * (x(0) - xc) / radius * std::exp(-0.5 * r2rad);
      const real_t vel2 = velX * velX + velY * velY;

      const real_t specific_heat =
         gas_constant * specific_heat_ratio * shrinv1;
      const real_t temp = temp_inf - 0.5 * (vel_inf * beta) *
                          (vel_inf * beta) / specific_heat *
                          std::exp(-r2rad);

      const real_t den = den_inf * std::pow(temp / temp_inf, shrinv1);
      const real_t pres = den * gas_constant * temp;
      const real_t energy = shrinv1 * pres / den + 0.5 * vel2;

      y(0) = den;
      y(1) = den * velX;
      y(2) = den * velY;
      y(3) = den * energy;
   };
}

Mesh EulerMesh(const int problem, const int requested_elements = 0)
{
   const auto count = [requested_elements](int fallback)
   {
      return requested_elements > 0 ? requested_elements : fallback;
   };
   switch (problem)
   {
      case 1:
      case 2:
      {
         Mesh mesh = Mesh::MakeCartesian2D(count(8), count(8), Element::QUADRILATERAL,
                                           true, 1.0, 1.0);
         Vector x_translation({1.0, 0.0});
         Vector y_translation({0.0, 1.0});
         std::vector<Vector> translations = {x_translation, y_translation};
         return Mesh::MakePeriodic(
                   mesh, mesh.CreatePeriodicVertexMapping(translations));
      }
      case 3:
      {
         Mesh mesh = Mesh::MakeCartesian2D(count(8), count(8), Element::QUADRILATERAL,
                                           true, 2.0, 2.0);
         Vector x_translation({2.0, 0.0});
         Vector y_translation({0.0, 2.0});
         std::vector<Vector> translations = {x_translation, y_translation};
         return Mesh::MakePeriodic(
                   mesh, mesh.CreatePeriodicVertexMapping(translations));
      }
      case 4:
      {
         Mesh mesh = Mesh::MakeCartesian1D(count(32), 2.0 * M_PI);
         std::vector<int> vertex_map(mesh.GetNV());
         for (int v = 0; v < mesh.GetNV(); ++v) { vertex_map[v] = v; }
         vertex_map.back() = 0;
         return Mesh::MakePeriodic(mesh, vertex_map);
      }
      case 5:
         return Mesh::MakeCartesian1D(count(200), 1.0);
      case 6:
         return Mesh::MakeCartesian2D(count(80), count(80), Element::QUADRILATERAL,
                                      true, 1.0, 1.0);
      case 7:
      case 8:
      {
         Mesh mesh = Mesh::MakeCartesian1D(
            count(problem == 7 ? 256 : 400), 10.0);
         for (int vertex = 0; vertex < mesh.GetNV(); ++vertex)
         {
            mesh.GetVertex(vertex)[0] -= 5.0;
         }
         return mesh;
      }
      default:
         MFEM_ABORT("Problem Undefined");
   }
}

// Initial condition
VectorFunctionCoefficient EulerInitialCondition(const int problem,
                                                const real_t specific_heat_ratio,
                                                const real_t gas_constant)
{
   switch (problem)
   {
      case 1: // fast moving vortex
         return VectorFunctionCoefficient(
                   4, GetMovingVortexInit(0.2, 0.5, 1. / 5., gas_constant,
                                          specific_heat_ratio));
      case 2: // slow moving vortex
         return VectorFunctionCoefficient(
                   4, GetMovingVortexInit(0.2, 0.05, 1. / 50., gas_constant,
                                          specific_heat_ratio));
      case 3: // moving sine wave
         return VectorFunctionCoefficient(4, [](const Vector &x, Vector &y)
         {
            MFEM_ASSERT(x.Size() == 2, "");
            const real_t density = 1.0 + 0.2 * std::sin(M_PI*(x(0) + x(1)));
            const real_t velocity_x = 0.7;
            const real_t velocity_y = 0.3;
            const real_t pressure = 1.0;
            const real_t energy =
               pressure / (1.4 - 1.0) +
               density * 0.5 * (velocity_x * velocity_x + velocity_y * velocity_y);

            y(0) = density;
            y(1) = density * velocity_x;
            y(2) = density * velocity_y;
            y(3) = energy;
         });
      case 4:
         return VectorFunctionCoefficient(3, [](const Vector &x, Vector &y)
         {
            MFEM_ASSERT(x.Size() == 1, "");
            const real_t wave = std::sin(x(0));
            const real_t density = 2.0 + 2.0 * wave * wave;
            const real_t velocity_x = 1.0;
            const real_t pressure = 2.0;
            const real_t energy =
               pressure / (1.4 - 1.0) + density * 0.5 * (velocity_x * velocity_x);

            y(0) = density;
            y(1) = density * velocity_x;
            y(2) = energy;
         });
      case 5: // Woodward--Colella interacting blast waves
         return VectorFunctionCoefficient(
                   3, [specific_heat_ratio](const Vector &x, Vector &y)
         {
            MFEM_ASSERT(x.Size() == 1, "");
            const real_t pressure = x(0) < 0.1 ? 1000.0
                                    : (x(0) < 0.9 ? 0.01 : 100.0);
            y(0) = 1.0;
            y(1) = 0.0;
            y(2) = pressure / (specific_heat_ratio - 1.0);
         });
      case 6: // first two-dimensional Riemann configuration
         return VectorFunctionCoefficient(
                   4, [specific_heat_ratio](const Vector &x, Vector &y)
         {
            MFEM_ASSERT(x.Size() == 2, "");

            real_t density, velocity_x, velocity_y, pressure;
            if (x(0) < 0.5 && x(1) < 0.5)
            {
               density = 0.8;
               velocity_x = velocity_y = 0.0;
               pressure = 1.0;
            }
            else if (x(0) < 0.5)
            {
               density = 1.0;
               velocity_x = 0.7276;
               velocity_y = 0.0;
               pressure = 1.0;
            }
            else if (x(1) < 0.5)
            {
               density = 1.0;
               velocity_x = 0.0;
               velocity_y = 0.7276;
               pressure = 1.0;
            }
            else
            {
               density = 0.5313;
               velocity_x = velocity_y = 0.0;
               pressure = 0.4;
            }

            y(0) = density;
            y(1) = density * velocity_x;
            y(2) = density * velocity_y;
            y(3) = pressure / (specific_heat_ratio - 1.0) +
                   0.5 * density *
                   (velocity_x * velocity_x + velocity_y * velocity_y);
         });
      case 7: // scaled Lax problem; scaling is applied by the driver
         return VectorFunctionCoefficient(
                   3, [specific_heat_ratio](const Vector &x, Vector &y)
         {
            const real_t density = x(0) < 0.0 ? 0.445 : 0.5;
            const real_t velocity = x(0) < 0.0 ? 0.698 : 0.0;
            const real_t pressure = x(0) < 0.0 ? 3.528 : 0.571;
            y(0) = density;
            y(1) = density * velocity;
            y(2) = pressure / (specific_heat_ratio - 1.0) +
                   0.5 * density * velocity * velocity;
         });
      case 8: // Shu--Osher shock--entropy interaction
         return VectorFunctionCoefficient(
                   3, [specific_heat_ratio](const Vector &x, Vector &y)
         {
            real_t density, velocity, pressure;
            if (x(0) < -4.0)
            {
               density = 3.857143;
               velocity = 2.629369;
               pressure = 10.33333;
            }
            else
            {
               density = 1.0 + 0.2 * std::sin(5.0 * x(0));
               velocity = 0.0;
               pressure = 1.0;
            }
            y(0) = density;
            y(1) = density * velocity;
            y(2) = pressure / (specific_heat_ratio - 1.0) +
                   0.5 * density * velocity * velocity;
         });
      default:
         MFEM_ABORT("Problem Undefined");
   }
}

VectorFunctionCoefficient EulerVortexExactCondition(
   const int problem, const real_t time,
   const real_t specific_heat_ratio, const real_t gas_constant)
{
   MFEM_VERIFY(problem >= 1 && problem <= 4,
               "Exact state is defined only for Euler problems 1--4.");
   if (problem == 3)
   {
      return VectorFunctionCoefficient(4, [time](const Vector &x, Vector &y)
      {
         const real_t density =
            1.0 + 0.2 * std::sin(M_PI * (x(0) + x(1) - time));
         const real_t velocity_x = 0.7;
         const real_t velocity_y = 0.3;
         const real_t pressure = 1.0;
         y(0) = density;
         y(1) = density * velocity_x;
         y(2) = density * velocity_y;
         y(3) = pressure / (1.4 - 1.0) + 0.5 * density *
                (velocity_x * velocity_x + velocity_y * velocity_y);
      });
   }
   if (problem == 4)
   {
      return VectorFunctionCoefficient(3, [time](const Vector &x, Vector &y)
      {
         const real_t wave = std::sin(x(0) - time);
         const real_t density = 2.0 + 2.0 * wave * wave;
         const real_t velocity = 1.0;
         const real_t pressure = 2.0;
         y(0) = density;
         y(1) = density * velocity;
         y(2) = pressure / (1.4 - 1.0) + 0.5 * density * velocity * velocity;
      });
   }
   const real_t mach = problem == 1 ? 0.5 : 0.05;
   const real_t beta = problem == 1 ? 1.0 / 5.0 : 1.0 / 50.0;
   auto initial = GetMovingVortexInit(0.2, mach, beta, gas_constant,
                                      specific_heat_ratio);
   return VectorFunctionCoefficient(4, [initial, time](const Vector &x,
                                                       Vector &state)
   {
      Vector foot(x);
      foot(0) -= time;
      foot(0) -= std::floor(foot(0));
      initial(foot, state);
   });
}

FunctionCoefficient EulerExactDensityCondition(
   const int problem, const real_t time,
   const real_t specific_heat_ratio, const real_t gas_constant)
{
   MFEM_VERIFY(problem >= 1 && problem <= 4,
               "Exact density is defined only for Euler problems 1--4.");
   if (problem == 1 || problem == 2)
   {
      const real_t mach = problem == 1 ? 0.5 : 0.05;
      const real_t beta = problem == 1 ? 1.0 / 5.0 : 1.0 / 50.0;
      auto initial = GetMovingVortexInit(0.2, mach, beta, gas_constant,
                                         specific_heat_ratio);
      return FunctionCoefficient([initial, time](const Vector &x)
      {
         Vector foot(x);
         foot(0) -= time;
         foot(0) -= std::floor(foot(0));
         Vector state(4);
         initial(foot, state);
         return state(0);
      });
   }
   if (problem == 3)
   {
      return FunctionCoefficient([time](const Vector &x)
      {
         return 1.0 + 0.2 * std::sin(M_PI * (x(0) + x(1) - time));
      });
   }
   return FunctionCoefficient([time](const Vector &x)
   {
      const real_t wave = std::sin(x(0) - time);
      return 2.0 + 2.0 * wave * wave;
   });
}

} // namespace mfem
