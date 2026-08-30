#pragma once

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

namespace mfem
{

/** Face-local propagation information used by OFDG and KXRCF.
 *
 * beta1 and beta2 are non-negative normal spectral radii for the two traces.
 * normal_transport is signed in the outward-normal direction of element 1 and
 * is used only to decide which cell receives inflow in KXRCF.
 */
struct FacePhysicsSample
{
   real_t beta1 = 0.0;
   real_t beta2 = 0.0;
   real_t normal_transport = 0.0;
};

class FacePhysics
{
public:
   virtual ~FacePhysics() = default;

   virtual FacePhysicsSample Evaluate(
      const Vector &state1,
      const Vector &state2,
      const Vector &unit_normal,
      FaceElementTransformations &transformations) const = 0;
};

/** Compatibility policy used by the historical OFDG constructor. */
class UnitFacePhysics final : public FacePhysics
{
public:
   FacePhysicsSample Evaluate(const Vector &, const Vector &state2,
                              const Vector &,
                              FaceElementTransformations &) const override
   {
      return {1.0, state2.Size() > 0 ? 1.0 : 0.0, 1.0};
   }
};

/** Face physics for a prescribed advection velocity coefficient. */
class AdvectionFacePhysics final : public FacePhysics
{
private:
   std::shared_ptr<VectorCoefficient> owned_velocity;
   VectorCoefficient *velocity = nullptr;

public:
   explicit AdvectionFacePhysics(VectorCoefficient *velocity_)
      : velocity(velocity_)
   {
      MFEM_VERIFY(velocity != nullptr,
                  "Advection face physics requires a velocity coefficient.");
   }

   explicit AdvectionFacePhysics(
      std::shared_ptr<VectorCoefficient> velocity_)
      : owned_velocity(std::move(velocity_)), velocity(owned_velocity.get())
   {
      MFEM_VERIFY(velocity != nullptr,
                  "Advection face physics requires a velocity coefficient.");
   }

   FacePhysicsSample Evaluate(const Vector &, const Vector &state2,
                              const Vector &unit_normal,
                              FaceElementTransformations &transformations) const override
   {
      Vector velocity1(unit_normal.Size());
      Vector velocity2(unit_normal.Size());

      const IntegrationPoint &ip1 = transformations.GetElement1IntPoint();
      velocity->Eval(velocity1, *transformations.Elem1, ip1);

      const real_t normal1 = velocity1 * unit_normal;
      real_t normal2 = normal1;

      if (transformations.Elem2 && state2.Size() > 0)
      {
         const IntegrationPoint &ip2 = transformations.GetElement2IntPoint();
         velocity->Eval(velocity2, *transformations.Elem2, ip2);
         normal2 = velocity2 * unit_normal;
      }

      return {std::abs(normal1),
              state2.Size() > 0 ? std::abs(normal2) : 0.0,
              state2.Size() > 0 ? 0.5 * (normal1 + normal2) : normal1};
   }
};

/** Face physics for MFEM's multidimensional Burgers flux F(u)=u^2/2 (1,...,1). */
class BurgersFacePhysics final : public FacePhysics
{
public:
   FacePhysicsSample Evaluate(const Vector &state1, const Vector &state2,
                              const Vector &unit_normal,
                              FaceElementTransformations &) const override
   {
      MFEM_VERIFY(state1.Size() == 1,
                  "Burgers face physics expects one state component.");

      const real_t direction = unit_normal.Sum();
      const real_t speed1 = state1(0) * direction;

      if (state2.Size() == 0)
      {
         return {std::abs(speed1), 0.0, speed1};
      }

      MFEM_VERIFY(state2.Size() == 1,
                  "Burgers face physics expects one state component.");
      const real_t speed2 = state2(0) * direction;
      return {std::abs(speed1), std::abs(speed2),
              0.5 * (speed1 + speed2)};
   }
};

struct EulerPrimitiveState
{
   real_t density = 0.0;
   Vector velocity;
   real_t pressure = 0.0;
   real_t sound_speed = 0.0;
};

inline bool DecodeEulerState(const Vector &conservative, int dim,
                             real_t gamma, EulerPrimitiveState &primitive,
                             std::string *reason = nullptr)
{
   if (conservative.Size() != dim + 2)
   {
      if (reason) { *reason = "unexpected number of conservative variables"; }
      return false;
   }

   primitive.density = conservative(0);
   if (!(primitive.density > 0.0) || !std::isfinite(primitive.density))
   {
      if (reason) { *reason = "density is not positive and finite"; }
      return false;
   }

   primitive.velocity.SetSize(dim);
   real_t momentum_squared = 0.0;
   for (int d = 0; d < dim; ++d)
   {
      const real_t momentum = conservative(1 + d);
      primitive.velocity(d) = momentum / primitive.density;
      momentum_squared += momentum * momentum;
   }

   const real_t energy = conservative(dim + 1);
   const real_t kinetic_energy =
      0.5 * momentum_squared / primitive.density;
   primitive.pressure = (gamma - 1.0) * (energy - kinetic_energy);

   if (!(primitive.pressure > 0.0) || !std::isfinite(primitive.pressure))
   {
      if (reason) { *reason = "pressure is not positive and finite"; }
      return false;
   }

   primitive.sound_speed =
      std::sqrt(gamma * primitive.pressure / primitive.density);
   return std::isfinite(primitive.sound_speed);
}

/** Reflect only the normal Euler momentum, as required by a slip wall. */
inline void ReflectEulerState(const Vector &interior,
                              const Vector &unit_normal, int dim,
                              Vector &exterior)
{
   MFEM_VERIFY(interior.Size() == dim + 2 && unit_normal.Size() == dim,
               "Euler reflection received inconsistent dimensions.");
   exterior = interior;
   real_t normal_momentum = 0.0;
   for (int d = 0; d < dim; ++d)
   {
      normal_momentum += interior(1 + d) * unit_normal(d);
   }
   for (int d = 0; d < dim; ++d)
   {
      exterior(1 + d) -= 2.0 * normal_momentum * unit_normal(d);
   }
}

/** State-dependent normal characteristic speeds for compressible Euler. */
class EulerFacePhysics final : public FacePhysics
{
private:
   int dim;
   real_t gamma;

public:
   EulerFacePhysics(int dim_, real_t gamma_) : dim(dim_), gamma(gamma_)
   {
      MFEM_VERIFY(dim >= 1 && dim <= 3,
                  "Euler face physics supports dimensions 1, 2, and 3.");
      MFEM_VERIFY(gamma > 1.0,
                  "Euler face physics requires gamma > 1.");
   }

   FacePhysicsSample Evaluate(const Vector &state1, const Vector &state2,
                              const Vector &unit_normal,
                              FaceElementTransformations &transformations) const override
   {
      EulerPrimitiveState primitive1;
      EulerPrimitiveState primitive2;
      std::string reason;

      MFEM_VERIFY(DecodeEulerState(state1, dim, gamma, primitive1, &reason),
                  "Nonphysical Euler state on element " << transformations.Elem1No
                  << ": " << reason);

      const real_t normal1 = primitive1.velocity * unit_normal;
      if (state2.Size() == 0)
      {
         return {std::abs(normal1) + primitive1.sound_speed,
                 0.0, normal1};
      }

      MFEM_VERIFY(DecodeEulerState(state2, dim, gamma, primitive2, &reason),
                  "Nonphysical Euler state on element " << transformations.Elem2No
                  << ": " << reason);
      const real_t normal2 = primitive2.velocity * unit_normal;

      return {std::abs(normal1) + primitive1.sound_speed,
              std::abs(normal2) + primitive2.sound_speed,
              0.5 * (normal1 + normal2)};
   }
};

} // namespace mfem
