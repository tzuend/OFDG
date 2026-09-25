#pragma once

#include "face_physics.hpp"
#include "mfem.hpp"

#include <memory>

namespace ofdg
{

namespace detail
{

struct OFDGSensorOptions
{
   enum class FaceRule
   {
      HighOrder,
      Trapezoidal
   };

   FaceRule face_rule = FaceRule::HighOrder;
   bool pool_components = false;
   bool use_face_height = false;
   bool use_element_mean_speed = false;
   bool relative_constant_tolerance = false;
};

} // namespace detail

class StabilizationCore
{
public:
   StabilizationCore(const mfem::FiniteElementSpace *fes, int basis_type,
                     std::shared_ptr<const FacePhysics> face_physics);
   StabilizationCore(const mfem::FiniteElementSpace *fes, int basis_type,
                     std::shared_ptr<const FacePhysics> face_physics,
                     detail::OFDGSensorOptions options);

   ~StabilizationCore();
   StabilizationCore(StabilizationCore &&) noexcept;
   StabilizationCore &operator=(StabilizationCore &&) noexcept;
   StabilizationCore(const StabilizationCore &) = delete;
   StabilizationCore &operator=(const StabilizationCore &) = delete;

   void ComputeStabilization(
      const mfem::Vector &state, mfem::Vector &result,
      const mfem::Array<bool> *active = nullptr) const;

   void CompDecay(const mfem::Vector &state, mfem::Vector &result,
                  double decay_time,
                  const mfem::Array<bool> *active = nullptr) const;

private:
   class Implementation;
   std::unique_ptr<Implementation> implementation;
};

} // namespace ofdg
