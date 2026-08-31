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

class FilterCore
{
public:
   FilterCore(const mfem::FiniteElementSpace *fes, int basis_type,
              std::shared_ptr<const FacePhysics> face_physics);
   FilterCore(const mfem::FiniteElementSpace *fes, int basis_type,
              std::shared_ptr<const FacePhysics> face_physics,
              detail::OFDGSensorOptions options);

   ~FilterCore();
   FilterCore(FilterCore &&) noexcept;
   FilterCore &operator=(FilterCore &&) noexcept;
   FilterCore(const FilterCore &) = delete;
   FilterCore &operator=(const FilterCore &) = delete;

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
