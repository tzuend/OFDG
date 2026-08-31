#include "oedg_2024.hpp"

#include "ofdg_core.hpp"

#include <utility>

namespace ofdg
{

namespace
{

detail::OFDGSensorOptions PaperOptions()
{
   detail::OFDGSensorOptions options;
   options.face_rule = detail::OFDGSensorOptions::FaceRule::Trapezoidal;
   options.pool_components = true;
   options.use_face_height = true;
   options.use_element_mean_speed = true;
   options.relative_constant_tolerance = true;
   return options;
}

} // namespace

class OEDG2024::Implementation
{
public:
   StabilizationCore stabilization;

   Implementation(const mfem::FiniteElementSpace *fes, int basis_type,
                  std::shared_ptr<const FacePhysics> face_physics)
      : stabilization(fes, basis_type, std::move(face_physics), PaperOptions())
   {
   }
};

OEDG2024::OEDG2024(const mfem::FiniteElementSpace *fes, int basis_type,
                   std::shared_ptr<const FacePhysics> face_physics)
   : implementation(std::make_unique<Implementation>(
        fes, basis_type, std::move(face_physics)))
{
}

OEDG2024::~OEDG2024() = default;
OEDG2024::OEDG2024(OEDG2024 &&) noexcept = default;
OEDG2024 &OEDG2024::operator=(OEDG2024 &&) noexcept = default;

void OEDG2024::ComputeStabilization(
   const mfem::Vector &state, mfem::Vector &result,
   const mfem::Array<bool> *active) const
{
   implementation->stabilization.ComputeStabilization(state, result, active);
}

void OEDG2024::CompDecay(const mfem::Vector &state, mfem::Vector &result,
                         double decay_time,
                         const mfem::Array<bool> *active) const
{
   implementation->stabilization.CompDecay(state, result, decay_time, active);
}

} // namespace ofdg
