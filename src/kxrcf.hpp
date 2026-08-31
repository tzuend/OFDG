#pragma once

#include "face_physics.hpp"
#include "mfem.hpp"

#include <memory>

namespace ofdg
{

/** KXRCF troubled-cell indicator for scalar fields and systems. */
class KXRCFIndicator
{
public:
   KXRCFIndicator(mfem::FiniteElementSpace *fes,
                  mfem::VectorCoefficient *velocity,
                  mfem::real_t threshold = 1.0,
                  mfem::real_t relative_scale_floor = 1e-14);

   KXRCFIndicator(mfem::FiniteElementSpace *fes,
                  std::shared_ptr<const FacePhysics> face_physics,
                  mfem::real_t threshold = 1.0,
                  mfem::real_t relative_scale_floor = 1e-14);

   ~KXRCFIndicator();
   KXRCFIndicator(KXRCFIndicator &&) noexcept;
   KXRCFIndicator &operator=(KXRCFIndicator &&) noexcept;
   KXRCFIndicator(const KXRCFIndicator &) = delete;
   KXRCFIndicator &operator=(const KXRCFIndicator &) = delete;

   void Compute(
      const mfem::Vector &state, mfem::Array<bool> &active,
      mfem::Vector *indicator_values = nullptr,
      mfem::DenseMatrix *component_indicator_values = nullptr) const;

private:
   class Implementation;
   std::unique_ptr<Implementation> implementation;
};

} // namespace ofdg
