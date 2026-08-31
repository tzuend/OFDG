#pragma once

#include "face_physics.hpp"
#include "mfem.hpp"

#include <memory>

namespace ofdg
{

/** Oscillation-eliminating DG method of Peng, Sun, and Wu (2024). */
class OEDG2024
{
public:
   OEDG2024(const mfem::FiniteElementSpace *fes, int basis_type,
            std::shared_ptr<const FacePhysics> face_physics);

   ~OEDG2024();
   OEDG2024(OEDG2024 &&) noexcept;
   OEDG2024 &operator=(OEDG2024 &&) noexcept;
   OEDG2024(const OEDG2024 &) = delete;
   OEDG2024 &operator=(const OEDG2024 &) = delete;

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
