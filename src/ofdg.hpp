#pragma once

#include "face_physics.hpp"
#include "mfem.hpp"

#include <memory>

namespace ofdg
{

/** Oscillation-free DG stabilization and exact polynomial-shell decay. */
class OFDG
{
public:
   OFDG(const mfem::FiniteElementSpace *fes, int basis_type);
   OFDG(const mfem::FiniteElementSpace *fes, int basis_type,
        mfem::VectorFunctionCoefficient velocity);
   OFDG(const mfem::FiniteElementSpace *fes, int basis_type,
        std::shared_ptr<const FacePhysics> face_physics);

   ~OFDG();
   OFDG(OFDG &&) noexcept;
   OFDG &operator=(OFDG &&) noexcept;
   OFDG(const OFDG &) = delete;
   OFDG &operator=(const OFDG &) = delete;

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
