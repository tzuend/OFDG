#include "ofdg.hpp"

#include "ofdg_core.hpp"

#include <utility>

namespace ofdg
{

class OFDG::Implementation
{
public:
   StabilizationCore stabilization;

   Implementation(const mfem::FiniteElementSpace *fes, int basis_type,
                  std::shared_ptr<const FacePhysics> face_physics)
      : stabilization(fes, basis_type, std::move(face_physics))
   {
   }
};

OFDG::OFDG(const mfem::FiniteElementSpace *fes, int basis_type)
   : implementation(std::make_unique<Implementation>(
        fes, basis_type, std::make_shared<UnitFacePhysics>()))
{
}

OFDG::OFDG(const mfem::FiniteElementSpace *fes, int basis_type,
           mfem::VectorFunctionCoefficient velocity)
   : implementation(std::make_unique<Implementation>(
        fes, basis_type,
        std::make_shared<AdvectionFacePhysics>(
           std::static_pointer_cast<mfem::VectorCoefficient>(
              std::make_shared<mfem::VectorFunctionCoefficient>(
                 std::move(velocity))))))
{
}

OFDG::OFDG(const mfem::FiniteElementSpace *fes, int basis_type,
           std::shared_ptr<const FacePhysics> face_physics)
   : implementation(std::make_unique<Implementation>(
        fes, basis_type, std::move(face_physics)))
{
}

OFDG::~OFDG() = default;
OFDG::OFDG(OFDG &&) noexcept = default;
OFDG &OFDG::operator=(OFDG &&) noexcept = default;

void OFDG::ComputeStabilization(const mfem::Vector &state,
                                mfem::Vector &result,
                                const mfem::Array<bool> *active) const
{
   implementation->stabilization.ComputeStabilization(state, result, active);
}

void OFDG::CompDecay(const mfem::Vector &state, mfem::Vector &result,
                     double decay_time,
                     const mfem::Array<bool> *active) const
{
   implementation->stabilization.CompDecay(state, result, decay_time, active);
}

} // namespace ofdg
