#pragma once

#include "ofdg.hpp"

#include <memory>

namespace ofdg
{

/**
 * Oscillation-eliminating DG filter of Peng, Sun, and Wu (2024).
 *
 * The implementation deliberately exposes a separate type even though it
 * reuses the verified projection, derivative, face-halo, and exact-decay
 * kernel of OFDG.  The mathematical choices which distinguish the 2024
 * comparator are fixed here and cannot silently drift with driver flags:
 *
 *  - scale-invariant normalization by the global deviation from the mean;
 *  - normal characteristic wave speeds;
 *  - interface L1 jumps, with the endpoint trapezoidal rule in two dimensions;
 *  - one pooled damping coefficient for the full system state.
 *
 * Paper-faithful use applies CompDecay after every Runge--Kutta stage.  The
 * class itself is independent of a time integrator so cadence experiments can
 * also apply it once per complete step.
 */
class OEDG2024 {
private:
    OFDG filter;

    static detail::OFDGSensorOptions PaperOptions()
    {
        detail::OFDGSensorOptions options;
        options.face_rule = detail::OFDGSensorOptions::FaceRule::Trapezoidal;
        options.pool_components = true;
        options.use_face_height = true;
        return options;
    }

public:
    OEDG2024(const mfem::FiniteElementSpace *fes, int basis_type,
             std::shared_ptr<const FacePhysics> face_physics)
        : filter(fes, basis_type, std::move(face_physics), PaperOptions())
    {
    }

    void ComputeStabilization(const mfem::Vector &state, mfem::Vector &result,
                              const mfem::Array<bool> *active = nullptr) const
    {
        filter.ComputeStabilization(state, result, active);
    }

    void CompDecay(const mfem::Vector &state, mfem::Vector &result,
                   double decay_time,
                   const mfem::Array<bool> *active = nullptr) const
    {
        filter.CompDecay(state, result, decay_time, active);
    }
};

} // namespace ofdg
