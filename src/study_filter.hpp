#pragma once

#include "kxrcf.hpp"
#include "oedg_2024.hpp"
#include "ofdg.hpp"
#include "study_method.hpp"

#include <memory>

namespace ofdg
{

struct StudyFilterStatistics {
    long applications = 0;
    long active_elements = 0;
    long indicator_evaluations = 0;
};

/** Uniform dispatch for DG, adapted OFDG, KXRCF-OFDG, and 2024 OEDG. */
class StudyFilter {
private:
    const mfem::FiniteElementSpace *space;
    StudyMethod method;
    FilterCadence cadence;
    std::unique_ptr<OFDG> ofdg;
    std::unique_ptr<OEDG2024> oedg;
    std::unique_ptr<KXRCFIndicator> indicator;
    mfem::Array<bool> active;
    mfem::Vector result;
    StudyFilterStatistics statistics;

public:
    StudyFilter(mfem::FiniteElementSpace *space, int basis_type,
                std::shared_ptr<const FacePhysics> face_physics,
                StudyMethod method_, FilterCadence cadence_,
                double kxrcf_threshold = 1.0)
        : space(space), method(method_),
          cadence(ResolveFilterCadence(method_, cadence_))
    {
        if (method == StudyMethod::OFDG ||
            method == StudyMethod::OFDGKXRCF) {
            ofdg = std::make_unique<OFDG>(space, basis_type, face_physics);
        }
        if (method == StudyMethod::OEDG2024) {
            oedg = std::make_unique<OEDG2024>(space, basis_type,
                                             std::move(face_physics));
        }
        if (method == StudyMethod::OFDGKXRCF) {
            indicator = std::make_unique<KXRCFIndicator>(
                space, std::move(face_physics), kxrcf_threshold);
        }
    }

    bool Enabled() const { return method != StudyMethod::DG; }
    FilterCadence Cadence() const { return cadence; }

    bool Apply(mfem::Vector &state, double step_size, bool final_stage)
    {
        if (!Enabled()) { return true; }
        if (cadence == FilterCadence::Step && !final_stage) { return true; }

        const mfem::Array<bool> *mask = nullptr;
        if (indicator) {
            indicator->Compute(state, active);
            mask = &active;
            ++statistics.indicator_evaluations;
            for (int e = 0; e < active.Size(); ++e) {
                statistics.active_elements += active[e] ? 1 : 0;
            }
        }

        if (oedg) {
            oedg->CompDecay(state, result, step_size);
        } else {
            ofdg->CompDecay(state, result, step_size, mask);
        }
        state = result;
        ++statistics.applications;
        return true;
    }

    StudyFilterStatistics Statistics() const
    {
        StudyFilterStatistics result = statistics;
#ifdef MFEM_USE_MPI
        if (const auto *parallel_space =
                dynamic_cast<const mfem::ParFiniteElementSpace *>(space)) {
            long counts[2] = {result.applications,
                              result.indicator_evaluations};
            MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_LONG, MPI_MAX,
                          parallel_space->GetComm());
            MPI_Allreduce(MPI_IN_PLACE, &result.active_elements, 1,
                          MPI_LONG, MPI_SUM, parallel_space->GetComm());
            result.applications = counts[0];
            result.indicator_evaluations = counts[1];
        }
#endif
        return result;
    }
};

} // namespace ofdg
