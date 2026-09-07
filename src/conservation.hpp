#pragma once

#include "mfem.hpp"
#include "curved_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

/** Globally integrate one scalar discontinuous GridFunction. */
inline mfem::real_t GlobalScalarIntegral(const mfem::GridFunction &field)
{
    const mfem::FiniteElementSpace *space = field.FESpace();
    mfem::Mesh *mesh = space->GetMesh();
    mfem::real_t integral = 0.0;
    for (int element = 0; element < space->GetNE(); ++element) {
        const mfem::FiniteElement *fe = space->GetFE(element);
        mfem::ElementTransformation *transformation =
            mesh->GetElementTransformation(element);
        const mfem::IntegrationRule &rule = ofdg::detail::VolumeRule(fe->GetOrder(), *transformation);
        for (int q = 0; q < rule.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &point = rule.IntPoint(q);
            transformation->SetIntPoint(&point);
            integral += point.weight * transformation->Weight() *
                        field.GetValue(element, point);
        }
    }
#ifdef MFEM_USE_MPI
    if (const auto *parallel_space =
            dynamic_cast<const mfem::ParFiniteElementSpace *>(space)) {
        MPI_Allreduce(MPI_IN_PLACE, &integral, 1,
                      mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM,
                      parallel_space->GetComm());
    }
#endif
    return integral;
}

inline mfem::real_t RelativeConservationDrift(mfem::real_t initial,
                                               mfem::real_t final)
{
    // Unit scaling keeps an initially zero conserved component meaningful:
    // its drift is then reported as an absolute nondimensional residual.
    const mfem::real_t scale = std::max(mfem::real_t(1.0), std::abs(initial));
    return std::abs(final - initial) / scale;
}
