#pragma once

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

/** Diagnostics returned by the conservative Euler positivity limiter. */
struct EulerPositivityDiagnostics {
    int limited_elements = 0;
    int inadmissible_means = 0;
    double minimum_density = std::numeric_limits<double>::infinity();
    double minimum_pressure = std::numeric_limits<double>::infinity();
};

/**
 * Cell-average-preserving scaling limiter for compressible Euler states.
 *
 * Every local polynomial is scaled towards its exact cell average.  The
 * admissibility check covers volume quadrature points and all face quadrature
 * points used by the DG operator.  Density is limited first; pressure then
 * uses a monotone bisection on the same line segment.  If a cell average is
 * itself inadmissible the limiter reports it and leaves that cell unchanged;
 * the time-step controller must reject the step.
 */
class EulerPositivityLimiter {
private:
    const mfem::FiniteElementSpace *fes;
    int dim;
    int ncomp;
    int order;
    double gamma;
    double density_floor;
    double pressure_floor;
    std::vector<std::vector<int>> local_faces;
#ifdef MFEM_USE_MPI
    std::vector<std::vector<int>> shared_faces;
#endif

    static double Pressure(const mfem::Vector &state, int dimension,
                           double heat_capacity_ratio)
    {
        const double density = state(0);
        if (!(density > 0.0) || !std::isfinite(density)) {
            return -std::numeric_limits<double>::infinity();
        }
        double momentum_squared = 0.0;
        for (int d = 0; d < dimension; ++d) {
            momentum_squared += state(1 + d) * state(1 + d);
        }
        return (heat_capacity_ratio - 1.0) *
               (state(dimension + 1) - 0.5 * momentum_squared / density);
    }

    template <typename Function>
    void ForEachEvaluationPoint(int element, Function &&function) const
    {
        const mfem::FiniteElement *fe = fes->GetFE(element);
        const mfem::IntegrationRule &volume_rule =
            mfem::IntRules.Get(fe->GetGeomType(), 2 * order + 1);
        for (int q = 0; q < volume_rule.GetNPoints(); ++q) {
            function(volume_rule.IntPoint(q));
        }

        mfem::Mesh *mesh = fes->GetMesh();
        for (int face : local_faces[element]) {
            mfem::FaceElementTransformations *transformations =
                mesh->GetFaceElementTransformations(face);
            if (!transformations) { continue; }
            const mfem::IntegrationRule &face_rule = mfem::IntRules.Get(
                transformations->GetGeometryType(), 2 * order + 1);
            for (int q = 0; q < face_rule.GetNPoints(); ++q) {
                transformations->SetAllIntPoints(&face_rule.IntPoint(q));
                if (transformations->Elem1No == element) {
                    function(transformations->GetElement1IntPoint());
                } else if (transformations->Elem2No == element) {
                    function(transformations->GetElement2IntPoint());
                }
            }
        }
#ifdef MFEM_USE_MPI
        if (auto *parallel_mesh = dynamic_cast<mfem::ParMesh *>(mesh)) {
            for (int face : shared_faces[element]) {
                mfem::FaceElementTransformations *transformations =
                    parallel_mesh->GetSharedFaceTransformations(face, false);
                const mfem::IntegrationRule &face_rule = mfem::IntRules.Get(
                    transformations->GetGeometryType(), 2 * order + 1);
                for (int q = 0; q < face_rule.GetNPoints(); ++q) {
                    transformations->SetAllIntPoints(&face_rule.IntPoint(q));
                    function(transformations->GetElement1IntPoint());
                }
            }
        }
#endif
    }

public:
    EulerPositivityLimiter(const mfem::FiniteElementSpace *space,
                           double heat_capacity_ratio,
                           double density_tolerance = 1e-12,
                           double pressure_tolerance = 1e-12)
        : fes(space), dim(space->GetMesh()->Dimension()),
          ncomp(space->GetVDim()), order(space->GetFE(0)->GetOrder()),
          gamma(heat_capacity_ratio), density_floor(density_tolerance),
          pressure_floor(pressure_tolerance)
    {
        MFEM_VERIFY(fes != nullptr, "Positivity limiter requires a finite element space.");
        MFEM_VERIFY(ncomp == dim + 2,
                    "Euler positivity limiter expects dim+2 conservative variables.");
        MFEM_VERIFY(gamma > 1.0, "Euler positivity limiter requires gamma > 1.");

        mfem::Mesh *mesh = fes->GetMesh();
        local_faces.resize(fes->GetNE());
        for (int face = 0; face < mesh->GetNumFaces(); ++face) {
            mfem::FaceElementTransformations *transformations =
                mesh->GetFaceElementTransformations(face);
            if (!transformations) { continue; }
            if (transformations->Elem1No >= 0) {
                local_faces[transformations->Elem1No].push_back(face);
            }
            if (transformations->Elem2No >= 0) {
                local_faces[transformations->Elem2No].push_back(face);
            }
        }
#ifdef MFEM_USE_MPI
        shared_faces.resize(fes->GetNE());
        if (auto *parallel_mesh = dynamic_cast<mfem::ParMesh *>(mesh)) {
            // Shared-face metadata is populated lazily by MFEM.  The limiter
            // may be the first parallel component constructed by a driver, so
            // it must not rely on another operator having initialized it.
            parallel_mesh->ExchangeFaceNbrData();
            for (int face = 0; face < parallel_mesh->GetNSharedFaces(); ++face) {
                mfem::FaceElementTransformations *transformations =
                    parallel_mesh->GetSharedFaceTransformations(face, false);
                shared_faces[transformations->Elem1No].push_back(face);
            }
        }
#endif
    }

    EulerPositivityDiagnostics Apply(mfem::Vector &state) const
    {
        MFEM_VERIFY(state.Size() == fes->GetVSize(),
                    "Positivity limiter state size does not match its space.");

        EulerPositivityDiagnostics diagnostics;
        mfem::Array<int> scalar_dofs;
        mfem::Array<int> component_dofs;
        mfem::DenseMatrix coefficients;
        mfem::Vector shape;
        mfem::Vector value(ncomp);
        mfem::Vector mean(ncomp);
        mfem::Vector density_limited(ncomp);

        for (int element = 0; element < fes->GetNE(); ++element) {
            const mfem::FiniteElement *fe = fes->GetFE(element);
            const int ndof = fe->GetDof();
            fes->GetElementDofs(element, scalar_dofs);
            coefficients.SetSize(ndof, ncomp);
            component_dofs.SetSize(ndof);
            for (int c = 0; c < ncomp; ++c) {
                for (int i = 0; i < ndof; ++i) {
                    component_dofs[i] = fes->DofToVDof(scalar_dofs[i], c);
                }
                mfem::Vector column(coefficients.GetColumn(c), ndof);
                state.GetSubVector(component_dofs, column);
            }

            mean = 0.0;
            double volume = 0.0;
            mfem::ElementTransformation *transformation =
                fes->GetMesh()->GetElementTransformation(element);
            const mfem::IntegrationRule &mean_rule =
                mfem::IntRules.Get(fe->GetGeomType(), 2 * order + 1);
            shape.SetSize(ndof);
            for (int q = 0; q < mean_rule.GetNPoints(); ++q) {
                const mfem::IntegrationPoint &ip = mean_rule.IntPoint(q);
                transformation->SetIntPoint(&ip);
                fe->CalcShape(ip, shape);
                const double weight = ip.weight * transformation->Weight();
                volume += weight;
                for (int c = 0; c < ncomp; ++c) {
                    const double *column = coefficients.GetColumn(c);
                    double point_value = 0.0;
                    for (int i = 0; i < ndof; ++i) {
                        point_value += shape(i) * column[i];
                    }
                    mean(c) += weight * point_value;
                }
            }
            mean /= volume;

            const double mean_pressure = Pressure(mean, dim, gamma);
            if (!(mean(0) > density_floor) ||
                !(mean_pressure > pressure_floor) ||
                !std::isfinite(mean_pressure)) {
                ++diagnostics.inadmissible_means;
                continue;
            }

            double minimum_density = std::numeric_limits<double>::infinity();
            auto evaluate = [&](const mfem::IntegrationPoint &ip,
                                const mfem::DenseMatrix &data,
                                mfem::Vector &output) {
                fe->CalcShape(ip, shape);
                for (int c = 0; c < ncomp; ++c) {
                    const double *column = data.GetColumn(c);
                    output(c) = 0.0;
                    for (int i = 0; i < ndof; ++i) {
                        output(c) += shape(i) * column[i];
                    }
                }
            };
            ForEachEvaluationPoint(element, [&](const mfem::IntegrationPoint &ip) {
                evaluate(ip, coefficients, value);
                minimum_density = std::min(minimum_density, value(0));
            });

            double theta_density = 1.0;
            if (minimum_density < density_floor) {
                theta_density = std::clamp(
                    (mean(0) - density_floor) /
                    (mean(0) - minimum_density), 0.0, 1.0);
            }

            double theta_pressure = 1.0;
            ForEachEvaluationPoint(element, [&](const mfem::IntegrationPoint &ip) {
                evaluate(ip, coefficients, value);
                for (int c = 0; c < ncomp; ++c) {
                    density_limited(c) =
                        mean(c) + theta_density * (value(c) - mean(c));
                }
                const double pressure = Pressure(density_limited, dim, gamma);
                diagnostics.minimum_density =
                    std::min(diagnostics.minimum_density, density_limited(0));
                diagnostics.minimum_pressure =
                    std::min(diagnostics.minimum_pressure, pressure);
                if (pressure >= pressure_floor) { return; }

                double lower = 0.0;
                double upper = 1.0;
                mfem::Vector trial(ncomp);
                for (int iteration = 0; iteration < 60; ++iteration) {
                    const double middle = 0.5 * (lower + upper);
                    for (int c = 0; c < ncomp; ++c) {
                        trial(c) = mean(c) + middle *
                                   (density_limited(c) - mean(c));
                    }
                    if (Pressure(trial, dim, gamma) >= pressure_floor) {
                        lower = middle;
                    } else {
                        upper = middle;
                    }
                }
                theta_pressure = std::min(theta_pressure, lower);
            });

            const double theta = theta_density * theta_pressure;
            if (theta < 1.0 - 1e-14) {
                ++diagnostics.limited_elements;
                for (int c = 0; c < ncomp; ++c) {
                    for (int i = 0; i < ndof; ++i) {
                        coefficients(i, c) =
                            mean(c) + theta * (coefficients(i, c) - mean(c));
                    }
                    for (int i = 0; i < ndof; ++i) {
                        component_dofs[i] = fes->DofToVDof(scalar_dofs[i], c);
                    }
                    mfem::Vector column(coefficients.GetColumn(c), ndof);
                    state.SetSubVector(component_dofs, column);
                }
            }
        }

#ifdef MFEM_USE_MPI
        if (const auto *parallel_space =
                dynamic_cast<const mfem::ParFiniteElementSpace *>(fes)) {
            MPI_Comm communicator = parallel_space->GetComm();
            MPI_Allreduce(MPI_IN_PLACE, &diagnostics.limited_elements, 1,
                          MPI_INT, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &diagnostics.inadmissible_means, 1,
                          MPI_INT, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &diagnostics.minimum_density, 1,
                          mfem::MPITypeMap<mfem::real_t>::mpi_type,
                          MPI_MIN, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &diagnostics.minimum_pressure, 1,
                          mfem::MPITypeMap<mfem::real_t>::mpi_type,
                          MPI_MIN, communicator);
        }
#endif
        return diagnostics;
    }
};
