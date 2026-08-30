#pragma once

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

using namespace mfem;

// =============================================================================
// Internal timing
// =============================================================================

struct OFDGInternalTiming {
    double compute_jumps = 0.0;
    double global_scaling = 0.0;
    double normal_speed = 0.0;
    double derivative_build = 0.0;
    double face_trace = 0.0;
    double sigma_accumulation = 0.0;
    double decay_application = 0.0;
    double stabilization_application = 0.0;
};

#ifdef OFDG_INTERNAL_TIMING

using OFDGTimingClock = std::chrono::steady_clock;

inline double OFDGSecondsSince(const OFDGTimingClock::time_point &begin)
{
    return std::chrono::duration<double>(OFDGTimingClock::now() - begin).count();
}

#endif

// =============================================================================
// Projection operators
//
// Builds the polynomial projection matrices used by OFDG:
//
//     Q_l = M_k - B_kl^T M_l^{-1} B_kl
//     S_l = M_k^{-1} Q_l
//
// This class contains only projection mathematics.
// =============================================================================

class OFDGProjection {
private:
    const FiniteElementSpace *fes;
    const int btype;

    int dim;
    int order;
    Geometry::Type geom;

public:
    OFDGProjection(const FiniteElementSpace *fes_, int btype_) : fes(fes_), btype(btype_)
    {
        const FiniteElement *fe = fes->GetFE(0);

        dim = fe->GetDim();
        order = fe->GetOrder();
        geom = fe->GetGeomType();
    }

    DenseMatrix AssembleMk(int k, int e = 0) const
    {
        DG_FECollection fec_k(k, dim, btype);

        const FiniteElement *fe_k = fec_k.FiniteElementForGeometry(geom);

        ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);

        const IntegrationRule &ir = IntRules.Get(fe_k->GetGeomType(), 2 * fe_k->GetOrder() + 1);

        DenseMatrix M_k(fe_k->GetDof());
        M_k = 0.0;

        Vector shape_k(fe_k->GetDof());

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const IntegrationPoint &ip = ir.IntPoint(q);

            T->SetIntPoint(&ip);
            fe_k->CalcShape(ip, shape_k);

            const double w = ip.weight * T->Weight();
            AddMult_a_VVt(w, shape_k, M_k);
        }

        return M_k;
    }

    DenseMatrix AssembleBkl(int k, int l, int e = 0) const
    {
        DG_FECollection fec_k(k, dim, btype);
        DG_FECollection fec_l(l, dim, btype);

        const FiniteElement *fe_k = fec_k.FiniteElementForGeometry(geom);

        const FiniteElement *fe_l = fec_l.FiniteElementForGeometry(geom);

        ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);

        const IntegrationRule &ir =
            IntRules.Get(fe_k->GetGeomType(), 2 * std::max(fe_k->GetOrder(), fe_l->GetOrder()) + 1);

        DenseMatrix B(fe_l->GetDof(), fe_k->GetDof());
        B = 0.0;

        Vector shape_k(fe_k->GetDof());
        Vector shape_l(fe_l->GetDof());

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const IntegrationPoint &ip = ir.IntPoint(q);

            T->SetIntPoint(&ip);

            fe_k->CalcShape(ip, shape_k);
            fe_l->CalcShape(ip, shape_l);

            const double w = ip.weight * T->Weight();
            AddMult_a_VWt(w, shape_l, shape_k, B);
        }

        return B;
    }

    DenseMatrix AssembleQl(int l, int e = 0) const
    {
        DenseMatrix M_k = AssembleMk(order, e);
        DenseMatrix M_l = AssembleMk(l, e);
        DenseMatrix B_l = AssembleBkl(order, l, e);

        const int N_k = M_k.Width();
        const int N_l = M_l.Width();

        DenseMatrix M_l_inv(M_l);
        M_l_inv.Invert();

        DenseMatrix MlInvB(N_l, N_k);
        Mult(M_l_inv, B_l, MlInvB);

        DenseMatrix BtMlInvB(N_k);
        MultAtB(B_l, MlInvB, BtMlInvB);

        DenseMatrix Q_l = M_k;
        Q_l -= BtMlInvB;

        return Q_l;
    }
};

// =============================================================================
// Reference derivative construction
//
// At interpolation points xi_q:
//
//     V(q,i)   = phi_i(xi_q)
//     G_d(q,i) = d phi_i / d xi_d (xi_q)
//
// G_d maps polynomial coefficients to derivative VALUES. To apply derivatives
// repeatedly we require derivative COEFFICIENTS:
//
//     V u_d = G_d u
//
// therefore
//
//     D_d = V^{-1} G_d
//
// and D_d maps coefficients -> derivative coefficients.
// =============================================================================

class OFDGReferenceDerivative {
private:
    const FiniteElementSpace *fes;

public:
    explicit OFDGReferenceDerivative(const FiniteElementSpace *fes_) : fes(fes_)
    {
    }

    DenseMatrix AssembleVandermonde(int e = 0, const IntegrationRule *ir = nullptr) const
    {
        const FiniteElement *fe = fes->GetFE(e);
        const int ndof = fe->GetDof();

        if (ir == nullptr) {
            ir = &fe->GetNodes();
        }

        DenseMatrix V(ndof);
        Vector shape(ndof);

        for (int q = 0; q < ndof; ++q) {
            const IntegrationPoint &ip = ir->IntPoint(q);

            fe->CalcShape(ip, shape);

            for (int i = 0; i < ndof; ++i) {
                V(q, i) = shape(i);
            }
        }

        return V;
    }

    std::vector<DenseMatrix> Assemble(int e = 0, const IntegrationRule *ir = nullptr) const
    {
        const FiniteElement *fe = fes->GetFE(e);

        const int ndof = fe->GetDof();
        const int dim = fe->GetDim();

        if (ir == nullptr) {
            ir = &fe->GetNodes();
        }

        std::vector<DenseMatrix> derivative_values;
        derivative_values.reserve(dim);

        for (int d = 0; d < dim; ++d) {
            derivative_values.emplace_back(ndof);
            derivative_values.back() = 0.0;
        }

        DenseMatrix dshape(ndof, dim);

        for (int q = 0; q < ndof; ++q) {
            const IntegrationPoint &ip = ir->IntPoint(q);

            fe->CalcDShape(ip, dshape);

            for (int i = 0; i < ndof; ++i) {
                for (int d = 0; d < dim; ++d) {
                    derivative_values[d](q, i) = dshape(i, d);
                }
            }
        }

        DenseMatrix V_inv = AssembleVandermonde(e, ir);
        V_inv.Invert();

        std::vector<DenseMatrix> D_ref;
        D_ref.reserve(dim);

        // D_xi = V^{-1} G_xi
        for (int d = 0; d < dim; ++d) {
            D_ref.emplace_back(ndof);
            Mult(V_inv, derivative_values[d], D_ref.back());
        }

        return D_ref;
    }
};

// =============================================================================
// OFDG mathematical operators
//
// Contains mesh-independent / reference-element mathematics:
//
//   - S_l projection operators
//   - reference derivative matrices
//   - derivative DAG
//   - sensor constants
//   - volume evaluation matrix
//
// No derivative cache or cache policy lives here.
// =============================================================================

struct OFDGDerivativeNode {
    int degree = 0;
    int parent = -1;
    int physical_direction = -1;
};

class OFDGOperators {
private:
    const FiniteElementSpace *fes;

    int dim;
    int order;
    int ndof;
    Geometry::Type geom;

    OFDGProjection projection;

    std::vector<DenseMatrix> S_l;
    std::vector<DenseMatrix> reference_derivatives;
    std::vector<OFDGDerivativeNode> derivative_nodes;

    Vector sensor_common;

    DenseMatrix volume_evaluation;
    Vector volume_reference_weights;
    double reference_volume = 0.0;

    static void GenerateMultiIndices(int dim, int pos, int remaining, std::vector<int> &alpha,
                                     std::vector<std::vector<int>> &out)
    {
        if (pos == dim - 1) {
            alpha[pos] = remaining;
            out.push_back(alpha);
            return;
        }

        for (int i = 0; i <= remaining; ++i) {
            alpha[pos] = i;
            GenerateMultiIndices(dim, pos + 1, remaining - i, alpha, out);
        }
    }

    std::size_t EncodeMultiIndex(const std::vector<int> &alpha) const
    {
        const std::size_t base = static_cast<std::size_t>(order + 1);

        std::size_t key = 0;
        std::size_t stride = 1;

        for (int d = 0; d < dim; ++d) {
            key += static_cast<std::size_t>(alpha[d]) * stride;
            stride *= base;
        }

        return key;
    }

    void BuildProjectionOperators()
    {
        DenseMatrix M_k = projection.AssembleMk(order, 0);
        DenseMatrix M_k_inv = M_k;

        M_k_inv.Invert();

        S_l.reserve(order + 1);

        for (int l = 0; l <= order; ++l) {
            DenseMatrix S(ndof);
            DenseMatrix Q_l = projection.AssembleQl(l, 0);

            Mult(M_k_inv, Q_l, S);
            S_l.push_back(S);
        }
    }

    void BuildReferenceDerivatives()
    {
        OFDGReferenceDerivative derivative_builder(fes);
        reference_derivatives = derivative_builder.Assemble(0);
    }

    void BuildDerivativeDAG()
    {
        derivative_nodes.clear();

        std::unordered_map<std::size_t, int> node_for_alpha;

        std::vector<int> zero(dim, 0);

        derivative_nodes.push_back({0, -1, -1});
        node_for_alpha.emplace(EncodeMultiIndex(zero), 0);

        for (int l = 1; l <= order; ++l) {
            std::vector<std::vector<int>> indices;
            std::vector<int> alpha(dim, 0);

            GenerateMultiIndices(dim, 0, l, alpha, indices);

            for (const auto &beta : indices) {
                int direction = -1;

                // Preserve the derivative ordering of the previous code.
                for (int d = dim - 1; d >= 0; --d) {
                    if (beta[d] > 0) {
                        direction = d;
                        break;
                    }
                }

                std::vector<int> parent_alpha = beta;
                --parent_alpha[direction];

                const auto parent_it = node_for_alpha.find(EncodeMultiIndex(parent_alpha));

                const int node = static_cast<int>(derivative_nodes.size());

                derivative_nodes.push_back({l, parent_it->second, direction});

                node_for_alpha.emplace(EncodeMultiIndex(beta), node);
            }
        }
    }

    void BuildSensorConstants()
    {
        sensor_common.SetSize(order + 1);

        double factorial = 1.0;

        for (int l = 0; l <= order; ++l) {
            if (l > 0) {
                factorial *= static_cast<double>(l);
            }

            sensor_common(l) = (2.0 * l + 1.0) / (2.0 * (2.0 * order - 1.0) * factorial);
        }
    }

    void BuildVolumeEvaluation()
    {
        const FiniteElement *fe = fes->GetFE(0);

        const IntegrationRule &ir = IntRules.Get(geom, 2 * order + 1);

        const int nq = ir.GetNPoints();

        volume_evaluation.SetSize(nq, ndof);
        volume_reference_weights.SetSize(nq);

        Vector shape(ndof);

        reference_volume = 0.0;

        for (int q = 0; q < nq; ++q) {
            const IntegrationPoint &ip = ir.IntPoint(q);

            fe->CalcShape(ip, shape);

            for (int i = 0; i < ndof; ++i) {
                volume_evaluation(q, i) = shape(i);
            }

            volume_reference_weights(q) = ip.weight;
            reference_volume += ip.weight;
        }
    }

public:
    OFDGOperators(const FiniteElementSpace *fes_, int btype)
        : fes(fes_), dim(0), order(0), ndof(0), geom(Geometry::POINT), projection(fes_, btype)
    {
        const FiniteElement *fe = fes->GetFE(0);

        dim = fe->GetDim();
        order = fe->GetOrder();
        ndof = fe->GetDof();
        geom = fe->GetGeomType();

        BuildProjectionOperators();
        BuildReferenceDerivatives();
        BuildDerivativeDAG();
        BuildSensorConstants();
        BuildVolumeEvaluation();
    }

    int Dim() const
    {
        return dim;
    }

    int Order() const
    {
        return order;
    }

    int NDof() const
    {
        return ndof;
    }

    Geometry::Type Geometry() const
    {
        return geom;
    }

    int DerivativeCount() const
    {
        return static_cast<int>(derivative_nodes.size());
    }

    const DenseMatrix &S(int l) const
    {
        return S_l[l];
    }

    const std::vector<DenseMatrix> &ReferenceDerivatives() const
    {
        return reference_derivatives;
    }

    const std::vector<OFDGDerivativeNode> &DerivativeNodes() const
    {
        return derivative_nodes;
    }

    double SensorCommon(int l) const
    {
        return sensor_common(l);
    }

    const DenseMatrix &VolumeEvaluation() const
    {
        return volume_evaluation;
    }

    const Vector &VolumeReferenceWeights() const
    {
        return volume_reference_weights;
    }

    double ReferenceVolume() const
    {
        return reference_volume;
    }

    void TestCommutators() const
    {
        DenseMatrix AB(ndof);
        DenseMatrix BA(ndof);

        for (int i = 0; i <= order; ++i) {
            for (int j = 0; j <= order; ++j) {
                Mult(S_l[i], S_l[j], AB);
                Mult(S_l[j], S_l[i], BA);

                AB -= BA;

                std::cout << "||[Si,Sj]|| = " << AB.FNorm() << std::endl;
            }
        }
    }
};

// =============================================================================
// Static mesh data
//
// Preconditions:
//
//   - fes is a nonempty DG space of polynomial order >= 1;
//   - all elements share their geometry, order, dimension and dof layout;
//   - every physical element mapping is affine;
//   - the mesh geometry and topology do not change after construction.
//
// The implementation is generic over the MFEM element geometry. For affine
// mappings, J^{-1} and det(J) may differ between elements but are constant
// inside each element.
//
// Owns element geometry:
//
//   - vdofs
//   - h
//   - J^{-1}
//   - affine Jacobian weight
//   - volume
//   - sensor h-prefactors
//
// No derivative coefficients or runtime cache state lives here.
// =============================================================================

class OFDGMeshData {
public:
    struct ElementData {
        Array<int> vdofs;

        double h = 0.0;
        double jacobian_weight = 0.0;
        double volume = 0.0;

        // dim <= 3
        std::array<double, 9> Jinv{};
    };

private:
    const FiniteElementSpace *fes;
    const OFDGOperators &operators;

    int dim;
    int order;
    int ndof;
    Geometry::Type geom;

    std::vector<ElementData> elements;

    std::vector<double> sensor_prefactor;

    void BuildElements()
    {
        const int ne = fes->GetNE();

        elements.resize(ne);

        sensor_prefactor.resize(static_cast<std::size_t>(ne) * static_cast<std::size_t>(order + 1));

        for (int e = 0; e < ne; ++e) {
            ElementData &element = elements[e];

            fes->GetElementVDofs(e, element.vdofs);

            element.h = fes->GetMesh()->GetElementSize(e);

            ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);
            T->SetIntPoint(&Geometries.GetCenter(geom));
            const DenseMatrix &Jinv = T->InverseJacobian();

            for (int r = 0; r < dim; ++r) {
                for (int c = 0; c < dim; ++c) {
                    element.Jinv[r * dim + c] = Jinv(r, c);
                }
            }

            element.jacobian_weight = T->Weight();

            element.volume = operators.ReferenceVolume() * element.jacobian_weight;

            double h_power = 1.0 / element.h;

            for (int l = 0; l <= order; ++l) {
                sensor_prefactor[e * (order + 1) + l] = operators.SensorCommon(l) * h_power;

                h_power *= element.h;
            }
        }
    }

public:
    OFDGMeshData(const FiniteElementSpace *fes_, const OFDGOperators &operators_)
        : fes(fes_), operators(operators_), dim(operators_.Dim()), order(operators_.Order()),
          ndof(operators_.NDof()), geom(operators_.Geometry())
    {
        BuildElements();
    }

    const ElementData &Element(int e) const
    {
        return elements[e];
    }

    double ElementJinv(int e, int row, int col) const
    {
        return elements[e].Jinv[row * dim + col];
    }

    double SensorPrefactor(int e, int l) const
    {
        return sensor_prefactor[e * (order + 1) + l];
    }

    void GetComponentVDofs(int element, int component, Array<int> &vdofs) const
    {
        vdofs.SetSize(ndof);

        const Array<int> &all = elements[element].vdofs;
        const int offset = component * ndof;

        for (int i = 0; i < ndof; ++i) {
            vdofs[i] = all[offset + i];
        }
    }
};

// =============================================================================
// Runtime scratch storage
//
// Temporary memory is grouped here so the mathematical classes do not carry
// unrelated scratch vectors throughout their implementations.
// =============================================================================

struct OFDGVolumeScratch {
    Vector element_data;
    DenseMatrix element_matrix;

    DenseMatrix values;

    Vector mean_integrals;
    Vector minimum_values;
    Vector maximum_values;

    void Initialize(int ndof, int ncomp, int nq)
    {
        element_data.SetSize(ndof * ncomp);

        element_matrix.UseExternalData(element_data.GetData(), ndof, ncomp);

        values.SetSize(nq, ncomp);

        mean_integrals.SetSize(ncomp);
        minimum_values.SetSize(ncomp);
        maximum_values.SetSize(ncomp);
    }
};

struct OFDGJumpScratch {
    DenseMatrix traces1;
    DenseMatrix traces2;

    void SetSize(int nq, int trace_columns)
    {
        traces1.SetSize(nq, trace_columns);
        traces2.SetSize(nq, trace_columns);
    }
};

struct OFDGScratch {
    OFDGVolumeScratch volume;
    OFDGJumpScratch jump;

    DenseMatrix jumps;

    void Initialize(int ndof, int ncomp, int volume_nq, int order)
    {
        volume.Initialize(ndof, ncomp, volume_nq);
        jumps.SetSize(ncomp, order + 1);
    }
};

// =============================================================================
// Face evaluation
//
// Owns:
//
//   - face quadrature
//   - reference evaluation matrices E_face
//   - runtime trace evaluation
//   - derivative jump accumulation
//
// Connectivity remains owned by MFEM. Cached evaluation data is indexed by
// the MFEM face number and contains no duplicate element/face metadata.
// =============================================================================

class OFDGFaceEvaluation {
public:
    struct FaceEvaluationData {
        DenseMatrix evaluation1;
        DenseMatrix evaluation2;
        Vector normalized_weights;
    };

private:
    const FiniteElementSpace *fes;
    const OFDGOperators &operators;
    [[maybe_unused]] OFDGInternalTiming &timing;

    int order;
    int ndof;
    int ncomp;

    std::vector<FaceEvaluationData> face_data;
    int interior_face_count = 0;

    void BuildFaceEvaluationData(FaceElementTransformations *Tr, DenseMatrix &evaluation1,
                                 DenseMatrix &evaluation2, Vector &normalized_weights) const
    {
        const FiniteElement *fe1 = fes->GetFE(Tr->Elem1No);
        const FiniteElement *fe2 = fes->GetFE(Tr->Elem2No);

        const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * order + 1);

        const int nq = ir.GetNPoints();

        evaluation1.SetSize(nq, ndof);
        evaluation2.SetSize(nq, ndof);
        normalized_weights.SetSize(nq);

        Vector shape1(ndof);
        Vector shape2(ndof);

        double weight_sum = 0.0;

        for (int q = 0; q < nq; ++q) {
            weight_sum += ir.IntPoint(q).weight;
        }

        for (int q = 0; q < nq; ++q) {
            const IntegrationPoint &ip = ir.IntPoint(q);

            Tr->SetAllIntPoints(&ip);

            const IntegrationPoint &ip1 = Tr->GetElement1IntPoint();

            const IntegrationPoint &ip2 = Tr->GetElement2IntPoint();

            fe1->CalcShape(ip1, shape1);
            fe2->CalcShape(ip2, shape2);

            for (int i = 0; i < ndof; ++i) {
                evaluation1(q, i) = shape1(i);
                evaluation2(q, i) = shape2(i);
            }

            // For affine faces:
            //
            //     (1 / |f|) int_f jump^2 dS
            //
            // has the constant physical face Jacobian in numerator and
            // denominator, so it cancels.
            normalized_weights(q) = ip.weight / weight_sum;
        }
    }

    void BuildCachedFaces()
    {
        Mesh *mesh = fes->GetMesh();
        face_data.clear();
        face_data.resize(mesh->GetNumFaces());
        interior_face_count = 0;

        for (int f = 0; f < mesh->GetNumFaces(); ++f) {
            FaceElementTransformations *Tr = mesh->GetFaceElementTransformations(f);
            if (!Tr || Tr->Elem1No < 0 || Tr->Elem2No < 0) {
                continue;
            }

            FaceEvaluationData &data = face_data[f];
            BuildFaceEvaluationData(Tr, data.evaluation1, data.evaluation2,
                                    data.normalized_weights);
            ++interior_face_count;
        }
    }

    void Evaluate(const DenseMatrix &derivatives1, const DenseMatrix &derivatives2,
                  const FaceEvaluationData &data, DenseMatrix &jumps,
                  OFDGJumpScratch &scratch) const
    {
        const int nq = data.evaluation1.Height();
        const int trace_columns = operators.DerivativeCount() * ncomp;

        jumps.SetSize(ncomp, order + 1);
        jumps = 0.0;

        scratch.SetSize(nq, trace_columns);

#ifdef OFDG_INTERNAL_TIMING
        const auto timer_begin = OFDGTimingClock::now();
#endif

        // All derivative states and components are evaluated in one matrix
        // multiplication per side.
        Mult(data.evaluation1, derivatives1, scratch.traces1);

        Mult(data.evaluation2, derivatives2, scratch.traces2);

        const auto &nodes = operators.DerivativeNodes();

        for (int node = 0; node < static_cast<int>(nodes.size()); ++node) {
            const int l = nodes[node].degree;
            const int column_offset = node * ncomp;

            for (int q = 0; q < nq; ++q) {
                const double weight = data.normalized_weights(q);

                for (int c = 0; c < ncomp; ++c) {
                    const int column = column_offset + c;

                    const double jump = scratch.traces1(q, column) - scratch.traces2(q, column);

                    jumps(c, l) += weight * std::abs(jump);
                }
            }
        }

#ifdef OFDG_INTERNAL_TIMING
        timing.face_trace += OFDGSecondsSince(timer_begin);
#endif
    }

public:
    OFDGFaceEvaluation(const FiniteElementSpace *fes_, const OFDGOperators &operators_,
                       OFDGInternalTiming &timing_)
        : fes(fes_), operators(operators_), timing(timing_), order(operators_.Order()),
          ndof(operators_.NDof()), ncomp(fes_->GetVDim())
    {
        BuildCachedFaces();
    }

    int FaceCount() const
    {
        return interior_face_count;
    }

    std::size_t StorageBytes() const
    {
        std::size_t bytes = 0;

        for (const auto &data : face_data) {
            bytes += static_cast<std::size_t>(data.evaluation1.Height()) *
                     static_cast<std::size_t>(data.evaluation1.Width()) * sizeof(double);
            bytes += static_cast<std::size_t>(data.evaluation2.Height()) *
                     static_cast<std::size_t>(data.evaluation2.Width()) * sizeof(double);
            bytes += static_cast<std::size_t>(data.normalized_weights.Size()) * sizeof(double);
        }

        return bytes;
    }

    void EvaluateDerivativeJumps(const DenseMatrix &derivatives1, const DenseMatrix &derivatives2,
                                 int face, DenseMatrix &jumps, OFDGJumpScratch &scratch) const
    {
        Evaluate(derivatives1, derivatives2, face_data[face], jumps, scratch);
    }

    // Used by the public one-face testing API.
    void EvaluateDerivativeJumps(const DenseMatrix &derivatives1, const DenseMatrix &derivatives2,
                                 FaceElementTransformations *Tr, DenseMatrix &jumps,
                                 OFDGJumpScratch &scratch) const
    {
        FaceEvaluationData data;
        BuildFaceEvaluationData(Tr, data.evaluation1, data.evaluation2, data.normalized_weights);
        Evaluate(derivatives1, derivatives2, data, jumps, scratch);
    }
};

// =============================================================================
// Derivative state
//
// This is the DATA produced for one element.
//
// all_coefficients:
//
//     ndof x (derivative_count * ncomp)
//
// columns:
//
//     [node 0 components]
//     [node 1 components]
//     ...
//
// coefficient_views are zero-copy views into this matrix.
//
// Physical derivative operators are stored separately once per element.
// =============================================================================

struct OFDGDerivativeState {
    DenseMatrix all_coefficients;
    Vector root_vector;

    std::vector<DenseMatrix> coefficient_views;

    OFDGDerivativeState(int ndof, int ncomp, int derivative_count)
    {
        const int total_columns = derivative_count * ncomp;

        all_coefficients.SetSize(ndof, total_columns);

        all_coefficients = 0.0;

        // Node zero contains the original polynomial coefficients.
        root_vector.SetDataAndSize(all_coefficients.GetData(), ndof * ncomp);

        coefficient_views.resize(derivative_count);

        const std::size_t node_stride =
            static_cast<std::size_t>(ndof) * static_cast<std::size_t>(ncomp);

        for (int node = 0; node < derivative_count; ++node) {
            coefficient_views[node].UseExternalData(
                all_coefficients.GetData() + static_cast<std::size_t>(node) * node_stride, ndof,
                ncomp);
        }
    }
};

// =============================================================================
// Full derivative storage
//
// It owns the mathematics of constructing physical derivative coefficient
// states:
//
//     D_xj = sum_r (J^{-1})_{rj} D_xi_r
//
// and then walks the derivative DAG:
//
//     u_node = D_direction u_parent
//
// Every element has one permanent derivative state. BuildAll() rebuilds
// all states once for the current solution; face evaluation is then a direct
// element-index lookup with no cache branches or eviction policy.
// =============================================================================

class OFDGDerivativeProvider {
public:
    struct Pair {
        const OFDGDerivativeState *first = nullptr;
        const OFDGDerivativeState *second = nullptr;
    };

private:
    const FiniteElementSpace *fes;
    const OFDGOperators &operators;
    const OFDGMeshData &mesh_data;
    [[maybe_unused]] OFDGInternalTiming &timing;

    int dim;
    int ndof;
    int ncomp;

    std::vector<std::vector<DenseMatrix>> physical_derivatives;
    std::vector<OFDGDerivativeState> states;
    std::size_t storage_bytes = 0;

    void BuildPhysicalDerivativeMatrices(int element, std::vector<DenseMatrix> &D_phys) const
    {
        const auto &D_ref = operators.ReferenceDerivatives();

        D_phys.resize(dim);

        for (int x = 0; x < dim; ++x) {
            D_phys[x].SetSize(ndof);
            D_phys[x] = 0.0;

            // D_x = sum_xi (J^{-1})_{xi,x} D_xi
            for (int xi = 0; xi < dim; ++xi) {
                D_phys[x].Add(mesh_data.ElementJinv(element, xi, x), D_ref[xi]);
            }
        }
    }

    void BuildPhysicalDerivatives()
    {
        physical_derivatives.resize(fes->GetNE());
        for (int e = 0; e < fes->GetNE(); ++e) {
            BuildPhysicalDerivativeMatrices(e, physical_derivatives[e]);
        }
    }

    void BuildState(const Vector &x, int element, OFDGDerivativeState &state) const
    {
        x.GetSubVector(mesh_data.Element(element).vdofs, state.root_vector);

        const auto &nodes = operators.DerivativeNodes();
        const auto &D = physical_derivatives[element];

        for (int node = 1; node < static_cast<int>(nodes.size()); ++node) {
            const OFDGDerivativeNode &info = nodes[node];

            Mult(D[info.physical_direction], state.coefficient_views[info.parent],
                 state.coefficient_views[node]);
        }
    }

public:
    OFDGDerivativeProvider(const FiniteElementSpace *fes_, const OFDGOperators &operators_,
                           const OFDGMeshData &mesh_data_, OFDGInternalTiming &timing_)
        : fes(fes_), operators(operators_), mesh_data(mesh_data_), timing(timing_),
          dim(operators_.Dim()), ndof(operators_.NDof()), ncomp(fes_->GetVDim())
    {
        BuildPhysicalDerivatives();

        states.reserve(fes->GetNE());
        for (int e = 0; e < fes->GetNE(); ++e) {
            states.emplace_back(ndof, ncomp, operators.DerivativeCount());
        }

        const std::size_t coefficient_bytes =
            static_cast<std::size_t>(fes->GetNE()) *
            static_cast<std::size_t>(operators.DerivativeCount()) * static_cast<std::size_t>(ndof) *
            static_cast<std::size_t>(ncomp) * sizeof(double);

        const std::size_t operator_bytes =
            static_cast<std::size_t>(fes->GetNE()) * static_cast<std::size_t>(dim) *
            static_cast<std::size_t>(ndof) * static_cast<std::size_t>(ndof) * sizeof(double);

        storage_bytes = coefficient_bytes + operator_bytes;
    }

    void BuildAll(const Vector &x)
    {
#ifdef OFDG_INTERNAL_TIMING
        const auto timer_begin = OFDGTimingClock::now();
#endif
        for (int e = 0; e < fes->GetNE(); ++e) {
            BuildState(x, e, states[e]);
        }
#ifdef OFDG_INTERNAL_TIMING
        timing.derivative_build += OFDGSecondsSince(timer_begin);
#endif
    }

    Pair GetPair(int e1, int e2) const
    {
        return {&states[e1], &states[e2]};
    }

    std::size_t StorageBytes() const
    {
        return storage_bytes;
    }
};

// =============================================================================
// OFDG
//
// This is now primarily the mathematical algorithm:
//
//     scaling
//       -> face derivative jumps
//       -> sigma accumulation
//       -> stabilization / decay
//
// Mesh and face data are precomputed once. Derivatives for every element are
// rebuilt once per solution evaluation and then reused by all adjacent faces.
//
// Preconditions:
//
//   - fes is a nonempty DG space of polynomial order >= 1;
//   - all elements share geometry, order, dimension and dof layout;
//   - every element mapping is affine;
//   - the mesh is not modified after this object is constructed.
// =============================================================================

class OFDG {
private:
    const FiniteElementSpace *fes;

    mutable VectorFunctionCoefficient velocity;

    bool use_unit_wave_speed;

    OFDGOperators operators;

    const int dim;
    const int order;
    const int ndof;
    const int ncomp;

    OFDGMeshData mesh_data;

    mutable OFDGInternalTiming internal_timing;

    OFDGFaceEvaluation face_evaluation;

    mutable OFDGDerivativeProvider derivative_provider;

    mutable OFDGScratch scratch;

    // sigma_elem[c](e,l)
    mutable std::vector<DenseMatrix> sigma_elem;

    // =========================================================================
    // Global mean / scaling
    // =========================================================================

    void ComputeMeansAndScaling(const Vector &x, Vector *means, Vector *scaling) const
    {
        OFDGVolumeScratch &volume = scratch.volume;

        volume.mean_integrals = 0.0;

        const double infinity = std::numeric_limits<double>::infinity();

        for (int c = 0; c < ncomp; ++c) {
            volume.minimum_values(c) = infinity;
            volume.maximum_values(c) = -infinity;
        }

        double total_volume = 0.0;

        const DenseMatrix &E = operators.VolumeEvaluation();

        const Vector &reference_weights = operators.VolumeReferenceWeights();

        const int nq = E.Height();

        for (int e = 0; e < fes->GetNE(); ++e) {
            const auto &element = mesh_data.Element(e);

            x.GetSubVector(element.vdofs, volume.element_data);

            Mult(E, volume.element_matrix, volume.values);

            total_volume += element.volume;

            for (int q = 0; q < nq; ++q) {
                const double physical_weight = reference_weights(q) * element.jacobian_weight;

                for (int c = 0; c < ncomp; ++c) {
                    const double value = volume.values(q, c);

                    volume.mean_integrals(c) += physical_weight * value;

                    volume.minimum_values(c) = std::min(volume.minimum_values(c), value);

                    volume.maximum_values(c) = std::max(volume.maximum_values(c), value);
                }
            }
        }

        Vector local_means(ncomp);

        for (int c = 0; c < ncomp; ++c) {
            local_means(c) = volume.mean_integrals(c) / total_volume;
        }

        if (means) {
            means->SetSize(ncomp);
            *means = local_means;
        }

        if (scaling) {
            scaling->SetSize(ncomp);

            for (int c = 0; c < ncomp; ++c) {
                (*scaling)(c) = std::max(local_means(c) - volume.minimum_values(c),
                                         volume.maximum_values(c) - local_means(c));
            }
        }
    }

public:
    // =========================================================================
    // Constructors
    // =========================================================================

    OFDG(const FiniteElementSpace *fes_, int btype_)
        : OFDG(fes_, btype_,
               VectorFunctionCoefficient(fes_->GetMesh()->Dimension(),
                                         [](const Vector &, Vector &v) { v = 1.0; }))
    {
        use_unit_wave_speed = true;
    }

    OFDG(const FiniteElementSpace *fes_, int btype_, VectorFunctionCoefficient vel_)
        : fes(fes_), velocity(vel_), use_unit_wave_speed(false), operators(fes_, btype_),
          dim(operators.Dim()), order(operators.Order()), ndof(operators.NDof()),
          ncomp(fes_->GetVDim()), mesh_data(fes_, operators), internal_timing(),
          face_evaluation(fes_, operators, internal_timing),
          derivative_provider(fes_, operators, mesh_data, internal_timing)
    {
        sigma_elem.resize(ncomp);

        for (int c = 0; c < ncomp; ++c) {
            sigma_elem[c].SetSize(fes->GetNE(), order + 1);
        }

        scratch.Initialize(ndof, ncomp, operators.VolumeEvaluation().Height(), order);

        const double derivative_storage_mib =
            static_cast<double>(derivative_provider.StorageBytes()) / (1024.0 * 1024.0);

        const double face_storage_kib =
            static_cast<double>(face_evaluation.StorageBytes()) / 1024.0;

        std::cout << "Optimized serial OFDG:\n"
                  << "  elements: " << fes->GetNE() << '\n'
                  << "  interior faces: " << face_evaluation.FaceCount() << '\n'
                  << "  components: " << ncomp << '\n'
                  << "  order: " << order << '\n'
                  << "  dofs/element: " << ndof << '\n'
                  << "  derivative states/element: " << operators.DerivativeCount() << '\n'
                  << "  face evaluation storage: ~" << std::fixed << std::setprecision(2)
                  << face_storage_kib << " KiB" << std::defaultfloat << '\n'
                  << "  stored derivative elements: " << fes->GetNE() << '\n'
                  << "  derivative coefficient/operator storage: ~" << std::fixed
                  << std::setprecision(2) << derivative_storage_mib << " MiB" << std::defaultfloat
                  << '\n';

#ifdef OFDG_INTERNAL_TIMING
        std::cout << "  internal timing: enabled" << std::endl;
#else
        std::cout << "  internal timing: disabled" << std::endl;
#endif
    }

    // =========================================================================
    // Internal timing
    // =========================================================================

    void ResetInternalTimings() const
    {
#ifdef OFDG_INTERNAL_TIMING
        internal_timing = OFDGInternalTiming{};
#endif
    }

    void PrintInternalTimings(std::ostream &out = std::cout) const
    {
#ifdef OFDG_INTERNAL_TIMING
        out << "\nOFDG internal timings\n"
            << "---------------------\n"
            << "ComputeJumps total:        " << internal_timing.compute_jumps << " s\n"
            << "  global scaling:          " << internal_timing.global_scaling << " s\n"
            << "  normal speed:            " << internal_timing.normal_speed << " s\n"
            << "  derivative construction: " << internal_timing.derivative_build << " s\n"
            << "  face trace/jumps:        " << internal_timing.face_trace << " s\n"
            << "  sigma accumulation:      " << internal_timing.sigma_accumulation << " s\n"
            << "Decay application:         " << internal_timing.decay_application << " s\n"
            << "Stabilization application: " << internal_timing.stabilization_application << " s\n"
            << '\n';
#else
        out << "OFDG internal timing is disabled. "
            << "Compile with -DOFDG_INTERNAL_TIMING to enable it." << std::endl;
#endif
    }

    // =========================================================================
    // Existing commutator diagnostic
    // =========================================================================

    void Test()
    {
        operators.TestCommutators();
    }

    // =========================================================================
    // Mean / scaling
    // =========================================================================

    void ComputeMean(const Vector &x, Vector &means) const
    {
        MFEM_VERIFY(x.Size() == fes->GetVSize(),
                    "OFDG input size does not match the finite element space.");

        ComputeMeansAndScaling(x, &means, nullptr);
    }

    void ComputeMean(const Vector &x, double &mean) const
    {
        MFEM_VERIFY(ncomp == 1, "Use ComputeMean(x, Vector&) for multi-component fields.");

        Vector means;

        ComputeMean(x, means);
        mean = means(0);
    }

    void ComputeGlobalMeanScaling(const Vector &x, Vector &scaling) const
    {
        MFEM_VERIFY(x.Size() == fes->GetVSize(),
                    "OFDG input size does not match the finite element space.");

#ifdef OFDG_INTERNAL_TIMING
        const auto timer_begin = OFDGTimingClock::now();
#endif

        ComputeMeansAndScaling(x, nullptr, &scaling);

#ifdef OFDG_INTERNAL_TIMING
        internal_timing.global_scaling += OFDGSecondsSince(timer_begin);
#endif
    }

    double ComputeGlobalMeanScaling(const Vector &x) const
    {
        MFEM_VERIFY(ncomp == 1,
                    "Use ComputeGlobalMeanScaling(x, Vector&) for multi-component fields.");

        Vector scaling;
        ComputeGlobalMeanScaling(x, scaling);

        return scaling(0);
    }

    // =========================================================================
    // Normal wave speed
    // =========================================================================

    void ComputeNormalSpeed(FaceElementTransformations *Tr, double &beta_1, double &beta_2) const
    {
        if (use_unit_wave_speed) {
            beta_1 = 1.0;
            beta_2 = Tr->Elem2 ? 1.0 : 0.0;
            return;
        }

        Vector normal(dim);
        Vector vel1(dim);
        Vector vel2(dim);

        beta_1 = 0.0;
        beta_2 = 0.0;

        const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * order + 1);

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const IntegrationPoint &ip = ir.IntPoint(q);

            Tr->SetAllIntPoints(&ip);

            if (dim == 1) {
                normal(0) = 1.0;
            } else {
                CalcOrtho(Tr->Face->Jacobian(), normal);

                normal /= normal.Norml2();
            }

            const IntegrationPoint &ip1 = Tr->GetElement1IntPoint();

            velocity.Eval(vel1, *Tr->Elem1, ip1);

            double s1 = 0.0;

            for (int d = 0; d < dim; ++d) {
                s1 += normal(d) * vel1(d);
            }

            beta_1 = std::max(beta_1, std::abs(s1));

            if (Tr->Elem2) {
                const IntegrationPoint &ip2 = Tr->GetElement2IntPoint();

                velocity.Eval(vel2, *Tr->Elem2, ip2);

                double s2 = 0.0;

                for (int d = 0; d < dim; ++d) {
                    s2 -= normal(d) * vel2(d);
                }

                beta_2 = std::max(beta_2, std::abs(s2));
            }
        }
    }

    // =========================================================================
    // Public derivative-jump API
    // =========================================================================

    void ComputeDerivativeJumpsAllComponents(const Vector &x, FaceElementTransformations *Tr,
                                             DenseMatrix &jumps) const
    {
        MFEM_VERIFY(x.Size() == fes->GetVSize(),
                    "OFDG input size does not match the finite element space.");

        derivative_provider.BuildAll(x);

        jumps.SetSize(ncomp, order + 1);
        jumps = 0.0;

        if (!Tr || Tr->Elem1No < 0 || Tr->Elem2No < 0) {
            return;
        }

        const auto derivatives = derivative_provider.GetPair(Tr->Elem1No, Tr->Elem2No);

        face_evaluation.EvaluateDerivativeJumps(derivatives.first->all_coefficients,
                                                derivatives.second->all_coefficients, Tr, jumps,
                                                scratch.jump);
    }

    void ComputeDerivativeJumps(const Vector &x, FaceElementTransformations *Tr, int component,
                                Vector &jumps) const
    {
        MFEM_VERIFY(component >= 0 && component < ncomp, "OFDG component index is out of range.");

        DenseMatrix all_jumps;

        ComputeDerivativeJumpsAllComponents(x, Tr, all_jumps);

        jumps.SetSize(order + 1);

        for (int l = 0; l <= order; ++l) {
            jumps(l) = all_jumps(component, l);
        }
    }

    void ComputeDerivativeJumps(const Vector &x, FaceElementTransformations *Tr,
                                Vector &jumps) const
    {
        MFEM_VERIFY(ncomp == 1, "Specify a component for a multi-component field.");

        ComputeDerivativeJumps(x, Tr, 0, jumps);
    }

    // =========================================================================
    // Sensor
    // =========================================================================

    void ComputeJumps(const Vector &x, const Array<bool> *active = nullptr) const
    {
#ifdef OFDG_INTERNAL_TIMING
        const auto total_timer_begin = OFDGTimingClock::now();
#endif

        for (int c = 0; c < ncomp; ++c) {
            sigma_elem[c] = 0.0;
        }

        Vector scaling;
#ifdef OFDG_INTERNAL_TIMING
        const auto scaling_timer_begin = OFDGTimingClock::now();
#endif
        ComputeMeansAndScaling(x, nullptr, &scaling);
#ifdef OFDG_INTERNAL_TIMING
        internal_timing.global_scaling += OFDGSecondsSince(scaling_timer_begin);
#endif

        // Build every element derivative once for the current solution.
        derivative_provider.BuildAll(x);

        DenseMatrix &jumps = scratch.jumps;

        Mesh *mesh = fes->GetMesh();

        for (int f = 0; f < mesh->GetNumFaces(); ++f) {
            FaceElementTransformations *FTr = mesh->GetFaceElementTransformations(f);
            if (!FTr || FTr->Elem1No < 0 || FTr->Elem2No < 0) {
                continue;
            }

            const int e1 = FTr->Elem1No;
            const int e2 = FTr->Elem2No;

            if (active && !(*active)[e1] && !(*active)[e2]) {
                continue;
            }

            double beta_1 = 0.0;
            double beta_2 = 0.0;

#ifdef OFDG_INTERNAL_TIMING
            const auto beta_timer_begin = OFDGTimingClock::now();
#endif

            ComputeNormalSpeed(FTr, beta_1, beta_2);

#ifdef OFDG_INTERNAL_TIMING
            internal_timing.normal_speed += OFDGSecondsSince(beta_timer_begin);
#endif

            const auto derivatives = derivative_provider.GetPair(e1, e2);

            face_evaluation.EvaluateDerivativeJumps(derivatives.first->all_coefficients,
                                                    derivatives.second->all_coefficients, f, jumps,
                                                    scratch.jump);

#ifdef OFDG_INTERNAL_TIMING
            const auto sigma_timer_begin = OFDGTimingClock::now();
#endif

            for (int c = 0; c < ncomp; ++c) {
                if (scaling(c) <= 1e-14) {
                    continue;
                }

                const double inverse_scaling = 1.0 / scaling(c);

                for (int l = 0; l <= order; ++l) {
                    sigma_elem[c](e1, l) +=
                        beta_1 * mesh_data.SensorPrefactor(e1, l) * jumps(c, l) * inverse_scaling;

                    sigma_elem[c](e2, l) +=
                        beta_2 * mesh_data.SensorPrefactor(e2, l) * jumps(c, l) * inverse_scaling;
                }
            }

#ifdef OFDG_INTERNAL_TIMING
            internal_timing.sigma_accumulation += OFDGSecondsSince(sigma_timer_begin);
#endif
        }

#ifdef OFDG_INTERNAL_TIMING
        internal_timing.compute_jumps += OFDGSecondsSince(total_timer_begin);
#endif
    }

    // =========================================================================
    // Stabilization
    //
    // Kept mathematically identical to the current implementation.
    // =========================================================================

    void ComputeStabilization(const Vector &x, Vector &S, const Array<bool> *active = nullptr) const
    {
        S.SetSize(x.Size());
        S = 0.0;

        ComputeJumps(x, active);

#ifdef OFDG_INTERNAL_TIMING
        const auto timer_begin = OFDGTimingClock::now();
#endif

        Vector u_e(ndof);
        Vector tmp(ndof);

        Array<int> vdofs;

        for (int e = 0; e < fes->GetNE(); ++e) {
            if (active && !(*active)[e]) {
                continue;
            }

            for (int c = 0; c < ncomp; ++c) {
                mesh_data.GetComponentVDofs(e, c, vdofs);

                x.GetSubVector(vdofs, u_e);

                for (int l = 0; l <= order; ++l) {
                    const int S_index = (l == 0) ? 0 : l - 1;

                    operators.S(S_index).Mult(u_e, tmp);

                    tmp *= sigma_elem[c](e, l);

                    S.AddElementVector(vdofs, tmp);
                }
            }
        }

#ifdef OFDG_INTERNAL_TIMING
        internal_timing.stabilization_application += OFDGSecondsSince(timer_begin);
#endif
    }

    // =========================================================================
    // Decay
    //
    // Kept mathematically identical to the current implementation.
    // =========================================================================

    void CompDecay(const Vector &x, Vector &decay, double t_decay,
                   const Array<bool> *active = nullptr) const
    {
        decay = x;

        ComputeJumps(x, active);

#ifdef OFDG_INTERNAL_TIMING
        const auto timer_begin = OFDGTimingClock::now();
#endif

        Vector u_e(ndof);
        Vector result(ndof);
        Vector shell(ndof);
        Vector tmp_a(ndof);
        Vector tmp_b(ndof);

        Array<int> vdofs;

        for (int e = 0; e < fes->GetNE(); ++e) {
            if (active && !(*active)[e]) {
                continue;
            }

            for (int c = 0; c < ncomp; ++c) {
                mesh_data.GetComponentVDofs(e, c, vdofs);

                x.GetSubVector(vdofs, u_e);

                result = u_e;

                double cumulative_delta = sigma_elem[c](e, 0);

                for (int j = 1; j <= order; ++j) {
                    cumulative_delta += sigma_elem[c](e, j);

                    const double factor = std::exp(-t_decay * cumulative_delta);

                    if (j < order) {
                        operators.S(j - 1).Mult(u_e, tmp_a);

                        operators.S(j).Mult(u_e, tmp_b);

                        shell = tmp_a;
                        shell -= tmp_b;
                    } else {
                        operators.S(j - 1).Mult(u_e, shell);
                    }

                    result.Add(factor - 1.0, shell);
                }

                decay.SetSubVector(vdofs, result);
            }
        }

#ifdef OFDG_INTERNAL_TIMING
        internal_timing.decay_application += OFDGSecondsSince(timer_begin);
#endif
    }
};
