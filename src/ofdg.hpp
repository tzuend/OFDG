#pragma once

#include "mfem.hpp"
#include "face_physics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

namespace ofdg
{

class OEDG2024;

namespace detail
{

/** Configuration of the face sensor used by the modal decay kernel. */
struct OFDGSensorOptions {
    enum class FaceRule {
        HighOrder,
        Trapezoidal
    };

    FaceRule face_rule = FaceRule::HighOrder;
    bool pool_components = false;
    bool use_face_height = false;
};

// Projection operators
//
// Builds the polynomial projection matrices used by OFDG:
//
//     Q_l = M_k - B_kl^T M_l^{-1} B_kl
//     S_l = M_k^{-1} Q_l
//
// This class contains only projection mathematics.

class OFDGProjection {
private:
    const mfem::FiniteElementSpace *fes;
    const int btype;

    int dim;
    int order;
    mfem::Geometry::Type geom;

public:
    OFDGProjection(const mfem::FiniteElementSpace *fes_, int btype_) : fes(fes_), btype(btype_)
    {
        const mfem::FiniteElement *fe = fes->GetFE(0);

        dim = fe->GetDim();
        order = fe->GetOrder();
        geom = fe->GetGeomType();
    }

    mfem::DenseMatrix AssembleMk(int k, int e = 0) const
    {
        mfem::DG_FECollection fec_k(k, dim, btype);

        const mfem::FiniteElement *fe_k = fec_k.FiniteElementForGeometry(geom);

        mfem::ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);

        const mfem::IntegrationRule &ir = mfem::IntRules.Get(fe_k->GetGeomType(), 2 * fe_k->GetOrder() + 1);

        mfem::DenseMatrix M_k(fe_k->GetDof());
        M_k = 0.0;

        mfem::Vector shape_k(fe_k->GetDof());

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            T->SetIntPoint(&ip);
            fe_k->CalcShape(ip, shape_k);

            const double w = ip.weight * T->Weight();
            mfem::AddMult_a_VVt(w, shape_k, M_k);
        }

        return M_k;
    }

    mfem::DenseMatrix AssembleBkl(int k, int l, int e = 0) const
    {
        mfem::DG_FECollection fec_k(k, dim, btype);
        mfem::DG_FECollection fec_l(l, dim, btype);

        const mfem::FiniteElement *fe_k = fec_k.FiniteElementForGeometry(geom);

        const mfem::FiniteElement *fe_l = fec_l.FiniteElementForGeometry(geom);

        mfem::ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);

        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(fe_k->GetGeomType(), 2 * std::max(fe_k->GetOrder(), fe_l->GetOrder()) + 1);

        mfem::DenseMatrix B(fe_l->GetDof(), fe_k->GetDof());
        B = 0.0;

        mfem::Vector shape_k(fe_k->GetDof());
        mfem::Vector shape_l(fe_l->GetDof());

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            T->SetIntPoint(&ip);

            fe_k->CalcShape(ip, shape_k);
            fe_l->CalcShape(ip, shape_l);

            const double w = ip.weight * T->Weight();
            mfem::AddMult_a_VWt(w, shape_l, shape_k, B);
        }

        return B;
    }

    mfem::DenseMatrix AssembleQl(int l, int e = 0) const
    {
        mfem::DenseMatrix M_k = AssembleMk(order, e);
        mfem::DenseMatrix M_l = AssembleMk(l, e);
        mfem::DenseMatrix B_l = AssembleBkl(order, l, e);

        const int N_k = M_k.Width();
        const int N_l = M_l.Width();

        mfem::DenseMatrix M_l_inv(M_l);
        M_l_inv.Invert();

        mfem::DenseMatrix MlInvB(N_l, N_k);
        mfem::Mult(M_l_inv, B_l, MlInvB);

        mfem::DenseMatrix BtMlInvB(N_k);
        mfem::MultAtB(B_l, MlInvB, BtMlInvB);

        mfem::DenseMatrix Q_l = M_k;
        Q_l -= BtMlInvB;

        return Q_l;
    }
};

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

class OFDGReferenceDerivative {
private:
    const mfem::FiniteElementSpace *fes;

public:
    explicit OFDGReferenceDerivative(const mfem::FiniteElementSpace *fes_) : fes(fes_)
    {
    }

    mfem::DenseMatrix AssembleVandermonde(int e = 0, const mfem::IntegrationRule *ir = nullptr) const
    {
        const mfem::FiniteElement *fe = fes->GetFE(e);
        const int ndof = fe->GetDof();

        if (ir == nullptr) {
            ir = &fe->GetNodes();
        }

        mfem::DenseMatrix V(ndof);
        mfem::Vector shape(ndof);

        for (int q = 0; q < ndof; ++q) {
            const mfem::IntegrationPoint &ip = ir->IntPoint(q);

            fe->CalcShape(ip, shape);

            for (int i = 0; i < ndof; ++i) {
                V(q, i) = shape(i);
            }
        }

        return V;
    }

    std::vector<mfem::DenseMatrix> Assemble(int e = 0, const mfem::IntegrationRule *ir = nullptr) const
    {
        const mfem::FiniteElement *fe = fes->GetFE(e);

        const int ndof = fe->GetDof();
        const int dim = fe->GetDim();

        if (ir == nullptr) {
            ir = &fe->GetNodes();
        }

        std::vector<mfem::DenseMatrix> derivative_values;
        derivative_values.reserve(dim);

        for (int d = 0; d < dim; ++d) {
            derivative_values.emplace_back(ndof);
            derivative_values.back() = 0.0;
        }

        mfem::DenseMatrix dshape(ndof, dim);

        for (int q = 0; q < ndof; ++q) {
            const mfem::IntegrationPoint &ip = ir->IntPoint(q);

            fe->CalcDShape(ip, dshape);

            for (int i = 0; i < ndof; ++i) {
                for (int d = 0; d < dim; ++d) {
                    derivative_values[d](q, i) = dshape(i, d);
                }
            }
        }

        mfem::DenseMatrix V_inv = AssembleVandermonde(e, ir);
        V_inv.Invert();

        std::vector<mfem::DenseMatrix> D_ref;
        D_ref.reserve(dim);

        // D_xi = V^{-1} G_xi
        for (int d = 0; d < dim; ++d) {
            D_ref.emplace_back(ndof);
            mfem::Mult(V_inv, derivative_values[d], D_ref.back());
        }

        return D_ref;
    }
};

/** Oscillation-free DG stabilization and exact modal decay.
 *
 * The finite-element space is a nonempty, uniform-order DG space whose
 * elements share one geometry and DOF signature. Mappings are affine and the
 * mesh topology remains static for the lifetime of the filter.
 */

struct OFDGDerivativeNode {
    int degree = 0;
    int parent = -1;
    int physical_direction = -1;
};

class OFDGOperators {
private:
    const mfem::FiniteElementSpace *fes;

    int dim;
    int order;
    int ndof;
    mfem::Geometry::Type geom;

    OFDGProjection projection;

    std::vector<mfem::DenseMatrix> S_l;
    std::vector<mfem::DenseMatrix> reference_derivatives;
    std::vector<OFDGDerivativeNode> derivative_nodes;

    mfem::Vector sensor_common;

    mfem::DenseMatrix volume_evaluation;
    mfem::Vector volume_reference_weights;
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
        mfem::DenseMatrix M_k = projection.AssembleMk(order, 0);
        mfem::DenseMatrix M_k_inv = M_k;

        M_k_inv.Invert();

        S_l.reserve(order + 1);

        for (int l = 0; l <= order; ++l) {
            mfem::DenseMatrix S(ndof);
            mfem::DenseMatrix Q_l = projection.AssembleQl(l, 0);

            mfem::Mult(M_k_inv, Q_l, S);
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
        const mfem::FiniteElement *fe = fes->GetFE(0);

        const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, 2 * order + 1);

        const int nq = ir.GetNPoints();

        volume_evaluation.SetSize(nq, ndof);
        volume_reference_weights.SetSize(nq);

        mfem::Vector shape(ndof);

        reference_volume = 0.0;

        for (int q = 0; q < nq; ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            fe->CalcShape(ip, shape);

            for (int i = 0; i < ndof; ++i) {
                volume_evaluation(q, i) = shape(i);
            }

            volume_reference_weights(q) = ip.weight;
            reference_volume += ip.weight;
        }
    }

public:
    OFDGOperators(const mfem::FiniteElementSpace *fes_, int btype)
        : fes(fes_), dim(0), order(0), ndof(0), geom(mfem::Geometry::POINT), projection(fes_, btype)
    {
        const mfem::FiniteElement *fe = fes->GetFE(0);

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

    mfem::Geometry::Type Geometry() const
    {
        return geom;
    }

    int DerivativeCount() const
    {
        return static_cast<int>(derivative_nodes.size());
    }

    const mfem::DenseMatrix &S(int l) const
    {
        return S_l[l];
    }

    const std::vector<mfem::DenseMatrix> &ReferenceDerivatives() const
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

    const mfem::DenseMatrix &VolumeEvaluation() const
    {
        return volume_evaluation;
    }

    const mfem::Vector &VolumeReferenceWeights() const
    {
        return volume_reference_weights;
    }

    double ReferenceVolume() const
    {
        return reference_volume;
    }

};

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

class OFDGMeshData {
public:
    struct ElementData {
        mfem::Array<int> vdofs;

        double h = 0.0;
        double jacobian_weight = 0.0;
        double volume = 0.0;

        // dim <= 3
        std::array<double, 9> Jinv{};
    };

private:
    const mfem::FiniteElementSpace *fes;
    const OFDGOperators &operators;

    int dim;
    int order;
    int ndof;
    mfem::Geometry::Type geom;

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

            mfem::ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);
            T->SetIntPoint(&mfem::Geometries.GetCenter(geom));
            const mfem::DenseMatrix &Jinv = T->InverseJacobian();

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
    OFDGMeshData(const mfem::FiniteElementSpace *fes_, const OFDGOperators &operators_)
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

    void GetComponentVDofs(int element, int component, mfem::Array<int> &vdofs) const
    {
        vdofs.SetSize(ndof);

        const mfem::Array<int> &all = elements[element].vdofs;
        const int offset = component * ndof;

        for (int i = 0; i < ndof; ++i) {
            vdofs[i] = all[offset + i];
        }
    }
};

// Runtime scratch storage
//
// Temporary memory is grouped here so the mathematical classes do not carry
// unrelated scratch vectors throughout their implementations.

struct OFDGVolumeScratch {
    mfem::Vector element_data;
    mfem::DenseMatrix element_matrix;

    mfem::DenseMatrix values;

    mfem::Vector mean_integrals;
    mfem::Vector minimum_values;
    mfem::Vector maximum_values;

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
    mfem::DenseMatrix traces1;
    mfem::DenseMatrix traces2;

    void SetSize(int nq, int trace_columns)
    {
        traces1.SetSize(nq, trace_columns);
        traces2.SetSize(nq, trace_columns);
    }
};

struct OFDGScratch {
    OFDGVolumeScratch volume;
    OFDGJumpScratch jump;

    mfem::DenseMatrix jumps;
    mfem::Vector face_state1;
    mfem::Vector face_state2;

    void Initialize(int ndof, int ncomp, int volume_nq, int order)
    {
        volume.Initialize(ndof, ncomp, volume_nq);
        jumps.SetSize(ncomp, order + 1);
        face_state1.SetSize(ncomp);
        face_state2.SetSize(ncomp);
    }
};

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

class OFDGFaceEvaluation {
public:
    struct FaceEvaluationData {
        mfem::DenseMatrix evaluation1;
        mfem::DenseMatrix evaluation2;
        mfem::Vector normalized_weights;
        double height1 = 1.0;
        double height2 = 1.0;
    };

private:
    const mfem::FiniteElementSpace *fes;
    const OFDGOperators &operators;
    int order;
    int ndof;
    int ncomp;
    OFDGSensorOptions::FaceRule face_rule;

    std::vector<FaceEvaluationData> face_data;

    double ElementFaceHeight(mfem::FaceElementTransformations *transformations,
                             bool first) const
    {
        const mfem::IntegrationPoint &face_center_reference =
            mfem::Geometries.GetCenter(transformations->GetGeometryType());
        transformations->Face->SetIntPoint(&face_center_reference);
        mfem::Vector face_center;
        transformations->Face->Transform(face_center_reference, face_center);

        mfem::ElementTransformation *element_transformation =
            first ? transformations->Elem1 : transformations->Elem2;
        const mfem::Geometry::Type geometry = operators.Geometry();
        const mfem::IntegrationRule &vertices = *mfem::Geometries.GetVertices(geometry);
        double height = 0.0;
        mfem::Vector physical_vertex;
        mfem::Vector displacement(operators.Dim());
        mfem::Vector normal(operators.Dim());
        if (operators.Dim() > 1) {
            mfem::CalcOrtho(transformations->Face->Jacobian(), normal);
            normal /= normal.Norml2();
        } else {
            normal = 1.0;
        }
        for (int vertex = 0; vertex < vertices.GetNPoints(); ++vertex) {
            element_transformation->Transform(vertices.IntPoint(vertex),
                                              physical_vertex);
            subtract(physical_vertex, face_center, displacement);
            height = std::max(height, std::abs(displacement * normal));
        }
        return height;
    }

    void BuildFaceEvaluationData(mfem::FaceElementTransformations *Tr, mfem::DenseMatrix &evaluation1,
                                 mfem::DenseMatrix &evaluation2, mfem::Vector &normalized_weights,
                                 double *height1 = nullptr,
                                 double *height2 = nullptr) const
    {
        const mfem::FiniteElement *fe1 = fes->GetFE(Tr->Elem1No);
        const mfem::FiniteElement *fe2 = fes->GetFE(Tr->Elem2No);
        if (height1) { *height1 = ElementFaceHeight(Tr, true); }
        if (height2) { *height2 = ElementFaceHeight(Tr, false); }

        mfem::IntegrationRule trapezoidal_rule;
        const mfem::IntegrationRule *rule = nullptr;
        if (face_rule == OFDGSensorOptions::FaceRule::Trapezoidal &&
            Tr->GetGeometryType() == mfem::Geometry::SEGMENT) {
            trapezoidal_rule.SetSize(2);
            trapezoidal_rule.IntPoint(0).Set1w(0.0, 0.5);
            trapezoidal_rule.IntPoint(1).Set1w(1.0, 0.5);
            rule = &trapezoidal_rule;
        } else {
            // A one-dimensional mesh has point faces, and three-dimensional
            // faces retain a proper face rule until a published 3D OEDG rule
            // is available.
            rule = &mfem::IntRules.Get(Tr->GetGeometryType(), 2 * order + 1);
        }
        const mfem::IntegrationRule &ir = *rule;

        const int nq = ir.GetNPoints();

        evaluation1.SetSize(nq, ndof);
        evaluation2.SetSize(nq, ndof);
        normalized_weights.SetSize(nq);

        mfem::Vector shape1(ndof);
        mfem::Vector shape2(ndof);

        double weight_sum = 0.0;

        for (int q = 0; q < nq; ++q) {
            weight_sum += ir.IntPoint(q).weight;
        }

        for (int q = 0; q < nq; ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            Tr->SetAllIntPoints(&ip);

            const mfem::IntegrationPoint &ip1 = Tr->GetElement1IntPoint();

            const mfem::IntegrationPoint &ip2 = Tr->GetElement2IntPoint();

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
        mfem::Mesh *mesh = fes->GetMesh();
        face_data.clear();
        face_data.resize(mesh->GetNumFaces());

        for (int f = 0; f < mesh->GetNumFaces(); ++f) {
            mfem::FaceElementTransformations *Tr = mesh->GetFaceElementTransformations(f);
            if (!Tr || Tr->Elem1No < 0 || Tr->Elem2No < 0) {
                continue;
            }

            FaceEvaluationData &data = face_data[f];
            BuildFaceEvaluationData(Tr, data.evaluation1, data.evaluation2,
                                    data.normalized_weights, &data.height1,
                                    &data.height2);
        }
    }

    void Evaluate(const mfem::DenseMatrix &derivatives1, const mfem::DenseMatrix &derivatives2,
                  const FaceEvaluationData &data, mfem::DenseMatrix &jumps,
                  OFDGJumpScratch &scratch) const
    {
        const int nq = data.evaluation1.Height();
        const int trace_columns = operators.DerivativeCount() * ncomp;

        jumps.SetSize(ncomp, order + 1);
        jumps = 0.0;

        scratch.SetSize(nq, trace_columns);

        // All derivative states and components are evaluated in one matrix
        // multiplication per side.
        mfem::Mult(data.evaluation1, derivatives1, scratch.traces1);

        mfem::Mult(data.evaluation2, derivatives2, scratch.traces2);

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

    }

public:
    OFDGFaceEvaluation(const mfem::FiniteElementSpace *fes_, const OFDGOperators &operators_,
                       OFDGSensorOptions::FaceRule face_rule_ =
                           OFDGSensorOptions::FaceRule::HighOrder)
        : fes(fes_), operators(operators_), order(operators_.Order()),
          ndof(operators_.NDof()), ncomp(fes_->GetVDim()), face_rule(face_rule_)
    {
        BuildCachedFaces();
    }


    double FaceHeight(int face, bool first) const
    {
        return first ? face_data[face].height1 : face_data[face].height2;
    }

    void EvaluateDerivativeJumps(const mfem::DenseMatrix &derivatives1, const mfem::DenseMatrix &derivatives2,
                                 int face, mfem::DenseMatrix &jumps, OFDGJumpScratch &scratch) const
    {
        Evaluate(derivatives1, derivatives2, face_data[face], jumps, scratch);
    }

    // Used by the public one-face testing API.
    void EvaluateDerivativeJumps(const mfem::DenseMatrix &derivatives1, const mfem::DenseMatrix &derivatives2,
                                 mfem::FaceElementTransformations *Tr, mfem::DenseMatrix &jumps,
                                 OFDGJumpScratch &scratch,
                                 double *height1 = nullptr,
                                 double *height2 = nullptr) const
    {
        FaceEvaluationData data;
        BuildFaceEvaluationData(Tr, data.evaluation1, data.evaluation2,
                                data.normalized_weights, height1, height2);
        Evaluate(derivatives1, derivatives2, data, jumps, scratch);
    }
};

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

struct OFDGDerivativeState {
    mfem::DenseMatrix all_coefficients;
    mfem::Vector root_vector;

    std::vector<mfem::DenseMatrix> coefficient_views;

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

class OFDGDerivativeProvider {
public:
    struct Pair {
        const OFDGDerivativeState *first = nullptr;
        const OFDGDerivativeState *second = nullptr;
    };

private:
    const mfem::FiniteElementSpace *fes;
    const OFDGOperators &operators;
    const OFDGMeshData &mesh_data;
    int dim;
    int ndof;
    int ncomp;

    std::vector<std::vector<mfem::DenseMatrix>> physical_derivatives;
    std::vector<OFDGDerivativeState> states;

    void BuildPhysicalDerivativeMatrices(int element, std::vector<mfem::DenseMatrix> &D_phys) const
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

    void BuildState(const mfem::Vector &x, int element, OFDGDerivativeState &state) const
    {
        x.GetSubVector(mesh_data.Element(element).vdofs, state.root_vector);

        const auto &nodes = operators.DerivativeNodes();
        const auto &D = physical_derivatives[element];

        for (int node = 1; node < static_cast<int>(nodes.size()); ++node) {
            const OFDGDerivativeNode &info = nodes[node];

            mfem::Mult(D[info.physical_direction], state.coefficient_views[info.parent],
                 state.coefficient_views[node]);
        }
    }

    void BuildStateFromLocalData(const mfem::Vector &local_state,
                                 mfem::ElementTransformation &transformation,
                                 OFDGDerivativeState &state) const
    {
        state.root_vector = local_state;

        transformation.SetIntPoint(&mfem::Geometries.GetCenter(operators.Geometry()));
        const mfem::DenseMatrix &inverse_jacobian = transformation.InverseJacobian();
        const auto &reference_derivatives = operators.ReferenceDerivatives();
        std::vector<mfem::DenseMatrix> physical_derivatives(dim);

        for (int direction = 0; direction < dim; ++direction) {
            physical_derivatives[direction].SetSize(ndof);
            physical_derivatives[direction] = 0.0;
            for (int reference_direction = 0;
                 reference_direction < dim; ++reference_direction) {
                physical_derivatives[direction].Add(
                    inverse_jacobian(reference_direction, direction),
                    reference_derivatives[reference_direction]);
            }
        }

        const auto &nodes = operators.DerivativeNodes();
        for (int node = 1; node < static_cast<int>(nodes.size()); ++node) {
            const OFDGDerivativeNode &info = nodes[node];
            mfem::Mult(physical_derivatives[info.physical_direction],
                 state.coefficient_views[info.parent],
                 state.coefficient_views[node]);
        }
    }

public:
    OFDGDerivativeProvider(const mfem::FiniteElementSpace *fes_, const OFDGOperators &operators_,
                           const OFDGMeshData &mesh_data_)
        : fes(fes_), operators(operators_), mesh_data(mesh_data_),
          dim(operators_.Dim()), ndof(operators_.NDof()), ncomp(fes_->GetVDim())
    {
        BuildPhysicalDerivatives();

        states.reserve(fes->GetNE());
        for (int e = 0; e < fes->GetNE(); ++e) {
            states.emplace_back(ndof, ncomp, operators.DerivativeCount());
        }

    }

    void BuildAll(const mfem::Vector &x)
    {
        for (int e = 0; e < fes->GetNE(); ++e) {
            BuildState(x, e, states[e]);
        }
    }

    Pair GetPair(int e1, int e2) const
    {
        return {&states[e1], &states[e2]};
    }

    const OFDGDerivativeState &Get(int element) const
    {
        return states[element];
    }

    void BuildFaceNeighborState(const mfem::Vector &local_state,
                                mfem::ElementTransformation &transformation,
                                OFDGDerivativeState &state) const
    {
        BuildStateFromLocalData(local_state, transformation, state);
    }

};

} // namespace detail

// OFDG
//
// This is now primarily the mathematical algorithm:
//
//     scaling
//       -> face derivative jumps
//       -> sigma accumulation
//       -> stabilization / decay
//
// mfem::Mesh and face data are precomputed once. Derivatives for every element are
// rebuilt once per solution evaluation and then reused by all adjacent faces.
//
// Preconditions:
//
//   - fes is a nonempty DG space of polynomial order >= 1;
//   - all elements share geometry, order, dimension and dof layout;
//   - every element mapping is affine;
//   - the mesh is not modified after this object is constructed.

class OFDG {
private:
    const mfem::FiniteElementSpace *fes;
    std::shared_ptr<const FacePhysics> face_physics;
    detail::OFDGSensorOptions sensor_options;

    detail::OFDGOperators operators;

    const int dim;
    const int order;
    const int ndof;
    const int ncomp;

    detail::OFDGMeshData mesh_data;

    detail::OFDGFaceEvaluation face_evaluation;

    mutable detail::OFDGDerivativeProvider derivative_provider;

    mutable detail::OFDGScratch scratch;

    // sigma_elem[c](e,l)
    mutable std::vector<mfem::DenseMatrix> sigma_elem;

    void ComputeScaling(const mfem::Vector &x, mfem::Vector &scaling) const
    {
        detail::OFDGVolumeScratch &volume = scratch.volume;

        volume.mean_integrals = 0.0;

        const double infinity = std::numeric_limits<double>::infinity();

        for (int c = 0; c < ncomp; ++c) {
            volume.minimum_values(c) = infinity;
            volume.maximum_values(c) = -infinity;
        }

        double total_volume = 0.0;

        const mfem::DenseMatrix &E = operators.VolumeEvaluation();

        const mfem::Vector &reference_weights = operators.VolumeReferenceWeights();

        const int nq = E.Height();

        for (int e = 0; e < fes->GetNE(); ++e) {
            const auto &element = mesh_data.Element(e);

            x.GetSubVector(element.vdofs, volume.element_data);

            mfem::Mult(E, volume.element_matrix, volume.values);

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

#ifdef MFEM_USE_MPI
        if (const auto *parallel_space =
               dynamic_cast<const mfem::ParFiniteElementSpace *>(fes)) {
            MPI_Comm communicator = parallel_space->GetComm();
            MPI_Allreduce(MPI_IN_PLACE, volume.mean_integrals.GetData(), ncomp,
                          mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, volume.minimum_values.GetData(), ncomp,
                          mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_MIN, communicator);
            MPI_Allreduce(MPI_IN_PLACE, volume.maximum_values.GetData(), ncomp,
                          mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_MAX, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &total_volume, 1,
                          mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM, communicator);
        }
#endif

        scaling.SetSize(ncomp);
        for (int c = 0; c < ncomp; ++c) {
            const double mean = volume.mean_integrals(c) / total_volume;
            scaling(c) = std::max(mean - volume.minimum_values(c),
                                  volume.maximum_values(c) - mean);
        }
    }

public:

    OFDG(const mfem::FiniteElementSpace *fes_, int btype_)
        : OFDG(fes_, btype_, std::make_shared<UnitFacePhysics>())
    {
    }

    OFDG(const mfem::FiniteElementSpace *fes_, int btype_, mfem::VectorFunctionCoefficient vel_)
        : OFDG(
              fes_, btype_,
              std::make_shared<AdvectionFacePhysics>(
                  std::static_pointer_cast<mfem::VectorCoefficient>(
                      std::make_shared<mfem::VectorFunctionCoefficient>(std::move(vel_)))))
    {
    }

    OFDG(const mfem::FiniteElementSpace *fes_, int btype_,
         std::shared_ptr<const FacePhysics> face_physics_)
        : OFDG(fes_, btype_, std::move(face_physics_), detail::OFDGSensorOptions{})
    {
    }

private:
    friend class OEDG2024;

    OFDG(const mfem::FiniteElementSpace *fes_, int btype_,
         std::shared_ptr<const FacePhysics> face_physics_,
         detail::OFDGSensorOptions sensor_options_)
        : fes(fes_),
          face_physics(std::move(face_physics_)),
          sensor_options(sensor_options_), operators(fes_, btype_),
          dim(operators.Dim()), order(operators.Order()), ndof(operators.NDof()),
          ncomp(fes_->GetVDim()), mesh_data(fes_, operators),
          face_evaluation(fes_, operators, sensor_options.face_rule),
          derivative_provider(fes_, operators, mesh_data)
    {
        sigma_elem.resize(ncomp);

        for (int c = 0; c < ncomp; ++c) {
            sigma_elem[c].SetSize(fes->GetNE(), order + 1);
        }

        scratch.Initialize(ndof, ncomp, operators.VolumeEvaluation().Height(), order);

    }

private:

    void ComputeNormalSpeed(mfem::FaceElementTransformations *Tr,
                            const mfem::DenseMatrix &traces1,
                            const mfem::DenseMatrix &traces2,
                            double &beta_1, double &beta_2) const
    {

        mfem::Vector normal(dim);
        mfem::Vector &state1 = scratch.face_state1;
        mfem::Vector &state2 = scratch.face_state2;

        beta_1 = 0.0;
        beta_2 = 0.0;

        const mfem::IntegrationRule &ir = mfem::IntRules.Get(Tr->GetGeometryType(), 2 * order + 1);

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            Tr->SetAllIntPoints(&ip);

            if (dim == 1) {
                normal(0) = 1.0;
            } else {
                mfem::CalcOrtho(Tr->Face->Jacobian(), normal);

                normal /= normal.Norml2();
            }

            for (int c = 0; c < ncomp; ++c) {
                state1(c) = traces1(q, c);
                state2(c) = traces2(q, c);
            }

            const FacePhysicsSample physics =
                face_physics->Evaluate(state1, state2, normal, *Tr);
            beta_1 = std::max(beta_1, static_cast<double>(physics.beta1));
            beta_2 = std::max(beta_2, static_cast<double>(physics.beta2));
        }
    }

    void AccumulateLocalFaces(const mfem::Vector &scaling,
                              const mfem::Array<bool> *active) const
    {
        mfem::DenseMatrix &jumps = scratch.jumps;
        mfem::Mesh *mesh = fes->GetMesh();

        for (int face = 0; face < mesh->GetNumFaces(); ++face) {
            mfem::FaceElementTransformations *transformations =
                mesh->GetFaceElementTransformations(face);
            if (!transformations || transformations->Elem1No < 0 ||
                transformations->Elem2No < 0) {
                continue;
            }

            const int element1 = transformations->Elem1No;
            const int element2 = transformations->Elem2No;
            if (active && !(*active)[element1] && !(*active)[element2]) {
                continue;
            }

            const auto derivatives =
                derivative_provider.GetPair(element1, element2);
            face_evaluation.EvaluateDerivativeJumps(
                derivatives.first->all_coefficients,
                derivatives.second->all_coefficients, face, jumps,
                scratch.jump);

            double speed1 = 0.0;
            double speed2 = 0.0;
            ComputeNormalSpeed(transformations, scratch.jump.traces1,
                               scratch.jump.traces2, speed1, speed2);

            for (int c = 0; c < ncomp; ++c) {
                if (scaling(c) <= 1e-14) { continue; }
                const double inverse_scaling = 1.0 / scaling(c);

                for (int l = 0; l <= order; ++l) {
                    const double prefactor1 = sensor_options.use_face_height
                        ? operators.SensorCommon(l) *
                              std::pow(face_evaluation.FaceHeight(face, true),
                                       l - 1)
                        : mesh_data.SensorPrefactor(element1, l);
                    const double prefactor2 = sensor_options.use_face_height
                        ? operators.SensorCommon(l) *
                              std::pow(face_evaluation.FaceHeight(face, false),
                                       l - 1)
                        : mesh_data.SensorPrefactor(element2, l);
                    sigma_elem[c](element1, l) +=
                        speed1 * prefactor1 * jumps(c, l) * inverse_scaling;
                    sigma_elem[c](element2, l) +=
                        speed2 * prefactor2 * jumps(c, l) * inverse_scaling;
                }
            }
        }
    }

    void AccumulateSharedFaces(const mfem::Vector &state,
                               const mfem::Vector &scaling,
                               const mfem::Array<bool> *active) const
    {
#ifdef MFEM_USE_MPI
        auto *parallel_space = dynamic_cast<mfem::ParFiniteElementSpace *>(
            const_cast<mfem::FiniteElementSpace *>(fes));
        if (!parallel_space) { return; }

        mfem::ParMesh *parallel_mesh = parallel_space->GetParMesh();
        mfem::ParGridFunction parallel_state(parallel_space);
        parallel_state = state;
        parallel_state.ExchangeFaceNbrData();

        mfem::Array<int> neighbor_vdofs;
        mfem::Vector neighbor_state;
        detail::OFDGDerivativeState neighbor_derivatives(
            ndof, ncomp, operators.DerivativeCount());
        mfem::DenseMatrix &jumps = scratch.jumps;

        // Shared faces are accumulated only for the element owned by this rank.
        for (int face = 0; face < parallel_mesh->GetNSharedFaces(); ++face) {
            mfem::FaceElementTransformations *transformations =
                parallel_mesh->GetSharedFaceTransformations(face, true);
            const int local_element = transformations->Elem1No;
            const int neighbor_element =
                transformations->Elem2No - parallel_mesh->GetNE();
            if (active && !(*active)[local_element]) { continue; }

            parallel_space->GetFaceNbrElementVDofs(neighbor_element,
                                                    neighbor_vdofs);
            neighbor_state.SetSize(neighbor_vdofs.Size());
            parallel_state.FaceNbrData().GetSubVector(neighbor_vdofs,
                                                       neighbor_state);
            derivative_provider.BuildFaceNeighborState(
                neighbor_state, *transformations->Elem2, neighbor_derivatives);

            const detail::OFDGDerivativeState &local_derivatives =
                derivative_provider.Get(local_element);
            double local_height = 1.0;
            double neighbor_height = 1.0;
            face_evaluation.EvaluateDerivativeJumps(
                local_derivatives.all_coefficients,
                neighbor_derivatives.all_coefficients, transformations, jumps,
                scratch.jump, &local_height, &neighbor_height);

            double local_speed = 0.0;
            double neighbor_speed = 0.0;
            ComputeNormalSpeed(transformations, scratch.jump.traces1,
                               scratch.jump.traces2, local_speed,
                               neighbor_speed);

            for (int c = 0; c < ncomp; ++c) {
                if (scaling(c) <= 1e-14) { continue; }
                const double inverse_scaling = 1.0 / scaling(c);

                for (int l = 0; l <= order; ++l) {
                    const double prefactor = sensor_options.use_face_height
                        ? operators.SensorCommon(l) *
                              std::pow(local_height, l - 1)
                        : mesh_data.SensorPrefactor(local_element, l);
                    sigma_elem[c](local_element, l) +=
                        local_speed * prefactor * jumps(c, l) * inverse_scaling;
                }
            }
        }
#else
        (void) state;
        (void) scaling;
        (void) active;
#endif
    }

    void PoolComponentSensors() const
    {
        if (!sensor_options.pool_components || ncomp == 1) { return; }

        for (int e = 0; e < fes->GetNE(); ++e) {
            for (int l = 0; l <= order; ++l) {
                double pooled = 0.0;
                for (int c = 0; c < ncomp; ++c) {
                    pooled = std::max(pooled, sigma_elem[c](e, l));
                }
                for (int c = 0; c < ncomp; ++c) {
                    sigma_elem[c](e, l) = pooled;
                }
            }
        }
    }

    void ComputeJumps(const mfem::Vector &state,
                      const mfem::Array<bool> *active = nullptr) const
    {
        for (mfem::DenseMatrix &sigma : sigma_elem) { sigma = 0.0; }

        mfem::Vector scaling;
        ComputeScaling(state, scaling);
        derivative_provider.BuildAll(state);
        AccumulateLocalFaces(scaling, active);
        AccumulateSharedFaces(state, scaling, active);
        PoolComponentSensors();
    }

public:

    void ComputeStabilization(const mfem::Vector &x, mfem::Vector &S, const mfem::Array<bool> *active = nullptr) const
    {
        S.SetSize(x.Size());
        S = 0.0;

        ComputeJumps(x, active);

        mfem::Vector u_e(ndof);
        mfem::Vector tmp(ndof);

        mfem::Array<int> vdofs;

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

    }

    void CompDecay(const mfem::Vector &x, mfem::Vector &decay, double t_decay,
                   const mfem::Array<bool> *active = nullptr) const
    {
        decay = x;

        ComputeJumps(x, active);

        mfem::Vector u_e(ndof);
        mfem::Vector result(ndof);
        mfem::Vector shell(ndof);
        mfem::Vector tmp_a(ndof);
        mfem::Vector tmp_b(ndof);

        mfem::Array<int> vdofs;

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

    }
};

} // namespace ofdg
