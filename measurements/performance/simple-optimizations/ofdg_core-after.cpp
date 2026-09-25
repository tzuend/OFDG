#include "ofdg_core.hpp"
#include "curved_geometry.hpp"

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

namespace detail
{

// Polynomial-shell projection:
//
//     Q_l = M_k - B_kl^T M_l^{-1} B_kl
//     S_l = M_k^{-1} Q_l

class OFDGProjection {
private:
    const mfem::FiniteElement *finite_element;
    const int btype;
    mfem::ElementTransformation *transformation;
    int quadrature_order;

    int dim;
    int order;
    mfem::Geometry::Type geom;

public:
    OFDGProjection(const mfem::FiniteElement &finite_element_, int btype_,
                   mfem::ElementTransformation *transformation_ = nullptr)
        : finite_element(&finite_element_), btype(btype_), transformation(transformation_),
          quadrature_order(transformation_ ? CurvedQuadratureOrder(finite_element_.GetOrder(), *transformation_) : 2 * finite_element_.GetOrder() + 1)
    {
        dim = finite_element->GetDim();
        order = finite_element->GetOrder();
        geom = finite_element->GetGeomType();
    }

    mfem::DenseMatrix AssembleMk(int k) const
    {
        mfem::DG_FECollection fec_k(k, dim, btype);

        const mfem::FiniteElement *fe_k =
            k == order ? finite_element : fec_k.FiniteElementForGeometry(geom);

        const mfem::IntegrationRule &ir = transformation ?
            CurvedRule(geom, quadrature_order) : mfem::IntRules.Get(geom, quadrature_order);

        mfem::DenseMatrix M_k(fe_k->GetDof());
        M_k = 0.0;

        mfem::Vector shape_k(fe_k->GetDof());

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            fe_k->CalcShape(ip, shape_k);

            if (transformation) { transformation->SetIntPoint(&ip); }
            const double weight = ip.weight * (transformation ? transformation->Weight() : 1.0);
            mfem::AddMult_a_VVt(weight, shape_k, M_k);
        }

        return M_k;
    }

    mfem::DenseMatrix AssembleBkl(int k, int l) const
    {
        mfem::DG_FECollection fec_k(k, dim, btype);
        mfem::DG_FECollection fec_l(l, dim, btype);

        const mfem::FiniteElement *fe_k =
            k == order ? finite_element : fec_k.FiniteElementForGeometry(geom);

        const mfem::FiniteElement *fe_l = fec_l.FiniteElementForGeometry(geom);

        const mfem::IntegrationRule &ir =
            (transformation ? CurvedRule(fe_k->GetGeomType(), quadrature_order) : mfem::IntRules.Get(fe_k->GetGeomType(), quadrature_order));

        mfem::DenseMatrix B(fe_l->GetDof(), fe_k->GetDof());
        B = 0.0;

        mfem::Vector shape_k(fe_k->GetDof());
        mfem::Vector shape_l(fe_l->GetDof());

        for (int q = 0; q < ir.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            fe_k->CalcShape(ip, shape_k);
            fe_l->CalcShape(ip, shape_l);

            if (transformation) { transformation->SetIntPoint(&ip); }
            const double weight = ip.weight * (transformation ? transformation->Weight() : 1.0);
            mfem::AddMult_a_VWt(weight, shape_l, shape_k, B);
        }

        return B;
    }

    mfem::DenseMatrix AssembleQl(int l) const
    {
        mfem::DenseMatrix M_k = AssembleMk(order);
        mfem::DenseMatrix M_l = AssembleMk(l);
        mfem::DenseMatrix B_l = AssembleBkl(order, l);

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

// At interpolation points xi_q,
//
//     V(q,i)   = phi_i(xi_q)
//     G_d(q,i) = d phi_i / d xi_d (xi_q)
// so D_d = V^{-1} G_d maps coefficients to derivative coefficients.

class OFDGReferenceDerivative {
private:
    const mfem::FiniteElement *finite_element;

public:
    explicit OFDGReferenceDerivative(const mfem::FiniteElement &finite_element_)
        : finite_element(&finite_element_)
    {
    }

    mfem::DenseMatrix AssembleVandermonde(const mfem::IntegrationRule *ir = nullptr) const
    {
        const mfem::FiniteElement *fe = finite_element;
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

    std::vector<mfem::DenseMatrix> Assemble(const mfem::IntegrationRule *ir = nullptr) const
    {
        const mfem::FiniteElement *fe = finite_element;

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

        mfem::DenseMatrix V_inv = AssembleVandermonde(ir);
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

/** Shared reference operators and element-specific curved operators. */

struct OFDGDerivativeNode {
    int degree = 0;
    int parent = -1;
    int physical_direction = -1;
};

class OFDGOperators {
private:
    const mfem::FiniteElement *finite_element;

    int dim;
    int order;
    int ndof;
    mfem::Geometry::Type geom;

    const OFDGOperators *reference;

    std::vector<mfem::DenseMatrix> S_l;
    std::vector<mfem::DenseMatrix> reference_derivatives;
    std::vector<OFDGDerivativeNode> derivative_nodes;

    std::vector<mfem::DenseMatrix> projected_derivatives;
    int quadrature_order;
    mfem::Vector sensor_common;

    mfem::DenseMatrix volume_evaluation;
    mfem::Vector volume_weights;
    double integration_volume = 0.0;

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

    void BuildProjectionOperators(const OFDGProjection &projection)
    {
        mfem::DenseMatrix M_k = projection.AssembleMk(order);
        mfem::DenseMatrix M_k_inv = M_k;

        M_k_inv.Invert();

        S_l.reserve(order + 1);

        for (int l = 0; l <= order; ++l) {
            mfem::DenseMatrix S(ndof);
            mfem::DenseMatrix Q_l = projection.AssembleQl(l);

            mfem::Mult(M_k_inv, Q_l, S);
            S_l.push_back(S);
        }
    }

    void BuildReferenceDerivatives()
    {
        OFDGReferenceDerivative derivative_builder(*finite_element);
        reference_derivatives = derivative_builder.Assemble();
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

    void BuildVolumeEvaluation(mfem::ElementTransformation *transformation)
    {
        const mfem::FiniteElement *fe = finite_element;

        const mfem::IntegrationRule &ir = transformation ?
            CurvedRule(geom, quadrature_order) : mfem::IntRules.Get(geom, quadrature_order);

        const int nq = ir.GetNPoints();

        volume_evaluation.SetSize(nq, ndof);
        volume_weights.SetSize(nq);

        mfem::Vector shape(ndof);

        integration_volume = 0.0;

        for (int q = 0; q < nq; ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            fe->CalcShape(ip, shape);

            for (int i = 0; i < ndof; ++i) {
                volume_evaluation(q, i) = shape(i);
            }

            if (transformation) { transformation->SetIntPoint(&ip); }
            const double weight = ip.weight * (transformation ? transformation->Weight() : 1.0);
            volume_weights(q) = weight;
            integration_volume += weight;
        }
    }

public:
    OFDGOperators(const mfem::FiniteElement &finite_element_, int btype,
                  mfem::ElementTransformation *transformation = nullptr,
                  const OFDGOperators *reference_ = nullptr)
        : finite_element(&finite_element_), dim(finite_element_.GetDim()),
          order(finite_element_.GetOrder()), ndof(finite_element_.GetDof()),
          geom(finite_element_.GetGeomType()), reference(reference_),
          quadrature_order(transformation ? CurvedQuadratureOrder(order, *transformation) : 2 * order + 1)
    {
        OFDGProjection projection(finite_element_, btype, transformation);
        BuildProjectionOperators(projection);
        if (!reference) {
            BuildReferenceDerivatives();
            BuildDerivativeDAG();
            BuildSensorConstants();
        }
        BuildVolumeEvaluation(transformation);
        if (transformation) { projected_derivatives = ProjectedPhysicalDerivatives(
            *finite_element, *transformation, quadrature_order); }
    }

    bool Curved() const { return !projected_derivatives.empty(); }
    int QuadratureOrder() const { return quadrature_order; }
    const std::vector<mfem::DenseMatrix> &ProjectedDerivatives() const
    { return projected_derivatives; }

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
        return static_cast<int>(DerivativeNodes().size());
    }

    const mfem::DenseMatrix &S(int l) const
    {
        return S_l[l];
    }

    const std::vector<mfem::DenseMatrix> &ReferenceDerivatives() const
    {
        return reference ? reference->ReferenceDerivatives() : reference_derivatives;
    }

    const std::vector<OFDGDerivativeNode> &DerivativeNodes() const
    {
        return reference ? reference->DerivativeNodes() : derivative_nodes;
    }

    double SensorCommon(int l) const
    {
        return reference ? reference->SensorCommon(l) : sensor_common(l);
    }

    const mfem::DenseMatrix &VolumeEvaluation() const
    {
        return volume_evaluation;
    }

    const mfem::Vector &VolumeWeights() const
    {
        return volume_weights;
    }

    double IntegrationVolume() const
    {
        return integration_volume;
    }

};

struct OFDGElementSignature
{
    int dimension = 0;
    int geometry = 0;
    int order = 0;
    int dofs = 0;
    int map_type = 0;

    bool operator==(const OFDGElementSignature &other) const
    {
        return dimension == other.dimension && geometry == other.geometry &&
               order == other.order && dofs == other.dofs &&
               map_type == other.map_type;
    }
};

struct OFDGElementSignatureHash
{
    std::size_t operator()(const OFDGElementSignature &signature) const
    {
        std::size_t value = static_cast<std::size_t>(signature.dimension);
        value = value * 11 + static_cast<std::size_t>(signature.geometry);
        value = value * 37 + static_cast<std::size_t>(signature.order);
        value = value * 101 + static_cast<std::size_t>(signature.dofs);
        return value * 7 + static_cast<std::size_t>(signature.map_type);
    }
};

class OFDGOperatorRepository
{
private:
    int basis_type;
    std::unordered_map<int, std::unique_ptr<OFDGOperators>> curved_entries;
    std::unordered_map<OFDGElementSignature,
                       std::unique_ptr<OFDGOperators>,
                       OFDGElementSignatureHash> entries;

    static OFDGElementSignature Signature(
        const mfem::FiniteElement &finite_element)
    {
        return {finite_element.GetDim(),
                static_cast<int>(finite_element.GetGeomType()),
                finite_element.GetOrder(), finite_element.GetDof(),
                static_cast<int>(finite_element.GetMapType())};
    }

public:
    explicit OFDGOperatorRepository(int basis_type_)
        : basis_type(basis_type_)
    {
    }

    const OFDGOperators &Get(const mfem::FiniteElement &finite_element,
                             mfem::ElementTransformation *transformation = nullptr,
                             int element_id = -1)
    {
        if (transformation && !IsAffine(*transformation)) {
            auto &entry = curved_entries[element_id];
            if (!entry) {
                const auto &reference = Get(finite_element);
                entry = std::make_unique<OFDGOperators>(
                    finite_element, basis_type, transformation, &reference);
            }
            return *entry;
        }
        const OFDGElementSignature signature = Signature(finite_element);

        auto entry = entries.find(signature);
        if (entry == entries.end()) {
            entry = entries.emplace(
                signature,
                std::make_unique<OFDGOperators>(finite_element, basis_type)).first;
        }
        return *entry->second;
    }

    const OFDGOperators &Find(const mfem::FiniteElement &finite_element, int element_id) const
    {
        auto curved = curved_entries.find(element_id);
        if (curved != curved_entries.end()) { return *curved->second; }
        return *entries.find(Signature(finite_element))->second;
    }
};

// Geometry is immutable after construction. Affine maps use constant metrics;
// curved maps store physically weighted integration data in their operators.

class OFDGMeshData {
public:
    struct ElementData {
        const OFDGOperators *operators = nullptr;
        mfem::Array<int> vdofs;
        std::vector<mfem::Array<int>> component_vdofs;

        double h = 0.0;
        double jacobian_weight = 0.0;
        double volume = 0.0;

        std::array<double, 9> Jinv{};
        std::vector<double> sensor_prefactors;
    };

private:
    const mfem::FiniteElementSpace *fes;
    OFDGOperatorRepository &repository;
    int ncomp;
    int max_ndof = 0;
    int max_volume_points = 0;

    std::vector<ElementData> elements;

    void BuildElements()
    {
        const int ne = fes->GetNE();

        elements.resize(ne);

        for (int e = 0; e < ne; ++e) {
            ElementData &element = elements[e];
            const mfem::FiniteElement &finite_element = *fes->GetFE(e);
            element.operators = &repository.Get(finite_element, fes->GetMesh()->GetElementTransformation(e), e);
            const OFDGOperators &operators = *element.operators;
            const int dim = operators.Dim();
            const int order = operators.Order();
            const int ndof = operators.NDof();

            fes->GetElementVDofs(e, element.vdofs);
            element.component_vdofs.resize(ncomp);
            mfem::Array<int> scalar_dofs;
            fes->GetElementDofs(e, scalar_dofs);
            for (int c = 0; c < ncomp; ++c) {
                element.component_vdofs[c].SetSize(ndof);
                for (int i = 0; i < ndof; ++i) {
                    element.component_vdofs[c][i] =
                        fes->DofToVDof(scalar_dofs[i], c);
                }
            }

            element.h = fes->GetMesh()->GetElementSize(e);

            mfem::ElementTransformation *T = fes->GetMesh()->GetElementTransformation(e);
            T->SetIntPoint(&mfem::Geometries.GetCenter(operators.Geometry()));
            const mfem::DenseMatrix &Jinv = T->InverseJacobian();

            for (int r = 0; r < dim; ++r) {
                for (int c = 0; c < dim; ++c) {
                    element.Jinv[r * dim + c] = Jinv(r, c);
                }
            }

            element.jacobian_weight = operators.Curved() ? 1.0 : T->Weight();

            element.volume = operators.IntegrationVolume() * element.jacobian_weight;
            if (operators.Curved()) {
                element.h = std::pow(element.volume / mfem::Geometry::Volume[operators.Geometry()], 1.0 / dim);
            }

            element.sensor_prefactors.resize(order + 1);
            double h_power = 1.0 / element.h;

            for (int l = 0; l <= order; ++l) {
                element.sensor_prefactors[l] = operators.SensorCommon(l) * h_power;

                h_power *= element.h;
            }

            max_ndof = std::max(max_ndof, ndof);
            max_volume_points =
                std::max(max_volume_points, operators.VolumeEvaluation().Height());
        }
    }

public:
    OFDGMeshData(const mfem::FiniteElementSpace *fes_,
                 OFDGOperatorRepository &repository_)
        : fes(fes_), repository(repository_), ncomp(fes_->GetVDim())
    {
        BuildElements();
    }

    const ElementData &Element(int e) const
    {
        return elements[e];
    }

    double ElementJinv(int e, int row, int col) const
    {
        const int dim = elements[e].operators->Dim();
        return elements[e].Jinv[row * dim + col];
    }

    double SensorPrefactor(int e, int l) const
    {
        return elements[e].sensor_prefactors[l];
    }

    void GetComponentVDofs(int element, int component, mfem::Array<int> &vdofs) const
    {
        vdofs = elements[element].component_vdofs[component];
    }

    int MaxNDof() const
    {
        return max_ndof;
    }

    int MaxVolumePoints() const
    {
        return max_volume_points;
    }

    void GatherElement(const mfem::Vector &state, int element,
                       mfem::Vector &storage,
                       mfem::DenseMatrix &components) const
    {
        const ElementData &data = elements[element];
        const int ndof = data.operators->NDof();
        storage.SetSize(ndof * ncomp);
        components.UseExternalData(storage.GetData(), ndof, ncomp);
        mfem::Vector component;
        for (int c = 0; c < ncomp; ++c) {
            component.SetDataAndSize(storage.GetData() + c * ndof, ndof);
            state.GetSubVector(data.component_vdofs[c], component);
        }
    }
};

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

// Every side keeps its true local size; padding never enters trace evaluation.

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
    const OFDGMeshData &mesh_data;
    int ncomp;
    OFDGSensorOptions::FaceRule face_rule;

    std::vector<FaceEvaluationData> face_data;

    double ElementFaceHeight(mfem::FaceElementTransformations *transformations,
                             bool first,
                             const OFDGOperators &operators) const
    {
        if (operators.Curved()) {
            const auto &rule = CurvedRule(transformations->GetGeometryType(),
                FaceQuadratureOrder(operators.Order(), *transformations));
            double area = 0.0;
            for (int q = 0; q < rule.GetNPoints(); ++q) {
                transformations->Face->SetIntPoint(&rule.IntPoint(q));
                area += rule.IntPoint(q).weight * transformations->Face->Weight();
            }
            return operators.IntegrationVolume() / area *
                mfem::Geometry::Volume[transformations->GetGeometryType()] /
                mfem::Geometry::Volume[operators.Geometry()];
        }
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

    void BuildFaceEvaluationData(mfem::FaceElementTransformations *Tr,
                                 const mfem::FiniteElement &fe1,
                                 const mfem::FiniteElement &fe2,
                                 const OFDGOperators &operators1,
                                 const OFDGOperators &operators2,
                                 mfem::DenseMatrix &evaluation1,
                                 mfem::DenseMatrix &evaluation2, mfem::Vector &normalized_weights,
                                 double *height1 = nullptr,
                                 double *height2 = nullptr) const
    {
        if (height1) { *height1 = ElementFaceHeight(Tr, true, operators1); }
        if (height2) { *height2 = ElementFaceHeight(Tr, false, operators2); }

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
            const int order = std::max(operators1.Order(), operators2.Order());
            rule = &detail::FaceRule(order, *Tr);
        }
        const mfem::IntegrationRule &ir = *rule;

        const int nq = ir.GetNPoints();

        evaluation1.SetSize(nq, operators1.NDof());
        evaluation2.SetSize(nq, operators2.NDof());
        normalized_weights.SetSize(nq);

        mfem::Vector shape1(operators1.NDof());
        mfem::Vector shape2(operators2.NDof());

        normalized_weights = NormalizedPhysicalFaceWeights(*Tr, ir);

        for (int q = 0; q < nq; ++q) {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);

            Tr->SetAllIntPoints(&ip);

            const mfem::IntegrationPoint &ip1 = Tr->GetElement1IntPoint();

            const mfem::IntegrationPoint &ip2 = Tr->GetElement2IntPoint();

            fe1.CalcShape(ip1, shape1);
            fe2.CalcShape(ip2, shape2);

            for (int i = 0; i < operators1.NDof(); ++i) {
                evaluation1(q, i) = shape1(i);
            }
            for (int i = 0; i < operators2.NDof(); ++i) {
                evaluation2(q, i) = shape2(i);
            }

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
            const auto &element1 = mesh_data.Element(Tr->Elem1No);
            const auto &element2 = mesh_data.Element(Tr->Elem2No);
            BuildFaceEvaluationData(Tr, *fes->GetFE(Tr->Elem1No),
                                    *fes->GetFE(Tr->Elem2No),
                                    *element1.operators, *element2.operators,
                                    data.evaluation1, data.evaluation2,
                                    data.normalized_weights, &data.height1,
                                    &data.height2);
        }
    }

    void Evaluate(const mfem::DenseMatrix &derivatives1, const mfem::DenseMatrix &derivatives2,
                  const FaceEvaluationData &data, mfem::DenseMatrix &jumps,
                  OFDGJumpScratch &scratch) const
    {
        const int nq = data.evaluation1.Height();
        const int trace_columns =
            static_cast<int>(derivatives1.Width());

        const int order = mesh_data.Element(0).operators->Order();
        jumps.SetSize(ncomp, order + 1);
        jumps = 0.0;

        scratch.SetSize(nq, trace_columns);

        // All derivative states and components are evaluated in one matrix
        // multiplication per side.
        mfem::Mult(data.evaluation1, derivatives1, scratch.traces1);

        mfem::Mult(data.evaluation2, derivatives2, scratch.traces2);

        const auto &nodes = mesh_data.Element(0).operators->DerivativeNodes();

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
    OFDGFaceEvaluation(const mfem::FiniteElementSpace *fes_,
                       const OFDGMeshData &mesh_data_,
                       OFDGSensorOptions::FaceRule face_rule_ =
                           OFDGSensorOptions::FaceRule::HighOrder)
        : fes(fes_), mesh_data(mesh_data_), ncomp(fes_->GetVDim()),
          face_rule(face_rule_)
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
                                 mfem::FaceElementTransformations *Tr,
                                 const mfem::FiniteElement &finite_element1,
                                 const mfem::FiniteElement &finite_element2,
                                 const OFDGOperators &operators1,
                                 const OFDGOperators &operators2,
                                 mfem::DenseMatrix &jumps,
                                 OFDGJumpScratch &scratch,
                                 double *height1 = nullptr,
                                 double *height2 = nullptr) const
    {
        FaceEvaluationData data;
        BuildFaceEvaluationData(Tr, finite_element1, finite_element2,
                                operators1, operators2,
                                data.evaluation1, data.evaluation2,
                                data.normalized_weights, height1, height2);
        Evaluate(derivatives1, derivatives2, data, jumps, scratch);
    }
};

// Columns are grouped by derivative-DAG node and component. The views alias
// all_coefficients and do not own storage.

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

// Affine derivatives use the constant inverse Jacobian; curved derivatives
// use cached L2 projections. Each state is rebuilt once per solution.

class OFDGDerivativeProvider {
public:
    struct Pair {
        const OFDGDerivativeState *first = nullptr;
        const OFDGDerivativeState *second = nullptr;
    };

private:
    const mfem::FiniteElementSpace *fes;
    const OFDGMeshData &mesh_data;
    int ncomp;

    std::vector<std::vector<mfem::DenseMatrix>> physical_derivatives;
    std::vector<OFDGDerivativeState> states;

    void BuildPhysicalDerivativeMatrices(int element, std::vector<mfem::DenseMatrix> &D_phys) const
    {
        const OFDGOperators &operators =
            *mesh_data.Element(element).operators;
        if (operators.Curved()) { return; }
        const int dim = operators.Dim();
        const int ndof = operators.NDof();
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
        const auto &data = mesh_data.Element(element);
        const OFDGOperators &operators = *data.operators;
        const int ndof = operators.NDof();
        mfem::Vector component;
        for (int c = 0; c < ncomp; ++c) {
            component.SetDataAndSize(state.root_vector.GetData() + c * ndof,
                                     ndof);
            x.GetSubVector(data.component_vdofs[c], component);
        }

        const auto &nodes = operators.DerivativeNodes();
        const auto &D = operators.Curved() ? operators.ProjectedDerivatives() :
                        physical_derivatives[element];

        for (int node = 1; node < static_cast<int>(nodes.size()); ++node) {
            const OFDGDerivativeNode &info = nodes[node];

            mfem::Mult(D[info.physical_direction], state.coefficient_views[info.parent],
                 state.coefficient_views[node]);
        }
    }

    void BuildStateFromLocalData(const mfem::Vector &local_state,
                                 mfem::ElementTransformation &transformation,
                                 const OFDGOperators &operators,
                                 OFDGDerivativeState &state) const
    {
        state.root_vector = local_state;

        std::vector<mfem::DenseMatrix> affine_derivatives;
        if (!operators.Curved()) {
            transformation.SetIntPoint(&mfem::Geometries.GetCenter(operators.Geometry()));
            const mfem::DenseMatrix &inverse_jacobian = transformation.InverseJacobian();
            const auto &reference_derivatives = operators.ReferenceDerivatives();
            const int dim = operators.Dim();
            affine_derivatives.resize(dim);
            for (int direction = 0; direction < dim; ++direction) {
                affine_derivatives[direction].SetSize(operators.NDof());
                affine_derivatives[direction] = 0.0;
                for (int reference_direction = 0; reference_direction < dim; ++reference_direction) {
                    affine_derivatives[direction].Add(
                        inverse_jacobian(reference_direction, direction),
                        reference_derivatives[reference_direction]);
                }
            }
        }
        const auto &derivatives = operators.Curved() ? operators.ProjectedDerivatives() :
                                  affine_derivatives;
        const auto &nodes = operators.DerivativeNodes();
        for (int node = 1; node < static_cast<int>(nodes.size()); ++node) {
            const OFDGDerivativeNode &info = nodes[node];
            mfem::Mult(derivatives[info.physical_direction],
                 state.coefficient_views[info.parent],
                 state.coefficient_views[node]);
        }
    }

public:
    OFDGDerivativeProvider(const mfem::FiniteElementSpace *fes_,
                           const OFDGMeshData &mesh_data_)
        : fes(fes_), mesh_data(mesh_data_), ncomp(fes_->GetVDim())
    {
        BuildPhysicalDerivatives();

        states.reserve(fes->GetNE());
        for (int e = 0; e < fes->GetNE(); ++e) {
            const OFDGOperators &operators = *mesh_data.Element(e).operators;
            states.emplace_back(operators.NDof(), ncomp,
                                operators.DerivativeCount());
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
                                const OFDGOperators &operators,
                                OFDGDerivativeState &state) const
    {
        BuildStateFromLocalData(local_state, transformation, operators, state);
    }

};

} // namespace detail

class StabilizationCoreImplementation {
private:
    const mfem::FiniteElementSpace *fes;
    std::shared_ptr<const FacePhysics> face_physics;
    detail::OFDGSensorOptions sensor_options;

    detail::OFDGOperatorRepository operator_repository;

    const int dim;
    const int order;
    const int ncomp;

    detail::OFDGMeshData mesh_data;

    detail::OFDGFaceEvaluation face_evaluation;

    mutable detail::OFDGDerivativeProvider derivative_provider;

    mutable detail::OFDGScratch scratch;

    // sigma_elem[c](e,l)
    mutable std::vector<mfem::DenseMatrix> sigma_elem;
    mutable mfem::DenseMatrix element_means;
    mutable mfem::Vector global_means;

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
        element_means.SetSize(ncomp, fes->GetNE());
        element_means = 0.0;

        for (int e = 0; e < fes->GetNE(); ++e) {
            const auto &element = mesh_data.Element(e);
            const detail::OFDGOperators &operators = *element.operators;
            const mfem::DenseMatrix &E = operators.VolumeEvaluation();
            const mfem::Vector &integration_weights =
                operators.VolumeWeights();
            const int nq = E.Height();

            mesh_data.GatherElement(x, e, volume.element_data,
                                    volume.element_matrix);

            mfem::Mult(E, volume.element_matrix, volume.values);

            total_volume += element.volume;

            for (int q = 0; q < nq; ++q) {
                const double physical_weight = integration_weights(q) * element.jacobian_weight;

                for (int c = 0; c < ncomp; ++c) {
                    const double value = volume.values(q, c);

                    volume.mean_integrals(c) += physical_weight * value;
                    element_means(c, e) += physical_weight * value;

                    volume.minimum_values(c) = std::min(volume.minimum_values(c), value);

                    volume.maximum_values(c) = std::max(volume.maximum_values(c), value);
                }
            }

            for (int c = 0; c < ncomp; ++c) {
                element_means(c, e) /= element.volume;
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
        global_means.SetSize(ncomp);
        for (int c = 0; c < ncomp; ++c) {
            const double mean = volume.mean_integrals(c) / total_volume;
            global_means(c) = mean;
            scaling(c) = std::max(mean - volume.minimum_values(c),
                                  volume.maximum_values(c) - mean);
        }
    }

#ifdef MFEM_USE_MPI
    void ComputeFaceNeighborMean(const mfem::ParGridFunction &state,
                                 int element,
                                 const detail::OFDGOperators &operators,
                                 mfem::Vector &mean) const
    {
        mean.SetSize(ncomp);
        mean = 0.0;
        const mfem::Vector &weights = operators.VolumeWeights();
        const mfem::IntegrationRule &rule = operators.Curved() ? detail::CurvedRule(
            operators.Geometry(), operators.QuadratureOrder()) : mfem::IntRules.Get(
            operators.Geometry(), operators.QuadratureOrder());
        for (int q = 0; q < rule.GetNPoints(); ++q) {
            const mfem::IntegrationPoint &point = rule.IntPoint(q);
            for (int c = 0; c < ncomp; ++c) {
                mean(c) += weights(q) * state.GetValue(element, point, c + 1);
            }
        }
        mean /= operators.IntegrationVolume();
    }
#endif

    bool HasResolvedAmplitude(int component,
                              const mfem::Vector &scaling) const
    {
        if (!sensor_options.relative_constant_tolerance) {
            return scaling(component) > 1e-14;
        }
        const double tolerance = 1e-14 *
            std::max(std::abs(global_means(component)), 1.0);
        return scaling(component) > tolerance;
    }

public:

    StabilizationCoreImplementation(const mfem::FiniteElementSpace *fes_,
                                    int btype_)
        : StabilizationCoreImplementation(
              fes_, btype_, std::make_shared<UnitFacePhysics>())
    {
    }

    StabilizationCoreImplementation(
        const mfem::FiniteElementSpace *fes_, int btype_,
        mfem::VectorFunctionCoefficient vel_)
        : StabilizationCoreImplementation(
              fes_, btype_,
              std::make_shared<AdvectionFacePhysics>(
                  std::static_pointer_cast<mfem::VectorCoefficient>(
                      std::make_shared<mfem::VectorFunctionCoefficient>(std::move(vel_)))))
    {
    }

    StabilizationCoreImplementation(
        const mfem::FiniteElementSpace *fes_, int btype_,
        std::shared_ptr<const FacePhysics> face_physics_)
        : StabilizationCoreImplementation(
              fes_, btype_, std::move(face_physics_),
              detail::OFDGSensorOptions{})
    {
    }

    StabilizationCoreImplementation(
        const mfem::FiniteElementSpace *fes_, int btype_,
        std::shared_ptr<const FacePhysics> face_physics_,
        detail::OFDGSensorOptions sensor_options_)
        : fes(detail::RequireSupportedGeometry(fes_)),
          face_physics(std::move(face_physics_)),
          sensor_options(sensor_options_), operator_repository(btype_),
          dim(fes_->GetFE(0)->GetDim()), order(fes_->GetFE(0)->GetOrder()),
          ncomp(fes_->GetVDim()), mesh_data(fes_, operator_repository),
          face_evaluation(fes_, mesh_data, sensor_options.face_rule),
          derivative_provider(fes_, mesh_data)
    {
#ifdef MFEM_USE_MPI
        if (auto *parallel_space = dynamic_cast<mfem::ParFiniteElementSpace *>(
                const_cast<mfem::FiniteElementSpace *>(fes))) {
            parallel_space->ExchangeFaceNbrData();
            mfem::ParMesh *parallel_mesh = parallel_space->GetParMesh();
            for (int neighbor = 0;
                 neighbor < parallel_mesh->GetNFaceNeighborElements();
                 ++neighbor) {
                operator_repository.Get(*parallel_space->GetFaceNbrFE(neighbor),
                    parallel_mesh->GetFaceNbrElementTransformation(neighbor), fes->GetNE() + neighbor);
            }
        }
#endif

        sigma_elem.resize(ncomp);

        for (int c = 0; c < ncomp; ++c) {
            sigma_elem[c].SetSize(fes->GetNE(), order + 1);
        }

        scratch.Initialize(mesh_data.MaxNDof(), ncomp,
                           mesh_data.MaxVolumePoints(), order);

    }

private:

    void ComputeNormalSpeed(mfem::FaceElementTransformations *Tr,
                            const mfem::DenseMatrix &traces1,
                            const mfem::DenseMatrix &traces2,
                            const mfem::Vector *mean1,
                            const mfem::Vector *mean2,
                            double &beta_1, double &beta_2) const
    {

        mfem::Vector normal(dim);
        mfem::Vector &state1 = scratch.face_state1;
        mfem::Vector &state2 = scratch.face_state2;

        beta_1 = 0.0;
        beta_2 = 0.0;

        const mfem::IntegrationRule &ir = detail::FaceRule(order, *Tr);

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
                state1(c) = mean1 ? (*mean1)(c) : traces1(q, c);
                state2(c) = mean2 ? (*mean2)(c) : traces2(q, c);
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
            int element1, element2;
            mesh->GetFaceElements(face, &element1, &element2);
            if (element1 < 0 || element2 < 0) { continue; }
            // Test topology and activity before rebuilding physical geometry.
            if (active && !(*active)[element1] && !(*active)[element2]) {
                continue;
            }
            mfem::FaceElementTransformations *transformations =
                mesh->GetFaceElementTransformations(face);

            const auto derivatives =
                derivative_provider.GetPair(element1, element2);
            face_evaluation.EvaluateDerivativeJumps(
                derivatives.first->all_coefficients,
                derivatives.second->all_coefficients, face, jumps,
                scratch.jump);

            double speed1 = 0.0;
            double speed2 = 0.0;
            mfem::Vector mean1;
            mfem::Vector mean2;
            if (sensor_options.use_element_mean_speed) {
                mean1.SetSize(ncomp);
                mean2.SetSize(ncomp);
                for (int c = 0; c < ncomp; ++c) {
                    mean1(c) = element_means(c, element1);
                    mean2(c) = element_means(c, element2);
                }
            }
            ComputeNormalSpeed(transformations, scratch.jump.traces1,
                               scratch.jump.traces2,
                               sensor_options.use_element_mean_speed ? &mean1 : nullptr,
                               sensor_options.use_element_mean_speed ? &mean2 : nullptr,
                               speed1, speed2);

            const auto &operators1 = *mesh_data.Element(element1).operators;
            const auto &operators2 = *mesh_data.Element(element2).operators;
            for (int l = 0; l <= order; ++l) {
                const double prefactor1 = sensor_options.use_face_height
                    ? operators1.SensorCommon(l) *
                          std::pow(face_evaluation.FaceHeight(face, true), l - 1)
                    : mesh_data.SensorPrefactor(element1, l);
                const double prefactor2 = sensor_options.use_face_height
                    ? operators2.SensorCommon(l) *
                          std::pow(face_evaluation.FaceHeight(face, false), l - 1)
                    : mesh_data.SensorPrefactor(element2, l);

                if (sensor_options.pool_components) {
                    double pooled1 = 0.0;
                    double pooled2 = 0.0;
                    for (int c = 0; c < ncomp; ++c) {
                        if (!HasResolvedAmplitude(c, scaling)) { continue; }
                        pooled1 = std::max(
                            pooled1, prefactor1 * jumps(c, l) / scaling(c));
                        pooled2 = std::max(
                            pooled2, prefactor2 * jumps(c, l) / scaling(c));
                    }
                    for (int c = 0; c < ncomp; ++c) {
                        sigma_elem[c](element1, l) += speed1 * pooled1;
                        sigma_elem[c](element2, l) += speed2 * pooled2;
                    }
                } else {
                    for (int c = 0; c < ncomp; ++c) {
                        if (!HasResolvedAmplitude(c, scaling)) { continue; }
                        const double inverse_scaling = 1.0 / scaling(c);
                        sigma_elem[c](element1, l) +=
                            speed1 * prefactor1 * jumps(c, l) * inverse_scaling;
                        sigma_elem[c](element2, l) +=
                            speed2 * prefactor2 * jumps(c, l) * inverse_scaling;
                    }
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
            const mfem::FiniteElement &local_finite_element =
                *parallel_space->GetFE(local_element);
            const mfem::FiniteElement &neighbor_finite_element =
                *parallel_space->GetFaceNbrFE(neighbor_element);
            const detail::OFDGOperators &local_operators =
                *mesh_data.Element(local_element).operators;
            const detail::OFDGOperators &neighbor_operators =
                operator_repository.Find(neighbor_finite_element, fes->GetNE() + neighbor_element);
            detail::OFDGDerivativeState neighbor_derivatives(
                neighbor_operators.NDof(), ncomp,
                neighbor_operators.DerivativeCount());
            neighbor_state.SetSize(neighbor_vdofs.Size());
            parallel_state.FaceNbrData().GetSubVector(neighbor_vdofs,
                                                       neighbor_state);
            derivative_provider.BuildFaceNeighborState(
                neighbor_state, *transformations->Elem2, neighbor_operators,
                neighbor_derivatives);

            const detail::OFDGDerivativeState &local_derivatives =
                derivative_provider.Get(local_element);
            double local_height = 1.0;
            double neighbor_height = 1.0;
            face_evaluation.EvaluateDerivativeJumps(
                local_derivatives.all_coefficients,
                neighbor_derivatives.all_coefficients, transformations,
                local_finite_element, neighbor_finite_element,
                local_operators, neighbor_operators, jumps,
                scratch.jump, &local_height, &neighbor_height);

            double local_speed = 0.0;
            double neighbor_speed = 0.0;
            mfem::Vector local_mean;
            mfem::Vector neighbor_mean;
            if (sensor_options.use_element_mean_speed) {
                local_mean.SetSize(ncomp);
                for (int c = 0; c < ncomp; ++c) {
                    local_mean(c) = element_means(c, local_element);
                }
                ComputeFaceNeighborMean(
                    parallel_state, parallel_mesh->GetNE() + neighbor_element,
                    neighbor_operators, neighbor_mean);
            }
            ComputeNormalSpeed(transformations, scratch.jump.traces1,
                               scratch.jump.traces2,
                               sensor_options.use_element_mean_speed
                                   ? &local_mean : nullptr,
                               sensor_options.use_element_mean_speed
                                   ? &neighbor_mean : nullptr,
                               local_speed,
                               neighbor_speed);

            for (int l = 0; l <= order; ++l) {
                const double prefactor = sensor_options.use_face_height
                    ? local_operators.SensorCommon(l) *
                          std::pow(local_height, l - 1)
                    : mesh_data.SensorPrefactor(local_element, l);

                if (sensor_options.pool_components) {
                    double pooled = 0.0;
                    for (int c = 0; c < ncomp; ++c) {
                        if (!HasResolvedAmplitude(c, scaling)) { continue; }
                        pooled = std::max(
                            pooled, prefactor * jumps(c, l) / scaling(c));
                    }
                    for (int c = 0; c < ncomp; ++c) {
                        sigma_elem[c](local_element, l) +=
                            local_speed * pooled;
                    }
                } else {
                    for (int c = 0; c < ncomp; ++c) {
                        if (!HasResolvedAmplitude(c, scaling)) { continue; }
                        sigma_elem[c](local_element, l) +=
                            local_speed * prefactor * jumps(c, l) / scaling(c);
                    }
                }
            }
        }
#else
        (void) state;
        (void) scaling;
        (void) active;
#endif
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
    }

public:

    void ComputeStabilization(const mfem::Vector &x, mfem::Vector &S, const mfem::Array<bool> *active = nullptr) const
    {
        S.SetSize(x.Size());
        S = 0.0;

        ComputeJumps(x, active);

        mfem::Vector u_e;
        mfem::Vector tmp;

        mfem::Array<int> vdofs;

        for (int e = 0; e < fes->GetNE(); ++e) {
            if (active && !(*active)[e]) {
                continue;
            }

            const detail::OFDGOperators &operators =
                *mesh_data.Element(e).operators;
            const int ndof = operators.NDof();
            u_e.SetSize(ndof);
            tmp.SetSize(ndof);

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

        mfem::Vector u_e;
        mfem::Vector result;
        mfem::Vector shell;
        mfem::Vector tail_projection;
        mfem::Vector next_tail_projection;

        mfem::Array<int> vdofs;

        for (int e = 0; e < fes->GetNE(); ++e) {
            if (active && !(*active)[e]) {
                continue;
            }

            const detail::OFDGOperators &operators =
                *mesh_data.Element(e).operators;
            const int ndof = operators.NDof();
            u_e.SetSize(ndof);
            result.SetSize(ndof);
            shell.SetSize(ndof);
            tail_projection.SetSize(ndof);
            next_tail_projection.SetSize(ndof);

            for (int c = 0; c < ncomp; ++c) {
                mesh_data.GetComponentVDofs(e, c, vdofs);

                x.GetSubVector(vdofs, u_e);

                result = u_e;

                double cumulative_delta = sigma_elem[c](e, 0);
                // Adjacent shells share S(j) u: evaluate each tail once.
                if (order > 0) { operators.S(0).Mult(u_e, tail_projection); }

                for (int j = 1; j <= order; ++j) {
                    cumulative_delta += sigma_elem[c](e, j);

                    const double factor = std::exp(-t_decay * cumulative_delta);

                    if (j < order) {
                        operators.S(j).Mult(u_e, next_tail_projection);
                        shell = tail_projection;
                        shell -= next_tail_projection;
                        tail_projection = next_tail_projection;
                    } else {
                        shell = tail_projection;
                    }

                    result.Add(factor - 1.0, shell);
                }

                decay.SetSubVector(vdofs, result);
            }
        }

    }
};

class StabilizationCore::Implementation
{
public:
    StabilizationCoreImplementation stabilization;

    Implementation(const mfem::FiniteElementSpace *fes, int basis_type,
                   std::shared_ptr<const FacePhysics> face_physics,
                   detail::OFDGSensorOptions options)
        : stabilization(fes, basis_type, std::move(face_physics), options)
    {
    }
};

StabilizationCore::StabilizationCore(
    const mfem::FiniteElementSpace *fes, int basis_type,
    std::shared_ptr<const FacePhysics> face_physics)
    : StabilizationCore(fes, basis_type, std::move(face_physics),
                        detail::OFDGSensorOptions{})
{
}

StabilizationCore::StabilizationCore(
    const mfem::FiniteElementSpace *fes, int basis_type,
    std::shared_ptr<const FacePhysics> face_physics,
    detail::OFDGSensorOptions options)
    : implementation(std::make_unique<Implementation>(
         fes, basis_type, std::move(face_physics), options))
{
}

StabilizationCore::~StabilizationCore() = default;
StabilizationCore::StabilizationCore(StabilizationCore &&) noexcept = default;
StabilizationCore &StabilizationCore::operator=(StabilizationCore &&) noexcept =
    default;

void StabilizationCore::ComputeStabilization(
    const mfem::Vector &state, mfem::Vector &result,
    const mfem::Array<bool> *active) const
{
    implementation->stabilization.ComputeStabilization(state, result, active);
}

void StabilizationCore::CompDecay(
    const mfem::Vector &state, mfem::Vector &result, double decay_time,
    const mfem::Array<bool> *active) const
{
    implementation->stabilization.CompDecay(state, result, decay_time, active);
}

} // namespace ofdg
