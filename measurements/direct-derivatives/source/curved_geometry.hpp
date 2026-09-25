#pragma once

#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <stdexcept>
#include <limits>
#include <array>

namespace ofdg { namespace detail {

// The pinned MFEM cannot supply two-sided native NURBS face maps. Reject
// before touching its face cache; callers may explicitly project geometry.
template <typename Space>
inline Space *RequireSupportedGeometry(Space *space)
{
    if (space->GetMesh()->NURBSext) {
        throw std::invalid_argument(
            "Native NURBS faces are unsupported by the pinned MFEM. "
            "Project the mesh with SetCurvature before constructing the DG space "
            "(this approximates the rational geometry). See docs/curved-validation.md.");
    }
    return space;
}

// Static geometry only. A higher-order representation of an affine map may
// use the curved path; linear tensor-product maps are checked for warping.
inline bool IsAffine(mfem::ElementTransformation &transformation)
{
    auto *iso = dynamic_cast<mfem::IsoparametricTransformation *>(&transformation);
    if (!iso || iso->GetFE()->GetOrder() > 1) { return false; }
    const auto geometry = transformation.GetGeometryType();
    transformation.SetIntPoint(&mfem::Geometries.GetCenter(geometry));
    const mfem::DenseMatrix center(transformation.Jacobian());
    const auto &rule = mfem::IntRules.Get(geometry, 4);
    for (int q = 0; q < rule.GetNPoints(); ++q) {
        transformation.SetIntPoint(&rule.IntPoint(q));
        const auto &jacobian = transformation.Jacobian();
        for (int j = 0; j < center.Width(); ++j) {
            for (int i = 0; i < center.Height(); ++i) {
                if (std::abs(jacobian(i,j) - center(i,j)) >
                    64 * std::numeric_limits<double>::epsilon() * center.MaxMaxNorm()) { return false; }
            }
        }
    }
    return true;
}

// Overintegration also resolves rational geometry approximately. Convergence
// under increased quadrature must be checked for each new geometry family.
inline int CurvedQuadratureOrder(int solution_order,
                                 mfem::ElementTransformation &transformation)
{
    return 2 * solution_order + 2 * std::max(0, transformation.OrderW()) + 8;
}

inline int FaceQuadratureOrder(int solution_order,
                               mfem::FaceElementTransformations &face)
{
    if (IsAffine(*face.Elem1) && IsAffine(*face.Elem2)) {
        return 2 * solution_order + 1;
    }
    return std::max(CurvedQuadratureOrder(solution_order, *face.Elem1),
                    CurvedQuadratureOrder(solution_order, *face.Elem2));
}

// Positive tensor/Duffy quadrature avoids cancellation in high-order
// simplex rules. MFEM's affine rules are deliberately left unchanged.
inline const mfem::IntegrationRule &CurvedRule(mfem::Geometry::Type geometry, int order)
{
    if (geometry != mfem::Geometry::TRIANGLE &&
        geometry != mfem::Geometry::TETRAHEDRON &&
        geometry != mfem::Geometry::PRISM) {
        return mfem::IntRules.Get(geometry, order);
    }
    using Key = std::pair<int,int>;
    static std::map<Key,std::unique_ptr<mfem::IntegrationRule>> rules;
    const Key key{geometry,order};
    auto &cached = rules[key];
    if (cached) { return *cached; }
    const auto &line = mfem::IntRules.Get(mfem::Geometry::SEGMENT, order + 2);
    const int n = line.GetNPoints();
    const bool triangle = geometry == mfem::Geometry::TRIANGLE;
    cached = std::make_unique<mfem::IntegrationRule>(triangle ? 6*n*n : n*n*n);
    int q = 0;
    for (int i=0; i<n; ++i) {
        const auto &u = line.IntPoint(i);
        for (int j=0; j<n; ++j) {
            const auto &v = line.IntPoint(j);
            for (int k=0; k<(triangle ? 1 : n); ++k) {
                const auto &w = line.IntPoint(k);
                const double x = u.x, y = (1-u.x)*v.x;
                const double weight = u.weight*v.weight*(1-u.x);
                if (triangle) {
                    // All six barycentric permutations make face samples
                    // identical under serial/MPI face reorientation.
                    const std::array<double,3> barycentric{1-x-y, x, y};
                    std::array<int,3> permutation{0,1,2};
                    do {
                        auto &point = cached->IntPoint(q++);
                        point.x = barycentric[permutation[1]];
                        point.y = barycentric[permutation[2]];
                        point.weight = weight/6;
                    } while (std::next_permutation(permutation.begin(), permutation.end()));
                } else {
                    auto &point = cached->IntPoint(q++);
                    point.x = x; point.y = y;
                    const double scale = geometry == mfem::Geometry::PRISM ? 1.0 : (1-u.x)*(1-v.x);
                    point.z = scale*w.x;
                    point.weight = weight*w.weight*scale;
                }
            }
        }
    }
    cached->SetOrder(order);
    return *cached;
}

inline const mfem::IntegrationRule &VolumeRule(
    int order, mfem::ElementTransformation &map, int affine_extra = 1)
{
    return IsAffine(map) ? mfem::IntRules.Get(map.GetGeometryType(), 2*order+affine_extra) :
        CurvedRule(map.GetGeometryType(), CurvedQuadratureOrder(order, map));
}

inline const mfem::IntegrationRule &FaceRule(
    int order, mfem::FaceElementTransformations &face, int affine_extra = 1)
{
    const int quadrature = FaceQuadratureOrder(order, face);
    return quadrature == 2*order+1 ?
        mfem::IntRules.Get(face.GetGeometryType(), 2*order+affine_extra) :
        CurvedRule(face.GetGeometryType(), quadrature);
}

inline mfem::Vector NormalizedPhysicalFaceWeights(
    mfem::FaceElementTransformations &face, const mfem::IntegrationRule &rule)
{
    mfem::Vector weights(rule.GetNPoints());
    double measure = 0.0;
    for (int q = 0; q < rule.GetNPoints(); ++q) {
        face.SetAllIntPoints(&rule.IntPoint(q));
        weights(q) = rule.IntPoint(q).weight * face.Face->Weight();
        measure += weights(q);
    }
    weights /= measure;
    return weights;
}

// Coefficient matrices for L2-projected first physical derivatives. Keeping
// this assembly separate also permits independent manufactured validation.
inline std::vector<mfem::DenseMatrix> ProjectedPhysicalDerivatives(
    const mfem::FiniteElement &element,
    mfem::ElementTransformation &transformation, int quadrature_order)
{
    const int ndof = element.GetDof(), dim = element.GetDim();
    mfem::DenseMatrix mass(ndof);
    mass = 0.0;
    std::vector<mfem::DenseMatrix> gradient(dim), derivatives(dim);
    for (auto &matrix : gradient) { matrix.SetSize(ndof); matrix = 0.0; }
    mfem::Vector shape(ndof);
    mfem::DenseMatrix dshape(ndof, dim);
    const auto &rule = CurvedRule(element.GetGeomType(), quadrature_order);
    for (int q = 0; q < rule.GetNPoints(); ++q) {
        const auto &point = rule.IntPoint(q);
        transformation.SetIntPoint(&point);
        element.CalcShape(point, shape);
        element.CalcPhysDShape(transformation, dshape);
        const double weight = point.weight * transformation.Weight();
        mfem::AddMult_a_VVt(weight, shape, mass);
        for (int d = 0; d < dim; ++d) {
            for (int j = 0; j < ndof; ++j) {
                for (int i = 0; i < ndof; ++i) {
                    gradient[d](i,j) += weight * shape(i) * dshape(j,d);
                }
            }
        }
    }
    mfem::DenseMatrixInverse inverse(mass);
    for (int d = 0; d < dim; ++d) {
        derivatives[d].SetSize(ndof);
        inverse.Mult(gradient[d], derivatives[d]);
    }
    return derivatives;
}

} } // namespace ofdg::detail
