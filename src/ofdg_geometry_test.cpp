#include "mfem.hpp"

#include "ofdg_serial_optimized.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace mfem;

namespace {

void Require(bool condition, const std::string &message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

FaceElementTransformations *FirstInteriorFace(Mesh &mesh)
{
    for (int f = 0; f < mesh.GetNumFaces(); ++f) {
        FaceElementTransformations *Tr = mesh.GetFaceElementTransformations(f);
        if (Tr && Tr->Elem1No >= 0 && Tr->Elem2No >= 0) {
            return Tr;
        }
    }
    return nullptr;
}

FaceElementTransformations *FirstBoundaryFace(Mesh &mesh)
{
    for (int f = 0; f < mesh.GetNumFaces(); ++f) {
        FaceElementTransformations *Tr = mesh.GetFaceElementTransformations(f);
        if (Tr && Tr->Elem1No >= 0 && Tr->Elem2No < 0) {
            return Tr;
        }
    }
    return nullptr;
}

Mesh MakeNonuniformRectangles()
{
    Mesh mesh(2, 6, 2, 6);
    mesh.AddVertex(0.0, 0.0);
    mesh.AddVertex(1.0, 0.0);
    mesh.AddVertex(3.0, 0.0);
    mesh.AddVertex(0.0, 1.0);
    mesh.AddVertex(1.0, 1.0);
    mesh.AddVertex(3.0, 1.0);
    mesh.AddQuad(0, 1, 4, 3);
    mesh.AddQuad(1, 2, 5, 4);
    mesh.AddBdrSegment(0, 1);
    mesh.AddBdrSegment(1, 2);
    mesh.AddBdrSegment(2, 5);
    mesh.AddBdrSegment(5, 4);
    mesh.AddBdrSegment(4, 3);
    mesh.AddBdrSegment(3, 0);
    mesh.FinalizeQuadMesh(1, 0, true);
    return mesh;
}

Mesh MakeSkewedParallelograms()
{
    Mesh mesh(2, 6, 2, 6);
    mesh.AddVertex(0.0, 0.0);
    mesh.AddVertex(1.0, 0.2);
    mesh.AddVertex(3.0, 0.0);
    mesh.AddVertex(0.3, 1.0);
    mesh.AddVertex(1.3, 1.2);
    mesh.AddVertex(3.3, 1.0);
    mesh.AddQuad(0, 1, 4, 3);
    mesh.AddQuad(1, 2, 5, 4);
    mesh.AddBdrSegment(0, 1);
    mesh.AddBdrSegment(1, 2);
    mesh.AddBdrSegment(2, 5);
    mesh.AddBdrSegment(5, 4);
    mesh.AddBdrSegment(4, 3);
    mesh.AddBdrSegment(3, 0);
    mesh.FinalizeQuadMesh(1, 0, true);
    return mesh;
}

Mesh MakeAffineTriangles()
{
    Mesh mesh(2, 4, 2, 4);
    mesh.AddVertex(0.0, 0.0);
    mesh.AddVertex(1.0, 0.1);
    mesh.AddVertex(0.2, 1.3);
    mesh.AddVertex(2.5, 1.6);
    mesh.AddTriangle(0, 1, 2);
    mesh.AddTriangle(1, 3, 2);
    mesh.AddBdrSegment(0, 1);
    mesh.AddBdrSegment(1, 3);
    mesh.AddBdrSegment(3, 2);
    mesh.AddBdrSegment(2, 0);
    mesh.FinalizeTriMesh(1, 0, true);
    return mesh;
}

void Polynomial(const Vector &x, Vector &value)
{
    value.SetSize(2);
    value(0) = 1.0 + 2.0 * x(0) - 0.7 * x(1) + 0.3 * x(0) * x(0) + 0.2 * x(0) * x(1);
    value(1) = -0.4 + 0.5 * x(0) + 1.3 * x(1) - 0.2 * x(1) * x(1);
}

double ScalarPolynomial(const Vector &x)
{
    return 1.0 + 2.0 * x(0) - 0.7 * x(1) + 0.3 * x(0) * x(0) + 0.2 * x(0) * x(1);
}

double YCoordinate(const Vector &x)
{
    return x(1);
}

void ValidateMeanAbsoluteJump()
{
    Mesh mesh = MakeNonuniformRectangles();
    DG_FECollection fec(2, 2, BasisType::GaussLobatto);
    FiniteElementSpace fes(&mesh, &fec);
    GridFunction x(&fes);
    FunctionCoefficient y_coordinate(YCoordinate);
    x.ProjectCoefficient(y_coordinate);

    // Across the interior face, element zero contains u = 0 and element one u = y.
    Array<int> vdofs;
    fes.GetElementVDofs(0, vdofs);
    x.SetSubVector(vdofs, 0.0);

    OFDG ofdg(&fes, BasisType::GaussLobatto);
    FaceElementTransformations *Tr = FirstInteriorFace(mesh);
    Require(Tr != nullptr, "mean-absolute-jump mesh has no interior face");

    DenseMatrix jumps;
    ofdg.ComputeDerivativeJumpsAllComponents(x, Tr, jumps);

    // mean_{y in [0,1]} |0-y| = 1/2, while the old RMS value was sqrt(1/3).
    Require(std::abs(jumps(0, 0) - 0.5) < 1e-12, "solution jump is not the mean absolute jump");
    Require(std::abs(jumps(0, 1) - 1.0) < 1e-12,
            "first-derivative jump is not the mean absolute jump");
    Require(std::abs(jumps(0, 2)) < 1e-12, "second-derivative jump should vanish");
}

void ValidateUniformCase(int components)
{
    Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true, 2.0, 2.0);
    DG_FECollection fec(2, 2, BasisType::GaussLobatto);
    FiniteElementSpace fes(&mesh, &fec, components, Ordering::byVDIM);
    GridFunction x(&fes);

    if (components == 1) {
        FunctionCoefficient polynomial(ScalarPolynomial);
        x.ProjectCoefficient(polynomial);
    } else {
        VectorFunctionCoefficient polynomial(2, Polynomial);
        x.ProjectCoefficient(polynomial);
    }

    OFDG ofdg(&fes, BasisType::GaussLobatto);

    FaceElementTransformations *interior = FirstInteriorFace(mesh);
    Require(interior != nullptr, "uniform mesh has no interior face");

    DenseMatrix jumps;
    ofdg.ComputeDerivativeJumpsAllComponents(x, interior, jumps);
    Require(jumps.FNorm() < 1e-12, "uniform continuous-polynomial derivative jumps are nonzero");

    FaceElementTransformations *boundary = FirstBoundaryFace(mesh);
    Require(boundary != nullptr, "uniform mesh has no boundary face");
    ofdg.ComputeDerivativeJumpsAllComponents(x, boundary, jumps);
    Require(jumps.FNorm() == 0.0, "boundary face produced derivative jumps");

    Vector stabilization;
    Vector decay;
    ofdg.ComputeStabilization(x, stabilization);
    ofdg.CompDecay(x, decay, 0.17);
    Require(stabilization.Norml2() < 1e-10, "uniform continuous polynomial was stabilized");

    decay -= x;
    Require(decay.Norml2() < 1e-10, "uniform continuous polynomial decayed");

    Vector mean;
    Vector scaling;
    ofdg.ComputeMean(x, mean);
    ofdg.ComputeGlobalMeanScaling(x, scaling);
    Require(std::isfinite(mean.Norml2()) && std::isfinite(scaling.Norml2()),
            "uniform mesh produced invalid mean/scaling");

    Array<int> vdofs;
    fes.GetElementVDofs(0, vdofs);
    for (int i = 0; i < vdofs.Size(); ++i) {
        x(vdofs[i]) += 0.1 * static_cast<double>(i + 1);
    }

    Array<bool> active(mesh.GetNE());
    active = false;
    active[0] = true;
    ofdg.ComputeStabilization(x, stabilization, &active);
    ofdg.CompDecay(x, decay, 0.17, &active);
    Require(std::isfinite(stabilization.Norml2()) && std::isfinite(decay.Norml2()),
            "uniform mesh failed active-element filtering");
}

void ValidateAffineMesh(Mesh &mesh, const std::string &name)
{
    DG_FECollection fec(2, 2, BasisType::GaussLobatto);
    FiniteElementSpace fes(&mesh, &fec, 2, Ordering::byVDIM);
    GridFunction x(&fes);
    VectorFunctionCoefficient polynomial(2, Polynomial);
    x.ProjectCoefficient(polynomial);

    OFDG ofdg(&fes, BasisType::GaussLobatto);

    FaceElementTransformations *Tr = FirstInteriorFace(mesh);
    Require(Tr != nullptr, name + " has no interior face");

    DenseMatrix jumps;
    ofdg.ComputeDerivativeJumpsAllComponents(x, Tr, jumps);
    Require(jumps.FNorm() < 1e-12, name + " continuous-polynomial derivative jumps are nonzero");

    Vector mean;
    Vector scaling;
    ofdg.ComputeMean(x, mean);
    ofdg.ComputeGlobalMeanScaling(x, scaling);
    Require(std::isfinite(mean.Norml2()) && std::isfinite(scaling.Norml2()),
            name + " produced invalid mean/scaling");

    Array<int> vdofs;
    fes.GetElementVDofs(0, vdofs);
    for (int i = 0; i < vdofs.Size(); ++i) {
        x(vdofs[i]) += 0.1 * static_cast<double>(i + 1);
    }

    Vector stabilization;
    Vector decay;
    ofdg.ComputeStabilization(x, stabilization);
    ofdg.CompDecay(x, decay, 0.13);
    Require(std::isfinite(stabilization.Norml2()) && std::isfinite(decay.Norml2()),
            name + " produced invalid stabilization/decay");

    Array<bool> active(mesh.GetNE());
    active = false;
    active[0] = true;
    ofdg.ComputeStabilization(x, stabilization, &active);
    ofdg.CompDecay(x, decay, 0.13, &active);
    Require(std::isfinite(stabilization.Norml2()) && std::isfinite(decay.Norml2()),
            name + " failed active-element filtering");
}

} // namespace

int main()
{
    try {
        ValidateUniformCase(1);
        ValidateUniformCase(2);
        ValidateMeanAbsoluteJump();

        Mesh rectangles = MakeNonuniformRectangles();
        ValidateAffineMesh(rectangles, "nonuniform rectangles");

        Mesh parallelograms = MakeSkewedParallelograms();
        ValidateAffineMesh(parallelograms, "skewed parallelograms");

        Mesh triangles = MakeAffineTriangles();
        ValidateAffineMesh(triangles, "affine triangles");
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }

    std::cout << "OFDG geometry and behavior tests passed.\n";
    return 0;
}
