#pragma once

#include "mfem.hpp"
#include "curved_geometry.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

// Integrate in physical coordinates. Euler callers average conservative variables
// first; primitive variables derived from that mean are not primitive averages.
template <typename Evaluate>
inline mfem::Vector PhysicalCellAverage(mfem::ElementTransformation &map,
                                       int order, int components, Evaluate evaluate)
{
    mfem::Vector average(components), value(components);
    average = 0.0;
    mfem::real_t volume = 0.0;
    const auto &rule = ofdg::detail::VolumeRule(order, map);
    for (int q = 0; q < rule.GetNPoints(); ++q) {
        const auto &ip = rule.IntPoint(q);
        map.SetIntPoint(&ip);
        const auto weight = ip.weight * map.Weight();
        evaluate(ip, value);
        average.Add(weight, value);
        volume += weight;
    }
    MFEM_VERIFY(volume > 0.0, "Nonpositive cell volume in profile export.");
    average /= volume;
    return average;
}

inline void WriteScalarSamples(const mfem::GridFunction &solution,
                               const std::string &prefix, int rank)
{
    if (prefix.empty()) { return; }
    std::ostringstream name;
    name << prefix << ".rank" << std::setfill('0') << std::setw(6)
         << rank << ".csv";
    std::ofstream output(name.str());
    MFEM_VERIFY(output, "Unable to open profile output " << name.str());
    output << "rank,element,sample,x,y,value\n";
    output << std::setprecision(17);
    const mfem::FiniteElementSpace *space = solution.FESpace();
    const mfem::Mesh *mesh = space->GetMesh();
    const int dim = mesh->Dimension();
    for (int element = 0; element < mesh->GetNE(); ++element) {
        const mfem::FiniteElement *fe = space->GetFE(element);
        mfem::ElementTransformation *transformation =
            space->GetMesh()->GetElementTransformation(element);
        const mfem::IntegrationRule &nodes = fe->GetNodes();
        auto write_point = [&](const mfem::IntegrationPoint &ip,
                               const char *sample) {
            mfem::Vector point;
            transformation->Transform(ip, point);
            output << rank << ',' << element << ',' << sample << ','
                   << point(0) << ',' << (dim > 1 ? point(1) : 0.0)
                   << ',' << solution.GetValue(element, ip) << '\n';
        };
        if (dim == 1) {
            for (int q = 0; q <= 20; ++q) {
                mfem::IntegrationPoint ip;
                ip.Set1w(q / 20.0, 1.0);
                write_point(ip, "polynomial");
            }
        } else {
            for (int q = 0; q < nodes.GetNPoints(); ++q) {
                write_point(nodes.IntPoint(q), "polynomial");
            }
        }
        const auto average = PhysicalCellAverage(*transformation, fe->GetOrder(), 1,
            [&](const mfem::IntegrationPoint &ip, mfem::Vector &value) {
                value(0) = solution.GetValue(element, ip);
            });
        mfem::Vector point;
        transformation->Transform(mfem::Geometries.GetCenter(fe->GetGeomType()), point);
        output << rank << ',' << element << ",cell_average," << point(0) << ','
               << (dim > 1 ? point(1) : 0.0) << ',' << average(0) << '\n';
        write_point(mfem::Geometries.GetCenter(fe->GetGeomType()), "center");
    }
}

