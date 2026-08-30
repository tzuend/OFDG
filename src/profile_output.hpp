#pragma once

#include "mfem.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

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
        for (int q = 0; q < nodes.GetNPoints(); ++q) {
            write_point(nodes.IntPoint(q), "polynomial");
        }
        write_point(mfem::Geometries.GetCenter(fe->GetGeomType()), "center");
    }
}

