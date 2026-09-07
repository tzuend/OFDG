// Minimal upstream MFEM reproducer. Run separately: the pinned release build
// may segfault; a debug MFEM build explicitly rejects the two-sided map.
#include "mfem.hpp"
#include <string>
int main()
{
    mfem::Mesh mesh((std::string(OFDG_MFEM_DATA_DIR)+"/disc-nurbs.mesh").c_str(),1,1);
    for(int f=0;f<mesh.GetNumFaces();++f) {
        int first,second;mesh.GetFaceElements(f,&first,&second);
        if(second<0){continue;}
        auto *face=mesh.GetFaceElementTransformations(f);
        auto point=mfem::Geometries.GetCenter(face->GetGeometryType());
        face->SetAllIntPoints(&point);
        mfem::out << face->Face->Weight() << std::endl;
    }
}
