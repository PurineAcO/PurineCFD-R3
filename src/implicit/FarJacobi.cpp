#include "implicit/Implicit.h"

void impl::far_jacobi(cc::face_class* face, cc::cell_class* cell, double s){
    const double len = face->nor.norm();
    const double un = std::abs(face->phy.u*face->nor.x + face->phy.v*face->nor.y +
                               face->phy.w*face->nor.z)/len;
    const double lamc = (un + face->otphy->a)*len;
    cell->diag += 0.5*lamc*cc::eyes;
    cell->diag += (-s)*face->imats.Bf;
}
