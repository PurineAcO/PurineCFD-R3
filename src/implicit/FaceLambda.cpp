#include "implicit/Implicit.h"
#include "physics/Physics.h"

double impl::face_lambda(const cc::face_class* f, const cc::cell_class& cell){
    const double len = f->nor.norm();
    const double un = std::abs(f->phy.u*f->nor.x + f->phy.v*f->nor.y + f->phy.w*f->nor.z)/len;
    const double Cn = std::max(4.0/3.0, ph::gamma/ph::Pr);
    return (un + f->otphy->a)*len + Cn*f->tur->mueff/f->phy.rho*len*len/cell.vol;
}
