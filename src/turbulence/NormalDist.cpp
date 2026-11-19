#include "turbulence/Turbulence.h"

double SA::normal_dist(const cc::face_class* f){
    const double len = f->nor.norm();
    if(f->type == cc::INTER && f->nei[0] != nullptr && f->nei[1] != nullptr){
        return std::abs(cc::dot(f->nei[0]->center - f->nei[1]->center, f->nor))/len;
    }
    const cc::cell_class* c = cc::boundary_findcell(const_cast<cc::face_class*>(f));
    if(f->type == cc::WALL) {return c->tur->sad;}
    return std::abs(cc::dot(f->mid - c->center, f->nor))/len;
}
