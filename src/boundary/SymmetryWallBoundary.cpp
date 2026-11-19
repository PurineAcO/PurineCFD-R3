// 原 src/boundary/SYM.cpp
#include "boundary/Boundary.h"

#include "geometry/Geometry.h"
#include "physics/Physics.h"

void symmetry_wall_boundary(){
    for(cc::bound& b : cc::bounds){
        if(b.type != cc::SYM) {continue;}
        for(cc::face_class* sym : b.fs){
            cc::cell_class* c = cc::boundary_findcell(sym);
            const cc::vec3 n = outward_unit_nor(sym, c)*cc::unit(sym->nor);
            const cc::vec3 v = c->phy.vel();
            const cc::vec3 vs = v - cc::dot(v, n)*n;
            sym->phy.rho = c->phy.rho;
            sym->phy.setvel(vs);
            sym->phy.T = c->phy.T;
            sym->tur->miubl = c->tur->miubl;
            *sym->otphy = ph::form_otphy(sym->phy);
        }
    }
}
