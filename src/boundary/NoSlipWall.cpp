#include "data_structure/classconfig.h"
#include "data_structure/face_class.h"
#include "turbulence/SA.h"
#include "NoSlipWall.h"
#include "physic/physic.h"

#if SAchosen
#endif

void noslip_wall_boundary(){
    for(cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) {continue;}
        for(cc::face_class* wall : b.fs){
            cc::cell_class* c = cc::boundary_findcell(wall);
            wall->phy.u = 0.0;
            wall->phy.v = 0.0;
            wall->phy.w = 0.0;
            wall->phy.T = c->phy.T;
            wall->otphy->p = c->otphy->p;
            wall->phy.rho = wall->otphy->p/(ph::R*wall->phy.T);
            wall->tur->miubl = 0.0;
            *wall->otphy = ph::form_otphy(wall->phy);
        }
    }
}