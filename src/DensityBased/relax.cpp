#include "implict.h"
#include "turbulence/SA.h"

void SA::relaxer(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < cc::cell_num; i++){
        double delta_nu = dnu[i];
        if(!std::isfinite(dnu[i])) {delta_nu = 0.0;}
        cc::CellList[i].tur->miubl = std::max(0.0, cc::CellList[i].tur->miubl + delta_nu);
    }
}

static bool ifnotneedrelax(const cc::vec5& newed, const cc::vec5& old){
    const double oldein = old.e - 0.5*(old.x*old.x + old.y*old.y + old.z*old.z)/old.c;
    const double newein = newed.e - 0.5*(newed.x*newed.x + newed.y*newed.y + newed.z*newed.z)/newed.c;
    return (newed.c > 0.5*old.c && newein > 0.5*oldein);
}

int impl::relax(){
    int clip = 0;
#pragma omp parallel for schedule(static) reduction(+ : clip)
    for(int i = 0; i < cc::cell_num; i++){
        cc::cell_class& c = cc::CellList[i];
        double relaxer = 1.0;
        bool clipped = false;
        for(int it = 0; it < 6; it++){
            const cc::vec5 q = c.conser + relaxer*dQ[i];
            if(q.c > 0.0 && std::isfinite(q.c) && ifnotneedrelax(q, c.conser)) {break;}
            relaxer *= 0.5;
            clipped = true;
        }
        if(clipped) {clip++;}
        c.conser = c.conser + relaxer*dQ[i];
    }
    n_clip = clip;
    return clip;
}