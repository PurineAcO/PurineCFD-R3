#include "implicit/Implicit.h"

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
