// 原 src/grad/BarthJespersen.cpp
#include <cmath>

#include "gradient/Gradient.h"

cc::vecp BarthJespersen_2D(const cc::cell_class* cell, const cc::vec3& off){
    const cc::vecp& c = cell->phy;
    cc::vecp lo = c, hi = c;
    allface(*cell){
        const cc::cell_class* o = cell->nei[i];
        if(o == nullptr) {continue;}
        for(int k = 0; k < 5; k++){
            if(o->phy[k] < lo[k]) lo[k] = o->phy[k];
            if(o->phy[k] > hi[k]) hi[k] = o->phy[k];
        }
    }
    const cc::vecp d = cell->phgrad*off;
    cc::vecp r = c;
    for(int k = 0; k < 5; k++){
        const double dd = std::abs(d[k]);
        if(dd <= 0.0) {continue;}
        const double dm = d[k] > 0.0 ? hi[k] - c[k] : c[k] - lo[k];
        const double eps = 1e-12*(c[k]*c[k] + 1e-30);
        const double psi = (dm*dm + 2.0*dm*dd + eps)/(dm*dm + dm*dd + 2.0*dd*dd + eps);
        r[k] = c[k] + psi*d[k];
    }
    return r;
}
