// 原 src/grad/BarthJespersen.cpp
#include "gradient/Gradient.h"

cc::vecp BarthJespersen(const cc::cell_class* cell, const cc::vec3& off){
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
    for(int k = 0; k < 5; k++){r[k] = c[k] + shrink(c[k], lo[k], hi[k], d[k])*d[k];}
    return r;
}
