// 原 src/grad/LSCB_turbulence.cpp
#include "gradient/Gradient.h"

void LSCB_turbulence(cc::cell_class &cell){
    cc::vec3 bmiubl{};
    #if SAchosen
        allface(cell){
            cc::cell_class* nei = cell.nei[i];
            if(nei == nullptr) {continue;}
            bmiubl += (nei->tur->miubl-cell.tur->miubl)*(nei->center-cell.center)*cell.LSCB->wi[i];
        }
        cell.tur->miublgrad = cell.LSCB->solve(bmiubl);
    #elif SSTchosen
    #endif
}
