#include "LSCB.h"
#include "data_structure/classconfig.h"
#include "turbulence/SA.h"

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
        /* waiting for SST */
    #endif
}