#include "implicit/Implicit.h"

void SA::relaxer(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < cc::cell_num; i++){
        double delta_nu = dnu[i];
        if(!std::isfinite(dnu[i])) {delta_nu = 0.0;}
        cc::CellList[i].tur->miubl = std::max(0.0, cc::CellList[i].tur->miubl + delta_nu);
    }
}
