#include "implict.h"

void SA::nucopy(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)dnu.size(); i++) {dnuold[i] = dnu[i];}
}

void impl::Qcopy(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)dQ.size(); i++) {dQold[i] = dQ[i];}
}