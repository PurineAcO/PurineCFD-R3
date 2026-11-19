#include "implicit/Implicit.h"

void impl::Qcopy(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)dQ.size(); i++) {dQold[i] = dQ[i];}
}
