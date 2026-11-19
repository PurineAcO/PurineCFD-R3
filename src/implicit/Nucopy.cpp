#include "implicit/Implicit.h"

void SA::nucopy(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)dnu.size(); i++) {dnuold[i] = dnu[i];}
}
