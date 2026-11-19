#include "implicit/Implicit.h"

void impl::clear_mat(){
#pragma omp parallel for schedule(static)
    allcell icell(i).diag.clear();
}
