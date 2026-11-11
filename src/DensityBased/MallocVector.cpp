#include "implict.h"

void SA::malloc_vector(){
    if((int)dnu.size() == cc::cell_num) {return;}
    dnu.resize(cc::cell_num);
    dnuold.resize(cc::cell_num);
    nurhs.resize(cc::cell_num);
}

void impl::malloc_vector(){
    if((int)dQ.size() == cc::cell_num) {return;}
    dQ.resize(cc::cell_num);
    dQold.resize(cc::cell_num);
    rhs.resize(cc::cell_num);
}

void impl::clear_mat(){
#pragma omp parallel for schedule(static)
    allcell icell(i).diag.clear();
}
