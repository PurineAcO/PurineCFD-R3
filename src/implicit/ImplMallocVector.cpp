#include "implicit/Implicit.h"

void impl::malloc_vector(){
    if((int)dQ.size() == cc::cell_num) {return;}
    dQ.resize(cc::cell_num);
    dQold.resize(cc::cell_num);
    rhs.resize(cc::cell_num);
}
