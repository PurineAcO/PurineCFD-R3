#include "implicit/Implicit.h"

void SA::malloc_vector(){
    if((int)dnu.size() == cc::cell_num) {return;}
    dnu.resize(cc::cell_num);
    dnuold.resize(cc::cell_num);
    nurhs.resize(cc::cell_num);
}
