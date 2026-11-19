#include "data_structure/DataStructure.h"

void malloc_LSCB(cc::cell_class& c,bool mall){
    if(c.LSCB) return;
    if(mall) c.LSCB = new cc::LSCBmatrix();
    else {delete c.LSCB;c.LSCB = NULL;}
}
