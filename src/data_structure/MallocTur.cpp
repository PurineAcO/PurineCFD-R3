#include "data_structure/DataStructure.h"

void malloc_tur(cc::cell_class& c,bool mall){
    if(c.tur != NULL) return;
    #if SAchosen
        if (mall) c.tur = new SA::turbulence();
        else {delete c.tur;c.tur = NULL;}
    #elif SSTchosen
    #endif
}

void malloc_tur(cc::face_class& f,bool mall){
    if(f.tur != NULL) return;
    #if SAchosen
        if(mall) f.tur = new SA::turbulence();
        else {delete f.tur;f.tur = NULL;}
    #elif SSTchosen
    #endif
}
