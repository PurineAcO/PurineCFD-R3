#include "data_structure/DataStructure.h"

void malloc_otphy(cc::cell_class& c,bool mall){
    if(c.otphy != NULL) return;
    if(mall) c.otphy = new cc::otphy();
    else {delete c.otphy;c.otphy = NULL;}
}

void malloc_otphy(cc::face_class& f,bool mall){
    if(f.otphy != NULL) return;
    if(mall) {f.otphy = new cc::otphy();f.otnei[0] = new cc::otphy();f.otnei[1] = new cc::otphy();}
    else {delete f.otphy;delete f.otnei[0];delete f.otnei[1];f.otphy = f.otnei[0] = f.otnei[1] = NULL;}
}
