#include "../data_structure/cell_class.h"
#include "../data_structure/face_class.h"
#include "data_structure/classconfig.h"
#include "data_structure/vectors.h"
#include "turbulence/SA.h"
#include "../physic/physic.h"
#include "grad/LSCB.h"

/*
执行任何初始化之前必须执行allocate_fields,否则没有内存
*/


void malloc_tur(cc::cell_class& c,bool mall){
    if(c.tur != NULL) return;
    #if SAchosen
        if (mall) c.tur = new SA::turbulence();
        else {delete c.tur;c.tur = NULL;}
    #elif SSTchosen
        /* waiting for SST */
    #endif
}

void malloc_tur(cc::face_class& f,bool mall){
    if(f.tur != NULL) return;
    #if SAchosen
        if(mall) f.tur = new SA::turbulence();
        else {delete f.tur;f.tur = NULL;}
    #elif SSTchosen
        /* waiting for SST */
    #endif
}

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

void malloc_LSCB(cc::cell_class& c,bool mall){
    if(c.LSCB) return;
    if(mall) c.LSCB = new cc::LSCBmatrix();
    else {delete c.LSCB;c.LSCB = NULL;}
}

void allocate_fields(){
    allcell{malloc_otphy(icell(i),1);malloc_tur(icell(i),1);malloc_LSCB(icell(i),1);}
    allfac{malloc_otphy(iface(i),1);malloc_tur(iface(i),1);}
}

void release_fields(){
    allcell{malloc_otphy(icell(i),0);malloc_tur(icell(i),0);malloc_LSCB(icell(i),0);}
    allfac{malloc_otphy(iface(i),0);malloc_tur(iface(i),0);}
}
