#include "../data_structure/cell_class.h"
#include "../data_structure/face_class.h"
#include "../physic/physic.h"
#include "data_structure/classconfig.h"
#include "initialization.h"
#include "turbulence/SA.h"

/*
执行任何初始化之前必须执行allocate_fields,否则没有内存
*/

void allocate_fields(){
    allcell{
        if(icell(i).otphy == nullptr) icell(i).otphy = new cc::otphy();
        if(icell(i).tur == nullptr) icell(i).tur = new SA::turbulence();
    }
    allfac{
        if(iface(i).otphy == nullptr) iface(i).otphy = new cc::otphy();
        if(iface(i).tur == nullptr) iface(i).tur = new SA::turbulence();
        for(int k = 0; k < 2; k++) if(iface(i).otnei[k] == nullptr) iface(i).otnei[k] = new cc::otphy();
    }
}

void release_fields(){
    allcell{delete icell(i).otphy;delete icell(i).tur;icell(i).otphy = nullptr;icell(i).tur = nullptr;}
    allfac{
        delete iface(i).otphy; delete iface(i).tur;iface(i).otphy = nullptr;iface(i).tur = nullptr;
        for(int k = 0; k < 2; k++) {delete iface(i).otnei[k];iface(i).otnei[k] = nullptr;}
    }
}
