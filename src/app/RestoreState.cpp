#include "app/App.h"

void restore_state(){
    for(int i = 0; i < cc::cell_num; i++){
        cc::CellList[i].conser = snap_q[i];
        cc::CellList[i].tur->miubl = snap_nu[i];
        cc::CellList[i].prim();
    }
}
