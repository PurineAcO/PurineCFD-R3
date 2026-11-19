#include "app/App.h"

void save_state(){
    if((int)snap_q.size() != cc::cell_num){
        snap_q.resize(cc::cell_num);
        snap_nu.resize(cc::cell_num);
    }
    for(int i = 0; i < cc::cell_num; i++){
        snap_q[i] = cc::CellList[i].conser;
        snap_nu[i] = cc::CellList[i].tur->miubl;
    }
}
