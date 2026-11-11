#include "data_structure/cell_class.h"
#include "data_structure/node_class.h"
#include "euler.h"

void assemble_flux(cc::cell_class& cell){
    cell.convect.clear();
    for(int i=0;i<cell.fcnt;i++){
        cc::face_class* face = cell.faces[i];
        const bool islow = (face->nei[0] == &cell);
        cell.convect += (2*islow - 1)*face->convect;
    }
}