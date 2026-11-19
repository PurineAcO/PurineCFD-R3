#include "mesh/Mesh.h"
#include "data_structure/DataStructure.h"

bool linkmesh(){
    for(int i = 0; i < cc::cell_num; i++){
        cc::cell_class& c = cc::CellList[i];
        for(int k = 0; k < c.fcnt; k++){
            cc::face_class* f = c.faces[k];
            c.nei[k] = (f->nei[0] == &c) ? f->nei[1] : f->nei[0];
        }
    }
    return true;
}
