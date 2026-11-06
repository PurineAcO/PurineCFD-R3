#include "face_class.h"

cc::cell_class* cc::boundary_findcell(cc::face_class* face){
    return face->nei[0] ? face->nei[0] : face->nei[1];
}