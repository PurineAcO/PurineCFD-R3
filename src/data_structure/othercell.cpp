#include "node_class.h"
#include "cell_class.h"
#include "face_class.h"

cc::cell_class* cc::othercell(cc::cell_class* c,cc::face_class* f){
    return f->nei[0] == c ? f->nei[1] : f->nei[0];
}