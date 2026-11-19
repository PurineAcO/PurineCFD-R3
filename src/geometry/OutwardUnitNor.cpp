#include "geometry/Geometry.h"

double outward_unit_nor(cc::face_class *f, const cc::cell_class *c){
    return cc::dot(f->nor,f->mid - c->center) > 0 ? 1.0 : -1.0 ;
}
