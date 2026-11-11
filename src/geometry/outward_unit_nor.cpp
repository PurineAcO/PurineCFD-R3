#include "data_structure/vectors.h"
#include "geometry.h"

double outward_unit_nor(cc::face_class *f, const cc::cell_class *c){
    /* 判断外法向方向,由于不能快速从face转移到cell,所以需要本函数 */
    return cc::dot(f->nor,f->mid - c->center) > 0 ? 1.0 : -1.0 ;
}