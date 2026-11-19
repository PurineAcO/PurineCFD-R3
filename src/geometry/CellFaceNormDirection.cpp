#include "geometry/Geometry.h"

void cell_face_norm_direction(cc::cell_class &c){
    allface(c) c.fnorm[i] = cc::dot(c.faces[i]->mid-c.center, c.faces[i]->nor) > 0 ? true : false;
}
