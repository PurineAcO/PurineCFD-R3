#include "geometry.h"

static double triangle_area(cc::vec3 a,cc::vec3 b,cc::vec3 c){
    return 0.5*std::abs(cc::cross((a)-(b),(a)-(c)).norm());
}

double face_area(cc::face_class &f){
    if(f.pcnt == 3) return triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p);
    return triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p) + 
           triangle_area(f.node[0]->p,f.node[3]->p,f.node[2]->p);
}