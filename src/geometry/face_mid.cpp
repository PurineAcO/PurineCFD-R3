#include "geometry.h"

static double ta(cc::vec3 a,cc::vec3 b,cc::vec3 c){
    return 0.5*std::abs(cc::cross((a)-(b),(a)-(c)).norm());
}

cc::vec3 face_mid(cc::face_class &f){
    if(f.pcnt == 3) return 1.0/3*(f.node[0]->p+f.node[1]->p+f.node[2]->p);
    else{return 1.0/f.area*
             (1.0/3*(f.node[0]->p+f.node[1]->p+f.node[2]->p)*ta(f.node[0]->p,f.node[1]->p,f.node[2]->p) + 
               1.0/3*(f.node[0]->p+f.node[3]->p+f.node[2]->p)*ta(f.node[0]->p,f.node[3]->p,f.node[2]->p));
    }
}