#include "geometry.h"

#if !chose2D
static double ta(cc::vec3 a,cc::vec3 b,cc::vec3 c){
    return 0.5*std::abs(cc::cross((a)-(b),(a)-(c)).norm());
}
#endif

cc::vec3 face_mid(cc::face_class &f){
#if chose2D
    return 0.5*(f.node[0]->p + f.node[1]->p);
#else
    if(f.pcnt == 3) return 1.0/3*(f.node[0]->p+f.node[1]->p+f.node[2]->p);
    else{return 1.0/f.area*
             (1.0/3*(f.node[0]->p+f.node[1]->p+f.node[2]->p)*ta(f.node[0]->p,f.node[1]->p,f.node[2]->p) + 
               1.0/3*(f.node[0]->p+f.node[3]->p+f.node[2]->p)*ta(f.node[0]->p,f.node[3]->p,f.node[2]->p));
    }
#endif
}