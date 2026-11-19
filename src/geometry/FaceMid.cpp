#include "geometry/Geometry.h"
#include <cmath>

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
