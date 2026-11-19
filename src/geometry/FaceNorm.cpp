#include "geometry/Geometry.h"

cc::vec3 face_norm(cc::face_class &f){
#if chose2D
    return cc::vec3(f.node[0]->p.y - f.node[1]->p.y, f.node[1]->p.x - f.node[0]->p.x, 0.0);
#else
    if(f.pcnt == 3) return 0.5*cc::cross(f.node[1]->p-f.node[0]->p, f.node[2]->p-f.node[0]->p);
    return 0.5*cc::cross(f.node[2]->p-f.node[0]->p, f.node[3]->p-f.node[1]->p);
#endif
}
