#include "geometry/Geometry.h"
#include <cmath>

double face_area(cc::face_class &f){
#if chose2D
    return f.nor.norm();
#else
    if(f.pcnt == 3) return triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p);
    return triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p) + 
           triangle_area(f.node[0]->p,f.node[3]->p,f.node[2]->p);
#endif
}
