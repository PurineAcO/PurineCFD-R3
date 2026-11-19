#include "geometry/Geometry.h"
#include <cmath>

#if chose2D
std::pair<double, cc::vec3> cube_volumn_cen_2D(cc::cell_class& c){
    cc::vec3 nd[4];
    int n = 0;
    for(int i = 0; i < c.fcnt && n < 4; i++){
        const cc::face_class* f = c.faces[i];
        for(int k = 0; k < 2 && n < 4; k++){
            const cc::vec3& p = f->node[k]->p;
            bool dup = false;
            for(int m = 0; m < n; m++){
                if(nd[m].x == p.x && nd[m].y == p.y) {dup = true; break;}
            }
            if(!dup) {nd[n++] = p;}
        }
    }
    if(n < 3) {return std::pair<double, cc::vec3>(0.0, cc::vec3{});}
    if(n == 3){
        const double s = tri_area_2d(nd[0], nd[1], nd[2]);
        return std::pair<double, cc::vec3>(s, tri_cen_2d(nd[0], nd[1], nd[2]));
    }
    const double s[4] = {tri_area_2d(nd[0], nd[1], nd[2]), tri_area_2d(nd[0], nd[1], nd[3]),
                         tri_area_2d(nd[0], nd[2], nd[3]), tri_area_2d(nd[1], nd[2], nd[3])};
    const cc::vec3 g[4] = {tri_cen_2d(nd[0], nd[1], nd[2]), tri_cen_2d(nd[0], nd[1], nd[3]),
                           tri_cen_2d(nd[0], nd[2], nd[3]), tri_cen_2d(nd[1], nd[2], nd[3])};
    const double sum = s[0] + s[1] + s[2] + s[3];
    cc::vec3 cen;
    for(int i = 0; i < 4; i++) {cen += s[i]*g[i];}
    cen *= 1.0/sum;
    cen.z = 0.0;
    return std::pair<double, cc::vec3>(0.5*sum, cen);
}
#endif
