#include "geometry/Geometry.h"
#include <cmath>

std::pair<double, cc::vec3> cube_volumn_cen(cc::cell_class& c){
#if chose2D
    return cube_volumn_cen_2D(c);
#else
    const std::vector<cc::node_class*> nd = sort_cube_node(c);
    static constexpr short ndi[6] = {1, 2, 3, 7, 4, 5};
    double vol = 0.0; cc::vec3 cen = {};
    for(int i=0;i<6;i++) {
        double dvol =  tetrahedron_volumn(nd[0]->p, nd[ndi[i]]->p, nd[ndi[(i+1)%6]]->p, nd[6]->p);
        cen += tetrahedron_cen(nd[0]->p, nd[ndi[i]]->p, nd[ndi[(i+1)%6]]->p, nd[6]->p)*dvol;
        vol += dvol;
    }
    return std::pair<double,cc::vec3>(vol,1.0/vol*cen);
#endif
}
