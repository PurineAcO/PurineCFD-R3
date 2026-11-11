#include "geometry.h"

void tangent_basis(const cc::vec3& n, cc::vec3& t1, cc::vec3& t2){
    const cc::vec3 seed = std::abs(n.x) < 0.9 ? cc::vec3(1.0,0.0,0.0) : cc::vec3(0.0,1.0,0.0);
    t1 = cc::cross(n, seed);
    t1 = (1.0/t1.norm())*t1;
    t2 = cc::cross(n, t1);
}