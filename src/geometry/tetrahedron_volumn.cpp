#include "geometry.h"
#include <cmath>

/*
return the volunm of a tetrahedron
*/

double tetrahedron_volumn(cc::vec3 a, cc::vec3 b, cc::vec3 c, cc::vec3 d){
    return std::abs(1.0/6*dot(b-a,cc::cross(c-a, d-a)));
}