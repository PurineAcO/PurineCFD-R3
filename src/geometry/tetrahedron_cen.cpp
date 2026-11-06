#include "geometry.h"

/*
return the center of a teteahedron 
*/

cc::vec3 tetrahedron_cen(cc::vec3 a, cc::vec3 b, cc::vec3 c, cc::vec3 d){
    return 0.25 * (a+b+c+d);
}