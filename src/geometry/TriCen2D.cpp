#include "geometry/Geometry.h"
#include <cmath>

#if chose2D
cc::vec3 tri_cen_2d(const cc::vec3& a, const cc::vec3& b, const cc::vec3& c){
    return (1.0/3.0)*(a + b + c);
}
#endif
