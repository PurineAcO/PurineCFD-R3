#include "geometry/Geometry.h"
#include <cmath>

#if chose2D
double tri_area_2d(const cc::vec3& a, const cc::vec3& b, const cc::vec3& c){
    return 0.5*std::abs((b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y));
}
#endif
