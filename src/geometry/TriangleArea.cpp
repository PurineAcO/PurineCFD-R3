#include "geometry/Geometry.h"
#include <cmath>

#if !chose2D
double triangle_area(cc::vec3 a,cc::vec3 b,cc::vec3 c){
    return 0.5*std::abs(cc::cross((a)-(b),(a)-(c)).norm());
}
#endif
