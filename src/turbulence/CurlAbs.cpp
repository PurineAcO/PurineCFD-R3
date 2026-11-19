#include "turbulence/Turbulence.h"

double curl_abs(const cc::vecgrad& g){
    const double cx = g.wgrad.y - g.vgrad.z;
    const double cy = g.ugrad.z - g.wgrad.x;
    const double cz = g.vgrad.x - g.ugrad.y;
    return std::sqrt(cx*cx + cy*cy + cz*cz);
}
