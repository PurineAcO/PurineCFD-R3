#include "turbulence/Turbulence.h"

double strain2(const cc::vecgrad& g){
    return 2.0*(g.ugrad.x*g.ugrad.x + g.vgrad.y*g.vgrad.y + g.wgrad.z*g.wgrad.z) +
           (g.ugrad.y + g.vgrad.x)*(g.ugrad.y + g.vgrad.x) +
           (g.ugrad.z + g.wgrad.x)*(g.ugrad.z + g.wgrad.x) +
           (g.vgrad.z + g.wgrad.y)*(g.vgrad.z + g.wgrad.y);
}
