#include "physics/Physics.h"

double ph::dynamic_pressure(const cc::vecp& refp) {
    return 0.5 * refp.rho * (refp.u * refp.u + refp.v * refp.v + refp.w * refp.w);
}
