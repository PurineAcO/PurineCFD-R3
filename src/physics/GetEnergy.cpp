#include "physics/Physics.h"

double ph::get_energy(const cc::vecp& phy) {
    return ph::Cv * phy.T + 0.5 * (phy.u * phy.u + phy.v * phy.v + phy.w * phy.w);
}
