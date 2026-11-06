#include "physic.h"

double ph::get_sonic_velocity(double T) {
    return std::sqrt(ph::gamma * ph::R * T);
}

double ph::get_energy(const cc::vecp& phy) {
    return ph::Cv * phy.T + 0.5 * (phy.u * phy.u + phy.v * phy.v + phy.w * phy.w);
}

cc::otphy::otphy(const cc::vecp& phy){
    e = ph::get_energy(phy);
    mu = ph::dynamic_viscosity(phy.T);
    a = ph::get_sonic_velocity(phy.T);
    p = phy.rho * ph::R * phy.T;
}

cc::otphy ph::form_otphy(const cc::vecp &phy){
   return cc::otphy(phy);
}

