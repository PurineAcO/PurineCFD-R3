#include "physics/Physics.h"

cc::otphy::otphy(const cc::vecp& phy){
    e = ph::get_energy(phy);
    mu = ph::dynamic_viscosity(phy.T);
    a = ph::get_sonic_velocity(phy.T);
    p = phy.rho * ph::R * phy.T;
}

cc::otphy::otphy(double a_, double p_, double e_) : a(a_), p(p_), e(e_) {}
