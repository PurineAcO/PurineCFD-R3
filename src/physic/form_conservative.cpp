#include "physic.h"

void cc::cell_class::form_conservative() {
    conser.c = phy.rho;
    conser.x = phy.rho * phy.u;
    conser.y = phy.rho * phy.v;
    conser.z = phy.rho * phy.w;
    conser.e = phy.rho * ph::get_energy(phy);
    conserformer = conser;
    lastconser = conser;
}