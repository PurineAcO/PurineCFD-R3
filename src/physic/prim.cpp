#include "physic.h"

void cc::cell_class::prim() {
    phy.rho = conser.c;
    phy.u = conser.x / conser.c;
    phy.v = conser.y / conser.c;
    phy.w = conser.z / conser.c;
    const double e = conser.e / conser.c;
    phy.T = (e - 0.5 * (phy.u * phy.u + phy.v * phy.v + phy.w * phy.w)) / ph::Cv;
    *otphy = ph::form_otphy(phy);
}