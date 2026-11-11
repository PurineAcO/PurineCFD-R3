#include "euler.h"
#include "physic/physic.h"

cc::mat5tri form_pfgpw_convect(const cc::cell_class& cell){
    const double u = cell.phy.u, v = cell.phy.v, w = cell.phy.w;
    const double g = ph::gamma, g1 = g - 1.0;
    const double phi = 0.5*g1*(u*u + v*v + w*w);
    const double H = ph::Cp*cell.phy.T + 0.5*(u*u + v*v + w*w);
    cc::mat5 Jx(cc::vec5(0.0, 1.0, 0.0, 0.0, 0.0),
                cc::vec5(phi - u*u, (3.0 - g)*u, (1.0 - g)*v, (1.0 - g)*w, g1),
                cc::vec5(-u*v, v, u, 0.0, 0.0),
                cc::vec5(-u*w, w, 0.0, u, 0.0),
                cc::vec5(u*(phi - H), H - g1*u*u, -g1*u*v, -g1*u*w, g*u));
    cc::mat5 Jy(cc::vec5(0.0, 0.0, 1.0, 0.0, 0.0),
                cc::vec5(-u*v, v, u, 0.0, 0.0),
                cc::vec5(phi - v*v, (1.0 - g)*u, (3.0 - g)*v, (1.0 - g)*w, g1),
                cc::vec5(-v*w, 0.0, w, v, 0.0),
                cc::vec5(v*(phi - H), -g1*u*v, H - g1*v*v, -g1*v*w, g*v));
    cc::mat5 Jz(cc::vec5(0.0, 0.0, 0.0, 1.0, 0.0),
                cc::vec5(-u*w, w, 0.0, u, 0.0),
                cc::vec5(-v*w, 0.0, w, v, 0.0),
                cc::vec5(phi - w*w, (1.0 - g)*u, (1.0 - g)*v, (3.0 - g)*w, g1),
                cc::vec5(w*(phi - H), -g1*u*w, -g1*v*w, H - g1*w*w, g*w));
    return cc::mat5tri(Jx, Jy, Jz);
}