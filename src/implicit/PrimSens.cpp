#include "implicit/Implicit.h"
#include "physics/Physics.h"
#include "turbulence/Turbulence.h"

void prim_sens(const cc::face_class& face, cc::vec5 cs[5]){
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const double T = face.phy.T, nu = face.tur->miubl, Cv = ph::Cv;
    const double q2 = u*u + v*v + w*w;
    cs[0] = cc::vec5(-u/rho, 1.0/rho, 0.0, 0.0, 0.0);
    cs[1] = cc::vec5(-v/rho, 0.0, 1.0/rho, 0.0, 0.0);
    cs[2] = cc::vec5(-w/rho, 0.0, 0.0, 1.0/rho, 0.0);
    cs[3] = cc::vec5(q2/(2.0*Cv*rho) - T/rho, -u/(Cv*rho), -v/(Cv*rho), -w/(Cv*rho), 1.0/(Cv*rho));
    cs[4] = cc::vec5(-nu/rho, 0.0, 0.0, 0.0, 1.0/rho);
}
