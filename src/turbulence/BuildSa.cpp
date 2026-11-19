#include "physics/Physics.h"
#include "turbulence/Turbulence.h"

void SA::build_SA(cc::cell_class& c){
    c.otphy->mu = ph::dynamic_viscosity(c.phy.T);
    c.tur->chi = _chi(c.phy.rho, c.otphy->mu, c.tur->miubl);
    c.tur->fv1 = _fv1(c.tur->chi);
    c.tur->ft2 = _ft2(c.tur->chi);
    c.tur->fv2 = _fv2(c.tur->chi);
    c.tur->Omega = curl_abs(c.phgrad);
    const double d2 = c.tur->sad*c.tur->sad;
    const double sn = c.tur->miubl*(1.0/d2)/(kappa*kappa);
    c.tur->Sbl = std::max(c.tur->Omega + c.tur->fv2*sn, std::max(0.3*c.tur->Omega, 1e-20));
    c.tur->r = std::min(sn/c.tur->Sbl, rmax);
    c.tur->g = _g(c.tur->r);
    c.tur->fw = _fw(c.tur->g);
    c.tur->S2 = strain2(c.phgrad);
}

void SA::build_SA(cc::face_class& f){
    f.otphy->mu = ph::dynamic_viscosity(f.phy.T);
    f.tur->chi = _chi(f.phy.rho, f.otphy->mu, f.tur->miubl);
    f.tur->fv1 = _fv1(f.tur->chi);
}
