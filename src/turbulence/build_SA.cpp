#include "SA.h"
#include "physic/physic.h"



static double curl_abs(const cc::vecgrad& g){
    const double cx = g.wgrad.y - g.vgrad.z;
    const double cy = g.ugrad.z - g.wgrad.x;
    const double cz = g.vgrad.x - g.ugrad.y;
    return std::sqrt(cx*cx + cy*cy + cz*cz);
}

static double strain2(const cc::vecgrad& g){
    return 2.0*(g.ugrad.x*g.ugrad.x + g.vgrad.y*g.vgrad.y + g.wgrad.z*g.wgrad.z) +
           (g.ugrad.y + g.vgrad.x)*(g.ugrad.y + g.vgrad.x) +
           (g.ugrad.z + g.wgrad.x)*(g.ugrad.z + g.wgrad.x) +
           (g.vgrad.z + g.wgrad.y)*(g.vgrad.z + g.wgrad.y);
}

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