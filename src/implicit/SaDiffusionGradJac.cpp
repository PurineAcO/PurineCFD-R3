#include "implicit/Implicit.h"
#include "physics/Physics.h"
#include "turbulence/Turbulence.h"

cc::mat5 SA::diffusion_grad_jac(const cc::face_class& face, double dn){
    const double len = face.nor.norm();
    const double nx = face.nor.x/len, ny = face.nor.y/len, nz = face.nor.z/len;
    const double u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const double mu = face.otphy->mu;
    const double mut = std::max(face.tur->mueff - mu, 0.0);
    const double m = face.tur->mueff;
    const double k = ph::Cp*(mu/ph::Pr + mut/Prt);
    cc::vec5 cs[5];
    prim_sens(face, cs);
    cc::mat5 a;
    for(int i = 0; i < 5; i++){
        const double du = cs[0][i], dv = cs[1][i], dw = cs[2][i];
        const double divi = nx*du + ny*dv + nz*dw;
        const double txx = (2.0*nx*du - (2.0/3.0)*divi)/dn;
        const double tyy = (2.0*ny*dv - (2.0/3.0)*divi)/dn;
        const double tzz = (2.0*nz*dw - (2.0/3.0)*divi)/dn;
        const double txy = (ny*du + nx*dv)/dn;
        const double txz = (nz*du + nx*dw)/dn;
        const double tyz = (nz*dv + ny*dw)/dn;
        const double fx = m*(nx*txx + ny*txy + nz*txz);
        const double fy = m*(nx*txy + ny*tyy + nz*tyz);
        const double fz = m*(nx*txz + ny*tyz + nz*tzz);
        a.a[0][i] = 0.0;
        a.a[1][i] = len*fx;
        a.a[2][i] = len*fy;
        a.a[3][i] = len*fz;
        a.a[4][i] = len*(u*fx + v*fy + w*fz + k*cs[3][i]/dn);
    }
    return a;
}
