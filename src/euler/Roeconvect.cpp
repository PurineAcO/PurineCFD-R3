#include "euler.h"
#include "physic/physic.h"

void convect_ROE(cc::face_class& face){
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const cc::otphy& Lo = *face.otnei[0];
    const cc::otphy& Ro = *face.otnei[1];
    const roe_waves w = roe_decompose(face);
    const cc::vec3& n = w.n;
    const double unl = cc::dot(L.vel(), n);
    const double unr = cc::dot(R.vel(), n);
    const cc::vec5 FL(L.rho*unl, L.rho*L.u*unl + Lo.p*n.x, L.rho*L.v*unl + Lo.p*n.y,
                      L.rho*L.w*unl + Lo.p*n.z, unl*(L.rho*Lo.e + Lo.p));
    const cc::vec5 FR(R.rho*unr, R.rho*R.u*unr + Ro.p*n.x, R.rho*R.v*unr + Ro.p*n.y,
                      R.rho*R.w*unr + Ro.p*n.z, unr*(R.rho*Ro.e + Ro.p));
    cc::vec5 d{};
    for(int i=0;i<5;i++) {d += (std::abs(w.lambda[i])*w.alpha[i])*w.r[i];}
    face.convect = face.nor.norm()*(0.5*(FL + FR) - 0.5*d);
}