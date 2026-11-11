#include "euler.h"
#include "physic/physic.h"

roe_waves roe_decompose(const cc::face_class& face){
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const cc::otphy& Lo = *face.otnei[0];
    const cc::otphy& Ro = *face.otnei[1];
    roe_waves w;
    const double len0 = face.nor.norm();
    const double sgn = face.nei[0] != nullptr
        ? (cc::dot(face.nor, face.mid - face.nei[0]->center) > 0.0 ? 1.0 : -1.0)
        : (cc::dot(face.nor, face.nei[1]->center - face.mid) > 0.0 ? 1.0 : -1.0);
    w.n = (sgn/len0)*face.nor;    
    const double sl = std::sqrt(L.rho), sr = std::sqrt(R.rho);
    const double den = sl + sr;
    const double Hl = ph::Cp*L.T + 0.5*(L.u*L.u + L.v*L.v + L.w*L.w);
    const double Hr = ph::Cp*R.T + 0.5*(R.u*R.u + R.v*R.v + R.w*R.w);
    const double rhob = sl*sr;
    const cc::vec3 vl{L.u,L.v,L.w}, vr{R.u,R.v,R.w};
    const cc::vec3 vb = (1.0/den)*(sl*vl + sr*vr);
    const double Hb = (sl*Hl + sr*Hr)/den;
    const double qb2 = cc::dot(vb, vb);
    const double ab = std::sqrt((ph::gamma - 1.0)*(Hb - 0.5*qb2));
    const double unb = cc::dot(vb, w.n);
    cc::vec3 t1 = vb - unb*w.n;
    if(t1.norm() > 1e-300){
        t1 = (1.0/t1.norm())*t1;
    }else{
        const cc::vec3 seed = std::abs(w.n.x) < 0.9 ? cc::vec3(1.0,0.0,0.0) : cc::vec3(0.0,1.0,0.0);
        t1 = cc::cross(w.n, seed);
        t1 = (1.0/t1.norm())*t1;
    }
    const cc::vec3 t2 = cc::cross(w.n, t1);
    const double dp = Ro.p - Lo.p;
    const double dun = cc::dot(vr - vl, w.n);
    w.lambda[0] = unb - ab;
    w.lambda[1] = unb;
    w.lambda[2] = unb;
    w.lambda[3] = unb;
    w.lambda[4] = unb + ab;
    w.alpha[0] = (dp - rhob*ab*dun)/(2.0*ab*ab);
    w.alpha[1] = R.rho - L.rho - dp/(ab*ab);
    w.alpha[2] = rhob*cc::dot(vr - vl, t1);
    w.alpha[3] = rhob*cc::dot(vr - vl, t2);
    w.alpha[4] = (dp + rhob*ab*dun)/(2.0*ab*ab);
    const cc::vec3 vm = vb - ab*w.n, vp = vb + ab*w.n;
    w.r[0] = cc::vec5(1.0, vm.x, vm.y, vm.z, Hb - unb*ab);
    w.r[1] = cc::vec5(1.0, vb.x, vb.y, vb.z, 0.5*qb2);
    w.r[2] = cc::vec5(0.0, t1.x, t1.y, t1.z, cc::dot(t1, vb));
    w.r[3] = cc::vec5(0.0, t2.x, t2.y, t2.z, cc::dot(t2, vb));
    w.r[4] = cc::vec5(1.0, vp.x, vp.y, vp.z, Hb + unb*ab);
    return w;
}