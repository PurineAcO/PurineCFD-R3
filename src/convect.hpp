#pragma once

#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "physic/physic.h"
#include "dissipation.hpp"
#include "udf/udf.h"
#include <algorithm>
#include <cmath>

// Roe波的分解量
struct roe_waves {
    cc::vec3 n;
    double lambda[5] = {};
    double alpha[5] = {};
    cc::vec5 r[5];
};

// JST无粘通量
void convect_JST(cc::face_class& face);
// 汇总单元各面的无粘和黏性通量
void assemble_flux(cc::cell_class& cell, char fluxtype = 'R');
// Roe无粘通量
void convect_ROE(cc::face_class& face);
// Roe分解
roe_waves roe_decompose(const cc::face_class& face);

inline void convect_JST(cc::face_class& face){
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const double H = ph::Cp*face.phy.T + 0.5*(u*u + v*v + w*w);
    const double p = face.otphy->p;
    const cc::vec5 F(rho*u, rho*u*u + p, rho*u*v, rho*u*w, rho*u*H);
    const cc::vec5 G(rho*v, rho*u*v, rho*v*v + p, rho*v*w, rho*v*H);
    const cc::vec5 K(rho*w, rho*u*w, rho*v*w, rho*w*w + p, rho*w*H);
    face.convect = face.nor.x*F + face.nor.y*G + face.nor.z*K;
    jst::dissipation(face);
}

inline void assemble_flux(cc::cell_class& cell, char fluxtype){
    cell.convect.clear();
    for(int i=0;i<cell.fcnt;i++){
        cc::face_class* face = cell.faces[i];
        if(fluxtype == 'J') {
            cell.convect += (2*cell.fnorm[i] - 1)*face->convect;
        }else{
            const bool islow = (face->nei[0] == &cell);
            cell.convect += (2*islow - 1)*face->convect;
        }
    }
}

inline roe_waves roe_decompose(const cc::face_class& face){
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const cc::otphy& Lo = *face.otnei[0];
    const cc::otphy& Ro = *face.otnei[1];
    roe_waves w;
    const double len0 = face.nor.norm();
    const double sgn = face.nei[0] != nullptr
        ? (cc::dot(face.nor, face.mid - face.nei[0]->center) > 0.0 ? 1.0 : -1.0)
        : (cc::dot(face.nor, face.nei[1]->center - face.mid) > 0.0 ? 1.0 : -1.0);
    w.n = (sgn/len0)*face.nor;    const double sl = std::sqrt(L.rho), sr = std::sqrt(R.rho);
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

inline void convect_ROE(cc::face_class& face){
    const cc::vecp& L = face.phynei[0];
    const cc::vecp& R = face.phynei[1];
    const cc::otphy& Lo = *face.otnei[0];
    const cc::otphy& Ro = *face.otnei[1];
    const roe_waves w = roe_decompose(face);
    const cc::vec3& n = w.n;
    const double unl = cc::dot(cc::vec3{L.u,L.v,L.w}, n);
    const double unr = cc::dot(cc::vec3{R.u,R.v,R.w}, n);
    const cc::vec5 FL(L.rho*unl, L.rho*L.u*unl + Lo.p*n.x, L.rho*L.v*unl + Lo.p*n.y,
                      L.rho*L.w*unl + Lo.p*n.z, unl*(L.rho*Lo.e + Lo.p));
    const cc::vec5 FR(R.rho*unr, R.rho*R.u*unr + Ro.p*n.x, R.rho*R.v*unr + Ro.p*n.y,
                      R.rho*R.w*unr + Ro.p*n.z, unr*(R.rho*Ro.e + Ro.p));
    cc::vec5 d{};
    for(int i=0;i<5;i++) {d += (std::abs(w.lambda[i])*w.alpha[i])*w.r[i];}
    face.convect = face.nor.norm()*(0.5*(FL + FR) - 0.5*d);
}
