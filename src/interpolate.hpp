#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// MUSCL左右插值
void muscl(cc::face_class* face);
// 面上梯度的处理
void grad_onface(cc::cell_class& cell);
// Venkatakrishnan限制器
cc::vecp limited_phy(const cc::cell_class* cell, const cc::vec3& off);
// 壁面法向修正
cc::vec3 wall_correct(cc::vec3 grad, double wall, double inner, const cc::vec3& n, double d);

inline void interpolate_mid(cc::face_class* face){
    if(face->type != cc::INTER) {return;}
    face->phy = 0.5*(face->nei[0]->phy + face->nei[1]->phy);
    face->tur.miubl = 0.5*(face->nei[0]->tur.miubl + face->nei[1]->tur.miubl);
    face->otphy = form_otphy(face->phy);
    const cc::vec3 v{face->phy.u, face->phy.v, face->phy.w};
    face->otphy.un = cc::dot(v, face->nor)/face->nor.norm();
}

inline cc::vec3 wall_correct(cc::vec3 grad, double wall, double inner, const cc::vec3& n, double d){
    const double dn = (wall - inner)/d;
    return grad + (dn - cc::dot(grad, n))*n;
}

inline cc::vecp limited_phy(const cc::cell_class* cell, const cc::vec3& off){
    const double c[5] = {cell->phy.rho, cell->phy.u, cell->phy.v, cell->phy.w, cell->phy.T};
    double lo[5], hi[5];
    for(int k=0;k<5;k++) {lo[k] = c[k]; hi[k] = c[k];}
    allface(*cell){
        const cc::face_class* f = cell->faces[i];
        const cc::cell_class* o = (f->nei[0] == cell) ? f->nei[1] : f->nei[0];
        if(o == nullptr) {continue;}
        const double n[5] = {o->phy.rho, o->phy.u, o->phy.v, o->phy.w, o->phy.T};
        for(int k=0;k<5;k++){
            if(n[k] < lo[k]) {lo[k] = n[k];}
            if(n[k] > hi[k]) {hi[k] = n[k];}
        }
    }
    const cc::vecp d = cell->phgrad*off;
    const double dk[5] = {d.rho, d.u, d.v, d.w, d.T};
    double r[5];
    for(int k=0;k<5;k++){
        r[k] = c[k];
        const double dd = std::abs(dk[k]);
        if(dd <= 0.0) {continue;}
        const double dm = dk[k] > 0.0 ? hi[k] - c[k] : c[k] - lo[k];
        const double eps = 1e-12*(c[k]*c[k] + 1e-30);
        const double psi = (dm*dm + 2.0*dm*dd + eps)/(dm*dm + dm*dd + 2.0*dd*dd + eps);
        r[k] = c[k] + psi*dk[k];
    }
    return cc::vecp(r[0], r[1], r[2], r[3], r[4]);
}

inline void muscl(cc::face_class* face){
    if(face->type == cc::INTER){
        face->phynei[0] = limited_phy(face->nei[0], face->mid - face->nei[0]->center);
        face->phynei[1] = limited_phy(face->nei[1], face->mid - face->nei[1]->center);
    }else{
        face->phynei[0] = face->phy;
        face->phynei[1] = face->phy;
    }
    face->otnei[0] = form_otphy(face->phynei[0]);
    face->otnei[1] = form_otphy(face->phynei[1]);
}

inline void grad_onface(cc::cell_class& cell){
    allface(cell){
        cc::face_class* face = cell.faces[i];
        if(face->type == cc::INTER){
            const cc::cell_class* lo = face->nei[0];
            const cc::cell_class* hi = face->nei[1];
            const double dlo = cc::dot(face->nor, face->mid - lo->center);
            const double dhi = cc::dot(face->nor, hi->center - face->mid);
            const double w = dlo/(dhi + dlo);
            face->phgrad = w*lo->phgrad + (1.0 - w)*hi->phgrad;
            face->tur.miublgrad = w*lo->tur.miublgrad + (1.0 - w)*hi->tur.miublgrad;
        }else if(face->type == cc::WALL){
            cc::vec3 n = (1.0/face->nor.norm())*face->nor;
            if(cc::dot(n, face->mid - cell.center) < 0.0) {n = -1.0*n;}
            const double d = cc::dot(face->mid - cell.center, n);
            face->phgrad = cc::vecgrad{};
            face->phgrad.rhograd = wall_correct(cell.phgrad.rhograd, cell.phy.rho, cell.phy.rho, n, d);
            face->phgrad.ugrad = wall_correct(cell.phgrad.ugrad, 0.0, cell.phy.u, n, d);
            face->phgrad.vgrad = wall_correct(cell.phgrad.vgrad, 0.0, cell.phy.v, n, d);
            face->phgrad.wgrad = wall_correct(cell.phgrad.wgrad, 0.0, cell.phy.w, n, d);
            face->phgrad.Tgrad = wall_correct(cell.phgrad.Tgrad, cell.phy.T, cell.phy.T, n, d);
            face->tur.miublgrad = wall_correct(cell.tur.miublgrad, 0.0, cell.tur.miubl, n, d);
        }else{
            face->phgrad = cell.phgrad;
            face->tur.miublgrad = cell.tur.miublgrad;
        }
    }
}
