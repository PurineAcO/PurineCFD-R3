#pragma once

#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "physic/physic.h"
#include "udf/udf.h"
#include <cmath>

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
    face->tur->miubl = 0.5*(face->nei[0]->tur->miubl + face->nei[1]->tur->miubl);
    *face->otphy = ph::form_otphy(face->phy);
    const cc::vec3 v{face->phy.u, face->phy.v, face->phy.w};
    face->otphy->un = cc::dot(v, face->nor)/face->nor.norm();
}

inline cc::vec3 wall_correct(cc::vec3 grad, double wall, double inner, const cc::vec3& n, double d){
    const double dn = (wall - inner)/d;
    return grad + (dn - cc::dot(grad, n))*n;
}

inline cc::vecp limited_phy(const cc::cell_class* cell, const cc::vec3& off){
    const double c[5] = {cell->phy.rho, cell->phy.u, cell->phy.v, cell->phy.w, cell->phy.T};
    double lo[5] = {c[0], c[1], c[2], c[3], c[4]};
    double hi[5] = {c[0], c[1], c[2], c[3], c[4]};
    allface(*cell){
        const cc::cell_class* o = cell->nei[i];
        if(o == nullptr) {continue;}
        const double n[5] = {o->phy.rho, o->phy.u, o->phy.v, o->phy.w, o->phy.T};
        for(int k=0;k<5;k++){
            if(n[k] < lo[k]) {lo[k] = n[k];}
            if(n[k] > hi[k]) {hi[k] = n[k];}
        }
    }
    const cc::vecp d = cell->phgrad*off;
    const double dk[5] = {d.rho, d.u, d.v, d.w, d.T};
    double r[5] = {c[0], c[1], c[2], c[3], c[4]};
    for(int k=0;k<5;k++){
        double psi = 1.0;
        if(dk[k] > 0.0){
            const double m = (hi[k] - c[k])/dk[k];
            psi = m < 1.0 ? m : 1.0;
        }else if(dk[k] < 0.0){
            const double m = (lo[k] - c[k])/dk[k];
            psi = m < 1.0 ? m : 1.0;
        }
        psi = psi > 0.0 ? psi : 0.0;
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
    *face->otnei[0] = ph::form_otphy(face->phynei[0]);
    *face->otnei[1] = ph::form_otphy(face->phynei[1]);
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
            face->tur->miublgrad = w*lo->tur->miublgrad + (1.0 - w)*hi->tur->miublgrad;
        }else if(face->type == cc::WALL){
            const double y1 = cell.tur->sad;
            cc::vec3 n = face->mid - cell.center;
            n = (1.0/n.norm())*n;
            const cc::vec3 seed = std::abs(n.x) < 0.9 ? cc::vec3(1.0,0.0,0.0) : cc::vec3(0.0,1.0,0.0);
            cc::vec3 t1 = cc::cross(n, seed);
            t1 = (1.0/t1.norm())*t1;
            const cc::vec3 t2 = cc::cross(n, t1);
            auto split = [&](const cc::vec3& g, double wall, double inner){
                const double dn = (wall - inner)/y1;
                return dn*n + cc::dot(g, t1)*t1 + cc::dot(g, t2)*t2;
            };
            face->phgrad.clear();
            face->phgrad.ugrad = split(cell.phgrad.ugrad, face->phy.u, cell.phy.u);
            face->phgrad.vgrad = split(cell.phgrad.vgrad, face->phy.v, cell.phy.v);
            face->phgrad.wgrad = split(cell.phgrad.wgrad, face->phy.w, cell.phy.w);
            face->phgrad.Tgrad = split(cell.phgrad.Tgrad, face->phy.T, cell.phy.T);
            face->tur->miublgrad = split(cell.tur->miublgrad, face->tur->miubl, cell.tur->miubl);
        }else{
            face->phgrad = cell.phgrad;
            face->tur->miublgrad = cell.tur->miublgrad;
        }
    }
}
