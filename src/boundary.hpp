#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "initialize.hpp"
#include "physic.hpp"
#include <cmath>

// 无滑移不穿透壁面
void noslip_wall_boundary();
// 对称滑移面
void symmetry_wall_boundary();
// 亚声速特征远场
void far_field_boundary();
// 面单位法向并由单元指向外
cc::vec3 outward_unit_nor(cc::face_class* f, const cc::cell_class* c);
// 由法向构造切向正交基
void tangent_basis(const cc::vec3& n, cc::vec3& t1, cc::vec3& t2);

inline cc::vec3 outward_unit_nor(cc::face_class* f, const cc::cell_class* c){
    cc::vec3 n = (1.0/f->nor.norm())*f->nor;
    if(cc::dot(n, f->mid - c->center) < 0.0) {n = -1.0*n;}
    return n;
}

inline void tangent_basis(const cc::vec3& n, cc::vec3& t1, cc::vec3& t2){
    const cc::vec3 seed = std::abs(n.x) < 0.9 ? cc::vec3(1.0,0.0,0.0) : cc::vec3(0.0,1.0,0.0);
    t1 = cc::cross(n, seed);
    t1 = (1.0/t1.norm())*t1;
    t2 = cc::cross(n, t1);
}

inline void noslip_wall_boundary(){
    for(cc::face_class* wall : cc::WallFaces){
        cc::cell_class* c = cc::boundary_findcell(wall);
        wall->phy.u = 0.0;
        wall->phy.v = 0.0;
        wall->phy.w = 0.0;
        wall->phy.T = c->phy.T;
        wall->otphy.p = c->otphy.p;
        wall->phy.rho = wall->otphy.p/(cc::R*wall->phy.T);
        wall->tur.miubl = 0.0;
        wall->otphy = form_otphy(wall->phy);
    }
}

inline void symmetry_wall_boundary(){
    for(cc::bound& b : cc::bounds){
        if(b.type != cc::SYM) {continue;}
        for(cc::face_class* sym : b.fs){
            cc::cell_class* c = cc::boundary_findcell(sym);
            const cc::vec3 n = outward_unit_nor(sym, c);
            const cc::vec3 v{c->phy.u, c->phy.v, c->phy.w};
            const cc::vec3 vs = v - (2.0*cc::dot(v, n))*n;
            sym->phy.rho = c->phy.rho;
            sym->phy.u = vs.x;
            sym->phy.v = vs.y;
            sym->phy.w = vs.z;
            sym->phy.T = c->phy.T;
            sym->tur.miubl = c->tur.miubl;
            sym->otphy = form_otphy(sym->phy);
        }
    }
}

inline void far_field_boundary(){
    const cc::vecp farp = from_boundary_to_phy();
    const cc::otphy faro = form_otphy(farp);
    const double miubl_inf = cc::FAR_config.chi*faro.mu/farp.rho;
    for(cc::face_class* far : cc::FarFaces){
        cc::cell_class* c = cc::boundary_findcell(far);
        const cc::vec3 n = outward_unit_nor(far, c);
        cc::vec3 t1{}, t2{};
        tangent_basis(n, t1, t2);
        const cc::vec3 v{c->phy.u, c->phy.v, c->phy.w};
        const double a = c->otphy.a;
        const double vn = cc::dot(v, n);
        const double vn_inf = cc::dot(cc::vec3{farp.u,farp.v,farp.w}, n);
        const double Rp = vn + 2.0*a/(cc::gamma - 1.0);
        const double Rm = vn_inf - 2.0*faro.a/(cc::gamma - 1.0);
        const double vn_star = 0.5*(Rp + Rm);
        const double a_star = 0.25*(cc::gamma - 1.0)*(Rp - Rm);
        const bool outflow = vn_star >= 0.0;
        const double s = outflow ? c->otphy.p/std::pow(c->phy.rho, cc::gamma)
                                 : cc::FAR_config.p/std::pow(farp.rho, cc::gamma);
        const cc::vec3 vt = outflow ? v : cc::vec3{farp.u,farp.v,farp.w};
        const cc::vec3 v_star = vn_star*n + cc::dot(vt,t1)*t1 + cc::dot(vt,t2)*t2;
        far->phy.rho = std::pow(a_star*a_star/(cc::gamma*s), 1.0/(cc::gamma - 1.0));
        far->phy.u = v_star.x;
        far->phy.v = v_star.y;
        far->phy.w = v_star.z;
        far->otphy.p = s*std::pow(far->phy.rho, cc::gamma);
        far->phy.T = far->otphy.p/(cc::R*far->phy.rho);
        far->tur.miubl = outflow ? c->tur.miubl : miubl_inf;
        far->otphy = form_otphy(far->phy);
    }
}
