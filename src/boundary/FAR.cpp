#include "FAR.h"
#include "initialize/initialization.h"
#include "physic/physic.h"
#include "geometry/geometry.h"

void far_field_boundary(){
    const cc::vecp farp = from_boundary_to_phy();
    const cc::otphy faro = ph::form_otphy(farp);
    const double miubl_inf = cc::FAR_config.chi*faro.mu/farp.rho;
    for(cc::bound& b : cc::bounds){
        if(b.type != cc::FAR) {continue;}
        for(cc::face_class* far : b.fs){
            cc::cell_class* c = cc::boundary_findcell(far);
            const cc::vec3 n = outward_unit_nor(far, c)*cc::unit(far->nor);
            cc::vec3 t1{}, t2{};
            tangent_basis(n, t1, t2);
            const cc::vec3 v = c->phy.vel();
            const cc::vec3 vinf = farp.vel();
            const double a = c->otphy->a;
            const double vn = cc::dot(v, n);
            const double vn_inf = cc::dot(vinf, n);
            const double Rp = vn + 2.0*a/(ph::gamma - 1.0);
            const double Rm = vn_inf - 2.0*faro.a/(ph::gamma - 1.0);
            const double vn_star = 0.5*(Rp + Rm);
            const double a_star = 0.25*(ph::gamma - 1.0)*(Rp - Rm);
            const bool outflow = vn_star >= 0.0;
            const double s = outflow ? c->otphy->p/std::pow(c->phy.rho, ph::gamma)
                                     : cc::FAR_config.p/std::pow(farp.rho, ph::gamma);
            const cc::vec3 vt = outflow ? v : vinf;
            const cc::vec3 v_star = vn_star*n + cc::dot(vt,t1)*t1 + cc::dot(vt,t2)*t2;
            far->phy.rho = std::pow(a_star*a_star/(ph::gamma*s), 1.0/(ph::gamma - 1.0));
            far->phy.setvel(v_star);
            far->otphy->p = s*std::pow(far->phy.rho, ph::gamma);
            far->phy.T = far->otphy->p/(ph::R*far->phy.rho);
            far->tur->miubl = outflow ? c->tur->miubl : miubl_inf;
            *far->otphy = ph::form_otphy(far->phy);
        }
    }
}