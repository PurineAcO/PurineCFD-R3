#include "physic.h"

double ph::dynamic_pressure(const cc::vecp& refp) {
    return 0.5 * refp.rho * (refp.u * refp.u + refp.v * refp.v + refp.w * refp.w);
}

void ph::wall_forces(double& cl, double& cd, const cc::vecp& refp) {
    const double q = dynamic_pressure(refp);
    const double a = std::atan2(refp.v, refp.u);
    const cc::vec3 dir{std::cos(a), std::sin(a), 0.0};
    double fx = 0.0, fy = 0.0, projected = 0.0;
    for(const cc::bound& bd : cc::bounds){
        if(bd.type != cc::WALL) {continue;}
        for(const cc::face_class* w : bd.fs){
            const cc::cell_class* c = cc::boundary_findcell(const_cast<cc::face_class*>(w));
            const double s = cc::dot(w->nor, w->mid - c->center) > 0 ? 1.0 : -1.0;
            fx += s * (w->otphy->p * w->nor.x - w->visflux.x);
            fy += s * (w->otphy->p * w->nor.y - w->visflux.y);
            projected += std::abs(cc::dot(w->nor, dir));
        }
    }
    const double aref = 0.5*projected;
    const double denom = q*aref > 0.0 ? q*aref : 1.0;
    cl = (-fx * std::sin(a) + fy * std::cos(a)) / denom;
    cd = (fx * std::cos(a) + fy * std::sin(a)) / denom;
}