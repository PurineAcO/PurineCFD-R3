#include "physic.h"

double ph::dynamic_pressure(const cc::vecp& refp) {
    return 0.5 * refp.rho * (refp.u * refp.u + refp.v * refp.v + refp.w * refp.w);
}

void ph::wall_forces(double& cl, double& cd, const cc::vecp& refp) {
    const double q = dynamic_pressure(refp);
    double fx = 0.0, fy = 0.0;
    cc::bound* wall;
    for(cc::bound bd : cc::bounds) if(bd.type == cc::WALL) {wall = &bd;break;}
    for (cc::face_class* w : wall->fs) {
        const cc::cell_class* c = cc::boundary_findcell(w);
        const double s = cc::dot(w->nor, w->mid - c->center) > 0 ? 1.0 : -1.0;
        fx += s * (w->otphy->p * w->nor.x - w->visflux.x);
        fy += s * (w->otphy->p * w->nor.y - w->visflux.y);
    }
    const double a = std::atan2(refp.v, refp.u);
    cl = (-fx * std::sin(a) + fy * std::cos(a)) / q;
    cd = (fx * std::cos(a) + fy * std::sin(a)) / q;
}