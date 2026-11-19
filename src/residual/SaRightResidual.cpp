#include "residual/Residual.h"

#include "turbulence/Turbulence.h"

double SA::right_residual(cc::cell_class& cell){
    double res = 0.0;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const cc::cell_class* nb = cc::othercell(&cell, f);
        const double nu_o = nb ? nb->tur->miubl : f->tur->miubl;
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0)*f->convect.c;
        res += mdot*(mdot > 0.0 ? cell.tur->miubl : nu_o);
        res += Gamma_eff(cell, f)*f->nor.norm()/normal_dist(f)*(cell.tur->miubl - nu_o);
    }
    res -= cell.vol*source_SA(cell);
    return res;
}
