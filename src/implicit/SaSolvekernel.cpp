#include "implicit/Implicit.h"
#include "turbulence/Turbulence.h"

void SA::solvekernel(int c, int lo, int hi){
    cc::cell_class& cell = cc::CellList[c];
    double acc = nurhs[c];
    allface(cell){
        cc::face_class* f = cell.faces[i];
        if(f->type != cc::INTER) {continue;}
        const cc::cell_class* nb = cc::othercell(&cell, f);
        const int j = nb->index - 1;
        const double dj = (j >= lo && j < hi) ? dnu[j] : dnuold[j];
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0)*f->convect.c;
        acc += (std::min(mdot, 0.0) - Gamma_eff(cell, f)*f->nor.norm()/normal_dist(f))*dj;
    }
    dnu[c] = -acc/cell.tur->diag;
}
