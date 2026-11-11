#include "implict.h"
#include "turbulence/SA.h"

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

void impl::solvekernel(int c, int lo, int hi){
    const cc::cell_class& cell = cc::CellList[c];
    cc::vec5 acc = rhs[c];
    allface(cell){
        const cc::face_class* f = cell.faces[i];
        if(f->type != cc::INTER) {continue;}
        const cc::cell_class* o = cc::othercell(const_cast<cc::cell_class*>(&cell),
                                                const_cast<cc::face_class*>(f));
        if(o == nullptr) {continue;}
        const int j = o->index - 1;
        if(j < 0 || j >= cc::cell_num) {continue;}
        const cc::vec5& dj = (j >= lo && j < hi) ? dQ[j] : dQold[j];
        const double s = 2*cell.fnorm[i] - 1;
        const cc::mat5 M = 0.5*s*f->toface_jacobi(o->conj) - 0.5*f->imats.Ab
                           - f->imats.Tf - (0.5*s)*f->imats.Bf;
        acc += M*dj;
    }
    acc *= -1.0;
    cell.diag.solve(acc);
    dQ[c] = acc;
}