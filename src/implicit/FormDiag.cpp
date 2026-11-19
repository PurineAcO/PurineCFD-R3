#include "implicit/Implicit.h"
#include "boundary/Boundary.h"

void impl::form_diag(cc::cell_class& cell, double inv_dt, double cfl){
    form_cell_lambda(cell);
    cell.diag.clear();
    cell.diag += (cell.vol*inv_dt + cell.lambda/cfl)*cc::eyes;
    cc::mat5 far;
    allface(cell){
        const double s = 2*cell.fnorm[i] - 1;
        cc::face_class* f = cell.faces[i];
        if(f->type == cc::INTER) {cell.diag += 0.5*s*f->toface_jacobi(cell.conj) + 0.5*f->imats.Ab
                                              + f->imats.Tf + (-0.5*s)*f->imats.Bf;}
        else if(f->type == cc::WALL){wall_jacobi(f, &cell, s);}
        else if(f->type == cc::SYM){sym_jacobi(f, &cell, s);}
        else if(f->type == cc::FAR){
            far_jacobi(f, &cell, s);
#if !chose2D
            far += far_convect_jacobi(f);
#endif
        }
    }
    cell.diag += far;
}
