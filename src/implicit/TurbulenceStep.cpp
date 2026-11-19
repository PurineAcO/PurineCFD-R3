#include "implicit/Implicit.h"
#include "residual/Residual.h"
#include "solver/Solver.h"

void SA::turbulence_step(double dt, double cfl){
    impl::right_residual();
#pragma omp parallel for schedule(static)
    allcell impl::form_cell_lambda(icell(i));
    SA::malloc_vector();
    const double inv_dt = dt > 0.0 ? 1.0/dt : 0.0;
#pragma omp parallel for schedule(static)
    for(int i = 0; i < cc::cell_num; i++){
        SA::nurhs[i] = SA::right_residual(cc::CellList[i]) +
                       cc::CellList[i].vol*inv_dt*
                       (cc::CellList[i].conser.c*cc::CellList[i].tur->miubl -
                        cc::CellList[i].lastconser.c*cc::CellList[i].tur->miubl_former);
    }
#pragma omp parallel for schedule(static)
    allcell SA::source_jacobi_diag_SA(icell(i), inv_dt, cfl);
#pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)SA::dnu.size(); i++) {SA::dnu[i] = 0.0;}
    for(int sw = 0; sw < solver::sweeps; sw++){SA::nucopy(); blockSgs(SA::solvekernel);}
    SA::relaxer();
}
