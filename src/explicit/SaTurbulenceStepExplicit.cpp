#include "explicit/Explicit.h"

void SA::turbulence_step_explicit(double dt){
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
    allcell SA::source_jacobi_diag_SA(icell(i), inv_dt, solver::exp_cfl);
#pragma omp parallel for schedule(static)
    for(int i = 0; i < cc::cell_num; i++){
        const double d = cc::CellList[i].tur->diag;
        SA::dnu[i] = d > 0.0 ? -SA::nurhs[i]/d : 0.0;
    }
    SA::relaxer();
}
