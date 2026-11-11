#include "implict.h"
#include "solver/LUSGS.h"
#include "turbulence/SA.h"

double impl::flow_step(double dt, double cfl){
    malloc_vector();
    clear_mat();
    right_residual();
    const double inv_dt = dt > 0.0 ? 1.0/dt : 0.0;
#pragma omp parallel for schedule(static)
    allcell form_euler_jacobi(icell(i));
#pragma omp parallel for schedule(static)
    allfac form_viscous_jacobi(iface(i));
#pragma omp parallel for schedule(static)
    allcell form_diag(icell(i), inv_dt, cfl);
    double res = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : res)
    for(int i = 0; i < cc::cell_num; i++){
        rhs[i] = cc::CellList[i].convect - cc::CellList[i].visflux +
                 (cc::CellList[i].vol*inv_dt)*(cc::CellList[i].conser - cc::CellList[i].lastconser);
        dQ[i].clear();
        res += rhs[i].c*rhs[i].c;
    }
    for(int sw = 0; sw < solver::sweeps; sw++){Qcopy(); block_sgs(solvekernel);}
    relax();
    return std::sqrt(res/cc::cell_num);
}

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
    for(int sw = 0; sw < solver::sweeps; sw++){SA::nucopy(); block_sgs(SA::solvekernel);}
    SA::relaxer();
}