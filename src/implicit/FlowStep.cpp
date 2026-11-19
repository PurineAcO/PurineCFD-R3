#include "implicit/Implicit.h"
#include "residual/Residual.h"
#include "solver/Solver.h"

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

    for(int sw = 0; sw < solver::sweeps; sw++){Qcopy(); blockSgs(impl::solvekernel);}
    relax();
    return std::sqrt(res/cc::cell_num);
}
