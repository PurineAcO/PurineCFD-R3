#include "explicit/Explicit.h"

double impl::flow_step_explicit(double dt){
    malloc_vector();
    right_residual();
    const double inv_dt = dt > 0.0 ? 1.0/dt : 0.0;
#pragma omp parallel for schedule(static)
    allcell form_cell_lambda(icell(i));
    double res = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : res)
    for(int i = 0; i < cc::cell_num; i++){
        const cc::cell_class& c = cc::CellList[i];
        const cc::vec5 R = c.convect - c.visflux + (c.vol*inv_dt)*(c.conser - c.lastconser);
        dQ[i] = (c.lambda > 0.0 ? -(solver::exp_cfl/c.lambda) : 0.0)*R;
        res += R.c*R.c;
    }
    relax();
    return std::sqrt(res/cc::cell_num);
}
