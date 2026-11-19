#include "solver/Solver.h"

void blockSgs(SolveKernel solve_kernel){
#pragma omp parallel
    {
#ifdef _OPENMP
        const int nt = omp_get_num_threads(), tid = omp_get_thread_num();
#else
        const int nt = 1, tid = 0;
#endif
        const int lo = (int)((long long)cc::cell_num*tid/nt);
        const int hi = (int)((long long)cc::cell_num*(tid + 1)/nt);
        for(int c = lo; c < hi; c++) {solve_kernel(c, lo, hi);}
        for(int c = hi - 1; c >= lo; c--) {solve_kernel(c, lo, hi);}
    }
}
