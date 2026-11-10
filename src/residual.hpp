#pragma once

#include "data_structure/cell_class.h"
#include "data_structure/classconfig.h"
#include <algorithm>
#include <cmath>

namespace res {
inline constexpr double drop_target = 1e-4;
inline int conv_step = 20;
inline double norm0[5] = {};
inline double current[5] = {};
}

// 计算当前残差
void calculate_residual();
// 残差下降量
double residual_drop();
// 输出残差
void report_residual(int step, double cl, double cd, double elapsed);

inline void calculate_residual(){
    double sum[5] = {};
#pragma omp parallel for schedule(static) reduction(+ : sum[:5])
    for(int i = 0; i < cc::cell_num; i++){
        const cc::cell_class& c = cc::CellList[i];
        const double inv = 1.0/c.vol;
        const double r[5] = {(c.convect.c - c.visflux.c)*inv, (c.convect.x - c.visflux.x)*inv,
                             (c.convect.y - c.visflux.y)*inv, (c.convect.z - c.visflux.z)*inv,
                             (c.convect.e - c.visflux.e)*inv};
        for(int k = 0; k < 5; k++) sum[k] += r[k]*r[k];
    }
    for(int k = 0; k < 5; k++){
        res::current[k] = std::sqrt(sum[k]/cc::cell_num);
        if(res::norm0[k] == 0.0) res::norm0[k] = res::current[k];
    }
}

inline double residual_drop(){
    double cur2 = 0.0, ref2 = 0.0;
    for(int k = 0; k < 5; k++){
        cur2 += res::current[k]*res::current[k];
        ref2 += res::norm0[k]*res::norm0[k];
    }
    if(ref2 <= 0.0) {return 0.0;}
    return std::sqrt(cur2/ref2);
}

inline void report_residual(int step, double cl, double cd, double elapsed){
    printf("step %6d  res %.4e  Cl %+.5f  Cd %+.5f  [%.0fs]\n",
           step, residual_drop(), cl, cd, elapsed);
    fflush(stdout);
}
