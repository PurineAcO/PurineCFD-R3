#pragma once

#include "classconfig.hpp"
#include <cmath>

// JST人工耗散
namespace jst {
inline constexpr double k2 = 0.5;
inline constexpr double k4 = 1.0/64;
// 激波(压力)检测器
void shockwave_recognize(cc::cell_class& cell);
// 单元的伪Laplace
void laplace_dissipation(cc::cell_class& cell);
// 面上耗散
void dissipation(cc::face_class& face);
}

inline void jst::shockwave_recognize(cc::cell_class& cell){
    double up = 0.0, down = 0.0;
    for(int i=0;i<cell.fcnt;i++){
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr) {continue;}
        up += std::abs(neighbor->otphy.p - cell.otphy.p);
        down += neighbor->otphy.p + cell.otphy.p;
    }
    cell.diss.Y = (down > 0.0) ? up/down : 0.0;
}

inline void jst::laplace_dissipation(cc::cell_class& cell){
    cc::vec5 L{};
    for(int i=0;i<cell.fcnt;i++){
        cc::cell_class* neighbor = cell.nei[i];
        if(neighbor == nullptr) {continue;}
        L += neighbor->conser - cell.conser;
    }
    cell.diss.L = L;
}

inline void jst::dissipation(cc::face_class& face){
    if(face.type != cc::INTER) {return;}
    const cc::cell_class* lo = face.nei[0];
    const cc::cell_class* hi = face.nei[1];
    const double lam = (std::abs(face.otphy.un) + face.otphy.a)*face.nor.norm();
    const double eps2 = k2*std::max(lo->diss.Y, hi->diss.Y);
    const double eps4 = std::max(0.0, k4 - eps2);
    const double k = lam*(2*face.outer - 1);
    face.convect += (-k*eps2)*(hi->conser - lo->conser);
    face.convect += (k*eps4)*(hi->diss.L - lo->diss.L);
}
