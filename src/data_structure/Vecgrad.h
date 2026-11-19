#pragma once

#include "data_structure/Vec3.h"
#include "data_structure/Vecp.h"

namespace cc {

struct vecgrad {
    vec3 rhograd, ugrad, vgrad, wgrad, Tgrad;
    vecgrad() = default;
    vecgrad(vec3 rho_, vec3 u_, vec3 v_, vec3 w_, vec3 t_): rhograd(rho_), ugrad(u_), vgrad(v_), wgrad(w_), Tgrad(t_) {}
    void clear();
};

// 梯度点乘位移 → 物理量增量
vecp operator*(const vecgrad& a, const vec3& r);
// 物理量乘面法向 → 梯度算子
vecgrad operator*(const vecp& phyf, const vec3& fnor);
vecgrad operator*(double a, const vecgrad& grad);
vecgrad operator+(const vecgrad& a, const vecgrad& b);

inline void vecgrad::clear() {rhograd.clear();ugrad.clear();vgrad.clear();wgrad.clear();Tgrad.clear();}

}
