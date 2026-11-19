#pragma once

#include "data_structure/Vec3.h"

namespace cc {

struct vecp {
    double rho = 0.0, u = 0.0, v = 0.0, w = 0.0, T = 0.0;
    vecp() = default;
    vecp(double rho_, double u_, double v_, double w_, double T_): rho(rho_), u(u_), v(v_), w(w_), T(T_) {}
    vecp& operator+=(const vecp& o) {rho += o.rho;u += o.u;v += o.v;w += o.w;T += o.T;return *this;}
    vecp& operator-=(const vecp& o) {rho -= o.rho;u -= o.u;v -= o.v;w -= o.w;T -= o.T;return *this;}
    vecp& operator*=(double k) {rho *= k;u *= k;v *= k;w *= k;T *= k;return *this;}
    double& operator[](int i);
    double operator[](int i) const;
    vec3 vel() const {return vec3{u, v, w};}
    void setvel(const vec3& q) {u = q.x;v = q.y;w = q.z;}
    void clear();
};

vecp operator+(vecp a, const vecp& b);
vecp operator-(vecp a, const vecp& b);
vecp operator*(vecp a, double s);
vecp operator*(double s, vecp a);

inline void vecp::clear() {rho = 0.0;u = 0.0;v = 0.0;w = 0.0;T = 0.0;}

}
