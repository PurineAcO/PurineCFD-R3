#pragma once

#include "data_structure/Vec5.h"

namespace cc {

struct mat5 {
    double a[5][5] = {};
    mat5() = default;
    mat5(const vec5& r0, const vec5& r1, const vec5& r2, const vec5& r3, const vec5& r4);
    mat5& operator+=(const mat5& o);
    mat5& operator*=(double k);
    double* operator[](int r);
    const double* operator[](int r) const;
    void clear();
    // 高斯消元求解 mat5 * x = b, 解回写 b; 主元下限 1e-14
    bool solve(vec5& b) const;
};

inline mat5::mat5(const vec5& r0, const vec5& r1, const vec5& r2, const vec5& r3, const vec5& r4){
    const vec5* r[5] = {&r0, &r1, &r2, &r3, &r4};
    for(int i = 0; i < 5; i++){
        a[i][0] = r[i]->c; a[i][1] = r[i]->x; a[i][2] = r[i]->y;
        a[i][3] = r[i]->z; a[i][4] = r[i]->e;
    }
}

// 矩阵加法
mat5 operator+(mat5 m, const mat5& o);
// 矩阵数乘
mat5 operator*(mat5 m, double k);
// 矩阵减法
mat5 operator-(mat5 m, const mat5& o);
// 矩阵数乘
mat5 operator*(double k, mat5 m);
// 矩阵乘向量
vec5 operator*(const mat5& m, const vec5& v);

inline mat5& mat5::operator+=(const mat5& o) {
    for (int r = 0; r < 5; r++)
        for (int c = 0; c < 5; c++) a[r][c] += o.a[r][c];
    return *this;
}

inline mat5& mat5::operator*=(double k) {
    for (int r = 0; r < 5; r++)
        for (int c = 0; c < 5; c++) a[r][c] *= k;
    return *this;
}

inline void mat5::clear() {
    for (int r = 0; r < 5; r++)
        for (int c = 0; c < 5; c++) a[r][c] = 0.0;
}

inline mat5 eyes = []{ mat5 I; for(int r = 0; r < 5; r++) I.a[r][r] = 1.0; return I; }();

}
