#pragma once
#include <cmath>

// 向量和算子
namespace cc {

// 三维向量
struct vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    vec3() = default;
    vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    vec3& operator+=(const vec3& o) {x += o.x;y += o.y;z += o.z;return *this;}
    vec3& operator-=(const vec3& o) {x -= o.x;y -= o.y;z -= o.z;return *this;}
    vec3& operator*=(double s) {x *= s;y *= s;z *= s;return *this;}
    void clear();
    double norm() const;
};

inline void vec3::clear() {x = 0.0;y = 0.0; z = 0.0;}
inline vec3 operator+(vec3 a, const vec3& b) {return a += b;}
inline vec3 operator-(vec3 a, const vec3& b) {return a -= b;}
inline vec3 operator*(vec3 a, double s) {return a *= s;}
inline vec3 operator*(double s, vec3 a) {return a *= s;}
inline vec3 unit(const vec3& a) {return (1.0/a.norm())*a;}
inline double dot(const vec3& a, const vec3& b) {return a.x * b.x + a.y * b.y + a.z * b.z;}
inline vec3 cross(const vec3& a, const vec3& b) {return vec3{a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
inline double vec3::norm() const {return std::sqrt(x * x + y * y + z * z);}

// 原始物理量
struct vecp {
    double rho = 0.0, u = 0.0, v = 0.0, w = 0.0, T = 0.0;
    vecp() = default;
    vecp(double rho_, double u_, double v_, double w_, double T_): rho(rho_), u(u_), v(v_), w(w_), T(T_) {}
    vecp& operator+=(const vecp& o) {rho += o.rho;u += o.u;v += o.v;w += o.w;T += o.T;return *this;}
    vecp& operator-=(const vecp& o) {rho -= o.rho;u -= o.u;v -= o.v;w -= o.w;T -= o.T;return *this;}
    vecp& operator*=(double k) {rho *= k;u *= k;v *= k;w *= k;T *= k;return *this;}
    double& operator[](int i) {return (&rho)[i];}
    double operator[](int i) const {return (&rho)[i];}
    vec3 vel() const {return vec3{u, v, w};}
    void setvel(const vec3& q) {u = q.x;v = q.y;w = q.z;}
    void clear();
};

inline void vecp::clear() {rho = 0.0;u = 0.0;v = 0.0;w = 0.0;T = 0.0;}
inline vecp operator+(vecp a, const vecp& b) {return a += b;}
inline vecp operator-(vecp a, const vecp& b) {return a -= b;}
inline vecp operator*(vecp a, double s) {return a *= s;}
inline vecp operator*(double s, vecp a) {return a *= s;}

// 无粘分量
struct vec5 {
    double c = 0.0,x = 0.0,y = 0.0,z =0.0 ,e =0.0;
    vec5() = default;
    vec5(double c_, double x_, double y_, double z_, double e_): c(c_), x(x_), y(y_), z(z_), e(e_) {}
    vec5& operator+=(const vec5& o) {c += o.c;x += o.x;y += o.y;z += o.z;e += o.e;return *this;}
    vec5& operator-=(const vec5& o) {c -= o.c;x -= o.x;y -= o.y;z -= o.z;e -= o.e;return *this;}
    vec5& operator*=(double k) {c *= k;x *= k; y *= k;z *= k;e *= k;return *this;}
    double& operator[](int i) {return (&c)[i];}
    double operator[](int i) const {return (&c)[i];}
    void clear();
};

inline void vec5::clear() {c = 0.0;x = 0.0;y = 0.0;z = 0.0;e = 0.0;}
inline vec5 operator+(vec5 a, const vec5& b) {return a += b;}
inline vec5 operator-(vec5 a, const vec5& b) {return a -= b;}
inline vec5 operator*(vec5 a, double s) {return a *= s;}
inline vec5 operator*(double s, vec5 a) {return a *= s;}

// 梯度向量
struct vecgrad {
    vec3 rhograd, ugrad, vgrad, wgrad, Tgrad;
    vecgrad() = default;
    vecgrad(vec3 rho_, vec3 u_, vec3 v_, vec3 w_, vec3 t_): rhograd(rho_), ugrad(u_), vgrad(v_), wgrad(w_), Tgrad(t_) {}
    void clear();
};

inline void vecgrad::clear() {rhograd.clear();ugrad.clear();vgrad.clear();wgrad.clear();Tgrad.clear();}
inline vecp operator*(const vecgrad& a, const vec3& r) {
    return vecp(dot(a.rhograd, r), dot(a.ugrad, r), dot(a.vgrad, r), dot(a.wgrad, r),
                dot(a.Tgrad, r));
}
inline vecgrad operator*(const vecp& phyf, const vec3& fnor) {
    return vecgrad(phyf.rho * fnor, phyf.u * fnor, phyf.v * fnor, phyf.w * fnor, phyf.T * fnor);
}
inline vecgrad operator*(double a, const vecgrad& grad) {
    return vecgrad(a * grad.rhograd, a * grad.ugrad, a * grad.vgrad, a * grad.wgrad,
                   a * grad.Tgrad);
}
inline vecgrad operator+(const vecgrad& a, const vecgrad& b) {
    return vecgrad(a.rhograd + b.rhograd, a.ugrad + b.ugrad, a.vgrad + b.vgrad,
                   a.wgrad + b.wgrad, a.Tgrad + b.Tgrad);
}

// 引申物理量
struct otphy;

// 五阶方阵
struct mat5 {
    double a[5][5] = {};
    mat5() = default;
    mat5(const vec5& r0, const vec5& r1, const vec5& r2, const vec5& r3, const vec5& r4);
    mat5& operator+=(const mat5& o);
    mat5& operator*=(double k);
    double* operator[](int r);
    const double* operator[](int r) const;
    void clear();
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
inline mat5 operator+(mat5 m, const mat5& o);
// 矩阵数乘
inline mat5 operator*(mat5 m, double k);
// 矩阵数乘
inline mat5 operator*(double k, mat5 m);
// 矩阵乘向量
inline vec5 operator*(const mat5& m, const vec5& v);

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

inline double* mat5::operator[](int r) { return a[r]; }

inline void mat5::clear() {
    for (int r = 0; r < 5; r++)
        for (int c = 0; c < 5; c++) a[r][c] = 0.0;
}

inline bool mat5::solve(vec5& b) const {
    double m[5][6];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) m[i][j] = a[i][j];
        m[i][5] = (&b.c)[i];
    }
    for (int k = 0; k < 5; k++) {
        int p = k;
        for (int i = k + 1; i < 5; i++) {
            if (std::abs(m[i][k]) > std::abs(m[p][k])) p = i;
        }
        if (p != k) {
            for (int j = 0; j < 6; j++) {
                const double t = m[k][j];
                m[k][j] = m[p][j];
                m[p][j] = t;
            }
        }
        if (std::abs(m[k][k]) < 1e-14) m[k][k] = m[k][k] >= 0.0 ? 1e-14 : -1e-14;
        for (int i = k + 1; i < 5; i++) {
            const double f = m[i][k] / m[k][k];
            for (int j = k; j < 6; j++) m[i][j] -= f * m[k][j];
        }
    }
    double x[5];
    for (int i = 4; i >= 0; i--) {
        double s = m[i][5];
        for (int j = i + 1; j < 5; j++) s -= m[i][j] * x[j];
        x[i] = s / m[i][i];
    }
    b = vec5(x[0], x[1], x[2], x[3], x[4]);
    return true;
}

inline const double* mat5::operator[](int r) const { return a[r]; }

inline mat5 operator+(mat5 m, const mat5& o) { return m += o; }inline mat5 operator*(mat5 m, double k) { return m *= k; }

inline mat5 operator-(mat5 m, const mat5& o) {
    for(int r = 0; r < 5; r++)
        for(int c = 0; c < 5; c++) m.a[r][c] -= o.a[r][c];
    return m;
}

inline mat5 operator*(double k, mat5 m) { return m *= k; }

inline vec5 operator*(const mat5& m, const vec5& v) {
    return vec5(m.a[0][0] * v.c + m.a[0][1] * v.x + m.a[0][2] * v.y + m.a[0][3] * v.z + m.a[0][4] * v.e,
                m.a[1][0] * v.c + m.a[1][1] * v.x + m.a[1][2] * v.y + m.a[1][3] * v.z + m.a[1][4] * v.e,
                m.a[2][0] * v.c + m.a[2][1] * v.x + m.a[2][2] * v.y + m.a[2][3] * v.z + m.a[2][4] * v.e,
                m.a[3][0] * v.c + m.a[3][1] * v.x + m.a[3][2] * v.y + m.a[3][3] * v.z + m.a[3][4] * v.e,
                m.a[4][0] * v.c + m.a[4][1] * v.x + m.a[4][2] * v.y + m.a[4][3] * v.z + m.a[4][4] * v.e);
}

// 单位阵
inline mat5 eyes = []{ mat5 I; for(int r = 0; r < 5; r++) I.a[r][r] = 1.0; return I; }();

// 三向无粘通量Jacobi
struct mat5tri {
    mat5 x, y, z;
    mat5tri() = default;
    mat5tri(mat5 x_, mat5 y_, mat5 z_) : x(x_), y(y_), z(z_) {}
};

}

// SA
namespace SA{

// 湍流变量
struct turbulence; 

} 

// SST
namespace SST {
struct turbulence;
}

// JST专用
namespace JST {
struct dissipation {
double Y = 0.0;   // 激波捕捉因子
cc::vec5 L;       // 伪Laplace
cc::vec5 Fd;      // JST耗散
};

}

// Roe专用
namespace ROE {
struct dissipation;
}

// LSCB的预留声明
namespace cc {
struct LSCBmatrix;
}

// 流场边界条件
namespace boundary_settings {
struct FARconfig{
    double p;
    double Ma;
    double T;
    double AOA;
    double chi;
};
}

// 支持的边界条件
namespace cc {
inline constexpr short VIL = 1;
inline constexpr short POL = 2;
inline constexpr short INTER = 3;
}

// 伪时间步参数
namespace fatime {
inline double CFL = 1.0;
}