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
struct otphy {
    double a = 0.0;  // 声速
    double p = 0.0;  // 压力
    double e = 0.0;  // 单位质量总能量
    double un = 0.0; // 法向速度
    double mu = 0.0; // 层流粘性系数
    otphy() = default;
    otphy(double a_, double p_, double e_) : a(a_), p(p_), e(e_) {}
    // 从原始物理量生成引申物理量
    void form_otphy(const vecp& phy);
};

} 

// SA
namespace SA{

// 湍流变量
struct turbulence {
    double miubl = 0.0;        // SA工作变量
    double miubl_former = 0.0; // 上一物理时间miubl
    double sad = 0.0;          // 到最近壁面中点的距离
    cc::vec3 miublgrad;            // miubl梯度
    double mueff = 0.0;        // 有效粘性系数
    double chi = 0.0;          // 湍流粘度比
    double fv1 = 0.0;          // 粘度阻尼函数
    double ft2 = 0.0;          // 转捩控制函数
    double fv2 = 0.0;          // 涡量辅助函数
    double Omega = 0.0;        // 涡量
    double Sbl = 0.0;          // 修正涡量
    double S2 = 0.0;           // 应变率张量模
    double r = 0.0;            // 无量纲距离I
    double g = 0.0;            // 无量纲距离II
    double fw = 0.0;           // 破坏项控制函数
    double diag = 0.0;         // SA方程隐式对角
};

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
struct LSCBmatrix {
    double Axx = 0.0;    // 对角xx
    double Axy = 0.0;    // 上三角xy
    double Axz = 0.0;    // 上三角xz
    double Ayy = 0.0;    // 对角yy
    double Ayz = 0.0;    // 上三角yz
    double Azz = 0.0;    // 对角zz
    double dxi[6] = {};  // delta x
    double dyi[6] = {};  // delta y
    double dzi[6] = {};  // delta z
    double wi[6] = {};   // 范数权
    // 解Ag=b
    bool solve(const double b[3], double g[3]) const;
};
}

// 设置区域
namespace cc {
// 支持的边界条件
inline constexpr short FAR = 0;
inline constexpr short VIL = 1;
inline constexpr short POL = 2;
inline constexpr short INTER = 3;
// 湍流模型开关
#define SSTchosen 1
#define SAchosen 0
// 耗散开关
#define ROEchosen 0
#define JSTchosen 1
// URANS选项
#define URANS 1;
}