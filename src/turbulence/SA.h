#pragma once
#include "../data_structure/vectors.h"
#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include <vector>

struct SA::turbulence
{
    double miubl = 0.0;        // SA工作变量
    double miubl_former = 0.0; // 上一物理时间miubl
    double sad = 0.0;          // 到最近壁面中点的距离
    cc::vec3 miublgrad;        // miubl梯度
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

namespace SA {
inline constexpr double Cb1 = 0.1355;
inline constexpr double Cb2 = 0.622;
inline constexpr double Cw1 = 3.2391;
inline constexpr double Cw2 = 0.3;
inline constexpr double Cw3 = 2.0;
inline constexpr double Cv1 = 7.1;
inline constexpr double Ct3 = 1.2;
inline constexpr double Ct4 = 0.5;
inline constexpr double kappa = 0.41;
inline constexpr double inv_sigma = 1.5;
inline constexpr double rmax = 10.0;
inline constexpr double relax = 0.8;
inline constexpr double Prt = 0.9;
inline constexpr double C5 = 3.5;

// 形成单元湍流中间变量
void build_SA(cc::cell_class& c);
// 形成面湍流中间变量
void build_SA(cc::face_class& f);
// 形成面上粘性通量
void diffusion_SA(cc::face_class& face);
// 组装单元粘性通量
void assemble_visflux(cc::cell_class& cell);
// 形成SA源项
double source_SA(cc::cell_class& cell);
// 面上有效扩散系数
double Gamma_eff(cc::cell_class& c, cc::face_class* f);
// 沿面法向的特征距离
double normal_dist(const cc::face_class* f);
// 粘性通量Jacobi(冻结状态, 梯度变化) - 对应R2的diffusion_grad_jac
cc::mat5 diffusion_grad_jac(const cc::face_class& face, double dn);
// 粘性通量Jacobi(冻结梯度, 状态变化) - 对应R2的form_pfpgq_diffusion+form_pfpw_from_pfpgq
cc::mat5 diffusion_state_jac(const cc::face_class& face);
// SA源项对nu的解析Jacobi对角
void source_jacobi_diag_SA(cc::cell_class& cell, double inv_dt, double cfl);
}

namespace SA {
inline double _chi(double rho, double mu, double miubl){return rho*(miubl > 0.0 ? miubl : 0.0)/mu;}
inline double _fv1(double chi){return chi*chi*chi/(chi*chi*chi + SA::Cv1*SA::Cv1*SA::Cv1);}
inline double _ft2(double chi){return SA::Ct3*std::exp(-SA::Ct4*chi*chi);}
inline double _fv2(double chi){return 1.0 - chi/(1.0 + chi*_fv1(chi));}
inline double _g(double r){return r + SA::Cw2*(r*r*r*r*r*r - r);}
inline double _fw(double g){
    constexpr double cw3_sixth = SA::Cw3*SA::Cw3*SA::Cw3*SA::Cw3*SA::Cw3*SA::Cw3;
    const double g_squared = g*g;
    const double g_sixth = g_squared*g_squared*g_squared;
    return g*std::pow((1.0 + cw3_sixth)/(cw3_sixth + g_sixth), 1.0/6.0);
}
}


