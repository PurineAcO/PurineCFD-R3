#pragma once

#include "data_structure/DataStructure.h"

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
// 形成SA源项
double source_SA(cc::cell_class& cell);
// 面上有效扩散系数
double Gamma_eff(cc::cell_class& c, cc::face_class* f);
// 沿面法向的特征距离
double normal_dist(const cc::face_class* f);
// 湍流粘度比
double _chi(double rho, double mu, double miubl);
// 粘度阻尼函数
double _fv1(double chi);
// 转捩控制函数
double _ft2(double chi);
// 涡量辅助函数
double _fv2(double chi);
// 无量纲距离II
double _g(double r);
// 破坏项控制函数
double _fw(double g);

}

// 涡量模
double curl_abs(const cc::vecgrad& g);
// 应变率张量模
double strain2(const cc::vecgrad& g);
