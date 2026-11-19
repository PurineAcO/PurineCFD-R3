#pragma once

#include "data_structure/DataStructure.h"

namespace ph {

inline constexpr double gamma = 1.4;                      // 气体绝热常数
inline constexpr double R = 287.05;                       // 气体常数R
inline constexpr double Cp = gamma*R/(gamma-1.0);         // 定压热容
inline constexpr double Cv = R/(gamma-1.0);               // 恒容热容
inline constexpr double Pr = 0.71;                        // 普朗特数
inline constexpr double Ts = 110.4;                       // Sutherland温度
inline constexpr double T0 = 273.15;                      // 参考温度
inline constexpr double mu0 = 1.716e-5;                   // 参考粘度

// Sutherland粘度
double dynamic_viscosity(double T);
// Sutherland粘度对温度的导数
double _MUT(double mu, double T);
// 形成声速
double get_sonic_velocity(double T);
// 形成能量
double get_energy(const cc::vecp& phy);
// 形成引申物理量
cc::otphy form_otphy(const cc::vecp& phy);
// 动压
double dynamic_pressure(const cc::vecp& refp);
// 壁面力系数
void wall_forces(double& cl, double& cd, const cc::vecp& refp);
// 检查坏单元
int bad_cells();

}
