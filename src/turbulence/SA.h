#pragma once
#include "../data_structure/vectors.h"

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