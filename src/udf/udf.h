#pragma once
#include "../data_structure/vectors.h"

namespace cc {
// 湍流模型开关
#define SSTchosen 0
#define SAchosen 1
// 耗散开关
#define ROEchosen 1
#define JSTchosen 0
// 梯度选项
#define LSCBchosen 1
// OpenMP线程数, 0表示自动选择
inline int threads = 32;
// 边界条件
inline boundary_settings::FARconfig FAR_config = {22.9,0.2,300,0.0,0.1};
}

// 求解器参数
namespace solver {
// 网格文件
inline const char* mesh = "cases/cyc/mesh.txt";
// 最大步数
inline int max_step = 100;
// CFL数
inline double cfl = 50.0;
// 湍流开关, 0为开启
inline int turboff = 0;
// 对称Gauss-Seidel扫描次数
inline int sweeps = 4;
// 残差与屏幕输出间隔
inline int conv_step = 10;
// 流场输出间隔
inline int field_interval = 100;
// 输出目录, 空则由程序按 3D-<n>-<tag> 自动编号
inline const char* outdir = "";
// 输出目录标签, 空则由网格名自动生成
inline const char* tag = "cyc-urans-rf";
}

// 双时间步/URANS参数
namespace dualtime {
// 物理时间步长, 不大于0表示定常
inline double dt = 1.5e-3;
// 内迭代最大次数
inline int inner = 20;
// 内迭代残差下降指标
inline double inner_tol = 1e-3;
// URANS之前的定常迭代次数
inline int steady_iters = 0;
// 定常迭代的CFL终值
inline double steady_cfl = 50.0;
// 初始扰动幅值, 占Uinf的比例
inline double seed = 0.02;
}