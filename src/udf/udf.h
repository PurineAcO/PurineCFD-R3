#pragma once
#include "../data_structure/vectors.h"

namespace cc {
// 2D还是3D
#define chose2D 1
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
#if chose2D
inline boundary_settings::FARconfig FAR_config = {63585.0302134262,0.72,300.0,4.5,3.0};
#else
inline boundary_settings::FARconfig FAR_config = {22.9,0.2,300,0.0,0.1};
#endif
}

// 求解器参数
namespace solver {
// 网格文件
#if chose2D
inline const char* mesh = "cases/naca0012/mesh.txt";
#else
inline const char* mesh = "cases/cyc/mesh.txt";
#endif
// 最大步数
#if chose2D
inline int max_step = 10000;
#else
inline int max_step = 100;
#endif
// CFL数
#if chose2D
inline double cfl = 100.0;
#else
inline double cfl = 50.0;
#endif
// 湍流开关, 0为开启
inline int turboff = 0;
// 对称Gauss-Seidel扫描次数
inline int sweeps = 4;
// 内迭代格式, 0为隐式LU-SGS, 1为显式伪时间推进
inline int inner_scheme = 0;
// 定常预热的内迭代格式, 0为隐式LU-SGS, 1为显式伪时间推进
inline int steady_scheme = 0;
// 显式伪时间推进的CFL数
inline double exp_cfl = 0.3;
// 残差与屏幕输出间隔
inline int conv_step = 20;
// 流场输出间隔
inline int field_interval = 2;
// 流场输出起始步, 小于该步不输出
inline int field_start_step = 100000000;
// 输出目录, 空则由程序按 3D-<n>-<tag> 自动编号
inline const char* outdir = "D:/PurineCFD-R3/run/NACA0012-Ma0.72-a4.5-Re1e7";
// 输出目录标签, 空则由网格名自动生成
#if chose2D
inline const char* tag = "naca0012";
#else
inline const char* tag = "cyc-urans-rf";
#endif
}

// 双时间步/URANS参数
namespace dualtime {
// 物理时间步长, 不大于0表示定常
inline double dt = 0.0;
// 内迭代最大次数
inline int inner = 40;
// 内迭代残差下降指标
inline double inner_tol = 1e-3;
// URANS之前的定常迭代次数
inline int steady_iters = 0;
// 定常迭代的CFL终值
inline double steady_cfl = 100.0;
// 初始扰动幅值, 占Uinf的比例
inline double seed = 0.0;
}
