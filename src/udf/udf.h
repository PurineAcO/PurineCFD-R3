#pragma once
#include "../data_structure/vectors.h"

namespace cc {
// 湍流模型开关
#define SSTchosen 0
#define SAchosen 1
// 耗散开关
#define ROEchosen 0
#define JSTchosen 1
// URANS选项
#define URANS 1;
// 梯度选项
#define LSCBchosen 1
#define GGCBchosen 0
// 无粘通量格式
inline char scheme = 'R';
// 是否瞬态
inline bool urans = false;
// OpenMP线程数, 0表示自动选择
inline int threads = 0;
// 边界条件
inline boundary_settings::FARconfig FAR_config = {22.9,0.2,300,0.0,0.1};
}