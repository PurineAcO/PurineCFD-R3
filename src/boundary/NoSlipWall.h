#pragma once
#include "data_structure/vectors.h"

namespace cc {
inline constexpr short WALL = 5;
}

// 无滑移不穿透壁面
void noslip_wall_boundary();
