#pragma once
#include "data_structure/face_class.h"
#include "data_structure/vectors.h"

namespace cc {
inline constexpr short FAR = 0;
}

// 远场边界更新
void far_field_boundary();
// 远场无粘通量Jacobian
cc::mat5 far_convect_jacobi(const cc::face_class* face);
