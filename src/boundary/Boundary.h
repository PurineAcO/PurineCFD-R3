#pragma once

#include "data_structure/DataStructure.h"

// 远场边界更新
void far_field_boundary();
// 无滑移不穿透壁面
void noslip_wall_boundary();
// 对称滑移面
void symmetry_wall_boundary();
// 远场无粘通量Jacobian（定义在 implicit/）
cc::mat5 far_convect_jacobi(const cc::face_class* face);
