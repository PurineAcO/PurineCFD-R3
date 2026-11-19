#pragma once

#include "data_structure/DataStructure.h"

// 标准初始化
// 标准初始化仅仅针对 cell，不针对任何 face 和边界面
void standard_initialization();
// 从流场初始化
// 当前无定义
void fromflowfield_initialization();
// 从边界条件形成物理量
cc::vecp from_boundary_to_phy();
