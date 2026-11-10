#pragma once
#include "data_structure/vectors.h"


// 标准初始化
void standard_initialization();
// 从流场初始化
void fromflowfield_initialization();
// 从边界条件形成物理量
cc::vecp from_boundary_to_phy();
