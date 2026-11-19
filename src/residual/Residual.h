#pragma once

#include "data_structure/DataStructure.h"

namespace impl {

// 形成右侧残差
void right_residual();

}

namespace SA {

// 单元SA残差
double right_residual(cc::cell_class& cell);

}
