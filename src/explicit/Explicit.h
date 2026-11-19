#pragma once

#include "implicit/Implicit.h"
#include "residual/Residual.h"

namespace impl {

// 进行一次显式伪时间推进
double flow_step_explicit(double dt);

}

namespace SA {

// 进行一次湍流模型显式伪时间推进
void turbulence_step_explicit(double dt);

}
