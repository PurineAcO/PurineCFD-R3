#pragma once

#include "data_structure/DataStructure.h"

namespace SA {

// 形成面上粘性通量
void diffusion_SA(cc::face_class& face);
// 组装单元粘性通量
void assemble_visflux(cc::cell_class& cell);

}
