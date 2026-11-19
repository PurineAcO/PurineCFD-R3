#pragma once

#include "data_structure/DataStructure.h"

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// MUSCL左右插值
void muscl(cc::face_class* face);
