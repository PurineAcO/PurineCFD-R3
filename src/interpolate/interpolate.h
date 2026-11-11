#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"

// 面上中心差分插值
void interpolate_mid(cc::face_class* face);
// MUSCL左右插值
void muscl(cc::face_class* face);
// 面上梯度的处理
void grad_onface(cc::cell_class& cell);