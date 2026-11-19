#pragma once

#include "data_structure/DataStructure.h"

#include "gradient/Limit.h"

// 网格的LSCB梯度
void least_square_cell_based(cc::cell_class& cell);
// 湍流变量的LSCB
void LSCB_turbulence(cc::cell_class& cell);
// 面上退化的梯度
void grad_onface(cc::cell_class& cell);
// LSCB梯度预处理
void LSCB_pre_process(cc::cell_class& c);
