#pragma once

#include "data_structure/DataStructure.h"

// 限制因子（q 在 [lo,hi] 内、增量为 dq 时返回 0..1 的缩放系数）
double shrink(double q, double lo, double hi, double dq);

// BarthJespersen 限制器
cc::vecp BarthJespersen(const cc::cell_class* cell, const cc::vec3& off);

// 2D路径的Venkatakrishnan限制器, 与R2保持一致
cc::vecp BarthJespersen_2D(const cc::cell_class* cell, const cc::vec3& off);
