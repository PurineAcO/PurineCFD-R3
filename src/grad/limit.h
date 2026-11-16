#include "data_structure/cell_class.h"
#include "data_structure/vectors.h"

// BarthJespersen 限制器
cc::vecp BarthJespersen(const cc::cell_class* cell, const cc::vec3& off);

// 2D路径的Venkatakrishnan限制器, 与R2保持一致
cc::vecp BarthJespersen_2D(const cc::cell_class* cell, const cc::vec3& off);