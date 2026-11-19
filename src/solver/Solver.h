#pragma once

#include <omp.h>

#include "data_structure/DataStructure.h"

using SolveKernel = void (*)(int c, int lo, int hi);

// 循环器GS
void blockSgs(SolveKernel solve_kernel);
