#pragma once

#include <cstdio>

#include "data_structure/DataStructure.h"
#include "physics/Physics.h"

// 输出流场
bool dump_field(const char* path);
// 输出壁面压力系数
bool write_wall_cp(const char* path, double pref, const cc::vecp& refp);
// 输出中截面壁面压力系数(逐步追加)
bool write_wall_cp_mid(FILE* fp, int step, double pref, const cc::vecp& refp);
